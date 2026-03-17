# 参考了https://blog.csdn.net/swangee/article/details/141646920
rm(list = ls())
setwd('./')

if(!dir.exists("./Processed"))  dir.create("./Processed")
if(!dir.exists("./Processed/RawMatrix"))  dir.create("./Processed/RawMatrix")
if(!dir.exists("./Processed/RawData"))  dir.create("./Processed/RawData")
if(!dir.exists("./Processed/RawData/csv"))  dir.create("./Processed/RawData/csv")
if(!file.exists("./Data/gdc_download_20260315_132413.267760.tar.gz")){
  warnings('源数据不存在，请下载https://github.com/Data708983/Difference_Analysis/releases/download/source/gdc_download_20260315_132413.267760.tar.gz到Data文件夹中\n')
}
tar_file = "./Data/gdc_download_20260315_132413.267760.tar.gz"
extract_dir = "./Processed/RawMatrix"
untar(tar_file, exdir = extract_dir)

# 整理
rm(list = ls())
setwd('./')
library(data.table)
library(dplyr)
sample_sheet <- fread("./Data/gdc_sample_sheet.2026-03-15.tsv")# 读取样本信息
sample_sheet$Barcode <- substr(sample_sheet$`Sample ID`,1,15) # 取ID前15字符作为barcode
## 有关TGCA_ID含义：
## > 关于数字（前15位）
## 1. https://zhuanlan.zhihu.com/p/572657767
## 2. https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
## > 关于字母（16位）
## 3. https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
## 可以看出TGCA的ID的14、15位决定了样本性质，而16位代表同一个患者的取样顺序
sample_sheet1 <- sample_sheet %>% filter(!duplicated(sample_sheet$Barcode)) # 去重
sample_sheet2 <- sample_sheet1 %>% filter(grepl("01$|11$|02$",sample_sheet1$Barcode)) # Barcode的最后两位：01表示肿瘤样本，11表示正常样本，02复发实体瘤


TCGA_GBM_Exp <- fread("./Processed/RawMatrix/00cce172-4b2b-4cb5-b6dd-47c9627972af/511393a1-da49-4cf1-80d6-7f50a6c3043d.rna_seq.augmented_star_gene_counts.tsv") # 任意读取一个文件
TCGA_GBM_Exp <- TCGA_GBM_Exp[-1:-4,c("gene_id","gene_name","gene_type")] # 删1~4行，保留包含"gene_id","gene_name","gene_type"的数据框，用于合并表达数据

## 将所有样本合并成一个数据框
## >这个算法直观但是太原始和低效了，建议使用现成的包和函数
for (i in 1:nrow(sample_sheet2)) {
  
  folder_name <- sample_sheet2$`File ID`[i]
  file_name <- sample_sheet2$`File Name`[i]
  sample_name <- sample_sheet2$Barcode[i]
  
  data1 <- fread(paste0("./Processed/RawMatrix/",folder_name,"/",file_name))
  #unstranded代表count值；如果要保存TPM，则改为tpm_unstranded
  data2 <- data1[-1:-4,c("gene_id","gene_name","gene_type","unstranded")] 
  colnames(data2)[4] <- sample_name
  
  TCGA_GBM_Exp <- inner_join(TCGA_GBM_Exp,data2)
}

# 根据需要的表达比例筛选满足条件的基因
zero_percentage <- rowMeans(TCGA_GBM_Exp[, 4:ncol(TCGA_GBM_Exp)] == 0) 
TCGA_GBM_Exp1 <- TCGA_GBM_Exp[zero_percentage < 0.6, ] # 筛选出表达超过60%的基因

install.packages("BiocManager")
library(BiocManager)
BiocManager::install("edgeR")
library(edgeR)
TCGA_GBM_Exp1 = avereps(TCGA_GBM_Exp1[,-c(1:3)],ID = TCGA_GBM_Exp$gene_name) # 对重复基因名取平均表达量，并将基因名作为行名
TCGA_GBM_Exp1 <- TCGA_GBM_Exp1[rowMeans(TCGA_GBM_Exp1)>100,] # 根据需要去除低表达基因，这里设置的平均表达量100为阈值

# 创建样本分组
library(stringr)
tumor <- colnames(TCGA_GBM_Exp1)[substr(colnames(TCGA_GBM_Exp1),14,15) == "01"]
normal <- colnames(TCGA_GBM_Exp1)[substr(colnames(TCGA_GBM_Exp1),14,15) == "11"]
tumor_sample <- TCGA_GBM_Exp1[,tumor]
normal_sample <- TCGA_GBM_Exp1[,normal]
exprSet_by_group <- cbind(tumor_sample,normal_sample)
gene_name <- rownames(exprSet_by_group)
exprSet <- cbind(gene_name, exprSet_by_group)  # 将gene_name列设置为数据框的行名，合并后又添加一列基因名


# 存储counts和TPM数据
fwrite(exprSet,"./Processed/RawData/TCGA_GBM_Count.txt") # txt格式
write.csv(exprSet, "./Processed/RawData/csv/TCGA_GBM_Count.csv", row.names = FALSE) # csv格式

