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


# 分别初始化 counts 和 TPM 的主数据框
TCGA_GBM_Exp_counts <- fread("./Processed/RawMatrix/00cce172-4b2b-4cb5-b6dd-47c9627972af/511393a1-da49-4cf1-80d6-7f50a6c3043d.rna_seq.augmented_star_gene_counts.tsv") # 任意读取一个文件
TCGA_GBM_Exp_counts <- TCGA_GBM_Exp_counts[-1:-4, c("gene_id","gene_name","gene_type")] # 删1~4行，保留基因注释列
TCGA_GBM_Exp_tpm <- TCGA_GBM_Exp_counts # 结构相同，复制一份用于TPM

## 将所有样本合并成两个数据框（counts 和 TPM）
## >这个算法直观但是太原始和低效了，建议使用现成的包和函数
for (i in 1:nrow(sample_sheet2)) {
  
  folder_name <- sample_sheet2$`File ID`[i]
  file_name <- sample_sheet2$`File Name`[i]
  sample_name <- sample_sheet2$Barcode[i]
  
  data1 <- fread(paste0("./Processed/RawMatrix/",folder_name,"/",file_name))
  
  # 提取 counts 列（unstranded）
  data_counts <- data1[-1:-4, c("gene_id","gene_name","gene_type","unstranded")] 
  colnames(data_counts)[4] <- sample_name
  TCGA_GBM_Exp_counts <- inner_join(TCGA_GBM_Exp_counts, data_counts)
  
  # 提取 TPM 列（tpm_unstranded）
  data_tpm <- data1[-1:-4, c("gene_id","gene_name","gene_type","tpm_unstranded")] 
  colnames(data_tpm)[4] <- sample_name
  TCGA_GBM_Exp_tpm <- inner_join(TCGA_GBM_Exp_tpm, data_tpm)
}

# 根据需要的表达比例筛选满足条件的基因（基于 counts）
zero_percentage <- rowMeans(TCGA_GBM_Exp_counts[, 4:ncol(TCGA_GBM_Exp_counts)] == 0) 
keep_genes <- zero_percentage < 0.6
TCGA_GBM_Exp_counts_filt <- TCGA_GBM_Exp_counts[keep_genes, ] # 筛选出表达超过60%的基因
TCGA_GBM_Exp_tpm_filt    <- TCGA_GBM_Exp_tpm[keep_genes, ]    # 同步过滤TPM

# 安装并加载 edgeR（用于 avereps）
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(BiocManager)
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
library(edgeR)

# 对重复基因名取平均表达量（分别处理 counts 和 TPM）
# 根据需要去除低表达基因，这里设置的平均表达量100为阈值（基于 counts）
# counts 矩阵
expr_counts = avereps(TCGA_GBM_Exp_counts_filt[,-c(1:3)],ID = TCGA_GBM_Exp_counts_filt$gene_name)
expr_counts_final <- expr_counts[rowMeans(expr_counts)>100,]

# TPM 矩阵
expr_tpm = avereps(TCGA_GBM_Exp_tpm_filt[,-c(1:3)],ID = TCGA_GBM_Exp_tpm_filt$gene_name)
expr_tpm_final <- expr_tpm[rowMeans(expr_tpm)>100,]

# 创建样本分组（使用 counts 的列名，TPM 列名相同）
library(stringr)
tumor <- colnames(expr_counts_final)[substr(colnames(expr_counts_final),14,15) == "01"]
normal <- colnames(expr_counts_final)[substr(colnames(expr_counts_final),14,15) == "11"]
tumor_sample_counts <- expr_counts_final[, tumor]
normal_sample_counts <- expr_counts_final[, normal]
tumor_sample_tpm <- expr_tpm_final[, tumor]
normal_sample_tpm <- expr_tpm_final[, normal]

# 合并肿瘤和正常样本（列顺序：先肿瘤后正常）
exprSet_counts_by_group <- cbind(tumor_sample_counts, normal_sample_counts)
exprSet_tpm_by_group    <- cbind(tumor_sample_tpm,    normal_sample_tpm)

# 添加基因名列
gene_name_count <- rownames(exprSet_counts_by_group)
gene_name_tpm <- rownames(exprSet_tpm_by_group)
exprSet_counts <- cbind(gene_name_count, as.data.frame(exprSet_counts_by_group))
exprSet_tpm    <- cbind(gene_name_tpm, as.data.frame(exprSet_tpm_by_group))

# 存储 counts 和 TPM 数据
fwrite(exprSet_counts, "./Processed/RawData/TCGA_GBM_Count.txt") # txt格式
write.csv(exprSet_counts, "./Processed/RawData/csv/TCGA_GBM_Count.csv", row.names = FALSE) # csv格式

fwrite(exprSet_tpm, "./Processed/RawData/TCGA_GBM_TPM.txt")
write.csv(exprSet_tpm, "./Processed/RawData/csv/TCGA_GBM_TPM.csv", row.names = FALSE)
