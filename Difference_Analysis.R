# 差异表达分析- edgeR包
# 参考了https://blog.csdn.net/swangee/article/details/141646920
rm(list = ls())
setwd('./')
load('./Processed/RawData/rda/Difference.rda')

library(edgeR)
library(data.table)

# 创建一个字符data.table# 创建一个字符向量 group_list，其长度等于总样本数
# 前 ncol(tumor_sample) 个元素为 "tumor"，后 ncol(normal_sample) 个元素为 "normal"
group_list <- c(rep('tumor',ncol(tumor_sample)),rep('normal',ncol(normal_sample)))
# 将字符向量转换为因子（factor），R 会自动按字母顺序设定水平（levels）
# 若需要指定顺序，可使用 factor(group_list, levels = c("normal", "tumor"))
group_list <- factor(group_list)
# 构建设计矩阵（design matrix），用于线性模型拟合
# ~0 + group_list 表示无截距项，为每个组别单独估计均值
# 结果是一个两列的矩阵，列名即为组的水平（如 group_listnormal 和 group_listtumor）
# 每行对应一个样本，相应组的位置为 1，否则为 0
design <- model.matrix(~0+group_list)
# 将设计矩阵的行名设为样本名（即表达矩阵的列名），确保样本与设计矩阵的行一一对应。
rownames(design) <- colnames(exprSet_by_group)
# 将设计矩阵的列名简化为组名本身（去掉前缀 group_list），便于后续代码阅读和对比设置。
# 即变为 "normal" 和 "tumor"
colnames(design) <- levels(group_list) # 创建分组

# 创建一个DGEList对象，用于存储差异表达分析所需的数据
DGElist <- DGEList(counts = exprSet_by_group, group = group_list)

# 使用cpm（counts per million，每百万计数）值对低表达量的基因进行过滤
keep_gene <- rowSums(cpm(DGElist) > 1 ) >= 2
# DGElist[keep_gene, ] 按行（基因）筛选，只保留 keep_gene 为 TRUE 的基因
# 参数 keep.lib.sizes = FALSE 表示不重新计算文库大小（即总计数）
# 默认情况下，edgeR 会在过滤后自动更新每个样本的文库大小（总计数），
# 但设置为 FALSE 则保留原始文库大小，
# 避免因过滤导致文库大小失真（尤其在差异分析后续的标准化中）
# 通常保持默认或设为 FALSE 都可以，具体取决于分析流程。
DGElist <- DGElist[keep_gene, keep.lib.sizes = FALSE] # 保留符合条件的基因

# 校正测序深度
DGElist <- calcNormFactors( DGElist )

# 估算离散度，该步骤会对DGElist进行添加或更新，一般看CommonDisp就行
DGElist <- estimateGLMCommonDisp(DGElist, design) # 共同离散度
DGElist <- estimateGLMTrendedDisp(DGElist, design) # 趋势离散度
DGElist <- estimateGLMTagwiseDisp(DGElist, design) # 基因特异的离散度

# 拟合广义线性模型
# 为每个基因拟合一个负二项分布的广义线性模型（GLM）。
# >输入：DGElist 包含原始计数和离散度信息；design 是之前构建的无截距设计矩阵（列对应 normal 和 tumor）
# >输出：fit 对象存储了每个基因的模型拟合结果，包括回归系数等，用于后续假设检验
fit <- glmFit(DGElist, design)
# 似然比检验，contrast = c(-1, 1)即对分组的两个条件检验，这里是tumor和normal
# c(-1, 1)正值表示肿瘤中表达上调，负值表示下调
results <- glmLRT(fit, contrast = c(-1, 1)) 


# 提取差异表达的top基因
GBM_nrDEG_edgeR <- topTags(results, n = nrow(DGElist)) # n = nrow(DGElist)即全部保存
GBM_nrDEG_edgeR <- as.data.frame(GBM_nrDEG_edgeR)

# 保存原始差异基因矩阵文件
fwrite(GBM_nrDEG_edgeR,"./Processed/RawData/GBM_nrDEG_edgeR.txt", row.names = TRUE) # csv格式
write.csv(GBM_nrDEG_edgeR, "./Processed/RawData/csv/GBM_nrDEG_edgeR.csv", row.names = T) # csv格式

# 筛选差异基因
library(data.table)
library(dplyr)
library(tibble)
GBM_Match_DEG <- GBM_nrDEG_edgeR # 或者从输出的文件里读取
GBM_Match_DEG$log10FDR <- -log10(GBM_Match_DEG$FDR)

GBM_Match_DEG <- rownames_to_column(GBM_Match_DEG, var = "gene_name")

GBM_Match_DEG <- GBM_Match_DEG %>% 
  mutate(DEG = case_when(logFC > 2 & FDR < 0.05 ~ "Up",
                         abs(logFC) < 2 | FDR > 0.05 ~ "None",
                         logFC < -2 & FDR < 0.05 ~ "Down")) # 打标签：logFC > 2 & FDR < 0.05：上调基因，logFC < -2 & FDR < 0.05：下调基因，其它认为无显著差异

# 保存添加标签后的基因
fwrite(GBM_Match_DEG,"./Processed/RawData/GBM_Match_DEG.txt") # txt格式
write.csv(GBM_Match_DEG, "./Processed/RawData/csv/GBM_Match_DEG.csv", row.names = F) # csv格式

# 可视化之火山图 (Volcano Plot) 
GBM_Match_DEG <- fread("./Processed/RawData/GBM_Match_DEG.txt") # 读取打完标签的基因列表
GBM_Match_DEG <- as.data.frame(GBM_Match_DEG)
rownames(GBM_Match_DEG) <- GBM_Match_DEG[,1]  # 将第一列作为行名

down_gene <- GBM_Match_DEG[GBM_Match_DEG$DEG == "Down", ]
up_gene <- GBM_Match_DEG[GBM_Match_DEG$DEG == "Up", ]


uptop <- rownames(up_gene)[1:10]  # 上调的前10基因
downtop <- rownames(down_gene)[1:10] # 下调的前10基因

GBM_Match_DEG$label <- ifelse(GBM_Match_DEG$gene_name %in% c(uptop,downtop), GBM_Match_DEG$gene_name, "") # 后面画图时用来突出显著表达的前10个基因


# 加载需要用到的程序包
library(data.table)
library(ggplot2)
if (!requireNamespace("ggprism", quietly = TRUE)) install.packages("ggprism")
library(ggprism)
library(ggrepel)

# 画图 volcano plot
p <- ggplot(GBM_Match_DEG, aes(x = logFC, y = log10FDR, colour = DEG)) +
  geom_point(alpha = 0.85, size = 1.5) + # 设置点的透明度和大小
  scale_color_manual(values = c('steelblue', 'gray', 'brown')) + # 调整点的颜色
  xlim(c(-11, 11)) +  # 调整x轴的范围
  geom_vline(xintercept = c(-2, 2), lty = 4, col = "black", lwd = 0.8) + # x轴辅助线
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) + # y轴辅助线
  labs(x = "logFC", y = "-log10FDR") + # x、y轴标签
  ggtitle("Difference Expression Genes of GBM") +  # 图表标题
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.title = element_blank()) +  # 设置图表标题和图例位置
  geom_label_repel(data = GBM_Match_DEG, aes(label = label),  # 添加标签
                   size = 3, box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black",
                   show.legend = FALSE, max.overlaps = 10000) +  # 标签设置
  theme_prism(border = TRUE)
pdf(file = './Processed/vocano.pdf', width = 10, height = 6)
p
dev.off()

# 热图
# 加载所需包
library(pheatmap)
library(RColorBrewer)

# --- 1. 准备数据 ---
# 提取显著差异的基因（这里沿你火山图中的阈值: |logFC| > 2 & FDR < 0.05）
sig_genes <- GBM_Match_DEG[GBM_Match_DEG$DEG %in% c("Up", "Down"), ]

# 可以根据需要选择展示的基因数量，例如：
# 1) 全部显著基因（如果数量不多）
# 2) 前50个上调 + 前50个下调（按logFC排序）
# 这里以选择前50上调、前50下调为例：
sig_genes_up <- sig_genes[sig_genes$DEG == "Up", ]
sig_genes_down <- sig_genes[sig_genes$DEG == "Down", ]

# 按logFC绝对值排序，选取前N个
N <- 10
top_up <- rownames(sig_genes_up[order(sig_genes_up$logFC, decreasing = TRUE), ])[1:min(N, nrow(sig_genes_up))]
top_down <- rownames(sig_genes_down[order(sig_genes_down$logFC, decreasing = FALSE), ])[1:min(N, nrow(sig_genes_down))]
selected_genes <- c(top_up, top_down)
selected_genes <- selected_genes[!is.na(selected_genes)] # 去除可能存在的NA

# 提取这些基因的表达矩阵（使用原始的、经过过滤和标准化的计数数据？还是logCPM？）
# 通常热图使用log2(CPM+1)或log2(TPM+1)来展示，使表达分布更接近正态，视觉效果更好。
# 我们可以从DGElist对象中获取logCPM值。
logCPM <- cpm(DGElist, log = TRUE, prior.count = 1) # 获取log2(CPM + prior.count)

# 筛选出我们选中的基因
heatmap_data <- logCPM[rownames(logCPM) %in% selected_genes, ]

# 确保基因顺序与selected_genes一致
heatmap_data <- heatmap_data[match(selected_genes, rownames(heatmap_data)), ]

# --- 2. 准备样本注释信息 ---
# 创建一个数据框，用于标注样本的组别（Tumor/Normal）
annotation_col <- data.frame(
  Group = group_list
)
rownames(annotation_col) <- colnames(heatmap_data)

# 定义注释的颜色
ann_colors <- list(
  Group = c(tumor = "brown", normal = "steelblue") # 与火山图颜色保持一致
)

# --- 3. 绘制热图并保存 ---
# 对表达矩阵进行行（基因）缩放，使每个基因在不同样本中的表达均值为0，标准差为1。
# 这样可以更好地展示基因在样本间的相对表达高低（上调/下调）。
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# 设置表达颜色的渐变
color_gradient <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# 绘制热图
pdf(file = './Processed/heatmap_topDEGs.pdf', width = 12, height = 10)
pheatmap(heatmap_data_scaled,
         cluster_rows = TRUE,        # 对基因聚类
         cluster_cols = TRUE,        # 对样本聚类
         show_rownames = FALSE,       # 基因太多时，行名显示为FALSE，若想显示可改为TRUE
         show_colnames = FALSE,       # 样本名显示
         annotation_col = annotation_col, # 添加样本分组注释条
         annotation_colors = ann_colors,
         color = color_gradient,
         main = "Top 10 Up-regulated and Top 10 Down-regulated Genes",
         fontsize = 8,
         border_color = NA)           # 去掉单元格边框，使图形更清爽
dev.off()
message("热图已保存至 ./Processed/heatmap_topDEGs.pdf")

