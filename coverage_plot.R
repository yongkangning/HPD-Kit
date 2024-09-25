#!/usr/bin/env Rscript
# Usage: 
# Rscript coverage_plot.R sampleId taxid seqid cov_file pathogen_type
# Rscript coverage_plot.R SRR12486989 1311 NZ_CP049938.1 SRR12486989.1311_NZ_CP049938.1_depth.txt virus

# rm(list = ls())

args = commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
taxid <- args[2]
seqid <- args[3]
cov_file <- args[4]
pathogen_type <- args[5]

# sample_id <- "SRR12486989"
# taxid <- "1311"
# seqid <- "NZ_CP049938.1"
# cov_file <- "QH6.1332244_NC_026428.1_depth.txt"
# pathogen_type <- "virus"

library(tidyverse)
library(data.table)

print(paste0("sample_id: ", sample_id))
print(paste0("taxid: ", taxid))
print(paste0("seqid: ", seqid))
print(paste0("cov_file: ", cov_file))

# 读入覆盖度结果
cov_data <- fread(cov_file, sep = "\t", col.names = c("Seqid", "Position", "Depth"))
# 对测序深度进行SymLog转换log1p(x) = log(1+x)
cov_data$coverage <- log1p(cov_data$Depth)
# 计算测序深度和覆盖度
coverage <- round(cov_data[, sum(Depth > 0) / .N] * 100, 4)
depth <- round(cov_data[, sum(Depth) / .N], 4)
text_infor <-  paste0("Avg Coverage:", coverage, "%  Avg Depth:", depth, "X")

# 待完成 怎么像文献一样计算区间内的覆盖度
# cov_data[, coverage := frollsum(Depth, n = 100, align = 'right', fill = 0) / 100]


p <- ggplot(cov_data, aes(x = Position, y = coverage)) +
  geom_line(color = "skyblue") +
  theme_classic() +
  labs(title = seqid, subtitle = text_infor, x = "Reference Position", y = "Coverage(SymLog)") +
  scale_x_continuous(expand = c(0,0), labels = scales::comma) + # 不用科学计数法
  scale_y_continuous(expand = c(0,0)) + # y轴和x轴从原点开始
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 1, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 1, size = 6),
    axis.title = element_text(face = "bold"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    axis.title.y = element_text(margin = margin(r = 10))
  )
ggsave(paste0(sample_id, ".", pathogen_type, ".", taxid, "_", seqid, "_cov.png"), width = 8, height = 4, dpi = 200)
