#!/usr/bin/env Rscript
# Usage: 
# Rscript blast_result_processed.R sampleId taxid  read1_file read2_file pathogen_type
# Rscript blast_result_processed.R SRR12486989 1311 SRR12486989.1311.R1.blast.txt  SRR12486989.1311.R2.blast.txt virus

# rm(list = ls())

args = commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
taxid <- args[2]
read1_file <- args[3]
read2_file <- args[4]
pathogen_type <- args[5]

# sample_id <- "SRR12486988"
# taxid <- "36330"
# read1_file <- "SRR12486988.36330.R1.blast.txt"
# read2_file <- "SRR12486988.36330.R2.blast.txt"
# pathogen_type <- "virus"

library(tidyverse)
library(data.table)

print(paste0("sample_id: ", sample_id))
print(paste0("taxid: ", taxid))
print(paste0("read1_file: ", read1_file))
print(paste0("read2_file: ", read2_file))

read1_input <- read_tsv(read1_file)
read2_input <- read_tsv(read2_file)

# Align to multiple locations, keep one of them
read1_res <- read1_input %>% group_by(qseqid, sseqid) %>% arrange(desc(bitscore), desc(pident)) %>% 
  distinct(qseqid, .keep_all = T) %>% dplyr::select(qseqid, sseqid, pident_r1 = pident) %>% ungroup()
  
read2_res <- read2_input %>% group_by(qseqid, sseqid) %>% arrange(desc(bitscore), desc(pident)) %>% 
  distinct(qseqid, .keep_all = T) %>% dplyr::select(qseqid, sseqid, pident_r2 = pident) %>% ungroup()

read_res <- read1_res %>% inner_join(read2_res, by = c("qseqid","sseqid")) %>%
  mutate(avg_pident_reads = (pident_r1 + pident_r2) / 2)

## 如何匹配不上，就输结果为0
if (nrow(read_res) == 0 ) {
  blast_out <- tibble(sample_id = sample_id,
                      taxid = taxid,
                      blast_reads = 0,
                      avg_pident = 0,
                      min_pident = 0,
                      max_pident = 0)
  write_tsv(blast_out, paste0(sample_id, ".", taxid, ".blast.res.txt"))
  q(save = 'no')
  quit(save = "no", status = 0)
  }

blast_out1 <- read_res %>% dplyr::select(avg_pident_reads) %>%
  summarise(blast_reads = n(), avg_pident = round(mean(avg_pident_reads), 2)) %>%
  mutate(sample_id = sample_id, taxid = taxid)
  
blast_out2<- read_res %>% summarise(max_pident = round(max(pident_r1,pident_r2), 2),
                                    min_pident = round(min(pident_r1,pident_r2), 2)) %>% 
  mutate(sample_id = sample_id, taxid = taxid)

blast_out <- blast_out1 %>% inner_join(blast_out2, by = c("sample_id", "taxid")) %>% 
  select(sample_id, taxid, blast_reads, avg_pident, min_pident, max_pident)

write_tsv(blast_out, paste0(sample_id, ".", pathogen_type, ".", taxid, ".blast.res.txt"))

