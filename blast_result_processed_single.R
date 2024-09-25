#!/usr/bin/env Rscript
# Usage: 
# Rscript blast_result_processed.R sampleId taxid  read1_file pathogen_type
# Rscript blast_result_processed.R SRR12486989 1311 SRR12486989.1311.R1.blast.txt virus

# rm(list = ls())

args = commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
taxid <- args[2]
read1_file <- args[3]
pathogen_type <- args[4]

# sample_id <- "SRR12486971"
# taxid <- "10298"
# read1_file <- "SRR12486971.virus.10298.R1.blast.txt"
# pathogen_type <- "virus"

library(tidyverse)
library(data.table)

print(paste0("sample_id: ", sample_id))
print(paste0("taxid: ", taxid))
print(paste0("read1_file: ", read1_file))

read1_input <- read_tsv(read1_file)

# Align to multiple locations, keep one of them
read1_res <- read1_input %>% group_by(qseqid, sseqid) %>% arrange(desc(bitscore), desc(pident)) %>% 
  distinct(qseqid, .keep_all = T) %>% dplyr::select(qseqid, sseqid, pident_r1 = pident) %>% ungroup()


blast_out <- read1_res %>% summarise(blast_reads = n(),
                                     avg_pident = round(mean(pident_r1), 2),
                                     min_pident = round(min(pident_r1), 2),
                                     max_pident = round(max(pident_r1), 2),
                                     ) %>%
  mutate(sample_id = sample_id, taxid = taxid) %>% 
  select(sample_id, taxid, blast_reads, avg_pident, min_pident, max_pident)


write_tsv(blast_out, paste0(sample_id, ".", pathogen_type, ".", taxid, ".blast.res.txt"))

