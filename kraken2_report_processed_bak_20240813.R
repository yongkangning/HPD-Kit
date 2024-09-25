#!/usr/bin/env Rscript
# Usage: 
# Rscript kraken2_report_processed.R sample_id k2_report_file min_k2_reads_count min_relative_abundance min_unique_kmers parasite genome_infor_file
# Rscript kraken2_report_processed.R SRR12486988 SRR12486988.k2.report.txt 15 0.005 1500 parasite bacteria.genome.infor.xls

# rm(list = ls())

args = commandArgs(trailingOnly = TRUE)
library(tidyverse)

sample_id <- args[1]
k2_report_file <- args[2]
min_k2_reads_count <- as.numeric(args[3])
min_relative_abundance <-  as.numeric(args[4])
min_unique_kmers <-  as.numeric(args[5])
min_unique_kmers_fold  <-  as.numeric(args[6])  # unique_kmers / k2_reads_count
pathogen_type <- args[7]   # bacteria virus fungi parasite
genome_infor_file <- args[8]

##########################  剔除假阳性，过滤条件 #############################################################
# 1. k2_reads_count >= 10
# 2. relative_abundance >= 0.005(%)
# 3. unique_kmers >= 800
##############################################################################################################
# sample_id <- "SRR12486989"
# k2_report_file <- "SRR12486989.k2.report.txt"
# min_k2_reads_count <- 10
# min_relative_abundance <-  0.001
# min_unique_kmers <-  800
# min_unique_kmers_fold <-  10
# pathogen_type <- "bacteria"
# genome_infor_file <- "bacteria.genome.infor.xls"

print(paste0("Working dir: ", getwd()))

print(paste0("sample_id: ", sample_id))
print(paste0("k2_report_file: ", k2_report_file))
print(paste0("min_k2_reads_count: ", min_k2_reads_count))
print(paste0("min_relative_abundance: ", min_relative_abundance, "%"))
print(paste0("min_unique_kmers: ", min_unique_kmers))
print(paste0("min_unique_kmers_fold: ", min_unique_kmers_fold))
print(paste0("pathogen_type: ", pathogen_type))
print(paste0("genome_infor_file: ", genome_infor_file))

k2_report <- read_tsv(k2_report_file, col_names = c("Percentage", "accumulate_reads_count", "k2_reads_count",
                                                     "kmers", "unique_kmers","rank", "taxid", "scientific_name"))

genome_infor <- read_tsv(genome_infor_file) %>% dplyr::select(taxid, assembly_accession, scaffold_count, genome_size)

# 界（Kingdom）门（Phylum）纲（Class）目（Order）科（Family）属（Genus）种（Species）
# 0 1 2 3 4 5 6 7 8 9
# U R D K P C O F G S

k2_raw <- k2_report %>%
  mutate(
    sample_id = sample_id,
    pathogen_type = pathogen_type,
    relative_abundance = k2_reads_count / sum(k2_reads_count) * 100,
    scientific_name = str_trim(scientific_name),
    species_level =
      case_when(
        str_starts(rank, "U") ~ "Unclassified",
        str_starts(rank, "R") ~ "Root",
        str_starts(rank, "D") ~ "Kingdom",
        str_starts(rank, "K") ~ "Kingdom",
        str_starts(rank, "P") ~ "Phylum",
        str_starts(rank, "C") ~ "Class",
        str_starts(rank, "O") ~ "Order",
        str_starts(rank, "F") ~ "Family",
        str_starts(rank, "G") ~ "Genus",
        rank == "S" ~ "Species",
        str_starts(rank, "S") & length(rank) > 1 ~ "Strain",
      ),
    rank_num =
      case_when(
        str_starts(rank, "U") ~ str_replace(rank, "U", "0"),
        str_starts(rank, "R") ~ str_replace(rank, "R", "1"),
        str_starts(rank, "D") ~ str_replace(rank, "D", "2"),
        str_starts(rank, "K") ~ str_replace(rank, "K", "3"),
        str_starts(rank, "P") ~ str_replace(rank, "P", "4"),
        str_starts(rank, "C") ~ str_replace(rank, "C", "5"),
        str_starts(rank, "O") ~ str_replace(rank, "O", "6"),
        str_starts(rank, "F") ~ str_replace(rank, "F", "7"),
        str_starts(rank, "G") ~ str_replace(rank, "G", "8"),
        str_starts(rank, "S") ~ str_replace(rank, "S", "9")
      ),
    has_child = rank_num < lead(rank_num),
    has_child = if_else(is.na(has_child), FALSE, TRUE)
    ) %>% 
  dplyr::select(sample_id, taxid, scientific_name, pathogen_type, species_level, rank, rank_num, has_child,
                accumulate_reads_count, k2_reads_count, kmers, unique_kmers, relative_abundance) %>%
  
  mutate(is_passed = if_else(k2_reads_count >= min_k2_reads_count &
                               relative_abundance >= min_relative_abundance &
                               (unique_kmers >= min_unique_kmers | unique_kmers/k2_reads_count >= min_unique_kmers_fold), "passed", "unpassed"))


write_tsv(k2_raw, paste0(sample_id, ".k2.", pathogen_type, ".raw.xls"))

# 只验证有序列的级别

# 只处理科、属、种级别
k2_filtered <- k2_raw %>% filter(is_passed == "passed" & str_sub(rank, 1, 1) %in% c("F", "G", "S")) %>%
  inner_join(genome_infor, by = "taxid") %>% arrange(desc(unique_kmers), desc(k2_reads_count)) %>%
  dplyr::select(-rank, -rank_num, -has_child, -accumulate_reads_count, -is_passed)

k2_result_count <- k2_filtered %>% dplyr::count() %>% pull(n)
print(paste0("[INFO] There are a total of " , k2_result_count, " pathogens that require further verification."))

write_tsv(k2_filtered, paste0(sample_id, ".k2.", pathogen_type, ".filtered.xls"))


# 输出未通过过滤条件的病原体，方便后期检查是否漏检
k2_unpassed <- k2_raw %>% filter(is_passed == "unpassed" & str_sub(rank, 1, 1) %in% c("F", "G", "S")) %>%
  inner_join(genome_infor, by = "taxid") %>% arrange(desc(unique_kmers), desc(k2_reads_count)) %>%
  dplyr::select(-rank, -rank_num, -has_child, -accumulate_reads_count, -is_passed)

write_tsv(k2_unpassed, paste0(sample_id, ".k2.", pathogen_type, ".unpassed.xls"))