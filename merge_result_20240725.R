#!/usr/bin/env Rscript
# Usage: 
# Rscript merge_result.R sampleId pathogen_type k2_resut_file bwt2_result_file blast_result_file min_unique_reads min_unique_reads min_unique_reads
# Rscript merge_result.R SRR12486989 bacteria SRR12486989.k2.bacteria.filtered.xls SRR12486989.bacteria.bwt2.xls SRR12486989.bacteria.blast.xls 50 20 1

# rm(list = ls())

args = commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
pathogen_type <- args[2]
k2_resut_file <- args[3]
bwt2_result_file <- args[4]
blast_result_file <- args[5]
genome_infor_file <- args[6]
min_unique_reads <- as.numeric(args[7])
min_size_abundance <- as.numeric(args[8])
min_coverage <- as.numeric(args[9])
clinical_infor_file <- args[10]


# sample_id <- "SRR12486971"
# pathogen_type <- "fungi"
# k2_resut_file <- "SRR12486971.k2.fungi.filtered.xls"
# bwt2_result_file <- "SRR12486971.fungi.bwt2.xls"
# blast_result_file <- "SRR12486971.fungi.blast.xls"
# genome_infor_file <- "fungi.genome.infor.xls"
# min_unique_reads <- 5
# min_size_abundance <- 20
# min_coverage <- 0.05
# clinical_infor_file <- "clinical.infor.xls"

library(tidyverse)
library(data.table)

output_file_all <- paste0(sample_id, ".", pathogen_type, ".final.all.xls")
output_file_filtered <- paste0(sample_id, ".", pathogen_type, ".final.filtered.xls")

print(paste0("sample_id: ", sample_id))
print(paste0("pathogen_type: ", pathogen_type))
print(paste0("k2_resut_file: ", k2_resut_file))
print(paste0("bwt2_result_file: ", bwt2_result_file))
print(paste0("blast_result_file: ", blast_result_file))
print(paste0("genome_infor_file: ", genome_infor_file))
print(paste0("min_unique_reads: ", min_unique_reads))
print(paste0("min_size_abundance: ", min_size_abundance))
print(paste0("min_coverage: ", min_coverage))
print(paste0("output_file_all: ", output_file_all))
print(paste0("output_file_filtered: ", output_file_filtered))
print(paste0("clinical_infor_file: ", clinical_infor_file))

k2_input <- read_tsv(k2_resut_file) %>% dplyr::select(-kmers)
bwt2_input <- read_tsv(bwt2_result_file)
blast_result_input <- read_tsv(blast_result_file)
genome_infor <- read_tsv(genome_infor_file)
clinical_infor <- read_tsv(clinical_infor_file) %>% dplyr::select(-organism_name)

merge_res <- k2_input %>% inner_join(bwt2_input, by = c("sample_id", "taxid")) %>%
  inner_join(blast_result_input, by = c("sample_id", "taxid")) %>%
  mutate(coverage_chrom = round(mapped_chrom_count / scaffold_count * 100, 2),
         size_abundance = round(blast_reads / genome_size * 1000000),
         score = round(size_abundance * (coverage_chrom / 100) * (avg_coverage / 100), 2)
         ) %>%
  dplyr::select(taxid,
                scientific_name,
                accession_version = accession,
                assembly_accession,
                # pathogen_type,
                species_level,
                k2_reads_count,
                bwt2_reads_count = bwt2_reads,
                blast_reads_count = blast_reads,
                k2_unique_kmers = unique_kmers,
                size_abundance,
                score,
                base_coverage = avg_coverage,
                base_depth = avg_depth,
                relative_abundance,
                genome_size,
                sequences_count = scaffold_count,
                mapped_sequences_count = mapped_chrom_count,
                sequences_coverage = coverage_chrom,
                avg_pident,
                min_pident,
                max_pident
                ) %>%
  arrange(desc(score), desc(size_abundance), desc(base_coverage), desc(blast_reads_count), desc(k2_unique_kmers))


# 同一个种内，只保留分数最高的那个
genome_infor2 <- genome_infor %>% dplyr::select(taxid, species_taxid)
merge_res_distinct <- merge_res %>% inner_join(genome_infor2, by = "taxid") %>%
  distinct(species_taxid, .keep_all = T) %>% dplyr::select(-species_taxid)

#add clinical infor
merge_res_distinct2 <- merge_res_distinct %>% left_join(clinical_infor, by = "taxid")


# add fiter
output_all <- merge_res_distinct2 %>% mutate(filter = if_else(blast_reads_count >= min_unique_reads &
                                                          size_abundance >= min_size_abundance &
                                                          base_coverage >= min_coverage,
                                                          "passed",
                                                          "unpassed"))

output_filtered <- output_all %>% filter(filter == "passed")
write_tsv(output_all, output_file_all, na = "")
write_tsv(output_filtered, output_file_filtered, na = "")

