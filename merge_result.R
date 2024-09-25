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
min_NPA <- as.numeric(args[8])
min_base_coverage <- as.numeric(args[9])
min_sequence_coverage <- as.numeric(args[10])
sample_reads <- as.numeric(args[11])


# sample_id <- "SRR12486971"
# pathogen_type <- "fungi"
# k2_resut_file <- "SRR12486971.k2.fungi.filtered.xls"
# bwt2_result_file <- "SRR12486971.fungi.bwt2.xls"
# blast_result_file <- "SRR12486971.fungi.blast.xls"
# genome_infor_file <- "fungi.genome.infor.xls"
# min_unique_reads <- 5
# min_NPA <- 10
# min_base_coverage <- 0.05
# min_sequence_coverage <- 1
# sample_reads <- 2127108

library(tidyverse)
library(data.table)

output_file_all <- paste0(sample_id, ".", pathogen_type, ".final.xls")
# output_file_filtered <- paste0(sample_id, ".", pathogen_type, ".final.filtered.xls")

print(paste0("sample_id: ", sample_id))
print(paste0("pathogen_type: ", pathogen_type))
print(paste0("k2_resut_file: ", k2_resut_file))
print(paste0("bwt2_result_file: ", bwt2_result_file))
print(paste0("blast_result_file: ", blast_result_file))
print(paste0("genome_infor_file: ", genome_infor_file))
print(paste0("min_unique_reads: ", min_unique_reads))
print(paste0("min_NPA: ", min_NPA))
print(paste0("min_base_coverage: ", min_base_coverage))
print(paste0("min_sequence_coverage: ", min_sequence_coverage))
print(paste0("output_file_all: ", output_file_all))
# print(paste0("output_file_filtered: ", output_file_filtered))
print(paste0("sample_reads: ", sample_reads))

k2_input <- read_tsv(k2_resut_file)
bwt2_input <- read_tsv(bwt2_result_file)
blast_result_input <- read_tsv(blast_result_file)
genome_infor <- read_tsv(genome_infor_file)
# clinical_infor <- read_tsv(clinical_infor_file) %>% dplyr::select(-organism_name)

merge_res <- k2_input %>% inner_join(bwt2_input, by = c("sample_id", "taxid")) %>%
  inner_join(blast_result_input, by = c("sample_id", "taxid")) %>%
  mutate(coverage_chrom = round(mapped_chrom_count / scaffold_count * 100, 2),
         sample_total_reads = sample_reads,
         NPA = round(1e9 * unique_kmers * blast_reads / (genome_size * sample_total_reads), 4),
         NPAS = round(log2(NPA * (coverage_chrom / 100) * (avg_coverage / 100) + 1), 4)
         ) %>%
  dplyr::select(taxid,
                scientific_name,
                accession_version = accession,
                assembly_accession,
                pathogen_type,
                # species_level,
                sample_total_reads,
                unique_reads = blast_reads,
                bwt2_reads_count = bwt2_reads,
                k2_reads_count,
                k2_kmers = kmers,
                k2_unique_kmers = unique_kmers,
                k2_relative_abundance = relative_abundance,
                NPA,
                NPAS,
                `base_coverage(%)` = avg_coverage,
                base_depth = avg_depth,
                genome_size,
                sequences_count = scaffold_count,
                mapped_sequences_count = mapped_chrom_count,
                `sequence_coverage(%)` = coverage_chrom,
                `avg_similarity(%)` = avg_pident,
                `min_similarity(%)` = min_pident,
                `max_similarity(%)` = max_pident
                ) %>%
  arrange(desc(NPAS), desc(NPA), desc(`base_coverage(%)`), desc(unique_reads), desc(k2_unique_kmers))


# Keep only the highest-scoring one within the same species
genome_infor2 <- genome_infor %>% dplyr::select(taxid, species_taxid)
merge_res_distinct <- merge_res %>% inner_join(genome_infor2, by = "taxid") %>%
  distinct(species_taxid, .keep_all = T) %>% dplyr::select(-species_taxid)

# add clinical infor
# merge_res_distinct2 <- merge_res_distinct %>% left_join(clinical_infor, by = "taxid")


# add fiter
output_all <- merge_res_distinct %>% mutate(passed_filter = if_else(unique_reads >= min_unique_reads &
                                                          NPA >= min_NPA &
                                                          `base_coverage(%)` >= min_base_coverage &
                                                          `sequence_coverage(%)` >= min_sequence_coverage,
                                                          "positive",
                                                          "negative"))

write_tsv(output_all, output_file_all, na = "")
# output_filtered <- output_all %>% filter(passed_filter == "passed")
# write_tsv(output_filtered, output_file_filtered, na = "")

