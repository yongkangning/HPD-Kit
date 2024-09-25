#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
treat_file <- args[1]
control_file <- args[2]
control_npa_method <- args[3]

# treat_file <- "SRR12486971.pathogens.final.xls"
# control_file <- "SRR12486971_vs_control.tsv"

sample_id = str_split(treat_file, '\\.')[[1]][1]

print(paste0("treat_file: ", treat_file))
print(paste0("control_file: ", control_file))
print(paste0("sample_id: ", sample_id))

treat <- read_tsv(treat_file)
control <- read_tsv(control_file)

# max or mean
control_mean <- control %>% group_by(taxid) %>% summarise(NPA_control = if_else(control_npa_method == "max", max(NPA_control), mean(NPA_control)))


treat_ajusted <- treat %>% left_join(control_mean, by = "taxid") %>% 
  mutate(NPA_control = if_else(is.na(NPA_control), 0, NPA_control),
         # NPA_adjusted = if_else(NPA - NPA_control > 0, NPA - NPA_control, 0),
         NPA_adjusted = NPA / (NPA_control + 1),
         NPAS_adjusted = round(log2(NPA_adjusted * (`base_coverage(%)` / 100) * (`sequence_coverage(%)` / 100) + 1), 4)
         ) %>% dplyr::select(-NPA_control) %>% relocate(ends_with("_adjusted"), .after = NPAS) %>% 
  arrange(desc(NPAS_adjusted))

write_tsv(treat_ajusted, paste0(sample_id, "_vs_control.pathogens.final.xls"))

# split by pathogen type
path_db <- treat_ajusted$pathogen_type %>% unique()
for (db in path_db) {
  path_res <- treat_ajusted %>% filter(pathogen_type == db) %>% arrange(desc(NPAS_adjusted))
  write_tsv(path_res, paste0(sample_id, "_vs_control.", db, ".final.xls"))
}
