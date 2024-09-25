#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]
sort_column <- args[3]
sort_type <- args[4]

# in_file <- "SRR8476435.pathogens.final.all_temp.xls"
# out_file <- "SRR8476435.pathogens.final.all.xls"
# sort_column <- "score"
# sort_type <- "reverse"

print(paste0("in_file: ", in_file))
print(paste0("out_file: ", out_file))
print(paste0("sort_column: ", sort_column))
print(paste0("sort_type: ", sort_type))



if (sort_type == "reverse") {
  read_tsv(in_file) %>% arrange(desc(!!sym(sort_column))) %>% write_tsv(out_file)
} else {
  read_tsv(in_file) %>% arrange(!!sym(sort_column)) %>% write_tsv(out_file)
}



