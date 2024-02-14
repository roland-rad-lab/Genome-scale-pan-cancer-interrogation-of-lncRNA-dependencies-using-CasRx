library(data.table)
library(tidyverse)

#script to run MAGeCK RRA algorithm (mageck test) alongg all sample/library pairs
`%notin%` <- Negate(`%in%`)
args = commandArgs(trailingOnly = T)
#get names of samples
results_dir <- args[1]
output_name <- args[2]
table_file <- args[3]
dir.create(results_dir, showWarnings = F)
setwd(results_dir)
table <- fread(table_file)
table <- table[, .SD, .SDcols = colnames(table) %notin% c("sgRNA", "Gene", "Library_rep1")]
pool <- "Library_rep1"
samples <- colnames(table)
samples <- gsub("*_rep[1-9]$", "", samples)
samples <- unique(samples)

 for (sample in samples) {
   cmd <- paste0("mageck test -k ", table_file, " -t ", sample, "_rep1,", sample, "_rep2",  " -c ", pool, " -n ", output_name, "_", sample, " --norm-method none")
   print(cmd)
   system(cmd)
 }

gene_table <- data.table("target" = character(), "RRA_score" = numeric(), "drop_pvalue" = numeric(), "enrich_pvalue" = numeric(), "drop_fdr" = numeric(), "enrich_fdr" = numeric(), "lfc" = numeric(), "sample" = character())
###merge results in one table
for (sample in samples) {
  rra_table_path <- paste0(results_dir, "/", output_name, "_", sample, ".gene_summary.txt")
  RRA_table <- fread(rra_table_path)
  RRA_table <- RRA_table[, c("id", "neg|score", "neg|p-value", "pos|p-value", "neg|fdr", "pos|fdr", "neg|lfc")]
  setnames(RRA_table, old = c("id", "neg|score", "neg|fdr", "pos|fdr", "neg|lfc", "neg|p-value", "pos|p-value"), new = c("target", "RRA_score", "drop_fdr", "enrich_fdr", "lfc", "drop_pvalue", "enrich_pvalue"))
  RRA_table[, sample := sample]
  gene_table <- rbind(gene_table, RRA_table)
}
write_tsv(gene_table, paste0(results_dir, "/gene_table.tsv"))


sgRNA_table <- data.table("sgrna" = character(), "Gene" = character(), "lfc" = numeric(), p_value = numeric(), "fdr" = numeric(), "sample" = character())
###merge sgrna results in one table
for (sample in samples) {
  rra_table_path <- paste0(results_dir, "/", output_name, "_", sample, ".sgrna_summary.txt")
  RRA_table <- fread(rra_table_path)
  RRA_table <- RRA_table[, c("sgrna", "Gene", "LFC", "p.low", "FDR")]
  setnames(RRA_table, old = c("LFC", "FDR", "p.low"), new = c("lfc", "fdr", "p_value"))
  RRA_table[, sample := sample]
  sgRNA_table <- rbind(sgRNA_table, RRA_table)
}
write_tsv(sgRNA_table, paste0(results_dir, "/sgRNA_table.tsv"))
