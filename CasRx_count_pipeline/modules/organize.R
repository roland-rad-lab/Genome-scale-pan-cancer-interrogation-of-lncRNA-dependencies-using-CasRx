library(data.table)
library(tidyverse)
#args section
args = commandArgs(trailingOnly=T)
 sample_table_file <- args[1]
 combinations_table_file <- args[2]
 out_dir <- args[3]
counts_dir <- paste0(out_dir, "/counts")
mapped_dir <- paste0(out_dir, "/mapped")



sample_table <- fread(sample_table_file)
combinations_table <- fread(combinations_table_file)
#merge different count tables together in one table
#get path of count tables
count_tables <- lapply(sample_table[,1], function(x){paste0(counts_dir, "/", x, ".txt")}) %>% unlist() %>% as.list()
#generate the initial table with combinations but no counts
full_table <- fread(combinations_table_file) %>% separate(combinations, into = c("reference_1", "reference_2"), sep = "~")
#add each count as a column of the full table
for (count_table_file in count_tables) {
   count_table <- fread(count_table_file) 
   sample_name <- gsub(pattern = "\\.txt$", "", count_table_file) %>% basename()
   setnames(count_table, old = "count", new = sample_name)
   full_table <- merge(full_table, count_table, all = T, by = c("reference_1", "reference_2"))
}
output_path <- paste0(out_dir, "/raw_counts.txt")
write_tsv(full_table, output_path)

#create mapping log table
log_full_table <- data.table(info = c("reads_processed_R1", 
                                      "reads_with_at_least_one_reported_alignment_R1", 
                                      "reads_that_failed_to_align_R1",
                                      "reads_processed_R2",
                                      "reads_with_at_least_one_reported_alignment_R2",
                                      "reads_that_failed_to_align_R2",
                                      "mapping_to_the_same_guide"
                                      ))
#get path of log files
log_tables <- lapply(sample_table[,1], function(x){paste0(counts_dir, "/", x, ".log")}) %>% unlist() %>% as.list()
#add each count as a column of the full table
for (log_table_file in log_tables) {
   log_table <- fread(log_table_file) 
   sample_name <- gsub(pattern = "\\.log$", "", log_table_file) %>% basename()
   setnames(log_table, old = "number", new = sample_name)
   log_full_table <- merge(log_full_table, log_table, all = T, by = "info")
}
output_log <- paste0(out_dir, "/log.txt")
write_tsv(log_full_table, output_log)

#############