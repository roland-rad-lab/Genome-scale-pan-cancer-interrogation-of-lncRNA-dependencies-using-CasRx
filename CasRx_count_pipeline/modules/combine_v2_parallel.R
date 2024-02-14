library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
######PARALLELIZATION ARGUMENTS########
#get number of cores
n_cores <- floor(parallel::detectCores()/2)
#create the cluster
my_cluster <- parallel::makeCluster(
   n_cores,
   #use PSOCK because of portability on windows
   type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my_cluster)
#create %notchin%
`%notchin%` <- Negate(`%chin%`)
#args section
args = commandArgs(trailingOnly=T)
sample_table_file <- args[1]
combinations_table_file <- args[2]
out_dir <- args[3]
counts_dir <- paste0(out_dir, "/counts")
mapped_dir <- paste0(out_dir, "/mapped")
dir.create(counts_dir, recursive = T)

sample_table <- fread(sample_table_file)
combinations_table <- fread(combinations_table_file)
combinations_table_char <- combinations_table[, c(combinations)]
column_names_R1 <- c("query", "flag", "reference_R1", "pos", "mapq", "cigar", "r_next", "p_next", "tlen", "seq", "qual", "na1", "na2", "na3", "na4")
column_names_R2 <- c("query", "flag", "reference_R2", "pos", "mapq", "cigar", "r_next", "p_next", "tlen", "seq", "qual", "na1", "na2", "na3", "na4")

 rows <- rownames(sample_table)
 foreach(i=1:length(rows), .packages = c("data.table", "tidyverse")) %dopar% {
   row <- rows[i]
   row <- as.integer(row)
   fwd_name = basename(as.character(sample_table[row, 3]))
   mapped_fwd_name <- paste0(sub('\\.fastq.gz$', '', fwd_name), ".sam")
   #print(fwd_name)
   mapped_fwd_full_name <- paste0(mapped_dir, "/", mapped_fwd_name)
   mappings_R1 <- fread(mapped_fwd_full_name, col.names = column_names_R1)
   ##log generation
   log_fwd_name <- paste0(mapped_dir, "/", sub('\\.fastq.gz$', '', fwd_name), ".log")
   log_fwd <- read_tsv(log_fwd_name, col_names = "info") %>%
   head(n = 3) %>%
   separate(info, into = c("info", "number"), sep = ":") %>%
   mutate(info = gsub("# ", "", info)) %>%
   mutate(info = gsub(" ", "_", info)) %>%
   mutate(info = paste0(info, "_R1")) %>%
   mutate(number = trimws(number)) %>%
   mutate(number = gsub("\\s.*$", "", number))
   if("R2" %in% colnames(sample_table)){
      rev_name = basename(as.character(sample_table[row, 4]))
      mapped_rev_name <- paste0(sub('\\.fastq.gz$', '', rev_name), ".sam")
      mapped_rev_full_name <- paste0(mapped_dir, "/", mapped_rev_name)
      mappings_R2 <- fread(mapped_rev_full_name, col.names = column_names_R2)
      ##log generation
      log_rev_name <- paste0(mapped_dir, "/", sub('\\.fastq.gz$', '', rev_name), ".log")
      log_rev <- read_tsv(log_rev_name, col_names = "info") %>%
      head(n = 3) %>%
      separate(info, into = c("info", "number"), sep = ":") %>%
      mutate(info = gsub("# ", "", info)) %>%
      mutate(info = gsub(" ", "_", info)) %>%
      mutate(info = paste0(info, "_R2")) %>%
      mutate(number = trimws(number)) %>%
      mutate(number = gsub("\\s.*$", "", number))
   } else {
      print(paste0("ERROR: column R2 not found in file ", sample_table_file, ". Please make sure to have a correct sample table."))
      stop()
   }
   mappings <- merge(mappings_R1, mappings_R2, all = TRUE, by = "query")
   mappings[, combinations := paste(reference_R1, reference_R2, sep = "~")]
   length_mapped_all <- nrow(mappings)
   off_map <- mappings[combinations %notchin% combinations_table_char]
   mappings <- mappings[combinations %chin% combinations_table_char]
   mappings[, copies := .N, by = query]
   mappings <- mappings[copies == 1]
   mappings <- mappings[, -c("copies")]
   length_mapped_correct_comb <- nrow(mappings)
   new_col <- c("mapping_to_the_same_guide", length_mapped_correct_comb)
   log <- rbind(log_fwd, log_rev, new_col)
   log_file_name <- paste0(counts_dir, "/", sample_table[row, 1],".log")
   write_tsv(log, log_file_name)
   mappings <- mappings[, .(count = .N), by = combinations]
   mappings <- merge(mappings, combinations_table, all = T, by = "combinations")
   mappings[is.na(mappings)] <- 0
   mappings <- mappings %>% separate(combinations, into = c("reference_1", "reference_2"), sep = "~")
   output_file <- paste0(counts_dir, "/", sample_table[row, 1], ".txt")
   write_tsv(mappings, output_file)
   #arrange off-maps
   off_map <- off_map[, .(count = .N), by = combinations]
   off_map <- off_map %>% separate(combinations, into = c("reference_1", "reference_2"), sep = "~")
   off_map <- off_map[reference_1 != "NA" & reference_2 != "NA"]
   output_file_off <- paste0(counts_dir, "/", sample_table[row, 1], "_off-map.txt")
   write_tsv(off_map, output_file_off)
}
stopCluster(my_cluster)
