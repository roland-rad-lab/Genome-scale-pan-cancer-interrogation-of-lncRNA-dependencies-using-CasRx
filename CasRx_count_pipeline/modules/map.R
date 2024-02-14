library(tidyverse)
library(parallel)
#args section
args = commandArgs(trailingOnly=T)
sample_table_file <- args[1]
reverse_complement_1 <- args[2]
reverse_complement_2 <- args[3]
read_1_pos <- args[4]
read_2_pos <- args[5]
path_index_1 <- args[6]
path_index_2 <- args[7]
trimmed_dir <- args[8]
output_dir <- args[9]
mismatches <- args[10]
n_threads <- detectCores()-1

##########FUNCTIONS##########
get_index_from_pos <- function(position){
  if (position == "first") {
    index <- path_index_1
  } else if (position == "second") {
    index <- path_index_2
  } else {
    print("ERROR: position not recognised, allowed positions are 'first' and 'second'")
    stop()
  }
  return(index)
}
select_samtools_strand <- function(reverse_complement){
  if (reverse_complement == "T") {
    view_strand <- "-f 16"
  } else if (reverse_complement == "F") {
    view_strand <- "-F 16"
  } else {
    print("ERROR: reverse complement information not recognised, allowed options are 'T', 'TRUE', 'F', 'FALSE'")
    stop()
  }
  return(view_strand)
}
################################
sample_table <- read_tsv(sample_table_file)

for(row in rownames(sample_table)){
  row <- as.integer(row)
  fwd_name = basename(as.character(sample_table[row, 3]))
  fwd_index_path <- get_index_from_pos(read_1_pos)
  fwd_samtools_strand_opt <- select_samtools_strand(reverse_complement_1)
  trimmed_fwd_name <- paste0(sub('\\.fastq.gz$', '', fwd_name), ".trimmed.fastq.gz")
  trimmed_fwd_full_name <- paste0(trimmed_dir, "/", trimmed_fwd_name)
  output_file_fwd <- paste0(output_dir, "/", sub('\\.trimmed.fastq.gz$', '', trimmed_fwd_name), ".sam")
  map_info_file_name_fwd <- sub('\\.sam$', '.log', output_file_fwd) 
  if (file.exists(trimmed_fwd_full_name)) {
   map_fwd_cmd <- paste('bowtie -S -a -v', mismatches, '--threads', n_threads, fwd_index_path, trimmed_fwd_full_name, "2>", map_info_file_name_fwd, "| samtools view -F 4 -h | samtools view", fwd_samtools_strand_opt, ">", output_file_fwd, sep = " ")
   #report all alignments
   print(map_fwd_cmd)
   system(map_fwd_cmd)
  } else {
    print(paste0("ERROR: trimmed file ", trimmed_fwd_full_name, " not found, exiting the script"))
    stop()
  }
  if("R2" %in% colnames(sample_table)){
    rev_name = basename(as.character(sample_table[row, 4]))
    rev_index_path <- get_index_from_pos(read_2_pos)
    rev_samtools_strand_opt <- select_samtools_strand(reverse_complement_2)
    trimmed_rev_name <- paste0(sub('\\.fastq.gz$', '', rev_name), ".trimmed.fastq.gz")
    trimmed_rev_full_name <- paste0(trimmed_dir, "/", trimmed_rev_name)
    output_file_rev <- paste0(output_dir, "/", sub('\\.trimmed.fastq.gz$', '', trimmed_rev_name), ".sam")
    map_info_file_name_rev <- sub('\\.sam$', '.log', output_file_rev)
    if (file.exists(trimmed_rev_full_name)) {
      map_rev_cmd <- paste('bowtie -S -a -v', mismatches, '--threads', n_threads, rev_index_path, trimmed_rev_full_name, "2>", map_info_file_name_rev, "| samtools view -F 4 -h | samtools view", rev_samtools_strand_opt, ">", output_file_rev, sep = " ")
      print(map_rev_cmd)
      system(map_rev_cmd)
    } else {
      print(paste0("ERROR: trimmed file ", trimmed_rev_full_name, " not found, exiting the script"))
      stop()
    }
  }
}
