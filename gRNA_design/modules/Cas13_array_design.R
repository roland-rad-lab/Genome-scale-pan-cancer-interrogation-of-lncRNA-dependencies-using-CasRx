##########################################################
####CasRx gRNA Design tool to maximize target efficacy####
##########################################################

# Load R packages
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(data.table)
library(Biostrings)
library(biomaRt)
library(dplyr)
library(tidyr)
library(BiocManager)
library(stringr)
library(randomForest)
library(R.utils)
library(zip)
library(seqRFLP)
library(rlist)
library(foreach)
library(doParallel)
library(farver)

array_of_min_dist <- function(ids, dist, guide_length){
  array_dist <- ids %>% 
    select(sequence) %>% 
    mutate(length=str_length(sequence)) %>% 
    mutate(dist_bp=length*dist/100) %>% 
    mutate(dist_bp=dist_bp+23) %>% #add guide length to avoid overlapping when min dist is less than 23(we free the user of thinking about the overlapping)
    select(dist_bp)
  return(array_dist$dist_bp)
}

array_of_min_dist_global <- function(ids, guide_length, dist_global){
  array_dist_global <- ids %>% 
    select(sequence) %>% 
    mutate(length=str_length(sequence)) %>% 
    mutate(global_dist_bp=length*dist_global/100) %>% 
    mutate(global_dist_bp=global_dist_bp+23) %>% #add guide length to avoid overlappings when min dist is less than 23(we free the user of thinking about the overlapping)
    select(global_dist_bp)
  return(array_dist_global$global_dist_bp)
}
#input: 
#number of guides
#number of spacers
#guide length
#extension
#minimum distance
#minimum distance global
#junction flag(avoid overlapping junctions)
#input type (1=ENSEMBLE ID, 2=sequence, 3=regions)
#organism
#file input(must match with type of input)
#sequence range

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

input_file <- args[1]
organism <- "human" #either human or mouse, useful if you give ENSEMBL IDs or Ranges to retrieve correct sequence
n_guides <- args[2] %>% as.numeric() #number of guides to design
# n_guides <- 2
n_spacers <- 2 #number of spacers within each guide
guide_length <- args[3] %>% as.numeric() #depends on type of Cas, for Cas13 is 23
extension <- args[4] %>% as.numeric() #extension of the spacers, useful if there is maturation(e.g. cutting) of the spacers from the Cas protein
min_dist_type <- "p"
min_dist <- args[5] #minimum distance between different spacers within the same guide
min_dist <- as.numeric(min_dist)
min_dist_global <- args[6] #minimum distance between different spacers between guides
min_dist_global <- as.numeric(min_dist_global)
junction_flag <- TRUE #TRUE or FALSE, if yes avoid creating guides in junction gaps
input_type <- args[7] #"ENSEMBL ID" = "1", "Sequence" = "2", Genomic Coordinates" = "3")
sequence_lower_range <- args[8] %>% as.numeric()
sequence_upper_range <- args[9] %>% as.numeric()
quality_threshold <- args[10] %>% as.numeric()
restriction_sequences <- args[11]
threads <- args[12] %>% as.numeric()
output_folder <- args[13]
registerDoParallel(threads)
path_to_scripts <- args[14]
path_to_scripts <- paste0(path_to_scripts, "/")
source(paste(path_to_scripts,"/scripts/predict.R", sep = ""))
source(paste(path_to_scripts,"/scripts/trim_sequence.R", sep = ""))
source(paste(path_to_scripts,"/scripts/extend_spacer.R", sep = ""))
source(paste(path_to_scripts,"/scripts/remove_joints.R", sep = ""))
source(paste(path_to_scripts,"/scripts/remove_redundancy_array_v3.R", sep = ""))

###set default values if not in range
if (n_guides < 1) {
  paste0("You entered an out of range value for 'Number of Guides:' Default value will be used instead")
  n_guides <- as.integer(2)
}
if (n_spacers < 1) {
  paste0("You entered an out of range value for 'Number of Spacers:' Default value will be used instead")
  n_spacers <- as.integer(3)
}
if (guide_length < 20 || guide_length > 80) {
  paste0("You entered an out of range value for 'Guide Length(bp):' Default value will be used instead")
  guide_length <- as.integer(23)
}
if (extension < 0) {
  paste0("You entered an out of range value for 'Guide Extension(bp):' Default value will be used instead")
  extension <- as.integer(0)
}
if (min_dist < 0) {
  paste0("You entered an out of range value for 'Minimum Distance Between Spacers(bp):' Default value will be used instead")
  min_dist <- as.integer(1)
}

if (min_dist_global < 0) {
  paste0("You entered an out of range value for 'Global Minimum Distance Between Spacers(bp):' Default value will be used instead")
  min_dist <- as.integer(1)
}
print(paste("n_guides=", n_guides, "N_spacers=", n_spacers, "guide_length", guide_length, "extension", extension, "minimum_dinstance", min_dist, "global_minimum_distance", min_dist_global))

final_tab <- c(n_guides,n_spacers)
if (organism == "human") {
  genome_ref <- BSgenome.Hsapiens.NCBI.GRCh38
  dataset_ensembl <- "hsapiens_gene_ensembl"
} else if (organism == "mouse") {
  genome_ref <- BSgenome.Mmusculus.UCSC.mm10
  dataset_ensembl <- "mmusculus_gene_ensembl"
}
print(dataset_ensembl)
if (input_type == "1") {
  ids <- read.csv(input_file, header = F)
  names(ids) <- "IDs"
  mart <- useMart("ensembl", dataset=dataset_ensembl) #search the sequence corresponding to each the ENSEMBL ID using biomaRt
  ids$sequence <- apply(ids, 1, function(x){  #add the sequence as a column of the ids dataframe
    seq <- getSequence(id = x, 
                       type = "ensembl_transcript_id", 
                       seqType = "cdna", 
                       mart = mart)
    seq$cdna
  })
} else if (input_type == "2") {
  ids <- read.csv(input_file)
  names(ids) <- c("IDs", "sequence")
} else if (input_type == "3") {
  ids <- read.csv(input_file, header = T) 
  ids_ranges <- makeGRangesFromDataFrame(ids)
  ids$sequence <- as.vector(getSeq(genome_ref, ids_ranges))
  ids_by_exon <- as_tibble(ids)
  ids <- ids_by_exon %>% 
    group_by(IDs) %>%
    arrange(chr, start, .by_group = TRUE) %>% 
    mutate(exon_number = if_else(strand == "+", true = row_number(), false = rev(row_number()))) %>% 
    ungroup() %>%
    group_by(IDs) %>% 
    arrange(exon_number, .by_group = TRUE) %>% 
    select(IDs, sequence) %>%
    mutate(sequence = paste(sequence, collapse = "")) %>% 
    distinct() 
}
print("got sequence")
n_5 <- sequence_lower_range
n_3 <- sequence_upper_range
ids <- trim_sides(ids, n_5, n_3)
print("trimmed sequence according to parameters")

if (min_dist_type=="p") {
  print("is percentage")
  min_dist_array <- array_of_min_dist(ids, min_dist, guide_length)
  min_dist_array_global <- array_of_min_dist_global(ids, guide_length, min_dist_global)
} else if (min_dist_type=="n") {
  array <- ids %>% 
    mutate(dist_bp=min_dist) %>% 
    select(dist_bp)
  array_global <- ids %>% 
    mutate(dist_bp=min_dist_global) %>% 
    select(dist_bp)
  min_dist_array <- array$dist_bp
  min_dist_array_global <- array_global$dist_bp
} else {
  print("ERROR: wrong type of minimum distance, will be 1 bp by default")
  min_dist_type <- "n"
  min_dist <- 1
  min_dist_global <- 1
  array <- ids %>% 
    mutate(dist_bp=min_dist) %>% 
    select(dist_bp)
  min_dist_array <- array$dist_bp
  array_global <- ids %>% 
    mutate(dist_bp=min_dist_global) %>% 
    select(dist_bp)
  min_dist_array_global <- array_global$dist_bp
}
print("generated array of minimum distances")

#n <- nrow(ids)

##extract from each csv result the number of best guides wanted by the user
number_of_fastas <- as.numeric(final_tab[1])*as.numeric(final_tab[2])
clust <- makeCluster(40, type = "FORK")
#clusterExport is needed to export all the variables needed inside the parallel clusters(http://gradientdescending.com/simple-parallel-processing-in-r/)
#clusterExport(clust, varlist = c("ids", "organism", "number_of_fastas","final_tab", "n_guides", "n_spacers", "guide_length", "extension", "min_dist_type", "min_dist", "min_dist_global", "junction_flag", "input_type", "ids_by_exon", "remove_joints", "extend_spacers", "remove_redundancy_array"), envir = environment())
#clusterEvalQ is needed to use functions inside the parallel clusters, taken from: clusterEvalQ(cl, library(...))
#clusterEvalQ(clust, { 
#  library(tidyverse)
#  library(seqRFLP)
#  library(rlist)
#  })
predictions <- parApply(clust, ids, 1, function(x){
  file <- read.csv(paste("/media/rad/HDD3/Riccardo/projects/ALBAROSSA_v1_paper/code_availability/gRNA_design/test_data/prediction_filtered/", x["IDs"], "_CasRxguides_filtered.csv", sep = ""), header = T)#read the result table, where the guides are
  file <- file %>% 
    filter(quartiles >= quality_threshold) %>% #filter everything lower than the quality threshold
    arrange(desc(standardizedGuideScores), GuideName)
  
  number_of_guides <- length(x["sequence"])
  
  if (input_type == "3" && junction_flag == TRUE) { #go in only if input type is ranges
    exons <- ids_by_exon %>% filter(str_detect(IDs, paste0(as.character(x["IDs"]), "\\b")))
    spacers <- remove_joints(x["sequence"], file, exons) #exons is instead of ids_by_exons cause it contains only the exons from the current transcript!
  } else {
    spacers <- as.vector(file$GuideSeq)
  }
  n_guides <- as.numeric(n_guides)
  print(n_guides)
  restriction_sequences <- c("CGTCTC")
  spacers <- remove_redundancy_array_v3(x["sequence"], spacers, min_dist_array, guide_length, number_of_fastas, n_guides, x["IDs"], min_dist_array_global, extension, restriction_sequences)
  
  spacers <- extend_spacers(x["sequence"], spacers, extension)
  
  spacers <- unlist(spacers)

  if (length(spacers) < number_of_fastas) {
    matr <- paste(x["IDs"], "has failed")
  } else {
    final_tab_1 <- as.numeric(final_tab[1])
    final_tab_2 <- as.numeric(final_tab[2])
    matr <- matrix(spacers, ncol = final_tab_2, nrow = final_tab_1)
    base::colnames(matr) <- base::paste("spacer", 1:final_tab_2, sep = "")
    matr <- as.data.frame(matr) %>% 
      mutate(IDs = x["IDs"]) %>% 
      mutate(guide = row_number()) %>% 
      unite("spacers", spacer1:paste0("spacer", final_tab_2)) %>% 
      select(IDs, guide, spacers) %>% 
      separate_rows(spacers, sep = "_") %>% 
      mutate(GuideSeq = str_sub(spacers, end=-(extension+1))) %>% 
      left_join(file, by = "GuideSeq") %>% 
      select(IDs, guide, spacers, standardizedGuideScores,quartiles) %>%
      group_by(guide) %>% 
      mutate(spacers_name = paste0("spacer", row_number())) %>% 
      mutate(standardizedGuideScores = paste0(standardizedGuideScores, collapse = ",")) %>% 
      mutate(quartiles = paste0(quartiles, collapse = ",")) %>% 
      pivot_wider(names_from = spacers_name, values_from = spacers) 
    matr
  }
})
stopCluster(clust)

dir.create("./score_tables_parallel")
dir.create("./score_images_parallel")
output_file <- paste0(output_folder, "/prediction.csv")
file.create(output_file)
names_of_table <- c("target", "guide", "score", "quartile", paste("spacer", 1:as.numeric(final_tab[2]), sep = "")) #col.names = names_of_table
matr_name <- matrix(names_of_table, ncol = length(names_of_table), nrow = 1)
write.table(matr_name, output_file, sep = "\t", append = TRUE, quote = FALSE, col.names = FALSE,row.names = FALSE)

for (df in predictions) {
  print("df")
  print(df)
  write.table(df, output_file, sep = "\t", append = TRUE, quote = FALSE, col.names = FALSE,row.names = FALSE)
}


