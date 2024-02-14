
# Load R packages
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(dplyr)
library(tidyr)
library(stringr)
library(randomForest)
library(parallel)
library(doParallel)

args = commandArgs(trailingOnly=T)
input_file <- args[1] #input file in csv format, according to type of input
organism <- args[2] #either human or mouse, useful if you give ENSEMBL IDs or Ranges to retrieve correct sequence
guide_length <- args[3] %>% as.numeric() #depends on type of Cas, for Cas13 is 23
input_type <- args[4] #"ENSEMBL ID" = "1", "Sequence" = "2", Genomic Coordinates" = "3")
threads <- args[5] %>% as.numeric() #number of threads for parallel prediction
sequence_lower_range <- args[6] %>% as.numeric()
sequence_upper_range <- args[7] %>% as.numeric()
prediction_dir <- args[8]
dir.create(prediction_dir)
path_to_scripts <- args[9]

path_to_scripts <- paste0(path_to_scripts, "/")

source(paste(path_to_scripts,"/scripts/predict.R", sep = ""))
source(paste(path_to_scripts,"/scripts/trim_sequence.R", sep = ""))

predictor_input <- paste0(path_to_scripts, "/data/Cas13designGuidePredictorInput.csv")

if (guide_length < 20 || guide_length > 80) {
  paste0("You entered an out of range value for 'Guide Length(bp):' Default value will be used instead")
  guide_length <- as.integer(23)
}

if (sequence_lower_range > 1 & sequence_lower_range < 0) {
  paste0("You entered an out of range value for 'Sequence Lower Range:' Default value will be used instead")
  sequence_lower_range <- as.integer(0)
}

if (sequence_upper_range > 1 & sequence_upper_range < 0) {
  paste0("You entered an out of range value for 'Sequence Upper Range:' Default value will be used instead")
  sequence_upper_range <- as.integer(1)
}

registerDoParallel(threads)

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
##trimming of sequence sides
n_5 <- sequence_lower_range
n_3 <- sequence_upper_range
ids <- trim_sides(ids, n_5, n_3)
print("trimmed sequence according to parameters")

##run the script for each of the geneIDs
n <- nrow(ids)
foreach (i=1:n) %dopar% {
 id <- as.character(ids[[i,1]])
 sequence <- as.character(ids[[i,2]])
 print(id)
 print(sequence)
 prediction(sequence, predictor_input, 'true', id, guide_length, path_to_scripts, prediction_dir)
}


