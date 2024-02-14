library(tidyverse)
library(foreach)
library(doParallel)

######PARALLELIZATION ARGUMENTS########
#get number of cores
n_cores <- floor(parallel::detectCores()/2)
#create the cluster
my_cluster <- parallel::makeCluster(
  n_cores, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my_cluster)
#args section
args = commandArgs(trailingOnly=T)
#sample table needs to have the following format:
#SAMPLE GROUP FORWARD_READ_PATH REVERSE_READ_PATH
#where SAMPLE is the name of the sample (e.g. LN-229_rep1); GROUP is the group to which many replicate belong (e.g. LN-229);
#FORWARD_READ_PATH is the full path to the fastq file of the forward read for a specific sample;
#REVERSE_READ_PATH is the full path to the fastq file of the reverse read for a specific sample. 
sample_table_file <- args[1]
#trim_jar is the full path to the jar file for trimmomatic (e.g. /home/user/packages/trimmomatic-0.38/trimmomatic-0.38.jar)
trim_jar <- args[2]
#crop_length gives the information on how the reads should be cropped
#it comes in the format fi:fn,ri:rn where fi is the base where the forward read start should be cropped; fn is the base where the
#forward read end should cropped; ri is the base where the reverse read start should cropped; rn is the base where the reverse read end should cropped;
#EXAMPLE: the string 1:27,26:56 will crop the forward read at position 1 and 27 and the reverse read at position 26 and 56
crop_length <- args[3]
#out_dir is the path to the output directory of the project
out_dir <- args[4]
#trim dir is the trimming dir inside the output directory of the project
trim_dir <- paste0(out_dir, "/trimmed")

##########FUNCTIONS##########
#this functions provide the HEADCROP and CROP values for trimmomatic based on the crop_length string
build_crop <- function(cs){
  cl = strsplit(cs,':')[[1]]
  outl = c()
  if(cl[1] != '.'){
    outl[1] = paste0(' HEADCROP:',cl[1])
  } else {outl[1]=''}
  if(cl[2] != '.'){
    outl[2] = paste0(' CROP:',cl[2])
  } else {outl[2]=''}
  return(paste0(outl, collapse= ' '))
}
################################
sample_table <- read_tsv(sample_table_file)

#cropping and launching trimmomatic taken from Olga Baranov's in-house Mageck pipeline and slightly adapted
cropin = strsplit(crop_length,',')[[1]]

cropfwd = build_crop(cropin[1])
if(("R2" %in% colnames(sample_table)) & (length(cropin) == 2)){
  if(cropin[2] == '.:.'){
    croprev = ''
  } else { croprev = build_crop(cropin[2]) }
  
} else {croprev = ''}


rows <- rownames(sample_table)
foreach(i=1:length(rows)) %dopar% {
  #row <- as.integer(row)
  row <- rows[i]
  fwd_name = basename(as.character(sample_table[row, 3]))
  trimmed_fwd_name <- paste0(sub('\\.fastq.gz$', '', fwd_name), ".trimmed.fastq.gz")
  fwd_cmd = paste0('java -jar ', trim_jar ,' SE  -phred33 ', sample_table[row, 3],' ', trim_dir, '/', trimmed_fwd_name, cropfwd)
  #print(fwd_cmd)
  system(fwd_cmd)
  if("R2" %in% colnames(sample_table)){
    rev_name = basename(as.character(sample_table[row, 4]))
    trimmed_rev_name <- paste0(sub('\\.fastq.gz$', '', rev_name), ".trimmed.fastq.gz")
    if (croprev != ''){
      rev_cmd = paste0('java -jar ', trim_jar ,' SE -phred33 ', sample_table[row, 4],' ', trim_dir, '/', trimmed_rev_name, croprev)
    } else { rev_cmd = paste0('scp ', sample_table[row, 3],' ', trim_dir, '/', trimmed_rev_name) }
    #print(rev_cmd)
    system(rev_cmd)
  }
}

