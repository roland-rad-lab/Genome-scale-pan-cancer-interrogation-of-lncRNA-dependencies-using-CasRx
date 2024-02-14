#!/bin/bash

#this pipeline is made for counting sgRNAs with two spacers, one on the forward read and one on the reverse read
#it works by mapping the reads to a custom index of the library in order to count the CRISPR arrays

#let's start by creating the index with spacer1 and spacer2
#we have to make two indexes, one for position 1 and the second for position 2, we cannot put them together cause it could create mapping artifacts overlapping the two of them
#to create the two indexes, the following command lines are used
#bowtie-build -f tables/ALBAROSA_lib_pos1_first-27.fa indexes/lib_spacer_1/lib_pos1 for position 1
#bowtie-build -f tables/ALBAROSA_lib_pos2_last-27.fa indexes/lib_spacer_2/lib_pos2 for position 2

######################################################################
###################GET ARGUMENTS FROM COMMAND LINE####################
######################################################################
usage()
{
echo "Usage: $0 "
echo '--sample_table                sample list'
echo '--combinations_table          combinations table'
echo '--lib                         library table'
echo '--out_dir                     output directory'
echo '--name                        name for the project'
echo '--trm_jar                     trimmomatic jar file'
echo '--crop_len                    crop length; forward.left:forward.right,reverse.left:reverse:right; as "right" the length from start to end is to be given'
echo '                              if no cropping is to be performed, use dot (.) for the repective position'
echo '--reverse_complement_R1       T if the forward read has to be reverse complemented, F if not'
echo '--reverse_complement_R2       T if the reverse read has to be reverse complemented, F if not'
echo '--read_1_pos                  first if read one is the sequencing of the first spacer, second if it is the sequencing of the second spacer'
echo '--read_2_pos                  first if read two is the sequencing of the first spacer, second if it is the sequencing of the second spacer'
echo '--mismatches                  number of mismatches allowed'

exit 1
}

if [ $# == 0 ]; then
     usage
fi

echo 'cas13 count pipeline was called with the following parameters:'
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        echo $1 $2 #// Optional to see the parameter:value result
   fi

  shift
done

skriptdir=$(dirname "$0")

# if exampleval is not declared, it is assigned the default value provided
out_dir=${out_dir:-'.'}
mkdir -p ${out_dir}
mkdir -p ${out_dir}/trimmed
mkdir -p ${out_dir}/mapped

#####################################################################################
#######################################TRIMMING######################################
#####################################################################################
#once we have the indexes we need to trim the guides in order to obtain only the portion of the read representing the guides
#to trim the reads we use trimmomatic
if [ -v sample_table ] && [ -v trim_dir ] && [ -v crop_length ] && [ -v out_dir ] && [ -v path_to_pipeline ]; then 
    echo "Trimming"
    Rscript ${path_to_pipeline}/modules/trim.R ${sample_table} ${trim_dir} ${crop_length} ${out_dir}
else
    echo "cannot run trimming with missing arguments"
fi
#####################################################################################
#######################################MAPPING#######################################
#####################################################################################
#once the files are trimmed we can use them to map the reads on the corresponding positions
if [ -v sample_table ] && [ -v reverse_complement_R1 ] && [ -v reverse_complement_R2 ] && [ -v read_1_pos ] && [ -v read_2_pos ] && [ -v path_to_index_1 ] && [ -v path_to_index_2 ]  && [ -v out_dir ]  && [ -v mismatches ] && [ -v path_to_pipeline ]; then 
     echo "Mapping"
     Rscript ${path_to_pipeline}/modules/map.R ${sample_table} ${reverse_complement_R1} ${reverse_complement_R2} ${read_1_pos} ${read_2_pos} ${path_to_index_1} ${path_to_index_2} ${out_dir}/trimmed ${out_dir}/mapped ${mismatches}
else
     echo "cannot run mapping with missing arguments"
     exit 0
fi
#####################################################################################
#######################################COUNTING######################################
#####################################################################################
#now that the mapping is done we have to put together the counting based on the sgRNA combination file
if [ -v sample_table ] && [ -v combination_table ] && [ -v out_dir ]&& [ -v path_to_pipeline ]; then 
     echo "counting arrays"
     Rscript ${path_to_pipeline}/modules/combine_v2_parallel.R ${sample_table} ${combination_table} ${out_dir}
else
     echo "cannot run combine script with missing arguments"
     exit 0
fi
#now we need to put the counts all together in a table and give the number of mapped and unmapped
if [ -v sample_table ] && [ -v combination_table ] && [ -v out_dir ]&& [ -v path_to_pipeline ]; then 
     echo "creting output files"
     Rscript ${path_to_pipeline}/modules/organize.R ${sample_table} ${combination_table} ${out_dir}
else
     echo "cannot run organize script with missing arguments"
     exit 0
fi