# Array CRISPR/CasRx counting pipeline

## Description

This pipeline is able to count gRNA arrays in CRISPR pooled screenings where the construct consists of an array with two gRNAs that are separately read by read1 and read2 from the sequencing machine.

It trims the reads based on user input parameters to leave only the portions that cover the gRNAs at the respective position. It first deals with the forward and reverse reads separately by using bowtie to map read1(forward) and read2(reverse) to the to the respective array position (user defined). It then collects the separate counting sam outputs and merges them together based on the Illumina read ID. It then assignes a count if the sequence from both read1 and read2 have mapped to the same array.

## Dependencies

In order to run the pipeline the following dependencies have to be downloaded:

- bowtie v1.2.2
- trimmomatic (v0.38)
- R v3.6.1 with the following dependencies
- - doParallel (v1.0.17)
- - foreach (v1.5.2)
- - tidyverse (v1.3.0)
- - data.table (v1.13.6)

## Creating the Library Indexes

First, the user has to create two separate bowtie indexes: one for the gRNA in position 1 and the other for gRNA in position 2. This can be done by running the _bowtie-build_ command (bowtie v1.2.2) on a fasta file containing the gRNA sequences for the position that is going to be indexed.

##### _test code_

In the _test data_ we provide two fasta files, one for the gRNA in position 1 and the other for gRNA in position 2. We can index this files using the following commands:

```
bowtie-build test_data/library/ALBAROSA_lib_pos1_first-27.fa test_data/library/indexes/lib_spacer_1/lib_pos1
bowtie-build test_data/library/ALBAROSA_lib_pos2_last-27.fa test_data/library/indexes/lib_spacer_2/lib_pos2 
```

## Running the Pipeline

#### Input Parameters

In order to run the pipeline, the user have to run the script _cas13_count_pipe_run_command.sh_. In the script, the following parameters need to be provided:

* path_to_pipeline
* sample_table
* combination_table
* trim_dir
* crop_length
* out_dir
* path_to_index_1
* path_to_index_2
* reverse_complement_R1
* reverse_complement_R2
* read_1_pos
* read_2_pos
* mismatches

##### path_to_pipeline

The path where the main pipeline and the modules are located

##### sample_table

The path to the _sample table_. The sample table contains: the sample information; the group to which the sample belongs to (e.g. two replicate samples belong to the same cell line group) and the locations for the forward and the reverse read

| sample | group | R1 | R2 |
| --- | --- | --- | --- |
| cell_rep1 | cell_line | path_to_R1 | path_to_R2 |
| ... | ... | ... | ... |

The fastq files specified in the _sample table_ have to come from amplicon sequencing (as described in the methods) and the first base pair has to start always in the same position of the sequenced amplicons in order for the trimming to be possible.

##### combination_table

The combination table contains the two gRNA IDs for each array separated by a "~". In our case the IDs are identical (it corresponds to the ID of the array), but they can differ in case we're targeting different genes (e.g. combinatorial CRISPR screens).

| combinations |
| --- |
| human_lncrna_fused_53558_1~human_lncrna_fused_53558_1 |
| human_lncrna_fused_53558_2~human_lncrna_fused_53558_2 |
|  ... |

##### trim_dir

The full path of the local trimmomatic jar file

##### crop_length

crop_length gives the information on how the reads should be cropped. It has the format: fi:fl,ri:rl, where fi is the base where the forward read start should be cropped; fl is the length of the fragment that has to be extracted from the forward read; ri is the base where the reverse read start should cropped; rn is the length of the fragment that has to be extracted from the reverse read.

**EXAMPLE**: the string _"1:27,26:30"_ will extract 27 bp starting from position one from the forward read and 30 bp starting from position 26 from the reverse

##### out_dir

The directory where the output files will be saved

##### path_to_index_1

The path to the index of the first position in the gRNA array

##### path_to_index_2

The path to the index of the second position in the gRNA array

##### reverse_complement_R1

T/TRUE if the read1 has to be reverse complemented before mapping, otherwise F/FALSE

##### reverse_complement_R2

T/TRUE if the read2 has to be reverse complemented before mapping, otherwise F/FALSE

##### read_1_pos

first if it's mapping to the index_1 and second if it's mapping to the index_2

##### read_2_pos

first if it's mapping to the index_1 and second if it's mapping to the index_2

##### mismatches

The number of mismatches to allow during the mapping procedure

##### _test code_

In the _test data_ folder, all the neccessary files to test the pipeline are provided. In detail, there are two small test fastq files; the fasta files representing position 1 and position 2 of the array and their respective index; the combinations table and the sample table; a folder with the expected output from the pipeline.

In order to run the test code, the user must run in the command line the script _test_cas13_count_pipe_run_command.sh_ using the following command:

```
bash CasRx_count_pipe_run_command.sh
```
expected runtime ~2 min