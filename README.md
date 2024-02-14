This is the Github repository for the publication "Genome-scale pan-cancer interrogation of lncRNA dependencies using CasRx" from Montero, Trozzo, et al. Nature Methods, 2024

Here you can find pipelines and scripts for CasRx gRNA arrays design, gRNA arrays quantification from amplicon-seq data and the computational analysis of the CasRx lncRNAs screens

#### Abstract
Although lncRNAs dominate the transcriptome, their functions are largely unexplored. lncRNA
characteristics, such as extensive overlap with coding and regulatory sequence restrict their systematic
interrogation by DNA-directed perturbation. Here, we developed genome-scale lncRNA-transcriptome
screening using Cas13d/CasRx. We show that RNA-targeting overcomes limitations inherent to other
screening methods, thereby considerably expanding the explorable space of the lncRNAome. By
evolving the screening system towards pan-cancer applicability, it supports molecular and phenotypic
data integration to contextualize screening hits or infer lncRNA function. We thereby addressed
challenges posed by the enormous transcriptome size and tissue-specificity through a size-reduced
multiplexed gRNA-library targeting 24,171 lncRNA-families. Its rational design incorporates target
prioritization based on expression, evolutionary conservation, and tissue-specificity, thereby reconciling
high discovery-power and pan-cancer representation with scalable experimental throughput. Applied
across entities, the screening platform identified numerous context-specific and common-essential
lncRNAs. Our work sets the stage for systematic exploration of lncRNA biology in health and disease.


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


# CRISPR/CasRx screen analysis

## Description
This script takes as input a median normalized count table and runs MAGeCK RRA algorithm (_test_ command) for each of the samples in the table against the library pool with default parameters. It then collects the output generated for each cell line and generates two common result files, one at guide level and one at gene level, where all the results are stored.

#### Input file
The input file must be a tsv file with the following format:
| sgRNA | Gene | Library_rep1 | CellLine1_rep1 | CellLine1_rep2 | CellLine2_rep1 | ... |
| --- | --- | --- | --- | --- | --- | --- |
| ENST00000005284_1 | ENST00000005284 | 102.77 | 180.59 | 174.31 | 88.24 | ... |
| ENST00000005284_2 | ENST00000005284 | 225.20 | 114.31 | 334.77 | 356.60 | ... |
| ... | ... | ... | ... | ... | ... | ... | 

Where ***sgRNA*** is the column containing a unique sgRNA identifier; ***Gene*** is the column containing the Gene name corresponding to the sgRNA; ***Library_rep1*** is the column where the sequenced library normalized counts are represented (There can be more than one replicate for the Library, e.g. Library_rep2) and the rest of the columns contain the normalized counts for each replicate of each cell line in the form _CellLineName_rep[n]_.

## Dependencies
In order to run the pipeline the following dependencies have to be downloaded:
- MAGeCK (v0.5.9.4)
- R v3.6.1 with the following dependencies
- - data.table (v1.13.6)
- - tidyverse (v1.3.0)

## Running the Pipeline
The pipeline can be run by launching the script _run_prediction.sh_ on the command line. All the parameters can be tweaked inside the script in order to design arrays for the desired targets, following the specifications given above.
In order for the off-target filtering step to work, a bowtie2 index of a reference containing the targeted transcripts and other transcripts the users would consider as off-targets (e.g. protein coding genes) should be provided. All the gRNAs mapping to more than a transcript in the reference will be removed.

##### _test code_
In the _test data_ we provide an example input file, which contains the median normalized counts for three different cell lines (KP4, MIAPACA2, MIAPACA2 without CasRx) and the Library pool counts, and a _test_run_analysis.sh_ file which can be launched from the command line in order to run the test analysis.
To run the test analysis, please specify the location where you saved the _scripts_ foler containing the _mageck_RRA_run.R_ script, paste the following command on the shell and press enter:

```sh
bash test_run_analysis.sh
```
estimated run time: ~2 min
 An _expected_output_ folder with the expected output is provided in the _test_data_ folder.
 
 
 # gRNA array design pipeline

## Description
This pipeline is is able to design CasRx gRNA arrays based on a number of user inputs that can be chosen to improve considerably the targeting efficacy. In particular, using an adaptation of the Wessel et al. (Nature Biotechnology, 2020) Cas13 guide prediction algorithm, all possible gRNA sequences (based on sequence length defined by the user) are generated and their efficiency is predicted based on several parameters defined by Wessel et al. (Nature Biotechnology, 2020). These sequences undergo a step of off-target filtering in order to remove sequences that map to multiple transcripts or to coding transcripts. Once the sequences are filtered, they are used to design the complex CRISPR/CasRx array based on the following parameters that the user can choose:

-  input type
-  organism
-  number of arrays for each transcript
-  number of gRNAs for each array
-  gRNA length
-  gRNA extension
-  minimum distance between gRNAs within each array (local minimum distance)
-  minimum distance between all gRNAs targeting a transcript (global minimum distance)
-  distance by percentage or by base pairs
-  skip exon junctions
-  avoid specific restriction enzyme sequences
-  restrict design to a portion of the transcript

##### input_type
The user can chose between three types of input: Ensembl Transcript ID (1), sequence (2) or Genomic coordinates (3). 
For Ensembl Transcript ID (1), a txt file with a list of Ensembl Transcript IDs and no header has to be provided (see example below).

|   |
| --- |
| ENST00000380152 |
| ENST00000614259 |
| ... |

For sequence (2), a csv file with two columns (IDs and sequence) has to be provided (see example below).

| IDs | sequence |
| --- | --- |
| target1 | ATCGACGATCGACTGGGGCTATCAGTGGCCCC |
| target2 | GCATGCGGCCTTTAATGCACGGATTACTGAGCAGCGTTAA |
| ... | ... |

For Genomic Coordinates (3), a csv file with the transcript ID, chr, start, end and strand for each exon forming the transcript (common id) has to be provided (see example below).

| IDs | chr | start | end | strand |
| --- | --- | --- | --- | --- |
| transcript1 | 1 | 1732257 | 1732523 | - |
| transcript1 | 1 | 1734689 | 1734835 | - |
| transcript2 | 2 | 2834342 | 2834598 | - |
| ... | ... | ... | ... | ... |

##### organism
It is possible to design gRNA arrays for human (hg38) and mouse (mm10).
##### number of arrays for each transcript
It's the number of gRNA arrays that will be designed
##### number of gRNAs for each array
It's the number of gRNAs that will form each array
##### gRNA length
It's the length of the gRNAs. The optimal gRNA length for CasRx is 23bp. Should be set to 23.
##### gRNA extension
It's the base pairs length of extension of the gRNA along the targeted transcript. For CasRx in array configuration (more than 1 gRNA) there is a maturation step in which a 30bp gRNA is maturated into a 23bp gRNA by cutting 7bp at the 3' end of the gRNA. Should be set to 7.
##### minimum local distance
It's the minum distance between each of the gRNAs that are part of the same array. It is used maximise the span of the transcript covered by gRNAs.
##### minimum global distance
It's the minum distance between all gRNAs that are targeting the same trasncript. It is used to avoid overlap between different arrays targeting the same transcript.
##### distance by percentage of by base pairs
The two distance parameters explained above (minimum local distance and minimum global distance) can be expressed as plain base pairs (n) or as percentage of the targeted transcript length (p). Using "p" is recommended and it's default.
##### skip exon junction
Logical. If the input type is Genomic Coordinates (3), the user can decide, by setting this parameter to T/TRUE, to avoid targeting exon-exon junctions. Default is T/TRUE.
##### avoid specific restriction enzyme sequences
Since the library has to undergo a process of cloning involving the use of restriction enzymes, it's possible to specify to the pipeline a restriction enzyme recognition sequence that should be avoided in the CRISPR/CasRX array design. The default is BsmBI ("CGTCTC")
##### restrict design to a portion of the transcript
The user can decide to restrict the portion of the transcript that will be targeted by gRNAs by selecting to parameter:
- sequence lower range
- sequence upper range

The default is 0 (beginning of transcript) for the lower range and 1 (end of transcript) for the upper range. Any range from 0  to 1 where the lower range is < the upper range will be accepted. 

#### How the Pipeline Works
The pipeline reads the input and retrieves the corresponding sequence (unless input is already a sequence). It then runs the prediction algorithm (Wessel et al. 2020) to predict targeting efficiency for each possible sequence of length defined by the user (guide_length parameter). The off-target filtering step gets in input the output of the prediction algorithm, feeds the sequences to bowtie to map each of these sequences to the off-target reference. It then return the prediction file, but without the entries that map to more than one location (potential off-target effect), retaining only unique mappers. This output is taken by the array design module. This input file is then taken by the array design module. There are then 4 main design steps:
1) quality filtering:
only the gRNA sequences that are below a user defined quality threshold (4th quartile for our design) are removed.
2) removal of exon-exon junctions:
Only gRNA sequences that are falling completely within an exon of the original transcript are retained (unless differently specified by the user), filtering out sequences that are overlapping exon-exon junctions
3) array design:
To design the array according to the distance parameters and the restriction enzyme sequence to avoid, all the gRNA sequences ordered by predicted efficiency score are fed to a design algorithm. This algorithm takes the first gRNA that fulfills the parameters (not containing restriction enzyme sequences) with the highest quality, this becames the gRNA in position 1 of array 1. It then scans through the remaining sequences ranked by quality until it finds the first one to fulfill the design parameters (minimum global distance with all the other gRNAs and doesn't contain the restriction enzyme sequence), until it filled the position 1 for all the arrays. It then starts to fill the gRNAs in position 2 of the arrays by scanning through the remaining sequences until it finds the first one to fulfill the design parameters (minimum global distance with all the other gRNAs; minimum local distance with the gRNA in the same array and doesn't contain the restriction enzyme sequence), until it filled the position 2 for all the arrays. If it finds enough gRNAs fulfilling all the parameters, it will output the designed arrays, otherwise, it will inform the user that the parameters were to strict for that specific target.
4) extension:
all the gRNAs are extended at their 3' end based on the extension length parameter defined by the user and following the sequence of the targeted transcript.

A tsv file is then given as output containing the target name, the guide number, the gRNA scores for that array separated by a comma, the quartile scores for that array separated by a comma and the sequences (called spacers) that form the array:
| target | guide | score | quartile | spacer1 | spacer2 | ... |
| --- | --- | --- | --- | --- | --- | --- |
| human_lncrna_fused_174 | 1 | 1,0.80 | 4,4 | TATCGATAAGAAGCTCACAGCCAAGGCTGT | TGAGCCAATACTGAGCCACTGAACTCCAGC | 
| human_lncrna_fused_174 | 2 | 0.96,0.77 | 4,4 | GAACTGATTATTACAGCAGCGAGGGAAACT | AATCACCTAAACGTGTGTGCGGGTCTCAAG | 
| ... | ... | ... | ... | ... | ... |


## Dependencies
In order to run the pipeline the following dependencies have to be downloaded:
- bowtie2 (v2.3.5.1)
- RNAfold (ViennaRNA Package 2.0) in /usr/local/bin
- RNAplfold (ViennaRNA Package 2.0) in /usr/local/bin
- R v3.6.1 with the following dependencies
- - GenomicRanges(v1.38.0)
- - BSgenome.Hsapiens.NCBI.GRCh38 (v1.3.1000)
- - BSgenome.Mmusculus.UCSC.mm10 (v1.4.0)
- - data.table (v1.13.6)
- - Biostrings (v2.54.0)
- - biomaRt (v2.42.1)
- - dplyr (v1.0.2)
- - tidyr (v1.1.2)
- - BiocManager (v1.30.19)
- - stringr (v1.4.0)
- - randomForest (v4.6-10)
- - R.utils (v2.10.1)
- - zip (v2.1.1)
- - seqRFLP (v1.0.1)
- - rlist (v0.4.6.2)
- - foreach (v1.5.2)
- - doParallel (v1.0.17)
- - farver (v2.0.3)

## Running the Pipeline
The pipeline can be run by launching the script _run_prediction.sh_ on the command line. All the parameters can be tweaked inside the script in order to design arrays for the desired targets, following the specifications given above.
In order for the off-target filtering step to work, a bowtie2 index of a reference containing the targeted transcripts and other transcripts the users would consider as off-targets (e.g. protein coding genes) should be provided. All the gRNAs mapping to more than a transcript in the reference will be removed.

```sh
cd ${pipeline_directory}
bash run_prediction.sh
```
