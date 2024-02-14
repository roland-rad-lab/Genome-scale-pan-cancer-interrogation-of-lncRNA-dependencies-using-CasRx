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

