#!/bin/bash

path_to_pipeline="./" #change this to the path your saving the pipeline in
path_to_trimmomatic="" # Insert here path to trimmomatic-0.38.jar

bash ${path_to_pipeline}/CasRx_count_pipe.sh \
    --path_to_pipeline ${path_to_pipeline} \
    --sample_table ${path_to_pipeline}/test_data/tables/sample_table.tsv \
    --trim_dir  ${path_to_trimmomatic}\
    --crop_length .:27,23:27 \
    --out_dir ${path_to_pipeline}/test_data/test_map \
    --reverse_complement_R1 T \
    --reverse_complement_R2 F \
    --read_1_pos second \
    --read_2_pos first \
    --path_to_index_1 "${path_to_pipeline}/test_data/library/indexes/lib_spacer_1/lib_pos_first-27_1" \
    --path_to_index_2 "${path_to_pipeline}/test_data/library/indexes/lib_spacer_2/lib_pos_last-27_2" \
    --mismatches 2 \
    --combination_table ${path_to_pipeline}/test_data/tables/combinations_table.tsv
