#!/bin/bash

working_directory="./" #change this to the path your saved the pipeline in
result_dir="${working_directory}/test_data/results"
output_name="test_analysis"
input_file="${working_directory}/test_data/median_normalized_counts.tsv"
#run analysis
Rscript ${working_directory}/scripts/mageck_RRA_run.R \
    ${result_dir} \
    ${output_name} \
    ${input_file}