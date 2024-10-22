#!/bin/bash

path_to_gRNA_design="./" #change this to the path your saved the pipeline in
scripts_dir="${path_to_gRNA_design}/modules"
input_file="${path_to_gRNA_design}/test_run/data/exons_for_prediction.csv"
organism="human"
guide_length="23"
extension="7"
n_guides="3"
input_type="3"
min_dist="15"
min_dist_global="10"
quality_threshold="4"
threads="40"
sequence_lower_range="0"
sequence_upper_range="1"
restriction_sequences="CGTCTC"
prediction_dir="${path_to_gRNA_design}/test_run/results/prediction"
filtered_prediction_dir="${path_to_gRNA_design}/test_run/results/prediction_filtered"
mkdir -p ${filtered_prediction_dir}
offtarget_index="${path_to_gRNA_design}/test_run/data/offtarget_index/bowtie2_index_22_09/off_targets_lnc_and_coding_22_09"
output_folder="${path_to_gRNA_design}/test_run/results"

#run predictions
Rscript ${scripts_dir}/gRNA_prediction.R \
    ${input_file} \
    ${organism} \
    ${guide_length} \
    ${input_type} \
    ${threads} \
    ${sequence_lower_range} \
    ${sequence_upper_range} \
    ${prediction_dir} \
    ${scripts_dir}
#run off-target filtering
bash ${scripts_dir}/off_targets_filtering.sh \
    ${prediction_dir} \
    ${filtered_prediction_dir} \
    ${offtarget_index}
#run gRNA design
Rscript ${scripts_dir}/Cas13_array_design.R \
    ${input_file} \
    ${n_guides} \
    ${guide_length} \
    ${extension} \
    ${min_dist} \
    ${min_dist_global} \
    ${input_type} \
    ${sequence_lower_range} \
    ${sequence_upper_range} \
    ${quality_threshold} \
    ${restriction_sequences} \
    ${threads} \
    ${output_folder} \
    ${scripts_dir} \
    ${filtered_prediction_dir}
