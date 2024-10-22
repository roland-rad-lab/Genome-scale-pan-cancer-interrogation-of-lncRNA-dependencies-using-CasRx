#!/bin/bash

#get paths
prediction_dir=$1
filtered_prediction_dir=$2
off_targets_index=$3

mkdir -p ${filtered_prediction_dir}

#run off-target filtering for every target
for prediction_file in ${prediction_dir}/*.csv; do
old_file=$(basename ${prediction_file})
new_file="${filtered_prediction_dir}/${old_file%.csv}_filtered.csv"
echo "filtering ${old_file} to ${new_file}"
head -n 1 ${prediction_file} > ${new_file}
tail +2 ${prediction_file} | awk -F "\"*,\"*" '{print">"$1,$2}' | tr " " "\n" | bowtie2 --end-to-end -D 5 -R 1 -N 0 -L 23 -i S,0,2.50 --nofw -f -p 75 -x ${off_targets_index} -U - | samtools view -q 42 | awk '{print$1}' | grep -f - ${prediction_file} >> ${new_file}
# awk -F "\"*,\"*" '{print">"$1,$2}' ${prediction_file} | tr " " "\n" | bowtie2 --end-to-end -D 5 -R 1 -N 0 -L 23 -i S,0,2.50 --nofw -f -p 75 -x ${off_targets_index} -U - | samtools view -q 42 | awk '{print$1}' | grep -f - ${prediction_file} >> ${new_file}
done