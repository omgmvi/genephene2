#!/usr/bin/env bash

# script to prepare a run of models
folder_input="/home/ubuntu/Models/GenePhene2/test.files"
folder_output="/home/ubuntu/Models/GenePhene2/test.results"

folder_model="/home/ubuntu/GenePhene2/Models/config.files"
file_model="model.glmnet_elasticnet.json"

files="${folder_input}/GenePhene2_*.dat"
for file in ${files}
do

    echo  $folder_input $(basename "${file}") $folder_model $file_model  $folder_output

    ./model_setup.R $folder_input $(basename "${file}") $folder_model $file_model $folder_output
#    mv "${file}" "${file}.complete"
done
