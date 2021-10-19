#!/usr/bin/env bash

# script to prepare a run of models
folder_input="/home/ubuntu/Models/GenePhene2/data.files"
folder_output="/home/ubuntu/Models/GenePhene2/model.results"
files="${folder_input}/GenePhene2_Input_*"
for file in ${files}
do

echo  $folder_input $(basename "${file}") $folder_output
./model_setup.R $folder_input $(basename "${file}") $folder_output

done
