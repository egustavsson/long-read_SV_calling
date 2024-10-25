#!/bin/bash

#########################################################################################################
#### This script merges pod5 files to a single file using POD5 and performs basecalling using dorado ####
#########################################################################################################

# directories to dorado and to the model to use
# Refer to https://github.com/nanoporetech/dorado and https://github.com/nanoporetech/dorado/issues/214 to choose the correct model
# these can be downloaded by running "dorado download --model"
dorado_dir="/home/MinaRyten/Emil/Tools/dorado-0.4.3-linux-x64/bin/"
model="/home/MinaRyten/Emil/Tools/dorado-0.4.3-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0"

# help function
helpFunction()
{
    echo ""
    echo "Usage: $0 -i input_directory1,input_directory2 -o output_directory -s sample_name"
    echo -e "\t-i input directories (comma-separated)"
    echo -e "\t-o output directory"
    echo -e "\t-s sample name"
    exit 1
}

while getopts "i:o:s:" opt
do
    case "$opt" in
        i ) input_directories="$OPTARG" ;;
        o ) output_directory="$OPTARG" ;;
        s ) sample_name="$OPTARG" ;;
        ? ) helpFunction ;;
    esac
done

IFS=',' read -ra directories <<< "$input_directories"

if [ -z "$input_directories" ] || [ -z "$output_directory" ] || [ -z "$sample_name" ]
then
    echo "Some or all of the parameters are empty";
    helpFunction
fi

# Temporary pod5 output file
temp_output_file="$(mktemp)"

# Merge pod5 files from multiple directories to a temporary pod5 file
for dir in "${directories[@]}"; do
    pod5 merge "$dir"/*.pod5 >> "$temp_output_file"
done

# Run basecalling using dorado
# Additional notes:
# - Make sure to replace "your_model" with the actual model you want to use for basecalling.
# - If dorado requires any specific options or parameters, include them in the dorado command.
${dorado_dir}dorado basecaller -x cpu ${model} "$temp_output_file" --skip-model-compatibility-check --modified-bases 5mCG_5hmCG > "${output_directory}/${sample_name}_pod5.5mCG_5hmCG.bam"

# Remove the temporary output file
rm "$temp_output_file"
