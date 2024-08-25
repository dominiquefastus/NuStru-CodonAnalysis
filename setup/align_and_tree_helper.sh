#!/bin/bash

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fasta) fasta_file="$2"; shift ;;
        --mad) mad_script="$2"; shift ;;
        --dir) parent_dir="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$mad_script" || (-z "$fasta_file" && -z "$parent_dir") ]]; then
    echo "Usage: $0 --mad <path_to_mad_script> [--fasta <path_to_fasta_file> | --dir <parent_directory>]"
    exit 1
fi

process_fasta() {
    local fasta_file="$1"
    local job_id=$(uuidgen)
    local cwd=$(pwd)
    local job_path="$cwd/nustruTREE/$job_id"
    local msa_path="$job_path/MSA"
    local tree_path="$job_path/TREE"

    mkdir -p "$msa_path" "$tree_path"

    if cp "$fasta_file" "$msa_path"; then
        echo "File copied successfully!"
    else
        echo "Error: Failed to copy the fasta file. Please check the path and permissions."
        return 1
    fi

    local working_dir="$job_path"
    local working_name=$(basename "$fasta_file" .fasta)

    mafft --auto "$working_dir/MSA/$working_name.fasta" > "$working_dir/MSA/${working_name}_aligned.fasta"

    FastTree "$working_dir/MSA/${working_name}_aligned.fasta" > "$working_dir/TREE/${working_name}.tree"

    python "$mad_script" "$working_dir/TREE/${working_name}.tree""

    echo "Pipeline completed successfully for $fasta_file."
}

if [[ -n "$fasta_file" ]]; then
    process_fasta "$fasta_file"
elif [[ -n "$parent_dir" ]]; then
    for dir in "$parent_dir"/*/; do
        for fasta in "$dir"/*.fasta; do
            if [[ -f "$fasta" ]]; then
                process_fasta "$fasta"
            fi
        done
    done
fi
