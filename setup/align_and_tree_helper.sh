#!/bin/bash

check_and_install() {
    local prog=$1
    local install_cmd=$2
    if ! command -v "$prog" &> /dev/null; then
        echo "$prog not found, attempting to install..."
        eval "$install_cmd"
        if ! command -v "$prog" &> /dev/null; then
            echo "Was not able to install $prog automatically. Please install manually."
        else
            echo "$prog installed successfully."
        fi
    else
        echo "$prog already installed: $(command -v $prog)"
    fi
}

check_and_install "mmseqs" "conda install -c conda-forge -c bioconda mmseqs2"
check_and_install "mafft" "if [[ $(uname) == 'Darwin' ]]; then brew install mafft; else echo 'Please install mafft manually for Linux'; fi"
check_and_install "FastTree" "if [[ $(uname) == 'Darwin' ]]; then brew install fasttree; else sudo apt install -y fasttree; fi"

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

    python "$mad_script" "$working_dir/MSA/$working_name.fasta"

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
