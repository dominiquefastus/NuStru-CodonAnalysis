#!/bin/bash

# Initialize variables
input_file=""
output_file=""
reduce_redundant=false

# Function to display usage
usage() {
    echo "Usage: $0 -i <input_file> -o <output_file> [-rr|--reduce-redundant]"
    echo "  -i, --input      Input file"
    echo "  -o, --output     Output file"
    echo "  -rr, --reduce-redundant  Remove redundant sequences (optional)"
    exit 1
}

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) input_file="$2"; shift ;;
        -o|--output) output_file="$2"; shift ;;
        -rr|--reduce-redundant) reduce_redundant=true ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if required arguments are provided
if [ -z "$input_file" ] || [ -z "$output_file" ]; then
    echo "Error: Input and output files must be provided."
    usage
fi

# Run the AWK command to process the sequence
awk '/^>/ {print ">"substr($1, 2, 6); next} 1' "$input_file" > "$output_file"

echo "Sequence processing completed. Output saved as $output_file"

# Check if reduce redundant flag is set
if [ "$reduce_redundant" = true ]; then
    # Run additional AWK commands to remove redundant sequences
    awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' "$output_file" | awk '!seen[$1]++' | awk -v OFS="\n" '{print $1,$2}' > "no_dupl_$output_file"
    
    echo "Redundant sequence removal completed. Output saved as no_dupl_$output_file"
fi
