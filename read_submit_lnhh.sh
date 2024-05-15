#!/bin/bash

# Set my inputs
start_line=$1
end_line=$2
filename=$3

i=$start_line
# Read lines from the theta file
sed -n "${start_line},${end_line}p" "${filename}" | while IFS= read -r line
##while IFS= read -r line
do
    # Split the line into parameters
    IFS=' ' read -ra params <<< "$line"

    # Assign parameters to variables
    qs02=${params[0]}
    gamma=${params[1]}
    c2=${params[2]}

    # Processing
    c2_invlog=$(awk -v c2="$c2" 'BEGIN {print 10^c2}')
    Cdown=0.1
    Cup=5.0
    c2step=$(awk -v c2_invlog="$c2_invlog" -v Cdown="$Cdown" -v Cup="$Cup" 'BEGIN {print 240 * (log(c2_invlog)/log(10) - log(Cdown)/log(10)) / (log(Cup/Cdown)/log(10))}')
    qs02_1000=$(awk -v qs02="$qs02" 'BEGIN {print qs02 * 1000}')
    sigma02=$(awk -v param="${params[3]}" 'BEGIN {print param * 2.56819}')
    n=$((i-1))

    # Submit job to SLURM
    echo sbatch -J nlosigmar_hh submitresumbk_nlohh.sh "$qs02_1000" "$c2step" "$sigma02" "$gamma" "$n"
    
    ((i++))

done