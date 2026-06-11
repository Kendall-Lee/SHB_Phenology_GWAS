#!/bin/bash

# Directory containing your *_samtools_coverage.txt files

# Loop through each *_samtools_coverage.txt file
for file in ./*_samtools_coverage.txt; do
    echo "Processing $file:"
    average=$(tail -n +2 "$file" | awk '{ total += $6; count++ } END { if (count > 0) print total/count; else print "No data" }')
    echo "Average of column 6: $average"
    echo  # Add an empty line for separation
done

