#!/bin/bash
total_start=4500
total_end=6640
chunk_size=100  

for (( i=$total_start; i<=$total_end; i+=$chunk_size )); do
    end=$((i + chunk_size - 1))
    if (( end > total_end )); then
        end=$total_end 
    fi
    nohup Rscript code/mr_03_mr.r $i $end &
    echo "Started from range $i to $end"
done