#================================================================
# run_ldsc.sh
#================================================================
conda activate ldsc
total_start=2
total_end=922
chunk_size=40

for (( i=$total_start; i<=$total_end; i+=$chunk_size )); do
    end=$((i + chunk_size - 1))
    if (( end > total_end )); then
        end=$total_end 
    fi
    echo $i $end
    nohup bash ldsc_03.1_ldsc.sh $i $end &
    echo "Started from range $i to $end"
done

# ls ~/easymr/data/ldsc/ | wc -l
# pkill -f ldsc