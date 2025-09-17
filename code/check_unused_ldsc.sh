#!/bin/bash

# Set paths
CSV_FILE="/home/zhang/easymr/info/ldsc_pair.csv"
LDSC_DIR="/home/zhang/easymr/data/ldsc"
UNUSED_FILE="./tmp_unused_ldsc_files.txt"
TEMP_ALL_FILES="./tmp_all_ldsc_files.txt"
TEMP_VALID_KEYS="./tmp_valid_keys.txt"

# Get all .log files and strip the .log extension
echo "Getting all LDSC files..."
ls "$LDSC_DIR"/*.log | xargs -n 1 basename | sed 's/\.log$//' > "$TEMP_ALL_FILES"

# Extract valid keys from CSV (third column)
echo "Extracting valid keys from CSV..."
cut -d',' -f3 "$CSV_FILE" > "$TEMP_VALID_KEYS"

# Find files that are not in the valid keys list
echo "Finding unused files..."
grep -v -f "$TEMP_VALID_KEYS" "$TEMP_ALL_FILES" > "$UNUSED_FILE"

# Count unused files
unused_count=$(wc -l < "$UNUSED_FILE")
echo "Found $unused_count unused files"

if [ $unused_count -gt 0 ]; then
    echo "Unused files:"
    while IFS= read -r key; do
        echo "$LDSC_DIR/$key.log"
    done < "$UNUSED_FILE"
fi

# # Clean up temporary files
# rm "$TEMP_ALL_FILES" "$TEMP_VALID_KEYS"

# # Uncomment below to actually delete the files
# if [ $unused_count -gt 0 ]; then
#     echo "Deleting unused files..."
#     while IFS= read -r key; do
#         rm -v "$LDSC_DIR/$key.log"
#     done < "$UNUSED_FILE"
#     echo "Unused files have been deleted"
# else
#     echo "No unused files found"
# fi 