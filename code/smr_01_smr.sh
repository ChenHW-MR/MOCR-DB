#================================================================
# info
#================================================================
# eqtl, gtex v8
# wget -c https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_sqtl_summary.html

# mqtl 
# wget -c https://yanglab.westlake.edu.cn/data/SMR/EUR.tar.gz

# sqtl
# wget -c https://yanglab.westlake.edu.cn/data/SMR/BrainMeta_cis_sqtl_summary.tar.gz

# caqtl
# wget -c https://yanglab.westlake.edu.cn/data/SMR/Bryois_caQTL_summary.tar.gz

# #================================================================
# # make cojo format
# #================================================================
# cd ~/db/gwas

# for source in fin ukb swe other; do
#     path_gwas="/bigdat/db/gwas/$source/"
#     mkdir -p ${path_gwas}cojo
#     for trait in $(tail -n +2 ${path_gwas}info/info.csv | cut -d, -f1 | sed 's/^"//; s/"$//'); do
#         echo $trait
#         file_cojo="$path_gwas/cojo/${trait}.ma"
#         if [ ! -s $file_cojo ]; then
#             zcat $path_gwas/clean/$trait.txt.gz | awk '{print $1, $4, $5, $6, $7, $8, $9, $10}' > $file_cojo
#         fi
#     done
# done

# # check if all rows have 8 columns
# for source in fin ukb swe other; do
#     path_gwas="/bigdat/db/gwas/$source/"
#     for f in ${path_gwas}cojo/*.ma; do
#         if ! awk 'NF == 8 { next } { exit 1 }' "$f"; then
#             echo "Deleting $f (some rows do not have exactly 8 columns)"
#             rm "$f"
#         fi
#     done
# done

#================================================================
# smr.sh
#================================================================
#!/bin/bash
file_bfile=/home/zhang/db/sldsc_ref/1000G_EUR_Phase3_plink/eur_nomhc_maf0.01
path_eqtl="/home/zhang/db/smr/GTEx_V8_cis_eqtl_summary/"
path_mqtl="/home/zhang/db/smr/EUR/"
path_sqtl="/home/zhang/db/smr/BrainMeta_cis_sqtl_summary/"
path_caqtl="/home/zhang/db/smr/Bryois_caQTL_summary/"

ntasks=0
for source in ukb fin swe other ; do
    path_gwas="/bigdat/db/gwas/$source/"
    for trait in $(tail -n +2 ${path_gwas}info/info.csv | cut -d, -f1 | sed 's/^"//; s/"$//' | cat); do
        echo $trait
        file_in="$path_gwas/cojo/${trait}.ma"
        for expo_source in eqtl mqtl sqtl caqtl; do
            path=path_$expo_source
            path=${!path}
            expos=$(ls ${path} | grep .besd | sed 's/.besd//g')         
            for expo in $expos; do
                file_qtl=$path$expo
                file_out="/home/zhang/easymr/data/smr/$expo_source/${expo}_TO_${trait}"
                if [ -s $file_out.smr ] || [ -s $file_out.smr.gz ]; then
                    continue
                fi
                if [ $ntasks -ge 20 ]; then
                    wait -n
                    ntasks=$((ntasks - 1))
                fi
                nohup smr_Linux \
                    --bfile $file_bfile \
                    --gwas-summary $file_in \
                    --beqtl-summary $file_qtl \
                    --diff-freq-prop 0.2 \
                    --out $file_out &
                ntasks=$((ntasks + 1))
            done
        done
    done
done




#================================================================
# gzip
#================================================================
#!/bin/bash
path_eqtl="/home/zhang/db/smr/GTEx_V8_cis_eqtl_summary/"
path_mqtl="/home/zhang/db/smr/EUR/"
path_sqtl="/home/zhang/db/smr/BrainMeta_cis_sqtl_summary/"
path_caqtl="/home/zhang/db/smr/Bryois_caQTL_summary/"

for source in ukb fin swe other ; do
    path_gwas="/bigdat/db/gwas/$source/"
    for trait in $(tail -n +2 ${path_gwas}info/info.csv | cut -d, -f1 | sed 's/^"//; s/"$//' | cat); do
        echo $trait
        file_in="$path_gwas/cojo/${trait}.ma"
        for expo_source in eqtl mqtl sqtl caqtl; do
            path=path_$expo_source
            path=${!path}
            expos=$(ls ${path} | grep .besd | sed 's/.besd//g')         
            for expo in $expos; do
                file_out="/home/zhang/easymr/data/smr/$expo_source/${expo}_TO_${trait}"
                if [ -f $file_out.smr ] && [ ! -f $file_out.smr.gz ]; then
                    gzip $file_out.smr
                fi
            done
        done
    done
done




# #================================================================
# # smr_data.sh make smr data for plot
# #================================================================
# file_smr="/home/zhang/easymr/result/smr_collect.csv"
# file_smrTmp="/home/zhang/easymr/result/smr_tmp.csv"

# # Extract the desired columns
# csvcut -c expo_source,expo,source,trait,probeID $file_smr > $file_smrTmp

# # setting
# file_bfile=/home/zhang/db/sldsc_ref/1000G_EUR_Phase3_plink/eur_nomhc_maf0.01
# file_glist="/home/zhang/soft/gcta_1.93.2beta//glist_hg19_refseq.txt"
# path_eqtl="/home/zhang/db/smr/GTEx_V8_cis_eqtl_summary/"
# path_mqtl="/home/zhang/db/smr/EUR/"
# path_sqtl="/home/zhang/db/smr/BrainMeta_cis_sqtl_summary/"
# path_caqtl="/home/zhang/db/smr/Bryois_caQTL_summary/"


# ntasks=0

# for line in $(tail -n +2 "$file_smrTmp"); do
#     # Split the line into fields based on the comma delimiter
#     IFS=',' read -r expo_source expo source trait probeID <<< "$line"
#     # Remove potential surrounding quotes
#     expo_source="${expo_source//\"/}"
#     expo="${expo//\"/}"
#     source="${source//\"/}"
#     trait="${trait//\"/}"
#     probeID="${probeID//\"/}"
#     # generate smr data
#     file_gwas="/bigdat/db/gwas/${source}/cojo/${trait}.ma"
#     path_qtl=path_$expo_source
#     file_qtl=${!path_qtl}$expo
#     file_out="/home/zhang/easymr/data/smr/$expo_source/${expo}_TO_${trait}"
#     file_out1="/home/zhang/easymr/data/smr/$expo_source/plot/${expo}_TO_${trait}.${probeID}.txt" # output will save in a plot folder
#     if [ -s $file_out1 ]; then
#         continue
#     fi
#     if [ $ntasks -ge 10 ]; then
#         wait -n
#         ntasks=$((ntasks - 1))
#     fi
#     nohup smr_Linux \
#         --bfile $file_bfile \
#         --gwas-summary $file_gwas \
#         --beqtl-summary $file_qtl \
#         --probe $probeID \
#         --out $file_out --plot --probe-wind 500 --gene-list $file_glist &
#     ntasks=$((ntasks + 1))
# done

