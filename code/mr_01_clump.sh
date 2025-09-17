
#================================================================
# clump
#================================================================
file_bfile="/bigdat/db/sldsc_ref/1000G_EUR_Phase3_plink/eur_nomhc_maf0.01"
file_nsnp="info/nsnp_5e-8"

cd /bigdat/db/gwas/
for source in other; do 
  # move to source path
  cd $source
  mkdir -p clump/input
  mkdir -p clump/output
  echo 'trait nsnp_sig nsnp_clump' > $file_nsnp
  for trait in $(tail -n +2 info/info.csv | cut -d, -f1 | sed 's/^"//; s/"$//'); do
    echo "$trait"
    file_sig="clump/input/$trait"
    file_clump="clump/output/$trait"
    # snp < 5e-8
    zcat clean/$trait.txt.gz | awk -v OFS='\t' 'NR==1 || $9 < 5e-8' > $file_sig
    # clump 
    nsnp_sig=$(( $(wc -l < "$file_sig") - 1 ))
    nsnp_clump=0
    if [ $nsnp_sig -ge 10 ] && [ ! -s "${file_clump}.clumped" ]; then
        plink \
        --clump $file_sig \
        --clump-kb 10000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.001 \
        --bfile $file_bfile   --out $file_clump
        nsnp_clump=$(($(grep -c -v '^$' "${file_clump}.clumped") - 1))
    fi
    echo $trait $nsnp_sig $nsnp_clump >> $file_nsnp
  done
  cd ..
done
