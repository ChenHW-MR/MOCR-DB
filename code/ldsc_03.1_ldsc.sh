#================================================================
# ldsc.sh
#================================================================
#!/bin/bash
path_ld=/home/zhang/db/sldsc_ref/1000g_ldscore_hg19/eur_w_ld_chr/
start=$1
end=$2
for (( i=$start; i<=$end; i++ )); do
  line=$(sed -n "${i}p" ~/easymr/result/ldsc_pair_torun.csv)
  IFS=, read -r trait1 trait2 key <<< "$line"
  log_file=~/easymr/data/ldsc/${key}.log
  # check exist and complete
  if [ -f "$log_file" ] && [[ $(wc -l < "$log_file") -ge 65 ]]; then
    echo "Skipping: $key (log file exists and is complete)"
    continue
  fi
  ~/soft/ldsc/ldsc.py \
    --rg /home/zhang/db/gwas/"${trait1%%_*}"/munge/${trait1}.sumstats.gz,/home/zhang/db/gwas/"${trait2%%_*}"/munge/${trait2}.sumstats.gz \
    --ref-ld-chr "$path_ld" \
    --w-ld-chr "$path_ld" \
    --out ~/easymr/data/ldsc/${key}
done


data/ldsc/other_COVID19_HGI_B2_EUR_AND_swe_CD34+.log:

trait1=other_COVID19_HGI_B2_EUR
trait2=swe_CD34+
key=other_COVID19_HGI_B2_EUR_AND_swe_CD34+
log_file=~/easymr/data/ldsc/${key}.log

~/soft/ldsc/ldsc.py \
    --rg /home/zhang/db/gwas/"${trait1%%_*}"/munge/${trait1}.sumstats.gz,/home/zhang/db/gwas/"${trait2%%_*}"/munge/${trait2}.sumstats.gz \
    --ref-ld-chr "$path_ld" \
    --w-ld-chr "$path_ld" \
    --out ~/easymr/data/ldsc/${key}




