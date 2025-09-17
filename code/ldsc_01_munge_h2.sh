#================================================================
# munge and h2
#================================================================
conda activate ldsc
file_hm3="/home/zhang/soft/ldsc/w_hm3.snplist"
path_ld="/home/zhang/db/sldsc_ref/eur_w_ld_chr/"

cd /bigdat/db/gwas/

for source in other
do 
  # move to source path
  cd $source
  mkdir -p result/h2
  mkdir -p munge
  for trait in $(tail -n +2 info/info.csv | cut -d, -f1 | sed 's/^"//; s/"$//'); do
    echo "$trait"
    # snp < 5e-8
    if [ ! -s result/h2/$trait.log ]; then
    # munge
    ~/soft/ldsc/munge_sumstats.py \
    --sumstats clean/$trait.txt.gz \
    --merge-alleles $file_hm3 \
    --out munge/$trait
    # single trait heritability
    ~/soft/ldsc/ldsc.py \
    --h2 munge/$trait.sumstats.gz \
    --ref-ld-chr $path_ld \
    --w-ld-chr $path_ld \
    --out result/h2/$trait
    fi; done
    # move back
  cd ../ 
done 