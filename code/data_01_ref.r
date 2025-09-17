# ================================================================
# data
# ================================================================
# bfile
# https://zenodo.org/record/7768714

# mhc
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
# MHC -- chr6 (NC_000006.11):28477797-33448354

# ld r2 for proxy iv
# https://mrcieu.github.io/TwoSampleMR/articles/outcome.html?q=proxy#ld-proxies
# https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/query.R
# ================================================================
# clean bfile
# ================================================================
# # tgz can further compress tar
# wget https://zenodo.org/records/7768714/files/1000G_Phase3_EAS_plinkfiles.tgz
# wget https://zenodo.org/records/7768714/files/1000G_Phase3_plinkfiles.tgz
# tar -xvf 1000G_Phase3_plinkfiles.tgz 
# tar -xvf 1000G_Phase3_EAS_plinkfiles.tgz 

cd /bigdat/db/sldsc_ref
for pop in eur eas
do
    if [ "$pop" = "eur" ]; then path="1000G_EUR_Phase3_plink/"
    else path="1000G_Phase3_EAS_plinkfiles/"; fi
    # merge
    cd $path
    plink2 --bfile 1000G.${pop^^}.QC.1 --pmerge-list merge.txt --make-bed --out ${pop}
    plink --bfile $pop --freq --out $pop
    # exclude mhc and maf < 0.01 snp
    awk '$5 < 0.01 {print $2}' ${pop}.frq > maf0.01_exclude_snps.txt
    awk '$4 > 28477797 && $4 < 33448354 {print $2}' 1000G.${pop^^}.QC.6.bim > mhc_exclude_snps.txt
    cat maf0.01_exclude_snps.txt mhc_exclude_snps.txt > exclude_snps.txt
    plink --bfile $pop --exclude exclude_snps.txt --make-bed --out ${pop}_nomhc_maf0.01
    plink --bfile ${pop}_nomhc_maf0.01 --freq --out ${pop}_nomhc_maf0.01
    cd ..
done

# in frq, a1 is minor
# df1 = read.table('./eur_nomhc_maf0.01.afreq')
# df2 = read.table('./eur_nomhc_maf0.01.frq', header=1)
# names(df1) = c('CHR', 'POS', 'A2', 'A1', 'PROVISIONAL_REF', 'FRQ', 'NCHROBS')
# t=df1$FRQ-df2$MAF


# ================================================================
# add pos to frq, lift 
# ================================================================
# library(rtracklayer)
# chain <- import.chain('/home/zhang/db/chain/hg19ToHg38.over.chain')

# paths_ref = list('eas' = '~/db/sldsc_ref/1000G_Phase3_EAS_plinkfiles/eas_nomhc_maf0.01',
#     'eur' = '~/db/sldsc_ref/1000G_EUR_Phase3_plink/eur_nomhc_maf0.01')

# frqs = list()
# for (pop in c('eas', 'eur')){
#     bim = read.table(sprintf("%s.bim", paths_ref[[pop]]))%>%dplyr::rename(SNP=V2, POS=V4)%>%select(SNP, POS)
#     frq = read.table(sprintf("%s.frq", paths_ref[[pop]]), header=1)%>%dplyr::rename(FRQ=MAF)%>%select(SNP, CHR, A1, A2, FRQ)
#     frq = frq%>%merge(bim, by='SNP')
#     # lift
#     gr = with(frq, GRanges(seqnames=paste0('chr', CHR), ranges=IRanges(start=POS, end=POS), id=SNP))
#     map <- as.data.frame(liftOver(gr, chain))
#     map = map%>%dplyr::rename(POS_38=start, SNP=id)%>%select(SNP, POS_38)
#     frq = frq%>%merge(map, by='SNP', all.x=T)
#     frqs[[pop]] = frq
# }


# save(frqs, file = '~/db/sldsc_ref/clean/nomhc_maf0.01.rdata')
# load('~/db/sldsc_ref/clean/nomhc_maf0.01.rdata')

#================================================================
# lift to hg38
#================================================================
library(dplyr)
for (pop in c('eas', 'eur')){
    path_prefix = paste0('/bigdat/db/sldsc_ref/', 
        ifelse(pop=='eur', "1000G_EUR_Phase3_plink/", '1000G_Phase3_EAS_plinkfiles/'), pop)
    frq = read.table(paste0(path_prefix, '.frq'), header=1)%>%select(-NCHROBS)
    bim = read.table(paste0(path_prefix, '.bim'))%>%select(V2, V4)%>%rename(SNP=V2, POS=V4)
    info = frq%>%merge(bim, by='SNP')
    write.table(info, paste0(path_prefix, '.temp'), sep='\t', row.names=F, quote=F)
} 

for pop in eur eas
do
    if [ "$pop" = "eur" ]; then path="/bigdat/db/sldsc_ref/1000G_EUR_Phase3_plink/"
    else path="/bigdat/db/sldsc_ref/1000G_Phase3_EAS_plinkfiles/"; fi

    cd /bigdat/soft/easylift

    python ./code/easylift.py \
    --lift hg19tohg38 \
    --chr_col CHR --pos_col POS \
    --file_in $path${pop}.temp \
    --file_out $path${pop}_hg38.info

done