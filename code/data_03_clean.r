# ================================================================
# clean fin
# ================================================================
## clean_fin.r
setwd('/bigdat/db/gwas/fin')
source('~/easymr/code/function.r')
load('~/db/sldsc_ref/clean/nomhc_maf0.01.rdata')
file_info = 'info/info.csv'
info = read.csv(file_info)
frq = frqs[['eur']]
# arg
args = commandArgs(T)
start = as.numeric(args[1]); end = as.numeric(args[2])
range = c(start:end)
# range = 1:nrow(info)
for (i in range){
    print(i)
    id_raw = info[i, 'id_raw']; n = info[i, 'n']; id = info[i, 'id']
    file_in = paste0('./raw/finngen_R11_', id_raw, '.gz')
    # file.remove(file_in)
    file_out = paste0('./clean/', id, '.txt')
    if (file.exists(paste0(file_out, '.gz'))){next}
    gwas = read.delim(file_in, header=1, comment.char="@")
    gwas = gwas%>%rename(SNP= rsids, CHR=X.chrom, POS=pos, A2=ref, A1=alt, BETA=beta, SE=sebeta, P=pval)%>%mutate(N=n)
    gwas = filter_gwas(gwas, frq, use='SNP')
    write.table(gwas, file_out, row.names=F, quote=F, sep='\t')
    gzip(file_out)
    rm(gwas)
}

# # run
# for ((i=1; i<=700; i+=30)); do
#   nohup Rscript clean_fin.r $i $((i+29)) &
# done


# ================================================================
# clean ukb
# ================================================================
## clean_ukb.r
setwd('/bigdat/db/gwas/ukb')
source(paste0('~/easymr/code/function.r'))
load('~/db/sldsc_ref/clean/nomhc_maf0.01.rdata')
file_info = 'info/info.csv'
info = read.csv(file_info)
frq = frqs[['eur']]
# arg
args = commandArgs(T)
start = as.numeric(args[1]); end = as.numeric(args[2])
range = c(start:end)
# range = 1:nrow(info)
for (i in range){
    print(i)
    id_raw = info[i, 'id_raw']; type = info[i, 'Data_type']; id = info[i, 'id']; n=info[i, 'n']
    suffix = ifelse(type=='Binary', '.v1.0.fastGWA.gz', '.v1.1.fastGWA.gz')
    file_in = paste0('./raw/', id_raw, suffix)
    file_out = paste0('./clean/', id, '.txt')
    if (file.exists(paste0(file_out, '.gz'))){next}
    gwas = read.delim(file_in, header=1, comment.char="@")
    gwas = gwas%>%mutate(N=n)
    gwas = filter_gwas(gwas, frq, use='SNP')
    write.table(gwas, file_out, row.names=F, quote=F, sep='\t')
    gzip(file_out)
    rm(gwas)
}

# # run
# for ((i=1800; i<=2000; i+=20)); do
#   nohup Rscript clean_ukb.r $i $((i+19)) &
# done


# ================================================================
# clean swe
# ================================================================
setwd('/bigdat/db/gwas/swe')
source(paste0('~/easymr/code/function.r'))
load('~/db/sldsc_ref/clean/nomhc_maf0.01.rdata')
file_info = 'info/info.csv'
info = read.csv(file_info)
frq = frqs[['eur']]
for (i in 1){
    print(i)
    id_raw = info[i, 'id_raw']; id = info[i, 'id']; n=info[i, 'n']
    file_in = paste0('./raw/', id_raw)
    file_out = paste0('./clean/', id, '.txt')
    if (file.exists(paste0(file_out, '.gz'))){next}
    gwas = read.delim(file_in, header=1, comment.char="@")
    gwas = gwas%>%mutate(N=n)%>%rename(SNP=variant_ID, A1=effect_allele, A2=other_allele, BETA=beta, P=p_value)
    gwas = gwas%>%mutate(SE = sqrt(BETA^2 / qchisq(P, df=1, lower=F))) # imputate SE, https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS9.html
    gwas = filter_gwas(gwas, frq, use='SNP')
    write.table(gwas, file_out, row.names=F, quote=F, sep='\t')
    gzip(file_out)
    rm(gwas)
}

# ================================================================
# clean other
# ================================================================
setwd('/bigdat/db/gwas/other')
source(paste0('~/easymr/code/function.r'))
load('~/db/sldsc_ref/clean/nomhc_maf0.01.rdata')
frq = frqs[['eur']]
file_info = 'info/info.csv'
info = read.csv(file_info)

for (i in 1:nrow(info)){
    print(i)
    id_raw = info[i, 'id_raw']; id = info[i, 'id']; n=info[i, 'n']
    file_in = paste0('./raw/', id_raw)
    file_out = paste0('./clean/', id, '.txt')
    if (file.exists(paste0(file_out, '.gz'))){next}
    n = info[i, 'n']
    if (grepl(".vcf.gz", file_in)){
      lines <- readLines(file_in)
      lines = lines[!grepl("^##", lines)]
      data_str = paste0(lines, collapse = "\n")
      gwas = read.table(text=data_str, header=1, sep="\t", comment.char="@")
      info_cols = str_split_fixed(gwas[,10], ':', 5)
      gwas = gwas%>%rename(CHR = "X.CHROM", A1 = ALT, A2 = REF, SNP = ID)%>%
        mutate(
          BETA = as.numeric(info_cols[,1]), SE = as.numeric(info_cols[,2]), LP = as.numeric(info_cols[,3]),
          P = 10^(-LP)
      )
    } else {
        gwas = read.delim(file_in, header=1)
        gwas = gwas%>%rename(CHR=chr, POS=pos, A1=effect.allele, A2=other.allele, BETA=beta, SE=se,P=pval)
    }
    gwas = gwas%>%mutate(N=n)
    gwas = filter_gwas(gwas, frq, use='SNP')
    write.table(gwas, file_out, row.names=F, quote=F, sep='\t')
    gzip(file_out)
    rm(gwas)
}









# ================================================================
# clean gtex
# ================================================================
## clean_gtex.r
setwd('/bigdat/db/gwas/gtex/')
source(paste0('~/easymr/code/function.r'))
load('~/db/sldsc_ref/clean/nomhc_maf0.01.rdata')
file_info = 'info/info.csv'
info = read.csv(file_info)
frq = frqs[['eur']]

tissues = list.files('raw')
tissues = gsub('.allpairs.txt.gz', '', get('txt.gz', tissues))

# arg
args = commandArgs(T)
range = as.numeric(args[1])
# range = 1:nrow(info)

for (tissue in tissues[range]){
  df = fread(sprintf('raw/%s.allpairs.txt.gz', tissue))
  df$N = info%>%filter(tissue==!!tissue)%>%pull(n)
  genes = unlist(unique(df$gene_id))
  for (gene in unique(df$gene_id)){
    dir.create(sprintf('./clean/%s', tissue))
    file_out = sprintf('./clean/%s/%s.txt', tissue, gene)
    if (file.exists(paste0(file_out, '.gz'))){next}
    gwas = df%>%filter(gene_id==gene)
    cols = str_split_fixed(gwas$variant_id, '_', 5)[,1:4] # The variant IDs are in the format chromosome_position_ref_alt_build
    cols[,1] = gsub('chr', '', cols[,1])
    gwas = gwas%>%mutate(CHR=cols[,1], POS_38=cols[,2], A1=cols[,4], A2=cols[,3])%>%
      rename(BETA=slope, SE=slope_se, P=pval_nominal)%>%
      mutate_if(is_numeric, as.numeric)
    gwas = filter_gwas(gwas, frq, use='POS_38')
    write.table(gwas, file_out, row.names=F, quote=F, sep='\t')
    gzip(file_out)
    rm(gwas)
  }
} 

# # run
# for ((i=1; i<=49; i+=1)); do
#   nohup Rscript clean_gtex.r $i &
# done

# ### check
# for i in {2200..2400}; do echo $i; gzip -t fin_$i.txt.gz; done
# du -sh * | sort


#================================================================
# check empty
#================================================================
# find . -type f -empty