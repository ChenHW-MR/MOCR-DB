
#================================================================
# info: ukb
#================================================================
# https://yanglab.westlake.edu.cn/software/gcta/#UKBiobankGWASresults

# source
setwd('/bigdat/db/gwas/ukb/info')
source('~/easymr/code/function.r')
keep = "(?i)thyroi" # keep trait
# read
df1 = read.csv('./UKB_impute_v1.1.csv') # https://yanglab.westlake.edu.cn/software/gcta/res/UKB_impute_v1.1.csv
df2 = read.csv('./UKB_binary_v1.11.csv') # https://yanglab.westlake.edu.cn/software/gcta/res/UKB_binary_v1.11.csv
# drop df1 sex spc
df1 = df1%>%filter(is.na(Gender_specific)|Gender_specific=='')%>%select(-Gender_specific)
# drop df1 dup
df1 = as.data.frame(df1%>%group_by(Description)%>%filter(N==max(N))%>%ungroup())
df1 = df1%>%filter(!Description%in%df2$Description) # df1 binary duplicate with df2
# df2 have dummy trait with dup description, keep them
head(df2%>%filter(duplicated(Description)|duplicated(Description, fromLast=T))%>%arrange(Description)%>%select(Description))
# drop df2 case_prop < 0.01
df2 = df2%>%mutate(case_prop=N_case/(N_case+N_control))
# merge
df1 = df1%>%mutate(N_case=NA, N_control=NA, case_prop=NA)%>%select(-Method, -Ncase)
df2 = df2%>%mutate(Data_type='Binary')%>%select(-ratio)%>%
    filter((case_prop>0.01&case_prop<0.99) | grepl(keep, Description))
df = rbind(df1, df2)
# filter
df = df%>%
    filter(!grepl('(?i)Treatment|Job |Average acceleration|Types of transport used|Data quality|Workplace', Description))
dim(df) # 2040
# clean
df = df%>%rename(id_raw=ID)%>%mutate(id=paste0('ukb_', id_raw))%>%relocate(id)%>%
    rename(n=N, ncase=N_case, ncontrol=N_control, trait=Description)
# run_umls
df = run_umls(df)
dim(df) # 1980
# save
url = df$URL
df = df%>%select(-URL)
write.table(url, './url.txt', row.names=F, col.names=F, quote=F)
write.csv(df, './info.csv', row.names=F) 



#================================================================
# info: fin r11 
#================================================================
# https://finngen.gitbook.io/documentation/v/r10/data-download
# source
setwd('/bigdat/db/gwas/fin/info')
source('~/easymr/code/function.r')
keep = "(?i)thyroi" # keep trait
# read
df = read.delim('./finngen_R11_manifest.tsv')
# drop df2 case_prop < 0.01
df = df%>%rename(trait=phenotype, id_raw=phenocode, ncase=num_cases, ncontrol=num_controls)
df = df%>%mutate(n=ncase+ncontrol, case_prop=ncase/(ncase+ncontrol))%>%
    filter((case_prop>0.01&case_prop<0.99)|grepl(keep, trait))
# clean
url = df%>%select(path_https)
df = df%>%select(-path_https, -path_bucket, -category)%>%mutate(id=paste0('fin_', id_raw), n=ncase+ncontrol)%>%relocate(id)%>%relocate(n, .before=ncase)
# check
dim(df)  # 725 
# run umls
df = run_umls(df)

# save
write.csv(df, './info.csv', row.names=F, quote=T)
write.table(url, './url.txt', row.names=F, col.names=F, quote=F)


# # check
# files = list.files('../raw')
# keep = paste0(df$id_raw, ifelse(df$Data_type=='Binary', '.v1.0.fastGWA.gz', '.v1.1.fastGWA.gz'))
# files[!files%in%keep]


#================================================================
# info swe
#================================================================
# this fill info.csv by hand
cat /home/zhang/db/gwas/swe/info/info.csv

## CD 34+
# https://pubmed.ncbi.nlm.nih.gov/35007327/
# wget -c https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90102001-GCST90103000/GCST90102460/GCST90102460_buildGRCh38.tsv.gz
# cd /bigdat/soft/public_mm_lite/ ; echo \"Food weight\" | ./metamaplite.sh --pipe

#================================================================
# info other
#================================================================
## COVID-19
# https://www.covid19hg.org/results/
# cd /bigdat/soft/public_mm_lite/ ; echo \"COVID-19\" | ./metamaplite.sh --pipe
# fill info.csv by hand
cat /home/zhang/db/gwas/other/info/info.csv

# ieu-b-4903 and ieu-b-4904

#================================================================
# download
#================================================================
## fin
# https://finngen.gitbook.io/documentation/data-download
cd /bigdat/db/gwas/fin/raw/
cat ../info/url.txt | while read url; do
  wget -c "$url"
done

nohup cat ../info/url.txt| xargs -n 1 -P 10 axel -c -n 10 &
nohup tac ../info/url.txt| xargs -n 1 -P 10 axel -c -n 10 &

## ukb
# https://yanglab.westlake.edu.cn/software/gcta/#UKBiobankGWASresults
cd /bigdat/db/gwas/ukb/raw
cat ../info/url.txt | while read url; do
  wget -c "$url"
  sleep 3
done

ps -ef | grep wget

## check
# gunzip -t *.gz




