libs = c('dplyr', 'R.utils', 'stringr', 'TwoSampleMR', 'data.table')
lapply(libs, require, character.only = TRUE) 
options(stringsAsFactors=F)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

path_ref = '/bigdat/db/sldsc_ref/1000G_EUR_Phase3_plink/eur_nomhc_maf0.01'
path_mm = '/bigdat/soft/public_mm_lite/'

# mr data
methods <- c("ivw", "egger", "wmedian", "wmode") # combine Wald ratio with IVW, when nsnp=1
link_names <- c(
    "exposure", "outcome", "nsnp", paste0(rep(methods, each = 3), c("_b", "_se", "_p")),
    "egger_int", "egger_int_se", "egger_int_p", "ivw_q", "ivw_q_p", "egger_q", "egger_q_p"
)


is_numeric <- function(x) {
   !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}

## extract element with pattern from vector
# e.g., get('a', c('AS', 'Ba', 'C')) ; get('a', c('AS', 'Ba', 'C'), exact=F)
get <- function(key, vec, exact = F, v = F) {
    vec <- unlist(vec)
    if (exact == T) {
        out <- vec[grep(key, vec, invert = v)]
    } else {
        out <- vec[grep(toupper(key), toupper(vec), invert = v)]
    }
    return(out)
}



## replace element in df
# map(data.frame(t=c('a', 'b')), c('a', 'b'), c('xxx', 'yyy'))
map <- function(df, keys, values) {
    for (i in 1:length(keys)) {
        df[df == keys[i]] <- values[i]
    }
    return(df)
}

# umls
# get_umls('diabetes', path_mm)
get_umls = function(trait, path_mm){
    command = sprintf('cd %s ; echo "%s" | ./metamaplite.sh --pipe ', path_mm, trait)
    map = system(command, intern=T)
    if (length(map)==0){umls = NA}
    umls = strsplit(map[1], '\\|')[[1]][4]
    return(umls)
}

# get_umls in info df, save log
run_umls = function(df, id_col='id', trait_col='trait', path_out='../umls/'){
    df$umls = NA
    for (i in 1:nrow(df)){
        file_umls = sprintf('%s%s.rdata', path_out, df[i, id_col])
        if (!file.exists(file_umls)){
            command = sprintf('cd %s ; echo "%s" | ./metamaplite.sh --pipe ', path_mm, df[i, trait_col])
            map = system(command, intern=T)
            save(map, file=file_umls)
        } 
        load(file_umls)
        if (length(map)==0){umls = NA}
        umls = strsplit(map[1], '\\|')[[1]][4]
        df[i, 'umls'] = umls
        rm('umls', 'map')
    }
    return(df)
}

# read 1000g frq
get_1000gfrq = function(){
    bim = read.table(sprintf("%s.bim", path_ref))%>%rename(SNP=V2, POS=V4)%>%select(SNP, POS)
    frq = read.table(sprintf("%s.frq", path_ref), header=1)%>%rename(FRQ=MAF)%>%select(SNP, CHR, A1, A2, FRQ)
    frq = frq%>%merge(bim, by='SNP')
    return(frq)
}

# filter gwas to 1000g ref, remove mhc, snp with maf < 0.01, map a1 a2, lift to hg 19.
filter_gwas = function(gwas, frq, use='SNP'){
    if (use == 'SNP'){
        df = gwas%>%filter(SNP%in%frq$SNP)%>%select(SNP, A1, A2, BETA, SE, P, N)
        df1 = df%>%merge(frq, by=c('SNP', 'A1', 'A2'))
        df2 = df%>%rename(A1=A2, A2=A1)%>%mutate(BETA=-BETA)%>%merge(frq, by=c('SNP', 'A1', 'A2'))
    }
    if (use == 'POS'){
        df = gwas%>%select(CHR, POS, A1, A2, BETA, SE, P, N)
        df1 = df%>%merge(frq, by=c('CHR', 'POS', 'A1', 'A2'))
        df2 = df%>%rename(A1=A2, A2=A1)%>%mutate(BETA=-BETA)%>%merge(frq, by=c('CHR', 'POS', 'A1', 'A2'))
    }
    if (use == 'POS_38'){
        df = gwas%>%select(CHR, POS_38, A1, A2, BETA, SE, P, N)
        df1 = df%>%merge(frq, by=c('CHR', 'POS_38', 'A1', 'A2'))
        df2 = df%>%rename(A1=A2, A2=A1)%>%mutate(BETA=-BETA)%>%merge(frq, by=c('CHR', 'POS_38', 'A1', 'A2'))
    }
    df_out = rbind(df1, df2)%>%select(SNP, CHR, POS, A1, A2, FRQ, BETA, SE, P, N)%>%na.omit()
    return(df_out)
}



# read ieu eqtl vcf
read_vcf = function(file_in){
    lines = readLines(file_in)
    lines = lines[!grepl("^##", lines)]
    if (length(lines)>1e7){ # if lines have too many elements, the paste0 function will error
        split_index = length(lines) %/% 2
        data_str1 = paste0(lines[1:split_index], collapse = "\n")
        data_str2 = paste0(lines[c(1, (split_index + 1):length(lines))], collapse = "\n") # add header
        gwas1 = read.table(text=data_str1, header=1, sep="\t", comment.char="@")
        gwas2 = read.table(text=data_str2, header=1, sep="\t", comment.char="@")
        gwas = rbind(gwas1, gwas2)
    } else{
        data_str = paste0(lines, collapse = "\n")
        gwas = read.table(text=data_str, header=1, sep="\t", comment.char="@")
    }
    n_format = length(strsplit(gwas$FORMAT[1], ':')[[1]])
    cols = data.frame(str_split_fixed(gwas[,ncol(gwas)], ":", n_format)) # LP is -log10 P
    names(cols) =  strsplit(gwas$FORMAT[1], ':')[[1]]
    cols$ID = NULL
    gwas = cbind(gwas, cols)
    gwas = gwas%>%rename(CHR=X.CHROM, SNP=ID, A1=ALT, A2=REF, BETA=ES)%>%
        mutate_if(is_numeric, as.numeric)%>%mutate(P=exp(-1*LP))
    gwas$FRQ = NULL
    return(gwas)
}


# ================================================================
# mr functions
# ================================================================

get_mr_res <- function(df1, df2, snp) {
  mr_res <- list()
  mr_res[["stat"]] <- c(from, to, rep(NA, 20)) # if not dat then all N
  sub1 <- df1 %>%select(-CHR, - POS)%>%filter(SNP%in%snp)%>%
      mutate(id.exposure = from, exposure = from) %>%
      rename("pval.exposure" = "P", "effect_allele.exposure" = "A1", "other_allele.exposure" = "A2", "samplesize.exposure" = "N", "beta.exposure" = "BETA", "se.exposure" = "SE", "eaf.exposure" = "FRQ") 
  sub2 <- df2 %>%select(-CHR, - POS)%>%
      filter(SE>0) %>% # some ukb gwas have 0 SE and P
      mutate(id.outcome = to, outcome = to) %>%
      rename("pval.outcome" = "P", "effect_allele.outcome" = "A1", "other_allele.outcome" = "A2", "samplesize.outcome" = "N", "beta.outcome" = "BETA", "se.outcome" = "SE", "eaf.outcome" = "FRQ")
  dat <- harmonise_data(sub1, sub2) %>% filter(mr_keep == T) %>%
    select(exposure, outcome, id.exposure, id.outcome, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure, eaf.outcome,
    beta.exposure, beta.outcome, se.exposure, se.outcome, pval.exposure, pval.outcome, samplesize.exposure, samplesize.outcome, mr_keep)
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")) %>%
      select(method, nsnp, b, se, pval)
  pleio <- mr_pleiotropy_test(dat)
  het <- mr_heterogeneity(dat)
  res <- map(res, c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode"), c("ivw", "egger", "wmedian", "wmode"))
  # collect
  stat <- c(from, to, res$nsnp[1])
  for (method in methods) {
      temp <- unlist(res %>% filter(method == UQ(method)) %>% select(b, se, pval))
      stat <- c(stat, temp)
  }
  stat <- c(
      stat, unlist(pleio[, c("egger_intercept", "se", "pval")]),
      unlist(het %>% filter(method == "Inverse variance weighted") %>% select(Q, Q_pval)),
      unlist(het %>% filter(method == "MR Egger") %>% select(Q, Q_pval))
  )
  stat[3:length(stat)] <- sprintf("%.2e", as.numeric(stat[3:length(stat)]))
  mr_res[["stat"]] <- stat # update
  mr_res[["dat"]] <- dat
  names(mr_res[["stat"]]) <- link_names
  if (length(mr_res[["stat"]]) != 22) {
      print("miss stat, break!")
      break
  }
  return(mr_res)
}

# ================================================================
# ldsc functions
# ================================================================
remove_brackets_split = function(x){
  if (length(line) >0) {
    x = str_squish(x) # multiple spaces to single
    x = gsub(')', '', gsub('(', '', x, fixed=T), fixed=T)
    x = strsplit(x, ' ')[[1]]
  }
  return(x)
}



collect_ldsc = function(file_in){
    out = c()
    lines = readLines(file_in)
    # break if incompelte
    if (length(lines) < 62) {
        stop("Error: The file does not contain enough lines (minimum 62 required).")
    }
    line = get('with valid alleles', lines); items = remove_brackets_split(line)
    nsnp = items[1]
    line = get('scale h2', lines)[1]; items = remove_brackets_split(line)
    h2_1 = as.numeric(items[5]); h2_1_se = as.numeric(items[6]); h2_1_p=2*pnorm(abs(h2_1/h2_1_se), lower.tail=F)
    line = get('scale h2', lines)[2]; items = remove_brackets_split(line)
    h2_2 = as.numeric(items[5]); h2_2_se = as.numeric(items[6]); h2_2_p=2*pnorm(abs(h2_2/h2_2_se), lower.tail=F)
    line = get('sumstats.gz', lines)[4]; items = remove_brackets_split(line)
    rg = items[3]; rg_se = items[4]; rg_p = items[6]; gcov_int = items[11]; gcov_int_se = items[12]
    out = c(file_in, nsnp, h2_1, h2_1_se, h2_1_p, h2_2, h2_2_se, h2_2_p, gcov_int, gcov_int_se, rg, rg_se, rg_p)
    names(out) = c('file', 'nsnp', 'h2_1', 'h2_1_se', 'h2_1_p', 'h2_2', 'h2_2_se', 'h2_2_p', 'gcov_int', 'gcov_int_se', 'rg', 'rg_se', 'rg_p')
    return(out)
}


collect_ldsc_h2 = function(path, files){
    out = c()
    for (file in files){
        cat('find', file, 'h2 ... \n')
        item = scan(paste0(path, file), character(), sep='\t')
        item = item[grepl('h2:', item)]
        item = suppressWarnings(as.numeric(strsplit(item, ' |[()]')[[1]]))
        item = item[!is.na(item)]
        item[3] = 2*pnorm(abs(item[1]/item[2]), lower.tail=FALSE)
        out = c(out, gsub('.log', '', file),  item)
    }
    collect = data.frame(matrix(out, ncol=4, byrow=T))
    colnames(collect) = c('trait', 'h2', 'se', 'p')
    return(collect)
}
