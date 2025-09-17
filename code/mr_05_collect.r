#================================================================
# collect result
#================================================================
setwd('~/easymr')
source('~/easymr/code/function.r')

pair = read.csv('result/mr_pair.csv')

res = data.frame()
for (i in 1:nrow(pair)){
  if (i%%100==0){
    print(sprintf('%s/%s: nrow res = %.f', i, nrow(pair), nrow(res)))
  }
  from = pair[i, 'from']
  to = pair[i, 'to']
  file_mr = sprintf('data/harmo_dat/%s_TO_%s.rdata', from, to)
  load(file_mr)
  stat = mr_res[['stat']]
  condition1 = as.numeric(stat['ivw_p'])<0.05
  condition2 = length(unique(sign(as.numeric(stat[get('_b', names(stat))]))))==1
  condition3 = all(as.numeric(stat[c('egger_int_p', 'egger_q_p', 'ivw_q_p')])>0.05)
  if (all(condition1, condition2, condition3)){
    # Convert stat to data.frame and preserve column names
    stat_df = as.data.frame(t(stat))
    colnames(stat_df) = names(stat)
    res = rbind(res, stat_df)
  }
  rm('stat', 'mr_res')
}

# add info
map = data.frame()
for (source in c('ukb', 'fin', 'swe', 'other')){
    info = read.csv(sprintf('~/db/gwas/%s/info/info.csv', source))%>%select(id, trait, n, case_prop, umls)
    map = rbind(map, info)
}

res1 = res%>%merge(map%>%rename(exposure=id, exposure_name=trait, exposure_n=n, exposure_case_prop=case_prop, exposure_umls=umls), by='exposure')
res1 = res1%>%merge(map%>%rename(outcome=id, outcome_name=trait, outcome_n=n, outcome_case_prop=case_prop, outcome_umls=umls), by='outcome')
write.csv(res1, 'result/mr_collect.csv', row.names=F)

# ================================================================
# extract interest trait
# ================================================================
res1 = read.csv('result/mr_collect.csv')

# res2 = res1%>%filter(grepl('swe_', outcome, ignore.case = T)|grepl('swe_', exposure, ignore.case = T))
# write.csv(res2, 'result/mr_collect_swe.csv', row.names=F)


#================================================================
# collect data
#================================================================
setwd('~/easymr')
source('~/easymr/code/function.r')

res = read.csv('result/mr_collect.csv')
for (i in 1:nrow(res)){
  print(i)
  file_in = sprintf('data/harmo_dat/%s_TO_%s.rdata', res[i, 'exposure'], res[i, 'outcome'])
  file_out = sprintf('data/harmo_dat_sig/%s_TO_%s.rdata', res[i, 'exposure'], res[i, 'outcome'])
  if (file.exists(file_out)){
    next
  }
  file.copy(file_in, file_out)
}

# gzip
tar -cvzf harmo_dat_sig.tar.gz harmo_dat_sig/
tar -cvzf harmo_dat.tar.gz harmo_dat/