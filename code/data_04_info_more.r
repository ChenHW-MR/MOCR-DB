# ================================================================
# pair
# ================================================================
setwd("~/easymr")
source("~/easymr/code/function.r")

df <- data.frame()
for (source in c("ukb", "fin", "swe", "other")) {
  info <- read.csv(sprintf("~/db/gwas/%s/info/info.csv", source)) %>%
    select(id, trait, umls, n, case_prop) %>%
    mutate(source = !!source)
  path_h2 <- sprintf("~/db/gwas/%s/result/h2/", source)
  nsnp = read.delim(sprintf('~/db/gwas/%s/info/nsnp_5e-8', source), sep=' ')
  sub <- collect_ldsc_h2(path_h2, paste0(info$id, ".log"))
  info_add <- info %>% merge(sub, by.x = "id", by.y = "trait")%>%merge(nsnp, by.x='id', by.y='trait')
  if (nrow(info_add) != nrow(info)) {
    stop('check source:', source)
  }
  df <- rbind(df, info_add)
}
df = df%>%mutate_if(is_numeric, as.numeric)%>%rename(h2_se=se, h2_p=p)

# impute missing umls with id
df = df%>%mutate(umls = ifelse(!is.na(umls), umls, id))

nrow(df) # 2710

# save
write.csv(df, 'result/info.csv', row.names=F)