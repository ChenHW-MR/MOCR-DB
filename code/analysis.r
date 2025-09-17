
# ================================================================
# collect ldsc result for analysis paper
# ================================================================
setwd("~/easymr")
source("~/easymr/code/db_gwas/prepare.r")

# pair
pair <- read.csv("info/ldsc_pair.csv")
pair = pair%>%filter(!(grepl('swe_', trait1)|grepl('swe_', trait2)))
write.csv(pair, "result/analysis/ldsc_pair.csv", row.names = F)

# result
res1 <- read.csv("result/ldsc_collect.csv")
res2 <- res1 %>% filter(!grepl("swe_", trait1) & !grepl("swe_", trait2)) %>% filter(h2_1_p<0.05 & h2_2_p<0.05)
write.csv(res2, "result/analysis/ldsc_collect.csv", row.names = F)


# ================================================================
# collect mr result for analysis paper
# ================================================================
setwd("~/easymr")
source("~/easymr/code/db_gwas/prepare.r")

# pair
pair <- read.csv("info/mr_pair.csv")
pair1 = pair%>%filter(!(grepl('swe_', from)|grepl('swe_', to)))
write.csv(pair1, "result/analysis/mr_pair.csv", row.names = F)

# result
res1 = read.csv("result/mr_collect.csv")
res2 <- res1 %>% filter(!grepl("swe_", exposure) & !grepl("swe_", outcome)) %>% filter(nsnp>=10)
write.csv(res2, "result/analysis/mr_collect.csv", row.names = F)

# ================================================================
# zip 
# ================================================================
# cd ~/easymr/data/

# tar -cvzf harmo_dat.tar.gz harmo_dat/
# tar -cvzf ldsc.tar.gz ldsc/



# ================================================================
# tmp 
# ================================================================
h2 = read.csv("result/analysis/h2.csv")
ldsc_pair <- read.csv("result/analysis/ldsc_pair.csv")
mr_pair <- read.csv("result/analysis/mr_pair.csv")
mr_collect <- read.csv("result/analysis/mr_collect.csv")
# find mr_pair in ldsc_pair
mr_collect = mr_collect%>%mutate(key1 = paste0(exposure, "_ANDAND_", outcome), key2 = paste0(outcome, "_ANDAND_", exposure))
idx1 = mr_collect$key1%in%ldsc_pair$key
idx2 = mr_collect$key2%in%ldsc_pair$key

sum(idx1 | idx2) # 6553


# grouping
ukb_info = read.csv('/home/zhang/db/gwas/ukb/info/info.csv')
fin_info = read.csv('/home/zhang/db/gwas/fin/info/info.csv')



traits = h2$trait

traits1 = get('intake|eat', traits)

remain = traits[!traits%in%traits1]

traits2 = get('ICD|disorder|diseases|hypertension|Heart', traits)

remain = remain[!remain%in%c(traits2)]