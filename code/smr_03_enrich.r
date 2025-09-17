library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(data.table)
library(dplyr)
library(purrr)
setwd("~/easymr")
source("code/smr_00_functions.r")


# Read data
df = read.csv("result/smr_collect.csv.gz")
df = df%>%dplyr::select(expo_source, expo, trait, ENSEMBL)

# Get all unique combinations of exposure source, expo, and trait
combinations <- df %>%
  distinct(expo_source, expo, trait) %>%
  arrange(expo_source, expo, trait)
print(paste("Total combinations found:", nrow(combinations)))
# ===============================
# prepare term2gene mapping
# ===============================

# Prepare GO mappings
go_mappings <- prepare_go_mappings()
ont_list <- go_mappings$ont_list
ont_name_list <- go_mappings$ont_name_list 

# ===============================
# run separately
# ===============================
# Get valid combinations to process
valid_combinations <- get_valid_combinations(combinations, df, min_genes = 5)

# Run enrichment analysis
ranges = 1:nrow(valid_combinations)
run_enrichment_analysis(valid_combinations, df, ont_list, ont_name_list, ranges = ranges)


# ===============================
# collect
# ===============================
res = list()
files = list.files("data/smr_enrich", full.names = TRUE)
files = files[grepl("\\.csv", files)]
for (file in files) {
    sub = read.csv(file)
    res[[file]] = sub
}

res_df = do.call(rbind, res)
row.names(res_df) = NULL
write.csv(res_df, "result/smr_enrich.csv", row.names = FALSE)
save(res_df, file = "result/smr_enrich.rdata")
nrow(res_df)


head(res_df)
res_df%>%filter(expo == 'Lung' & grepl ('COVID', trait))

# # ===============================
# # test
# # ===============================

# df = read.csv("result/smr_collect.csv.gz")
# df = df%>%dplyr::select(expo_source, expo, trait, Gene, ENSEMBL)

# sub = df%>%filter(trait == 'other_COVID19_HGI_A2_EUR' & expo=='Lung')
# ENSEMBLs = unique(sub%>%pull(ENSEMBL))

# sub1 <- enrichGO(gene = ENSEMBLs, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1)
# sub1 = as.data.frame(sub1)

# sub1 <- enrichGO(
#   gene = ENSEMBLs,
#   OrgDb = org.Hs.eg.db,
#   keyType = "ENSEMBL",
#   ont = "ALL",
#   pvalueCutoff = 1,
#   qvalueCutoff = 1
# )
# sub1 <- as.data.frame(sub1) %>%
#   dplyr::filter(pvalue < 0.05)
