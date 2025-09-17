## prepare
source("./code/mr_server/prepare.r")
load(file_ao)
load(file_pairs)
links <- fread("./data/link.csv")
str(links)


# ================================================================
# filter by significance, add trait name, combine p
# ================================================================
# filter links by ivw_p and b sign
links <- links %>%
  filter(!is.na(nsnp)) %>%
  filter(ivw_p < 0.05)
links <- links %>% filter(sign(ivw_b) == sign(egger_b) & sign(ivw_b) == sign(wmedian_b) & sign(ivw_b) == sign(wmode_b))
links <- links %>% filter(egger_int_p > 0.05 & egger_q_p > 0.05 & ivw_q_p > 0.05)

# add trait name
ao <- ao %>%
  mutate(db = toupper(db), trait = paste0(db, ": ", trait))
map_trait <- ao %>% select(id, trait, db, umls)
links <- links %>%
  merge(map_trait, by.x = "exposure", by.y = "id") %>%
  merge(map_trait, by.x = "outcome", by.y = "id") %>%
  rename(exposure_id = exposure, outcome_id = outcome, exposure = trait.x, outcome = trait.y, exposure_db = db.x, 
    outcome_db = db.y, exposure_umls=umls.x, outcome_umls=umls.y) %>%
  relocate(exposure, outcome, exposure_db, outcome_db)

# filter with pairs
links = links%>%merge(pairs%>%rename(exposure_id=exposure, outcome_id=outcome), by=c('exposure_id', 'outcome_id'))

# strength of edge
links$mr_combine_p <- sapply(1:nrow(links), function(x) {
  pvals <- unlist(links[x, c("ivw_p", "egger_p", "wmedian_p", "wmode_p")])
  combine_p <- combinePvalVector(pvals)[2]
  return(combine_p)
})
links$weight <- (1 - links$mr_combine_p) * pmin(links$nsnp / 10, 1) # nsnp=5, weight=0.5, nsnp>=10, weight=1

# save
write.csv(links, file = "./data/db_mr.csv", row.names = F)

temp = links%>%select(exposure, outcome, weight, exposure_umls, outcome_umls)
write.table(temp, file = "./data/temp.txt", row.names = F, sep='\t', quote=F)


# ================================================================
# temp
# ================================================================
source("./code/mr_server/prepare.r")
db_mr = fread(file_mrdb)

## default value
input_default = list()
input_default$pthres_ivw = 0.05; input_default$pthres_egger = 0.05; input_default$pthres_wmedian = 0.05; input_default$pthres_wmode = 0.05
input_default$thres_nsnp = 10; input_default$db = dbs

# get links
links = select_links(input_default, db_mr)
links = links %>%
  select(exposure, outcome, weight) %>%
  rename(from = exposure, to = outcome)

# make graph 
graph <- graph.data.frame(d = links, directed = TRUE, vertices = NULL)
E(graph)$weight <- links$weight
length(V(graph))
length(E(graph))

# cluster
communities <- cluster_edge_betweenness(graph, directed = TRUE)
modularity(graph, membership = communities$membership) # modularity is based on the whole graph


tab = table(communities$membership)