#================================================================
# pair
#================================================================
setwd('~/easymr')
source('~/easymr/code/function.r')

df = read.csv('result/info.csv')

# set expo with more than 10 snp
df = df%>%mutate(expo=ifelse(nsnp_clump>=10, 1, 0))

# set by hand, customized expo with less than 10 snp!
df = df%>%mutate(expo=ifelse(source=='swe', 1, expo))

# make pair
pair = data.frame()
for (from in df%>%filter(expo==1)%>%pull(id)) {
  from_row <- df%>%filter(id==from)
  outcomes = df%>%filter(id!=from & source != from_row$source & umls != from_row$umls)%>%select(id)
  sub = cbind(from, outcomes)
  pair = rbind(pair, sub)
}

pair = data.frame(pair)
names(pair) = c('from', 'to')
dim(pair) # 634754

write.csv(pair, 'result/mr_pair.csv', row.names=F)

# ================================================================
# pair to run
# ================================================================
# nohup find /home/zhang/easymr/data/harmo_dat -type f > /home/zhang/easymr/mr_filelist.txt &
pair = read.csv('result/mr_pair.csv')
file_list <- "data/mr_filelist.txt"
existing_files <- readLines(file_list)

# Extract keys from existing files
existing_keys <- gsub(".*/([^/]+)_TO_([^/]+)\\.rdata", "\\1_TO_\\2", existing_files)

# Find missing pairs
missing_pairs <- pair %>%
  mutate(key = paste(from, to, sep="_TO_")) %>%
  filter(!key %in% existing_keys)
cat(sprintf("Number of missing pairs: %d\n", nrow(missing_pairs)))

write.csv(missing_pairs, "result/mr_pair_torun.csv", row.names = FALSE)

# Find files to delete (files that don't match any valid pair)
files_to_delete <- existing_files[!existing_keys %in% paste(pair$from, pair$to, sep="_TO_")]
cat(sprintf("Number of files to delete: %d\n", length(files_to_delete)))

# Delete files that don't match any valid pair
if(length(files_to_delete) > 0) {
    for(file in files_to_delete) {
        if(file.exists(file)) {
            unlink(file)
        }
    }
}







