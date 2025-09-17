# ================================================================
# pair
# ================================================================
setwd("~/easymr")
source("~/easymr/code/function.r")

df <- read.csv('result/info.csv')

df <- df %>%
  filter(h2_p < 0.05)

# make pair
choose(nrow(df), 2) # 2795430
pair <- data.frame()
for (i in 1:nrow(df)) {
  trait1 <- df[i, "id"]
  trait2s <- df %>%
    filter(source != df[i, "source"] & umls != df[i, "umls"]) %>%
    pull(id)
  sub <- cbind(trait1, trait2s)
  pair <- rbind(pair, sub)
}

# remove dup
pair$key <- apply(pair, 1, function(x) paste(sort(c(x[1], x[2])), collapse = "_AND_"))
pair <- pair %>% distinct(key, .keep_all = T)

names(pair) <- c("trait1", "trait2", "key")
dim(pair) # 1123830

write.csv(pair, "result/ldsc_pair.csv", quote = F, row.names = F)

# ================================================================
# pair to run
# ================================================================
# nohup find /home/zhang/easymr/data/ldsc -type f -printf "%p\t%s\n" > /home/zhang/easymr/data/ldsc_files_with_size.txt &

pair = read.csv('result/ldsc_pair.csv')

file_list <- "data/ldsc_files_with_size.txt"
existing_files <- read.delim(file_list, header = F, stringsAsFactors = F)
names(existing_files) <- c("path", "size")

# Extract keys from existing files
existing_keys <- gsub(".*/ldsc/(.*)\\.log$", "\\1", existing_files$path)

# Check for empty files (size = 0)
empty_files <- existing_files$path[existing_files$size == 0]
empty_keys <- gsub(".*/ldsc/(.*)\\.log$", "\\1", empty_files)

# Find missing pairs and pairs with empty files
missing_pairs <- pair[!pair$key %in% existing_keys | pair$key %in% empty_keys, ]
cat(sprintf("Number of pairs to run (missing + empty): %d\n", nrow(missing_pairs)))

# Save missing pairs to run
write.csv(missing_pairs, "result/ldsc_pair_torun.csv", quote = F, row.names = F)

# Find files to delete (files that don't match any valid pair)
files_to_delete <- existing_files$path[!existing_keys %in% pair$key]
files_to_delete <- unique(c(files_to_delete, empty_files))
cat(sprintf("Number of files to delete (including empty files): %d\n", length(files_to_delete)))

# Delete the files
if (length(files_to_delete) > 0) {
  file.remove(files_to_delete)
  cat(sprintf("Deleted %d files (including empty files)\n", length(files_to_delete)))
}
