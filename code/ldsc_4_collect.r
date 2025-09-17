# ================================================================
# collect result
# ================================================================
setwd("~/easymr")
source("~/easymr/code/function.r")
library(parallel)

# Create temp directory if it doesn't exist
dir.create("data/temp", showWarnings = FALSE)

# Delete existing chunk files
chunk_files <- list.files("data/temp", pattern = "ldsc_collect_chunk_\\d+\\.csv", full.names = TRUE)
if (length(chunk_files) > 0) {
  message(sprintf("Deleting %d existing chunk files...", length(chunk_files)))
  file.remove(chunk_files)
}

pair <- read.csv("result/ldsc_pair.csv")

# Pre-allocate vectors for faster access
trait1_vec <- pair$trait1
trait2_vec <- pair$trait2
key_vec    <- pair$key

# Function to process a single file
process_file <- function(i) {
  if (i %% 1000 == 0) {
    message(sprintf("Processing file %d/%d", i, length(key_vec)))
  }
  
  file_in <- file.path("data/ldsc", paste0(key_vec[i], ".log"))
  tryCatch({
    sub <- c(trait1_vec[i], trait2_vec[i], collect_ldsc(file_in))
    pval <- as.numeric(sub[["rg_p"]])
    rg   <- as.numeric(sub[["rg"]])
    
    if (!is.na(pval) && !is.na(rg) && pval < 0.05 && abs(rg) < 1) {
      return(sub)
    }
    return(NULL)
  }, error = function(e) {
    message(sprintf("Error processing file %s: %s", file_in, e$message))
    return(NULL)
  })
}

# Process in chunks
n_cores <- 20
chunk_size <- 50000  # Adjust this based on your memory constraints
n_chunks <- ceiling(length(key_vec) / chunk_size)

message(sprintf("Processing %d files in %d chunks of size %d", length(key_vec), n_chunks, chunk_size))

# Process each chunk and save intermediate results
for (chunk in 1:n_chunks) {
  chunk_file <- sprintf("data/temp/ldsc_collect_chunk_%d.csv", chunk)
  if (file.exists(chunk_file)) {
    message(sprintf("Skipping chunk %d because it already exists", chunk))
    next
  }
  start_idx <- (chunk - 1) * chunk_size + 1
  end_idx <- min(chunk * chunk_size, length(key_vec))
  
  message(sprintf("Processing chunk %d/%d (files %d-%d)", chunk, n_chunks, start_idx, end_idx))
  
  # Process current chunk
  chunk_indices <- start_idx:end_idx
  res_list <- mclapply(chunk_indices, process_file, mc.cores = n_cores)
  
  # Remove NULL elements and convert to data frame
  res_list <- res_list[!sapply(res_list, is.null)]
  res <- do.call(rbind, res_list)
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  names(res)[1:2] <- c("trait1", "trait2")
  
  # Save chunk results
  chunk_file <- sprintf("data/temp/ldsc_collect_chunk_%d.csv", chunk)
  write.csv(res, chunk_file, row.names = FALSE)
  rm(res_list, res)  # clean memory
  gc() # clean memory
}

# Merge all chunks
message("Merging all chunks...")
chunk_files <- list.files("data/temp", pattern = "ldsc_collect_chunk_\\d+\\.csv", full.names = TRUE)
res <- do.call(rbind, lapply(chunk_files, read.csv))
res <- distinct(res)  # Remove duplicate rows

# Add trait information
info <- read.csv("~/easymr/result/info.csv") %>% select(id, trait, n, case_prop, umls)

# Merge trait information
res1 <- res %>% merge(info %>% rename(trait1 = id, trait1_name = trait, trait1_n = n, trait1_case_prop = case_prop, trait1_umls = umls), by = "trait1")
res1 <- res1 %>% merge(info %>% rename(trait2 = id, trait2_name = trait, trait2_n = n, trait2_case_prop = case_prop, trait2_umls = umls), by = "trait2")
res1$file <- NULL

# Save final results
write.csv(res1, "result/ldsc_collect.csv", row.names = F)
message("Final results saved to result/ldsc_collect.csv")

# # ================================================================
# # Extract Swedish traits
# # ================================================================
# res1 <- read.csv("result/ldsc_collect.csv")

# # Filter for Swedish traits and save
# res2 <- res1 %>% filter(grepl("swe", trait1, ignore.case = T) | grepl("swe", trait2, ignore.case = T))
# write.csv(res2, "result/ldsc_collect_swe.csv", row.names = F)

# gzip
# nohup tar -cvzf ldsc.tar.gz ldsc/ &
# tar -tzf ldsc.tar.gz | tail -n 5   # check if done