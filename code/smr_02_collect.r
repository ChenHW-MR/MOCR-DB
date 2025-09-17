# ================================================================
# Collect significant SMR results in chunks to avoid memory issues
# ================================================================
setwd("~/easymr")
source("~/easymr/code/function.r")
library(parallel)

# Define paths for different QTL summary data
path_eqtl <- "/home/zhang/db/smr/GTEx_V8_cis_eqtl_summary/"
path_mqtl <- "/home/zhang/db/smr/EUR/"
path_sqtl <- "/home/zhang/db/smr/BrainMeta_cis_sqtl_summary/"
path_caqtl <- "/home/zhang/db/smr/Bryois_caQTL_summary/"
expo_sources <- c("eqtl", "mqtl", "sqtl", "caqtl")
sources <- c("ukb", "fin", 'swe', 'other')

# Read trait information
map = read.csv('result/info.csv')

# Create temp directory if it doesn't exist
if (!dir.exists("data/temp")) dir.create("data/temp", showWarnings = FALSE)

# 1. Generate all tasks (as a list, not a data.frame)
tasks <- list()
task_idx <- 0
for (expo_source in expo_sources) {
  expos <- gsub(".besd", "", get("besd", list.files(base::get(sprintf("path_%s", expo_source)))))
  message(sprintf("Processing expo_source: %s (%d expos)", expo_source, length(expos)))
  for (trait_idx in seq_along(map$id)) {
    trait <- map$id[trait_idx]
    for (expo in expos) {
      file_in <- sprintf("~/easymr/data/smr/%s/%s_TO_%s.smr.gz", expo_source, expo, trait)
      if (file.exists(file_in)) {
        tasks[[length(tasks) + 1]] <- list(expo_source=expo_source, expo=expo, trait=trait, file_in=file_in)
        task_idx <- task_idx + 1
        if (task_idx %% 1000 == 0) {
          message(sprintf("Generated %d tasks so far...", task_idx))
        }
      }
    }
    # Print progress every 100 traits
    if (trait_idx %% 100 == 0) {
      message(sprintf("expo_source %s: trait %d/%d", expo_source, trait_idx, length(map$id)))
    }
  }
}

chunk_size <- 5000  # Number of tasks per chunk
n_chunks <- ceiling(length(tasks) / chunk_size)
n_cores <- 20  # Number of CPU cores to use for parallel processing

# 2. Define the function to process a single file
def_process_file <- function(task) {
  sub <- tryCatch({
    # Read the SMR result file
    df <- read.delim(task$file_in)
    # Filter and annotate the results
    df <- df %>% filter(if_all(.cols = -Gene, ~ !is.na(.))) %>% filter(nsnp_HEIDI>=10 & p_HEIDI>=0.05) %>%
      mutate(expo_source=task$expo_source, expo=task$expo, trait=task$trait) %>%
      filter(p_SMR<1e-3)
    # Return the filtered data frame if not empty
    if (nrow(df) > 0) return(df)
    else return(NULL)
  }, error=function(e) NULL)
  return(sub)
}

# 3. Process tasks in chunks and save each chunk's result to a temporary file
for (chunk in 1:n_chunks) {
  chunk_file <- sprintf("data/temp/smr_collect_chunk_%d.csv", chunk)
  if (file.exists(chunk_file)) {
    message(sprintf("Skipping chunk %d because it already exists", chunk))
    next
  }
  start_idx <- (chunk - 1) * chunk_size + 1
  end_idx <- min(chunk * chunk_size, length(tasks))
  chunk_tasks <- tasks[start_idx:end_idx]
  message(sprintf("Processing chunk %d/%d (tasks %d-%d)", chunk, n_chunks, start_idx, end_idx))
  res_list <- mclapply(chunk_tasks, def_process_file, mc.cores = n_cores)
  res_list <- res_list[!sapply(res_list, is.null)]
  if (length(res_list) > 0) {
    res <- do.call(rbind, res_list)
    write.csv(res, chunk_file, row.names = FALSE)
    rm(res)
  }
  rm(res_list, chunk_tasks); gc()
}

# 4. Merge all chunk results and annotate
message("Merging all chunks...")
chunk_files <- list.files("data/temp", pattern = "smr_collect_chunk_\\d+\\.csv", full.names = TRUE)
res <- do.call(rbind, lapply(chunk_files, read.csv))

# ================================================================
# Map methylation to nearest gene 
# ================================================================
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# split
sub1 = res%>%filter(expo_source%in%c('mqtl', 'caqtl'))
sub2 = res%>%filter(!expo_source%in%c('mqtl', 'caqtl'))

# make granges
chr <- paste0('chr', sub1$ProbeChr)
pos <- sub1$Probe_bp
probes_gr <- GRanges(seqnames = chr, ranges   = IRanges(start = pos, width = 1))
genes_gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
nn <- distanceToNearest(probes_gr, genes_gr, ignore.strand = TRUE)

# initialize outputs
nearest_entrez <- rep(NA_character_, nrow(sub1))
nearest_dist   <- rep(NA_integer_,  nrow(sub1))

# fill back for kept rows
qh <- queryHits(nn); sh <- subjectHits(nn)
# Entrez IDs are stored in mcols(genes_gr)$gene_id (as integer/character)
entrez_vec <- as.character(mcols(genes_gr)$gene_id)
nearest_entrez[qh] <- entrez_vec[sh]
nearest_dist[qh]   <- mcols(nn)$distance

# map Entrez -> HGNC symbol
nearest_symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys     = unique(na.omit(nearest_entrez)),
  keytype  = "ENTREZID",
  column   = "SYMBOL",
  multiVals = "first"
)
nearest_symbol_full <- unname(nearest_symbol[nearest_entrez])

# bind to original, note there are na genes
sub1 = sub1 %>% mutate(Gene = nearest_symbol_full)

# merge back, add ensembl id
res1 = rbind(sub1, sub2)%>%filter(!is.na(Gene))
sym2ens <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                 keys = unique(res1$Gene),
                                 keytype = "SYMBOL",
                                 column = "ENSEMBL",
                                 multiVals = "first")
map_tbl = tibble(
  Gene = names(sym2ens),
  ENSEMBL  = unname(sym2ens)
)
res1 = res1%>%left_join(map_tbl, by = "Gene")%>%filter(!is.na(ENSEMBL))

# add trait info
res1 = res1%>%merge(map %>% dplyr::rename(trait = id, trait_name = trait, trait_n = n, trait_case_prop = case_prop, trait_umls = umls), by = "trait")

write.csv(res1, "result/smr_collect.csv", row.names = F)

# check
table(res1$expo_source)
#  caqtl    eqtl    mqtl    sqtl 
#    8888 1154739  985512  254557 