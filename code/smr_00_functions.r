# ===============================
# SMR Enrichment Analysis Functions
# ===============================

#' Prepare GO term-to-gene mappings for enrichment analysis
#' @param org_db Organism database (default: org.Hs.eg.db)
#' @param keytype Key type for gene identifiers (default: "ENSEMBL")
#' @return List containing ont_list and ont_name_list for enrichment analysis
prepare_go_mappings <- function(org_db = org.Hs.eg.db, keytype = "ENSEMBL") {
    print("Preparing GO term-to-gene mappings...")
    
    # Get annotation data
    ann <- AnnotationDbi::select(
        org_db,
        keys     = keys(org_db, keytype = keytype),
        keytype  = keytype,
        columns  = c("GOALL", "ONTOLOGYALL")
    )
    
    # Keep valid rows and de-duplicate TERM<->GENE pairs
    ann <- ann[!is.na(ann$GOALL) & !is.na(ann$ONTOLOGYALL), c("ENSEMBL", "GOALL", "ONTOLOGYALL")]
    ann <- unique(ann)
    
    # Add GO term names (may be NA for some obsolete terms)
    ann$TERM <- Term(ann$GOALL)
    
    # Helper function to split by ontology and prepare TERM2GENE (and TERM2NAME) tables
    to_dt <- function(df, onto) {
        dt <- as.data.table(df[df$ONTOLOGYALL == onto, c("GOALL","ENSEMBL","TERM")])
        # TERM2GENE
        T2G <- unique(dt[, .(term = GOALL, gene = ENSEMBL)])
        # TERM2NAME (optional but useful for nicer output)
        T2N <- unique(dt[, .(term = GOALL, name = TERM)])
        list(T2G = T2G, T2N = T2N)
    }
    
    # Process each ontology
    bp <- to_dt(ann, "BP")
    cc <- to_dt(ann, "CC")
    mf <- to_dt(ann, "MF")
    
    # Create final lists
    ont_list      <- list(BP = bp$T2G, CC = cc$T2G, MF = mf$T2G)     # for TERM2GENE
    ont_name_list <- list(BP = bp$T2N, CC = cc$T2N, MF = mf$T2N)
    
    print("GO mappings prepared successfully!")
    print(paste("BP terms:", nrow(bp$T2G)))
    print(paste("CC terms:", nrow(cc$T2G)))
    print(paste("MF terms:", nrow(mf$T2G)))
    
    return(list(ont_list = ont_list, ont_name_list = ont_name_list))
}

#' Get valid combinations that need to be processed
#' @param combinations Data frame with expo_source, expo, trait columns
#' @param df Main data frame with ENSEMBL genes
#' @param output_dir Directory where enrichment results are saved
#' @param min_genes Minimum number of genes required (default: 5)
#' @return Data frame of valid combinations to process
get_valid_combinations <- function(combinations, df, output_dir = "data/smr_enrich", min_genes = 5) {
    print("Checking which combinations have already been processed...")
    
    # Create directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Get list of already processed files - optimized
    existing_files <- list.files(output_dir, pattern = "\\.csv$", full.names = FALSE)
    
    # Create processed combinations data frame efficiently
    if (length(existing_files) > 0) {
        # Remove .csv extension and split efficiently
        file_names <- gsub("\\.csv$", "", existing_files)
        
        # Use strsplit with sapply for better performance
        split_parts <- strsplit(file_names, "__")
        
        # Filter valid splits and create data frame efficiently
        valid_splits <- split_parts[sapply(split_parts, length) >= 3]
        
        if (length(valid_splits) > 0) {
            processed_combinations <- data.frame(
                expo_source = sapply(valid_splits, function(x) x[1]),
                expo = sapply(valid_splits, function(x) x[2]),
                trait = sapply(valid_splits, function(x) paste(x[3:length(x)], collapse = "__")),
                stringsAsFactors = FALSE
            )
        } else {
            processed_combinations <- data.frame(
                expo_source = character(),
                expo = character(),
                trait = character(),
                stringsAsFactors = FALSE
            )
        }
    } else {
        processed_combinations <- data.frame(
            expo_source = character(),
            expo = character(),
            trait = character(),
            stringsAsFactors = FALSE
        )
    }
    
    # Find combinations that haven't been processed yet
    if (nrow(processed_combinations) > 0) {
        # Merge to find which combinations are missing
        missing_combinations <- combinations %>%
            anti_join(processed_combinations, by = c("expo_source", "expo", "trait"))
    } else {
        missing_combinations <- combinations
    }
    
    print(paste("Total combinations:", nrow(combinations)))
    print(paste("Already processed:", nrow(processed_combinations)))
    print(paste("Remaining to process:", nrow(missing_combinations)))
    
    # Optimized gene counting - pre-calculate gene counts for all combinations
    if (nrow(missing_combinations) > 0) {
        # Create a lookup table for gene counts
        gene_counts <- df %>%
            group_by(expo_source, expo, trait) %>%
            summarise(gene_count = n_distinct(ENSEMBL), .groups = 'drop')
        
        # Join with missing combinations to get gene counts
        valid_combinations <- missing_combinations %>%
            left_join(gene_counts, by = c("expo_source", "expo", "trait")) %>%
            filter(gene_count >= min_genes)
    } else {
        valid_combinations <- missing_combinations
    }
    
    print(paste("Valid combinations (>=", min_genes, " genes):", nrow(valid_combinations)))
    
    return(valid_combinations)
}

#' Run GO enrichment analysis for a list of combinations
#' @param valid_combinations Data frame of combinations to process
#' @param df Main data frame with ENSEMBL genes
#' @param ont_list List of ontology term-to-gene mappings
#' @param ont_name_list List of ontology term-to-name mappings
#' @param output_dir Directory to save results (default: "data/smr_enrich")
#' @param progress_interval Show progress every N iterations (default: 10)
#' @param ranges Vector of row indices to process (default: 1:nrow(valid_combinations))
run_enrichment_analysis <- function(valid_combinations, df, ont_list, ont_name_list, 
                                   output_dir = "data/smr_enrich", progress_interval = 10,
                                   ranges = 1:nrow(valid_combinations)) {
    
    if (nrow(valid_combinations) == 0) {
        print("No valid combinations to process.")
        return(NULL)
    }
    
    # Validate ranges
    if (max(ranges) > nrow(valid_combinations)) {
        stop("Ranges exceed the number of valid combinations")
    }
    
    # Filter combinations based on ranges
    combinations_to_process <- valid_combinations[ranges, ]
    
    if (nrow(combinations_to_process) == 0) {
        print("No combinations to process with the specified ranges.")
        return(NULL)
    }
    
    # Create directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Define output columns
    cols <- c("expo_source","expo","trait","ID","Description","GeneRatio","BgRatio",
              "pvalue","p.adjust","qvalue","geneID","Count","ONTOLOGY")
    empty_df <- setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)
    
    print(paste("Starting enrichment analysis for", nrow(combinations_to_process), "combinations..."))
    print(paste("Processing rows:", min(ranges), "to", max(ranges)))
    
    for (i in 1:nrow(combinations_to_process)) {
        if (i %% progress_interval == 0) {
            print(paste0("Processing: ", i, "/", nrow(combinations_to_process), " (", round(i/nrow(combinations_to_process)*100, 1), "%)"))
        }
        
        combo <- combinations_to_process[i, ]
        ENSEMBLs = df%>%filter(expo_source == combo$expo_source, expo == combo$expo, trait == combo$trait)%>%pull(ENSEMBL)%>%unique()
        file_out = paste0(output_dir, "/", combo$expo_source, "__", combo$expo, "__", combo$trait, ".csv")
        
        # Perform GO enrichment analysis directly with ENSEMBL
        res_list <- imap(ont_list, ~{
            res <- tryCatch(enricher(ENSEMBLs, pvalueCutoff = 1, qvalueCutoff = 1, TERM2GENE = .x, TERM2NAME = ont_name_list[[.y]]), error = function(e) NULL)
            if (is.null(res)) return(NULL)
            df <- as.data.frame(res)
            if (!nrow(df)) return(NULL)
            mutate(df, ONTOLOGY = .y)
        })

        sub <- bind_rows(compact(res_list))
        
        # Check if sub has the required columns and filter if it does
        if (nrow(sub) > 0 && "pvalue" %in% colnames(sub)) {
            sub <- sub %>% filter(pvalue < 0.05)
        }
        
        # original method, slow
        # sub1 <- enrichGO(gene = ENSEMBLs, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "ALL")
        # sub1 = as.data.frame(sub1)

        if (nrow(sub) == 0) {sub = empty_df}
        sub = sub%>%mutate(expo_source = combo$expo_source, expo = combo$expo, trait = combo$trait)%>%
            dplyr::select(expo_source, expo, trait, everything())
        write.csv(sub, file_out, row.names = FALSE)
    }
    
    print("Enrichment analysis completed!")
}
