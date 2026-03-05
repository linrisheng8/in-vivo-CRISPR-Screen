# ==========================================================
# Batch ORA Enrichment (GO BP/CC/MF + KEGG + Hallmark)
# for ALL / UP / DOWN DEG sets
# grouped by celltype_coarse
#
# Input DEG columns required:
# gene, celltype_coarse, logFC, FDR
#
# Your DEG file columns:
# gene celltype_coarse contrast logFC logCPM F PValue FDR
#
# Species: Mouse (org.Mm.eg.db)
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(msigdbr)
})

# ------------------------------
# 0) Paths
# ------------------------------
setwd("/Users/linrisheng/Library/CloudStorage/OneDrive-个人/invivo/10x/pRS050-perturb-seq-3targets/Fam50a/20260125/Enrichment")

infile <- "/Users/linrisheng/Library/CloudStorage/OneDrive-个人/invivo/10x/pRS050-perturb-seq-3targets/Fam50a/20260125/edgeR_QLF_by_celltype_blockReplicate.csv"

out_dir <- "."
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------
# 1) Thresholds
# ------------------------------
# DEG thresholds
fdr_cutoff   <- 0.05
logfc_cutoff <- 0.25

# enrichment cutoffs
p_cutoff <- 0.05
q_cutoff <- 0.05

# minimal gene number for enrichment
min_gene_n <- 10

# ------------------------------
# 2) Hallmark gene sets (mouse)
#    Robust for different msigdbr versions
# ------------------------------
msig_raw <- msigdbr(species = "Mus musculus", collection = "H")

# TERM2GENE must be: data.frame(term, gene)
# Here: (gs_name, ENTREZID)
if ("entrez_gene" %in% colnames(msig_raw)) {
  
  msig_h <- msig_raw %>%
    dplyr::select(gs_name, entrez_gene) %>%
    dplyr::rename(ENTREZID = entrez_gene) %>%
    distinct(gs_name, ENTREZID)
  
} else if ("entrez_gene_id" %in% colnames(msig_raw)) {
  
  msig_h <- msig_raw %>%
    dplyr::select(gs_name, entrez_gene_id) %>%
    dplyr::rename(ENTREZID = entrez_gene_id) %>%
    distinct(gs_name, ENTREZID)
  
} else if ("gene_symbol" %in% colnames(msig_raw)) {
  
  # only SYMBOL -> map to ENTREZID
  tmp <- msig_raw %>%
    dplyr::select(gs_name, gene_symbol) %>%
    distinct()
  
  map_h <- bitr(tmp$gene_symbol,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Mm.eg.db) %>%
    distinct()
  
  msig_h <- tmp %>%
    left_join(map_h, by = c("gene_symbol" = "SYMBOL")) %>%
    dplyr::select(gs_name, ENTREZID) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct()
  
} else {
  stop("❌ msigdbr output missing entrez columns and gene_symbol.\nColumns:\n",
       paste(colnames(msig_raw), collapse = ", "))
}

# ------------------------------
# 3) Read DEG
# ------------------------------
deg <- read.csv(infile, stringsAsFactors = FALSE)
deg$gene <- trimws(deg$gene)

required_cols <- c("gene", "celltype_coarse", "logFC", "FDR")
if (!all(required_cols %in% colnames(deg))) {
  stop("❌ Your CSV must contain columns: gene, celltype_coarse, logFC, FDR\n",
       "Current columns are:\n",
       paste(colnames(deg), collapse = ", "))
}

# ------------------------------
# 4) Helper: robust mapping to ENTREZ (SYMBOL + ENSEMBL)
# ------------------------------
map_gene_to_entrez <- function(gene_vec) {
  
  gene_vec <- unique(trimws(gene_vec))
  gene_vec <- gene_vec[!is.na(gene_vec) & gene_vec != ""]
  
  # treat ENSMUSG... as ENSEMBL id
  genes_symbol  <- gene_vec[!grepl("^ENSMUSG", gene_vec)]
  genes_ensembl <- gene_vec[ grepl("^ENSMUSG", gene_vec)]
  
  # strip version: ENSMUSG00000000001.1 -> ENSMUSG00000000001
  genes_ensembl <- sub("\\..*$", "", genes_ensembl)
  
  map_symbol  <- NULL
  map_ensembl <- NULL
  
  if (length(genes_symbol) > 0) {
    map_symbol <- tryCatch({
      bitr(genes_symbol,
           fromType = "SYMBOL",
           toType   = "ENTREZID",
           OrgDb    = org.Mm.eg.db) %>%
        distinct(SYMBOL, .keep_all = TRUE) %>%
        transmute(gene = SYMBOL, ENTREZID = ENTREZID)
    }, error = function(e) {
      message("⚠️ SYMBOL mapping failed: ", e$message)
      NULL
    })
  }
  
  if (length(genes_ensembl) > 0) {
    map_ensembl <- tryCatch({
      bitr(genes_ensembl,
           fromType = "ENSEMBL",
           toType   = "ENTREZID",
           OrgDb    = org.Mm.eg.db) %>%
        distinct(ENSEMBL, .keep_all = TRUE) %>%
        transmute(gene = ENSEMBL, ENTREZID = ENTREZID)
    }, error = function(e) {
      message("⚠️ ENSEMBL mapping failed: ", e$message)
      NULL
    })
  }
  
  bind_rows(map_symbol, map_ensembl) %>%
    distinct(gene, .keep_all = TRUE)
}

# ------------------------------
# ✅ NEW Helper: add SetSize + Overlap
# ------------------------------
add_overlap_cols <- function(df) {
  # Count: k (# hits)
  # BgRatio: "M/N" where M = term size in background
  if (!all(c("Count", "BgRatio") %in% colnames(df))) return(df)
  
  bg <- strsplit(df$BgRatio, "/", fixed = TRUE)
  
  # ✅ SetSize = M (numerator of BgRatio)
  bg_num <- suppressWarnings(as.integer(vapply(bg, function(x) x[1], "")))
  
  df$SetSize <- bg_num
  df$Overlap <- paste0(df$Count, "/", df$SetSize)
  
  df
}

# ------------------------------
# 5) Helper: run ORA (GO/KEGG/Hallmark)
# ------------------------------
run_ora_all <- function(sig_entrez, group_dir, prefix_tag) {
  
  # ---------------- GO: BP/CC/MF ----------------
  for (ont_type in c("BP", "CC", "MF")) {
    
    ego <- tryCatch({
      enrichGO(
        gene          = sig_entrez,
        OrgDb         = org.Mm.eg.db,
        keyType       = "ENTREZID",
        ont           = ont_type,
        pAdjustMethod = "BH",
        pvalueCutoff  = p_cutoff,
        qvalueCutoff  = q_cutoff,
        readable      = TRUE
      )
    }, error = function(e) {
      message("❌ ", prefix_tag, " GO-", ont_type, " failed: ", e$message)
      NULL
    })
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      
      df <- as.data.frame(ego)
      df <- add_overlap_cols(df)  # ✅ add SetSize + Overlap
      
      write.csv(df,
                file.path(group_dir, paste0(prefix_tag, "_GO_", ont_type, ".csv")),
                row.names = FALSE)
      
      pdf(file.path(group_dir, paste0(prefix_tag, "_GO_", ont_type, "_dotplot.pdf")),
          width = 9, height = 6)
      print(dotplot(ego, showCategory = 20))
      dev.off()
      
    } else {
      message("⚠️ ", prefix_tag, " GO-", ont_type, " empty.")
    }
  }
  
  # ---------------- KEGG ----------------
  ekk <- tryCatch({
    enrichKEGG(
      gene         = sig_entrez,
      organism     = "mmu",
      pvalueCutoff = p_cutoff
    )
  }, error = function(e) {
    message("❌ ", prefix_tag, " KEGG failed: ", e$message)
    NULL
  })
  
  if (!is.null(ekk) && nrow(as.data.frame(ekk)) > 0) {
    
    ekk_readable <- setReadable(ekk, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    df <- as.data.frame(ekk_readable)
    df <- add_overlap_cols(df)  # ✅ add SetSize + Overlap
    
    write.csv(df,
              file.path(group_dir, paste0(prefix_tag, "_KEGG.csv")),
              row.names = FALSE)
    
    pdf(file.path(group_dir, paste0(prefix_tag, "_KEGG_dotplot.pdf")),
        width = 9, height = 6)
    print(dotplot(ekk_readable, showCategory = 20))
    dev.off()
    
  } else {
    message("⚠️ ", prefix_tag, " KEGG empty.")
  }
  
  # ---------------- Hallmark ----------------
  eh <- tryCatch({
    enricher(
      gene          = sig_entrez,
      TERM2GENE     = msig_h,
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff
    )
  }, error = function(e) {
    message("❌ ", prefix_tag, " Hallmark failed: ", e$message)
    NULL
  })
  
  if (!is.null(eh) && nrow(as.data.frame(eh)) > 0) {
    
    eh_readable <- setReadable(eh, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    df <- as.data.frame(eh_readable)
    df <- add_overlap_cols(df)  # ✅ add SetSize + Overlap
    
    write.csv(df,
              file.path(group_dir, paste0(prefix_tag, "_HALLMARK.csv")),
              row.names = FALSE)
    
    pdf(file.path(group_dir, paste0(prefix_tag, "_HALLMARK_dotplot.pdf")),
        width = 9, height = 6)
    print(dotplot(eh_readable, showCategory = 20))
    dev.off()
    
  } else {
    message("⚠️ ", prefix_tag, " Hallmark empty.")
  }
}

# ------------------------------
# 6) Loop groups by celltype_coarse
# ------------------------------
groups <- deg %>%
  distinct(celltype_coarse) %>%
  arrange(celltype_coarse)

cat("Total celltypes:", nrow(groups), "\n")

for (i in seq_len(nrow(groups))) {
  
  ct <- groups$celltype_coarse[i]
  
  cat("\n==============================\n")
  cat("Celltype:", ct, "\n")
  
  sub_deg <- deg %>% filter(celltype_coarse == ct)
  
  safe_ct <- gsub("[^A-Za-z0-9_\\-]", "_", ct)
  group_dir <- file.path(out_dir, safe_ct)
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------- ALL ----------
  sig_all <- sub_deg %>%
    filter(!is.na(FDR)) %>%
    filter(FDR < fdr_cutoff, abs(logFC) > logfc_cutoff) %>%
    distinct(gene, .keep_all = TRUE)
  
  cat("ALL sig genes:", nrow(sig_all), "\n")
  
  if (nrow(sig_all) >= min_gene_n) {
    map_all <- map_gene_to_entrez(sig_all$gene)
    entrez_all <- unique(map_all$ENTREZID)
    entrez_all <- entrez_all[!is.na(entrez_all)]
    
    if (length(entrez_all) >= min_gene_n) {
      run_ora_all(entrez_all, group_dir, "ALL")
    } else {
      message("⚠️ ALL mapped Entrez too few, skip.")
    }
  } else {
    message("⚠️ ALL too few genes, skip.")
  }
  
  # ---------- UP ----------
  sig_up <- sub_deg %>%
    filter(!is.na(FDR)) %>%
    filter(FDR < fdr_cutoff, logFC > logfc_cutoff) %>%
    distinct(gene, .keep_all = TRUE)
  
  cat("UP sig genes:", nrow(sig_up), "\n")
  
  if (nrow(sig_up) >= min_gene_n) {
    map_up <- map_gene_to_entrez(sig_up$gene)
    entrez_up <- unique(map_up$ENTREZID)
    entrez_up <- entrez_up[!is.na(entrez_up)]
    
    if (length(entrez_up) >= min_gene_n) {
      run_ora_all(entrez_up, group_dir, "UP")
    } else {
      message("⚠️ UP mapped Entrez too few, skip.")
    }
  } else {
    message("⚠️ UP too few genes, skip.")
  }
  
  # ---------- DOWN ----------
  sig_down <- sub_deg %>%
    filter(!is.na(FDR)) %>%
    filter(FDR < fdr_cutoff, logFC < -logfc_cutoff) %>%
    distinct(gene, .keep_all = TRUE)
  
  cat("DOWN sig genes:", nrow(sig_down), "\n")
  
  if (nrow(sig_down) >= min_gene_n) {
    map_down <- map_gene_to_entrez(sig_down$gene)
    entrez_down <- unique(map_down$ENTREZID)
    entrez_down <- entrez_down[!is.na(entrez_down)]
    
    if (length(entrez_down) >= min_gene_n) {
      run_ora_all(entrez_down, group_dir, "DOWN")
    } else {
      message("⚠️ DOWN mapped Entrez too few, skip.")
    }
  } else {
    message("⚠️ DOWN too few genes, skip.")
  }
  
  cat("✅ Finished:", ct, "\n")
}

cat("\nALL DONE. Results saved under:\n", normalizePath(out_dir), "\n")
