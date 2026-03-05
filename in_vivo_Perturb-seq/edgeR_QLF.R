library(edgeR)
library(dplyr)

IN_DIR  <- "/Users/linrisheng/Library/CloudStorage/OneDrive-个人/invivo/10x/pRS050-perturb-seq-3targets/Fam50a/for_github"
OUT_DIR <- IN_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

counts <- read.csv(file.path(IN_DIR, "pseudobulk_counts.csv"),
                   row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(IN_DIR, "pseudobulk_metadata.csv"),
                 check.names = FALSE)

colnames(meta) <- trimws(colnames(meta))

stopifnot(nrow(counts) == nrow(meta))

if ("pb_group_id" %in% colnames(meta)) {
  rownames(meta) <- meta$pb_group_id
  counts <- counts[rownames(meta), , drop = FALSE]
}

meta$celltype_coarse <- factor(meta$celltype_coarse)
meta$replicate <- factor(meta$batch)                 # ✅ 用 batch 当 replicate
meta$target <- factor(meta$target, levels = c("NTC", "Fam50a"))

all_results <- list()

for (ct in levels(meta$celltype_coarse)) {
  
  message("\n==============================")
  message("Celltype: ", ct)
  message("==============================")
  
  keep <- meta$celltype_coarse == ct
  meta_ct   <- meta[keep, , drop = FALSE]
  counts_ct <- counts[keep, , drop = FALSE]
  
  if (length(unique(meta_ct$target)) < 2) {
    message("  Skip (not enough conditions)")
    next
  }
  
  y <- DGEList(counts = t(counts_ct), samples = meta_ct)
  
  keep_genes <- filterByExpr(y, group = meta_ct$target)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ replicate + target, data = meta_ct)
  
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  
  coef_name <- "targetFam50a"
  if (!coef_name %in% colnames(design)) {
    message("  Skip (coef not found: ", coef_name, ")")
    next
  }
  
  qlf <- glmQLFTest(fit, coef = coef_name)
  tab <- topTags(qlf, n = Inf)$table
  
  tab$gene <- rownames(tab)
  tab$celltype_coarse <- ct
  tab$contrast <- "Fam50a_vs_NTC"
  
  all_results[[ct]] <- tab
}

final_res <- do.call(rbind, all_results)

final_res <- final_res[, c("gene","celltype_coarse","contrast",
                           "logFC","logCPM","F","PValue","FDR")]

out_all <- file.path(OUT_DIR, "edgeR_QLF_by_celltype_blockReplicate.csv")  # ✅
write.csv(final_res, out_all, row.names = FALSE)

message("\n✅ Saved: ", out_all)