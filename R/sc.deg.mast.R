#' @title Systematic Cohort-Level Differential Expression Analysis (Full Parameters)
#'
#' @description
#' Performs Differential Expression (DE) analysis across all cell types in a Seurat object,
#' comparing two specific conditions. Fully integrated with Seurat's FindMarkers parameters.
#'
#' @author Hyundong Yoon
#' @export
sc.deg.mast <- function(seurat_obj,
                        ident.1,
                        ident.2,
                        celltype_col = "detailed.celltypes",
                        group_col = "condition",
                        # ── [Seurat FindMarkers Parameters] ──
                        assay = NULL,
                        features = NULL,
                        logfc.threshold = 0.1,
                        test.use = "wilcox",
                        slot = "data",
                        min.pct = 0.01,
                        min.diff.pct = -Inf,
                        only.pos = FALSE,
                        max.cells.per.ident = Inf,
                        random.seed = 1,
                        latent.vars = NULL,
                        min.cells.feature = 3,
                        min.cells.group = 3,
                        mean.fxn = NULL,
                        fc.name = NULL,
                        base = 2,
                        densify = FALSE,
                        # ──────────────────────────────────────
                        gc_interval = 5,
                        ...) { # ... 을 통해 명시되지 않은 파라미터도 통과시킴

  # Required libraries
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tibble)
  })

  # 1. Preparation
  cat("[1/3] Preparing aggregate metadata...\n")
  seurat_obj$tmp_aggregate <- paste(seurat_obj[[celltype_col, drop = TRUE]], 
                                    seurat_obj[[group_col, drop = TRUE]], sep = "_")
  
  all_celltypes <- unique(as.character(seurat_obj[[celltype_col, drop = TRUE]]))
  total_types <- length(all_celltypes)
  valid_groups <- unique(seurat_obj$tmp_aggregate)

  DEG_list <- list()
  start_time <- Sys.time()

  cat(sprintf("[2/3] Starting DE Analysis for %s cell types (Test: %s)\n", total_types, test.use))
  cat(sprintf("Comparison: %s vs %s\n", ident.1, ident.2))

  # 2. Main Loop
  for (i in seq_along(all_celltypes)) {
    celltype <- all_celltypes[i]
    group1_id <- paste0(celltype, "_", ident.1)
    group2_id <- paste0(celltype, "_", ident.2)

    cat(sprintf("\n[%d/%d] (%s) Processing: %s ... ", i, total_types, format(Sys.time(), "%H:%M:%S"), celltype))

    # Check if both comparison groups exist
    if (group1_id %in% valid_groups && group2_id %in% valid_groups) {
      
      sub_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj[[celltype_col, drop = TRUE]] == celltype])
      
      tryCatch({
        # ── 모든 파라미터를 내부 FindMarkers로 전달 ──
        deg <- FindMarkers(
          object = sub_obj,
          ident.1 = group1_id,
          ident.2 = group2_id,
          group.by = "tmp_aggregate",
          assay = assay,
          features = features,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          slot = slot,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          only.pos = only.pos,
          max.cells.per.ident = max.cells.per.ident,
          random.seed = random.seed,
          latent.vars = latent.vars,
          min.cells.feature = min.cells.feature,
          min.cells.group = min.cells.group,
          mean.fxn = mean.fxn,
          fc.name = fc.name,
          base = base,
          densify = densify,
          ...
        ) %>%
          as.data.frame() %>%
          rownames_to_column(var = "genes") %>%
          mutate(
            celltype = celltype,
            comparison = paste0(ident.1, "_vs_", ident.2),
            Expression = case_when(
              avg_log2FC > logfc.threshold & p_val < 0.05 ~ "Up",
              avg_log2FC < -logfc.threshold & p_val < 0.05 ~ "Down",
              TRUE ~ "Unchanged"
            ),
            color = case_when(
              Expression == "Up" ~ "#f2210a",
              Expression == "Down" ~ "#3158e8",
              TRUE ~ "#afb0b3"
            )
          )

        DEG_list[[celltype]] <- deg
        cat("Success!")

      }, error = function(e) {
        cat(sprintf("Failed: %s", e$message))
      })

      rm(sub_obj)
      if (i %% gc_interval == 0) gc(verbose = FALSE)

    } else {
      cat("Skipped (One or both groups missing)")
    }
  }

  # 3. Consolidation
  cat("\n\n[3/3] Merging results into final dataframe...\n")
  DEG_final_df <- dplyr::bind_rows(DEG_list)
  
  end_time <- Sys.time()
  processing_duration <- end_time - start_time
  
  cat("========================================================\n")
  cat(sprintf("Analysis Completed in: %s\n", format(processing_duration)))
  cat("========================================================\n")

  return(DEG_final_df)
}
