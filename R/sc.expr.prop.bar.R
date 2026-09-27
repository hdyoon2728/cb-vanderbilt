sc.expr.prop.bar <- function(seurat_obj,
                             target_gene,
                             group_var = "condition",
                             celltype_var = "detailed.celltypes",
                             celltypes = "all",
                             color_pal = NULL,
                             show_legend = TRUE,
                             show_errorbar = TRUE,
                             show_stat = TRUE,
                             p_adjust_method = "BH") {
  
  require(dplyr)
  require(ggplot2)
  require(Seurat)
  require(ggsignif)
  
  missing_genes <- target_gene[!target_gene %in% rownames(seurat_obj)]
  if (length(missing_genes) > 0) {
    stop(paste0("Genes not found in object: ", paste(missing_genes, collapse = ", ")))
  }
  
  meta_data <- seurat_obj@meta.data
  if (length(celltypes) == 1 && celltypes == "all") {
    valid_cells <- rownames(meta_data)
  } else {
    valid_cells <- rownames(meta_data[meta_data[[celltype_var]] %in% celltypes, ])
    if (length(valid_cells) == 0) stop("No matching cell types found.")
  }
  
  raw_counts <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")[target_gene, valid_cells, drop = FALSE]
  
  cell_list <- lapply(target_gene, function(g) {
    data.frame(
      group = meta_data[valid_cells, group_var],
      gene = g,
      expressed = as.numeric(as.numeric(raw_counts[g, ]) > 0)
    )
  })
  cell_df <- do.call(rbind, cell_list) %>% filter(!is.na(group))
  
  prop_list <- lapply(target_gene, function(g) {
    cell_df %>%
      filter(gene == g) %>%
      group_by(group, gene) %>%
      summarise(
        n_total = n(),
        n_pos   = sum(expressed),
        pct = sum(expressed) / n() * 100,
        .groups = "drop"
      )
  })
  prop_df <- do.call(rbind, prop_list)
  
  if (show_errorbar) {
    z <- qnorm(0.975)
    prop_df <- prop_df %>%
      rowwise() %>%
      mutate(
        p_hat  = n_pos / n_total,
        denom  = 1 + z^2 / n_total,
        center = (p_hat + z^2 / (2 * n_total)) / denom,
        margin = (z * sqrt(p_hat * (1 - p_hat) / n_total + z^2 / (4 * n_total^2))) / denom,
        ci_lower = pmax((center - margin) * 100, 0),
        ci_upper = pmin((center + margin) * 100, 100)
      ) %>%
      ungroup() %>%
      select(-p_hat, -denom, -center, -margin)
  } else {
    prop_df$ci_upper <- prop_df$pct
  }
  
  label_offset <- max(prop_df$pct, na.rm = TRUE) * 0.04
  
  ## pairwise.t.test 결과를 comparisons/annotations 리스트로 gene별로 준비 ----
  stat_list <- list()
  if (show_stat) {
    for (g in target_gene) {
      sub <- cell_df %>% filter(gene == g)
      if (length(unique(sub$group)) < 2) next
      
      pt <- pairwise.t.test(sub$expressed, sub$group, p.adjust.method = p_adjust_method)
      pmat <- pt$p.value
      
      comps <- list(); annos <- c()
      for (i in seq_len(nrow(pmat))) {
        for (j in seq_len(ncol(pmat))) {
          pval <- pmat[i, j]
          if (!is.na(pval)) {
            comps[[length(comps) + 1]] <- c(rownames(pmat)[i], colnames(pmat)[j])
            annos <- c(annos, if (pval < 0.001) "***" else if (pval < 0.01) "**"
                       else if (pval < 0.05) "*" else "ns")
          }
        }
      }
      stat_list[[g]] <- list(comparisons = comps, annotations = annos)
    }
  }
  
  ## x축 = group, facet = gene 구조로 변경 -------------------------------------
  p <- ggplot(prop_df, aes(x = group, y = pct, fill = group)) +
    geom_col(color = "black", linewidth = 0.5, width = 0.7) +
    
    {if (show_errorbar) geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                                        width = 0.15, linewidth = 0.4)} +
    
    geom_text(aes(label = sprintf("%.1f%%", pct), y = pct + label_offset), 
              vjust = 0, fontface = "bold", size = 3.5) +
    
    facet_wrap(~ gene, scales = "free_x") +
    
    theme_classic() +
    labs(title = NULL, subtitle = NULL,
         y = "Percent of Cells (%)", x = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    theme(
      axis.text.x = element_text(face = "bold", size = 10, angle = 0),
      axis.title.y = element_text(face = "bold", size = 12),
      legend.title = element_blank(),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  ## gene별로 geom_signif를 따로 추가 (facet 대응) ------------------------------
  if (show_stat) {
    base_y <- max(prop_df$ci_upper, na.rm = TRUE)
    for (g in names(stat_list)) {
      info <- stat_list[[g]]
      if (length(info$comparisons) == 0) next
      
      p <- p + geom_signif(
        data = data.frame(gene = g),
        comparisons = info$comparisons,
        annotations = info$annotations,
        y_position = base_y * (1.05 + 0.12 * seq_along(info$comparisons)),
        tip_length = 0.01,
        textsize = 4,
        inherit.aes = FALSE,
        mapping = aes(x = NULL, y = NULL)   # 메인 plot aes 상속 방지
      )
    }
  }
  
  if (!is.null(color_pal)) p <- p + scale_fill_manual(values = color_pal)
  if (!show_legend) p <- p + theme(legend.position = "none")
  
  return(p)
}
