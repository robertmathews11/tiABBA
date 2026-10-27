#' Render individual results panel
#'
#' @param sil_res Results from [tiABBA::scan_silhouette]
#' @param cor_mat A square correlation matrix
#' @param cor_cut Numeric value of correlation parameter
#' @param sil_cut Numeric value of silhouette parameter
#' @param outdir File path to group_condition
#' @param size_px Image size to be saved
#' @param max_dim Maximum image size to save
#'
#' @returns A saved heatmap and annotation png
#'
#' @export
#'
render_panel <- function(sil_res, cor_mat,
                         cor_cut, sil_cut,
                         outdir, size_px = 300, max_dim = size_px) {



  # Prepare data for plotting
  ord_names <- row.names(sil_res[["labels"]])

  caption_text <- sprintf(
    "n=%d  var=%.3f  sil=%.3f",
    sil_res[["metric"]][["n_gene"]],
    sil_res[["metric"]][["var"]],
    sil_res[["metric"]][["silhouette"]]
  )

  param_id <- make_param_id(cor_cut, sil_cut)

  # Create correlation heatmap
  cor_sub <- downsample_cor(cor_mat[ord_names, ord_names], max_dim = max_dim)
  heat_names <- rownames(cor_sub)
  cor_df <- cor_sub %>%
    as.data.frame() %>%
    rownames_to_column("row_gene") %>%
    pivot_longer(-.data$row_gene, names_to = "col_gene", values_to = "correlation") %>%
    mutate(
      row_gene = factor(.data$row_gene, levels = heat_names),
      col_gene = factor(.data$col_gene, levels = heat_names)
    )

  heat_plot <- ggplot(cor_df, aes(cor_df$col_gene, cor_df$row_gene, fill = cor_df$correlation)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-1, 1)
    ) +
    coord_cartesian(expand = FALSE) +
    labs(caption = caption_text) +
    guides(fill = "none") +
    theme_void() +
    theme(plot.caption = element_text(size = 4, hjust = 0.5, margin = margin(t = 1)))

  ggsave(file.path(outdir,"heatmap",paste0(param_id,".png")),
         plot = heat_plot, width = size_px, height = size_px, units = "px",
         device = ragg::agg_png)
  rm(cor_sub)
  rm(cor_df)
  rm(heat_plot)

  color_pal <- c("#88CCEEFF", "#CC6677FF", "#DDCC77FF", "#117733FF", "#332288FF", "#AA4499FF", "#44AA99FF", "#999933FF", "#882255FF", "#661100FF", "#6699CCFF", "#888888FF")
  # Create main cluster plot
  cluster_types <- colnames(sil_res[["labels"]])
  walk(cluster_types, \(test_i){

    pull_names <- sil_res[["labels"]] %>% select(test_i)
    annot_df <- tibble(
      col_gene = factor(ord_names, levels = ord_names),
      cluster = factor(pull_names[ord_names,])
    )

    annot_plot <- ggplot(annot_df, aes(annot_df$col_gene, 1, fill = annot_df$cluster)) +
      geom_tile() +
      coord_cartesian(expand = FALSE) +
      scale_fill_manual(values = setNames(color_pal[seq(1:length(unique(annot_df$cluster)))], unique(annot_df$cluster))) +
      guides(fill = "none") +
      theme_void()

    ggsave(file.path(outdir,test_i,paste0(param_id,".png")),
           plot = annot_plot, width = size_px, height = 50, units = "px",
           device = ragg::agg_png)
  })

  return()
}
