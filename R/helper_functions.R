#' Prepare directories
#'
#' This helper function prepares directories for saving data
#'
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#'
#' @return This function creates file directories and has no return
#'
#' @export

prepare_directory <- function(

  group_file

  ) {

  if(!dir.exists("app")){dir.create("app")}
  if(!dir.exists("app/data")){dir.create("app/data")}
  if(!dir.exists("app/images")){dir.create("app/images")}

  outdir <- file.path("app/data", group_file)
  filedir <- file.path("app/images", group_file)
  if(!dir.exists(outdir)){dir.create(outdir)}
  if(!dir.exists(filedir)){dir.create(filedir)}

  if(!dir.exists(file.path(filedir,"heatmap"))){dir.create(file.path(filedir,"heatmap"))}
  if(!dir.exists(file.path(filedir,"sil_k"))){dir.create(file.path(filedir,"sil_k"))}
  if(!dir.exists(file.path(filedir,"gap_0p5"))){dir.create(file.path(filedir,"gap_0p5"))}
  if(!dir.exists(file.path(filedir,"gap_0p6"))){dir.create(file.path(filedir,"gap_0p6"))}
  if(!dir.exists(file.path(filedir,"gap_0p7"))){dir.create(file.path(filedir,"gap_0p7"))}
  if(!dir.exists(file.path(filedir,"gap_0p8"))){dir.create(file.path(filedir,"gap_0p8"))}
  if(!dir.exists(file.path(filedir,"gap_0p9"))){dir.create(file.path(filedir,"gap_0p9"))}

  return("Directory setup complete")
}


#' Prepare expression data
#'
#' This helper function averages sample data by time point for use in expression plots
#'
#' @param expr A data frame or matrix object, with rows of gene expression readings and columns of sample-time names.
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#'
#' @return A data frame matching the input row length, and columns averaged by time point
#'
#' @export

prepare_expr <- function(

  expr,
  group_file

  ) {

  if(!"matrix" %in% class(expr)){stop("Error: Please input a matrix class object")}

  res <- normalize_rows(expr)

  # Save data
  saveRDS(res, file.path("app/data", group_file, "expression"))

  return(paste("Finished Time Normalization, saved to", file.path("app/data",group_file,"expression")))
}


#' Min-Max normalize data
#'
#' This helper function applied min-max normalization to the input data
#'
#' @param expr A data frame or matrix object to be normalized
#'
#' @return A matrix with rows min-max normalized between -1 and 1
#'
#' @export

normalize_rows <- function(

  expr

  ) {
  # Apply normalization to each row
  normalized_expr <- t(apply(expr, 1, function(row) {
    min_val <- min(row)
    max_val <- max(row)

    # Normalize the row values between -1 and 1
    normalized_row <- (row - min_val) / (max_val - min_val) * 2 - 1
    return(normalized_row)
  }))

  # Convert the result back to a data frame
  normalized_expr <- as.matrix(normalized_expr)

  return(normalized_expr)
}


#' Condense silhouette result
#'
#' @param sil_res A list object containing the results from [totoInferno::silhouette_filter]
#'
#' @returns A list object containing result metrics, clustered gene order, and two cluster assignments.
#'
#' @export

filter_res <- function(

  sil_res

  ) {

  eigen_metric <- sil_res[["metric"]][["eigen"]]
  eigen_metric <- eigen_metric %>%
    distinct(.data$cluster) %>%
    map(\(cluster_i){
      sub_cluster <- eigen_metric[eigen_metric$cluster == cluster_i,] %>%
        arrange(desc(eigen)) %>%
        select(eigen) %>%
        deframe()

      eigen_ratio <- sub_cluster[2] / sub_cluster[1]
      return(data.frame(
        cluster = cluster_i,
        eigen = eigen_ratio
      ))
    })

  res <- list()

  res[["metric"]] <- list(
    silhouette = sil_res[["metric"]][[1]],
    var = sil_res[["metric"]][["within_cluster_var_overall"]],
    eigen = eigen_metric[[1]],
    n_cluster = sil_res[["metric"]][["n_clusters_final"]],
    n_gene = sil_res[["metric"]][["n_genes_out"]]
  )

  order_labels <- sil_res[["data"]][["labels_final"]][sil_res[["data"]][["hclust_final"]][["order"]]]
  res[["labels"]] <- order_labels

  gap_seq <- seq(.5, .9, by = 0.1)
  dist_mat <- as.matrix(sil_res[["data"]][["dist_final"]])

  dynamic_res <- map(gap_seq, \(gap_i){

    dynamic_cut <- cutreeHybrid(sil_res[["data"]][["hclust_final"]], distM = dist_mat,
                                verbose = 0, deepSplit = 0, minClusterSize = 5,
                                minGap = gap_i)
    alt_res <- setNames(dynamic_cut$labels, names(sil_res[["data"]][["labels_final"]]))
    order_alt <- alt_res[sil_res[["data"]][["hclust_final"]][["order"]]]

    return(order_alt)

  }) %>%
    setNames(str_c("gap_", str_replace(as.character(gap_seq), "\\.", "p"))) %>%
    as.data.frame()


  res[["labels"]] <- cbind(
    data.frame(sil_k = order_labels, row.names = names(order_labels)),
    dynamic_res
  )

  return(res)
}


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

  # Create main cluster plot
  cluster_types <- colnames(sil_res[["labels"]])
  walk(cluster_types, \(test_i){

    pull_names <- sil_res[["labels"]] %>% select(test_i) %>% deframe()
    annot_df <- tibble(
      col_gene = factor(ord_names, levels = ord_names),
      cluster = factor(pull_names[ord_names])
    )

    annot_plot <- ggplot(annot_df, aes(annot_df$col_gene, 1, fill = annot_df$cluster)) +
      geom_tile() +
      coord_cartesian(expand = FALSE) +
      guides(fill = "none") +
      theme_void()

    ggsave(file.path(outdir,test_i,paste0(param_id,".png")),
           plot = annot_plot, width = size_px, height = 50, units = "px",
           device = ragg::agg_png)
  })

  return()
}


#' Make parameter file ID
#'
#' @param cor_cut Numeric value of correlation parameter
#' @param sil_cut Numeric value of silhouette parameter
#'
#' @returns Transformed, formatted parameter ID
#'
#' @export

make_param_id <- function(

  cor_cut,
  sil_cut

  ) {
  fmt <- function(x) {
    x <- formatC(as.numeric(x), format = "f", digits = 2)
    x <- gsub("-", "n", x)
    gsub("\\.", "p", x)
  }
  paste0("cor", fmt(cor_cut), "_sil", fmt(sil_cut))
}


#' Downsample large plot data
#'
#' @param sub_mat A subsetting matrix
#' @param max_dim Maximum dimensions for the matrix
#'
#' @returns A downsampled matrix for rasterization
#'
#' @export

downsample_cor <- function(

  sub_mat,
  max_dim

  ) {
  n <- nrow(sub_mat)
  if (n <= max_dim) return(sub_mat)

  bin <- pmin(findInterval(seq_len(n), seq(1, n, length.out = max_dim + 1),
                           rightmost.closed = TRUE), max_dim)
  bin_size <- tabulate(bin, nbins = max_dim)

  row_avg <- rowsum(sub_mat, bin) / bin_size
  binned <- rowsum(t(row_avg), bin) / bin_size

  block_names <- sprintf("b%03d", seq_len(max_dim))
  dimnames(binned) <- list(block_names, block_names)
  binned
}
