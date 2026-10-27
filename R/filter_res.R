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
