#' Create index of results
#'
#' @param sil_res A list containing the results from [tiABBA::scan_silhouette]
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#' @param max_count Numeric containing the initial number of genes, before filtering results.
#'
#' @returns A data.frame object with more accessible indices for the web app
#'
#' @export

scan_index <- function(

  sil_res,
  group_file,
  max_count

){

  index <- imap_dfr(sil_res, \(res_i, cor_i){
    imap_dfr(res_i, \(res_j, sil_j){
      res <- data.frame(
        cor = format(as.numeric(cor_i), nsmall = 2, justify = "none"),
        sil = format(as.numeric(sil_j), nsmall = 2, justify = "none"),
        m_var = res_j$metric$var,
        m_eigen = mean(res_j$metric$eigen$eigen),
        m_gene = (res_j$metric$n_gene / max_count) * 100,
        m_sil = res_j$metric$silhouette
      ) %>%
        mutate(
          file_id = str_c("cor",.data$cor,"_","sil",.data$sil),
          file_id = str_replace_all(.data$file_id, pattern = "\\.", replacement = "p"),
          file_id = str_replace_all(.data$file_id, pattern = "-", replacement = "n"),
          file_id = str_c(.data$file_id, ".png")
        )
    })
  })
  saveRDS(index, file = file.path("app/data", group_file, "index"))
  return(index)
}
