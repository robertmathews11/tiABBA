#' Iterate Correlation Filters
#'
#' This function implements [totoInferno::cor_filter] across a parameter sequence. Saves intermediate results
#'
#' @param expr A matrix of expression data. Rows as genes and columns as time points.
#' @param cor_seq A sequence of correlation filters to run on the data. Default value is the recommended range.
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#' @param verbose Optional TRUE/FALSE to print progress along iterations.
#'
#' @returns A list object containing correlation filters tested and the genes kept. Duplicated results are removed.
#'
#' @export

scan_correlation <- function(

  expr,
  cor_seq = seq(.1, .4, by = .01),
  group_file,
  verbose = TRUE

  ) {
  if(!"matrix" %in% class(expr)){stop("Error: Please input a matrix class object")}

  expr <- as.matrix(expr)
  cor_res <- list()
  cor_data <- expr
  if(verbose == TRUE){
    cat("Running correlation scan: ")
    progress_bar <- txtProgressBar(min = 0, max = length(cor_seq), initial = 0, style = 3)
    }
  for(i in seq_along(cor_seq)){
    res_i <- tryCatch(cor_filter(expr = cor_data,
                                 cor_threshold = cor_seq[i],
                                 return_output = "expr"),

                      error = function(e) NULL)

    if(verbose == TRUE){setTxtProgressBar(progress_bar, i)}
    if(!is.null(res_i)){
      cor_res[[str_c(cor_seq[i])]] <- row.names(res_i)
      cor_data <- cor_data[row.names(res_i), , drop = FALSE]
    }else{break}
  }
  if(verbose == TRUE){close(progress_bar)}
  rm(cor_data)


  # filter to exclude duplicates
  cor_keep <- lengths(cor_res) %>% diff()
  cor_keep <- cor_keep[cor_keep != 0]

  cor_res <- cor_res[c("0.1", names(cor_keep))]

  saveRDS(cor_res, file.path("app/data", group_file, "cor_results"))
  return(cor_res)
}


#' Iterate Silhouette Filters
#'
#' This function implements [totoInferno::silhouette_filter] across a parameter sequence. Saves intermediate results.
#'
#' @param expr A matrix of expression data. Rows as genes and columns as time points.
#' @param cor_res A list containing the results from [tiABBA::scan_correlation]
#' @param sil_seq A sequence of silhouette filters to run on the data. Default value is the recommended range.
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#' @param k_max The maximum cluster amount to test. Results greater than 8 clusters may break plotting, along with being must more expensive to run.
#' @param run_parallel FALSE or an integer setting the number of cores to use in parallel processing. Recommended to use [parallel::detectCores] and not exceed free memory divided by maximum object size.
#' @param verbose Optional TRUE/FALSE to print progress along iterations.
#'
#' @returns A list object with multiple levels: correlation parameter, silhouette parameter, then results
#'
#' @export

scan_silhouette <- function(

  expr,
  cor_res,
  sil_seq = seq(-.5, .7, by = .01),
  group_file,
  k_max = 8L,
  run_parallel = FALSE,
  verbose = TRUE

  ) {
  if(!"matrix" %in% class(expr)){stop("Error: Please input a matrix class object")}

  ## Run silhouette, output panels
  sil_seq = seq(-.5, .7, by = .01)

  if(verbose == TRUE){cat("Running silhouette scan: ")}
  if(run_parallel == FALSE){
  sil_res <- imap(cor_res, \(cor_i, cor_name) {
    cor_filter <- expr[cor_i, , drop = FALSE]
    pre_cor <- cor(t(cor_filter), method = "pearson")

    sil_res <- list()
    for(i in seq_along(sil_seq)){
      res_i <- tryCatch(silhouette_filter(expr = cor_filter,
                                          sil_threshold = sil_seq[i],
                                          cor_expr = pre_cor,
                                          k_max = k_max),

                        error = function(e) NULL)

      if(!is.null(res_i)){

        sil_res[[str_c(sil_seq[i])]] <- filter_res(res_i)

        cor_filter <- cor_filter[row.names(res_i[["expr"]]), , drop = FALSE]
      }else{break}
    }
    return(sil_res)
  }, .progress = {verbose})
  }else{
    plan(multisession, workers = run_parallel)
    sil_res <- future_imap(cor_res, \(cor_i, cor_name) {
      cor_filter <- expr[cor_i, , drop = FALSE]
      pre_cor <- cor(t(cor_filter), method = "pearson")

      sil_res <- list()
      for(i in seq_along(sil_seq)){
        res_i <- tryCatch(silhouette_filter(expr = cor_filter,
                                            sil_threshold = sil_seq[i],
                                            cor_expr = pre_cor,
                                            k_max = 8L),

                          error = function(e) NULL)

        if(!is.null(res_i)){

          sil_res[[str_c(sil_seq[i])]] <- filter_res(res_i)

          cor_filter <- cor_filter[row.names(res_i[["expr"]]), , drop = FALSE]
        }else{break}
      }
      return(sil_res)
    }, .progress = {verbose})
    plan(sequential)
  }

  sil_res <- map(sil_res, \(res_i){
    correct_names <- names(res_i) %>%
      as.numeric() %>%
      round(digits = 2)
    return(setNames(res_i, correct_names))
  })

  saveRDS(sil_res, file = file.path("app/data", group_file, "results"))
  return(sil_res)
}


#' Plot parameter scan results
#'
#' This function prepares plots the results generated by [tiABBA::scan_silhouette] for use in the web interface
#'
#' @param expr A matrix of expression data. Rows as genes and columns as time points.
#' @param sil_res A list containing the results from [tiABBA::scan_silhouette]
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#' @param set_workers An integer to set the amount of parallel workers used in plot generation.
#' @param verbose Optional TRUE/FALSE to print progress along iterations.
#'
#' @returns Saved PNG files based on each correlation and silhouette combination. Saved to images/group_condition/ folders
#'
#' @export

scan_plot <- function(

  expr,
  sil_res,
  group_file,
  set_workers = {detectCores() - 1},
  verbose = TRUE

  ) {

  filedir <- file.path("app/images", group_file)

  # Plot parameter results
  cor_input <- cor(t(expr), method = "pearson")

  if(verbose == TRUE){cat("Plotting results: ")}
  plan(multisession, workers = set_workers)
  iwalk(sil_res, \(cor_i, cor_name){

    cor_genes <- cor_i %>%
      map(\(res_j) row.names(res_j[["labels"]])) %>%
      reduce(union)

    cor_sub <- cor_input[cor_genes, cor_genes]

    future_imap(cor_i, \(res_j, sil_name){
      render_panel(res_j, cor_mat = cor_sub,
                   cor_cut = cor_name, sil_cut = sil_name,
                   outdir = filedir)
    }, .options = furrr_options(seed = TRUE))
  }, .progress = {verbose})

  plan(sequential)

  message("Saved plots to ", filedir)
  return(filedir)
}

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


#' Prepare app data
#'
#' This function is a wrapper for other functions to create app data based off initial data
#'
#' @param expr A matrix of expression data. Rows as genes and columns as time points.
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#' @param cor_seq A sequence of correlation filters to run on the data. Default value is the recommended range.
#' @param sil_seq A sequence of silhouette filters to run on the data. Default value is the recommended range.
#' @param k_max The maximum cluster amount to test. Results greater than 8 clusters may break plotting, along with being must more expensive to run.
#' @param parallel_silhouette FALSE or an integer setting the number of cores to use in parallel processing. Recommended to use [parallel::detectCores] and not exceed free memory divided by maximum object size.
#' @param set_workers An integer to set the amount of parallel workers used in plot generation.
#' @param verbose Optional TRUE/FALSE to print progress along iterations.
#'
#' @returns Locally saved results to data/group_condition and images to images/group_condition
#'
#' @export

wrap_app_inputs <- function(

  expr,
  group_file,
  cor_seq = seq(.1, .4, by = .01),
  sil_seq = seq(-.5, .7, by = .01),
  k_max = 8L,
  parallel_silhouette = FALSE,
  set_workers = {detectCores() - 1},
  verbose = TRUE

  ) {

  if(!"matrix" %in% class(expr)){stop("Error: Please input a matrix class object")}

  max_count <- nrow(expr)
  prepare_directory(group_file)

  cor_res <- scan_correlation(
    expr = expr,
    cor_seq = cor_seq,
    group_file = group_file,
    verbose = verbose
  )

  sil_res <- scan_silhouette(
    expr = expr,
    cor_res = cor_res,
    sil_seq = sil_seq,
    group_file = group_file,
    k_max = k_max,
    run_parallel = parallel_silhouette,
    verbose = verbose
  )

  scan_index(
    sil_res = sil_res,
    group_file = group_file,
    max_count = max_count
  )

  scan_plot(
    expr = expr,
    sil_res = sil_res,
    group_file = group_file,
    set_workers = set_workers,
    verbose = verbose
  )

  return(paste("Finished generating app data for", group_file))
}


#' Run shiny app
#'
#' This function sources and wraps app.R to display the web app
#'
#' @param results_dir A string file.path which connects the app to the saved results. Uses root directory by default.
#'
#' @returns A locally run instance of the shiny app
#'
#' @export

launch_tiABBA <- function(

  results_dir = "app"

  ) {

  set_dir <- file.path(getwd(), results_dir)

  addResourcePath("images", normalizePath(file.path(set_dir, "images")))

  app_dir <- system.file("app", package = "tiABBA")
  shinyOptions(results_dir = set_dir)
  runApp(app_dir)
}
