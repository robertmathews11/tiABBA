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
