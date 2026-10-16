#addResourcePath("images", normalizePath("images"))
library(tiABBA)
library(shiny)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(patchwork)
library(generics)

color_pal <- c("#88CCEEFF", "#CC6677FF", "#DDCC77FF", "#117733FF", "#332288FF", "#AA4499FF", "#44AA99FF", "#999933FF", "#882255FF", "#661100FF", "#6699CCFF", "#888888FF")

locf_lookup <- function(query, existing) {
  existing_num <- as.numeric(existing)
  ord <- order(existing_num)
  existing_sorted <- existing[ord]
  existing_num_sorted <- existing_num[ord]
  idx <- findInterval(as.numeric(query), existing_num_sorted)
  ifelse(idx == 0, NA_character_, existing_sorted[idx])
}

normalize_key <- function(x) as.character(as.numeric(x))

filter_index <- function(input_index){
  res <- input_index %>%
    mutate(
      r_var = percent_rank(desc(m_var)),
      r_gene = percent_rank(m_gene),
      r_sil = percent_rank(m_sil),
      r_eigen = percent_rank(m_eigen)
    ) %>%
      filter(r_var > 0.05 &
               r_gene > 0.05 &
               r_sil > 0.05 &
               r_eigen > 0.05) %>%
      mutate(r_var = percent_rank(desc(m_var)),
             r_gene = percent_rank(m_gene),
             r_sil = percent_rank(m_sil),
             r_eigen = percent_rank(m_eigen),
             r_total = (r_gene + sqrt(r_var + r_sil + r_eigen)))
  return(res)
}

make_image_cell <- function(match_row, cor_val, sil_val, label = NULL, is_top = FALSE,
                             is_selected = FALSE, cluster_annotation = "sil_k") {
  if (nrow(match_row) == 0) {
    return(div(style = "background: #eee; min-width: 150px; min-height: 150px;"))
  }

  annotation_url <- match_row[[str_c("url_",cluster_annotation)]][1]

  bg <- if (is_top) "#61FF76" else if (is_selected) "#CCCCCC" else "transparent"

  cell_style <- sprintf(
    paste(
      "display: flex; flex-direction: column; width: 100%%;",
      "min-width: 150px; min-height: 150px; cursor: pointer;",
      "background: %s; padding: 10px; box-sizing: border-box;"
    ),
    bg
  )

  div(
    style = cell_style,
    onclick = sprintf(
      "Shiny.setInputValue('clicked_plot', {cor: '%s', sil: '%s'}, {priority: 'event'})",
      cor_val, sil_val
    ),
    if (!is.null(label)) div(label, style = "font-weight: bold; text-align: center;"),
    tags$img(
      src = annotation_url, loading = "lazy",
      style = "width: 100%; aspect-ratio: 300 / 35; display: block;"
    ),
    tags$img(
      src = match_row$url_heatmap[1], loading = "lazy",
      style = "width: 100%; aspect-ratio: 300 / 300; display: block;"
    )
  )
}

group_option <- list.files(file.path(getShinyOption("results_dir"),"data")) %>% str_split_i(pattern = "_", 1) %>% unique()
condition_option <- list.files(file.path(getShinyOption("results_dir"),"data")) %>% str_split_i(pattern = "_", 2) %>% unique()
cluster_option <- c("sil_k", "gap_0p5", "gap_0p6", "gap_0p7", "gap_0p8", "gap_0p9")

ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        # Shared across both tabs so selection state stays in sync
        selectInput("group", "Group:", choices = group_option),
        selectInput("condition", "Condition:", choices = condition_option),
        selectInput("cluster_type", "Assignment:", choices = cluster_option),
        actionButton("clear_selection", "Clear selection"),

        # Only relevant to the Scan Results tab; hidden while on Visualize Group
        conditionalPanel(
          condition = "input.main_tabs == 'Scan Results'",

          ## Maybe have it as some scrolling window instead of a target pointer
          sliderInput("zoom_cor","Correlation Zoom:", min = 0.01, max = 0.1, step = 0.01,
                      value = 0.1, ticks = F),
          sliderInput("limit_cor", "Correlation Limit:", min = 0.1, max = 0.4, step = 0.01,
                      value = c(0.1, 0.4)),
          sliderInput("zoom_sil","Silhouette Zoom:", min = 0.01, max = 0.2, step = 0.01,
                      value = 0.2, ticks = F),
          sliderInput("limit_sil", "Silhouette Limit:", min = -0.50, max = 0.70, step = 0.01,
                      value = c(-0.50, 0.70))
        ),
        selectInput("download_pick", "Download results: ", choices = character(0)),
        downloadButton("download_result")

      ),
      mainPanel(
        tabsetPanel(
          id = "main_tabs",
          tabPanel("Scan Results",
                   textOutput("seq_cor"),
                   textOutput("seq_sil"),
                   uiOutput("plot_grid") ## Use this as the grouped plots
                   ),
          tabPanel("Visualize Group",
                   textOutput("n_selected"),
                   uiOutput("group_columns") ## One column per selected parameter:
                                              ## param text, expression plot, 4 histograms
                   ),
          tabPanel("Top Results",
                   uiOutput("top_grid") ## One row per top-ranked cor value:
                                         ## lower/top/higher sil cells + expression plot
                   )
        )
      )
    )
)

server <- function(input, output, session){
  ### Need:
  # Directory to the image + data cache
  setdir <- reactive({
    str_c(input$group,"_",input$condition)
  })
  # Call in the scan_data for results
  scan <- reactive({readRDS(file.path(getShinyOption("results_dir"),"data", setdir(), "results"))})
  # (later) call in expression data for later use
  expr <- reactive({readRDS(file.path(getShinyOption("results_dir"),"data", setdir(), "expression"))})
  # A reduced scan_data which can be filtered and used for indexing
  index <- reactive({readRDS(file.path(getShinyOption("results_dir"),"data", setdir(), "index"))})

  # Set rows/columns based off the sweep
  make_grid <- reactive({
    res <- index() %>%
      mutate(cor = as.numeric(cor),
             sil = as.numeric(sil))
    list(
      cor = format(seq(input$limit_cor[1], input$limit_cor[2], by = input$zoom_cor), nsmall = 2, trim = TRUE),
      sil = format(seq(input$limit_sil[1], input$limit_sil[2], by = input$zoom_sil), nsmall = 2, trim = TRUE)
    )
  })

  output$seq_cor <- renderText({
    grid <- make_grid()
    paste0("Cor Row: ", length(grid[[1]]), ", Sil Col: ", length(grid[[2]]))
  })

  all_files <- reactive({
    index() %>%
      mutate(
        file_id = str_replace(file_id, pattern = "\\s", replacement = ""),
        url_heatmap = file.path("images", setdir(), "heatmap", file_id),
        url_sil_k = file.path("images", setdir(), "sil_k", file_id),
        url_gap_0p5 = file.path("images", setdir(), "gap_0p5", file_id),
        url_gap_0p6 = file.path("images", setdir(), "gap_0p6", file_id),
        url_gap_0p7 = file.path("images", setdir(), "gap_0p7", file_id),
        url_gap_0p8 = file.path("images", setdir(), "gap_0p8", file_id),
        url_gap_0p9 = file.path("images", setdir(), "gap_0p9", file_id)
      )
  })

  plot_files <- reactive({
    all_files() %>%
      filter(cor %in% make_grid()[["cor"]] & sil %in% make_grid()[["sil"]])
  })

  top_res <- reactive({
    res <- index() %>%
      filter_index() %>%
      arrange(desc(r_total))

    res_filter <- res[!duplicated(res$cor),]


    res_filter <- res_filter %>%
      mutate(top_key = paste(cor, sil, sep = "|")) %>%
      select(top_key) %>%
      deframe()
    return(res_filter)
  })

  top_rows <- reactive({
    keys <- top_res()
    req(length(keys) > 0)

    idx <- index() %>% mutate(sil_num = as.numeric(sil))

    purrr::map_dfr(keys, \(key){
      parts <- strsplit(key, "|", fixed = TRUE)[[1]]
      cor_val <- parts[1]
      sil_val <- parts[2]
      sil_num <- as.numeric(sil_val)

      same_cor <- idx %>% filter(cor == cor_val)

      target_low <- max(sil_num - 0.2, -0.5)
      target_high <- min(sil_num + 0.2, 0.7)

      nearest_sil <- function(target) {
        same_cor$sil[which.min(abs(same_cor$sil_num - target))]
      }

      tibble(
        cor_val = cor_val,
        sil_low = nearest_sil(target_low),
        sil_top = sil_val,
        sil_high = nearest_sil(target_high),
        key = key
      )
    })
  })

  selected_keys <- reactiveVal(character(0))

  observeEvent(input$clicked_plot, {
    key <- paste(input$clicked_plot$cor, input$clicked_plot$sil, sep = "|")
    current <- selected_keys()
    if (key %in% current) {
      selected_keys(setdiff(current, key))
    } else {
      selected_keys(c(current, key))
    }
  })

  observeEvent(input$clear_selection, {
    selected_keys(character(0))
  })

  output$plot_grid <- renderUI({

    cors <- make_grid()[["cor"]]
    sils <- make_grid()[["sil"]]
    files_df <- plot_files()
    cluster_annotation <- input$cluster_type
    keys <- selected_keys()

    effective_cor <- locf_lookup(cors, unique(files_df$cor))

    label_style <- "font-weight: bold; text-align: center; align-self: center;"

    # Blank corner cell + one column-header label per sil value
    header_cells <- c(
      list(div()),
      map(sils, \(s) div(s, style = label_style))
    )

    # Each body row: a cor-value row label followed by one image cell per sil value
    body_cells <- map2(cors, effective_cor, \(cor_val, eff_cor){
      row_label <- div(cor_val, style = label_style)

      img_cells <- map(sils, \(sil_val){
        match_row <- if (is.na(eff_cor)) {
          files_df[0, ]
        } else {
          files_df %>% filter(cor == eff_cor, sil == sil_val)
        }
        cell_key <- paste(eff_cor, sil_val, sep = "|")

        make_image_cell(
          match_row, eff_cor, sil_val,
          is_top = cell_key %in% top_res(),
          is_selected = cell_key %in% keys,
          cluster_annotation = cluster_annotation
        )
      })

      c(list(row_label), img_cells)
    })

    div(
      style = sprintf(
        "display: grid; grid-template-columns: auto repeat(%d, minmax(150px, 1fr)); gap: 2px; align-items: center;",
        length(sils)
      ),
      header_cells,
      body_cells
    )

  })

  ### Top Results tab ----
  output$top_grid <- renderUI({
    rows <- top_rows()
    req(nrow(rows) > 0)

    files_df <- all_files()
    cluster_annotation <- input$cluster_type
    keys <- selected_keys()

    label_style <- "font-weight: bold; text-align: center; align-self: center;"

    header_cells <- list(
      div(),
      div("Lower (-0.2)", style = label_style),
      div("Top Result", style = label_style),
      div("Higher (+0.2)", style = label_style),
      div("Expression", style = label_style)
    )

    body_cells <- pmap(rows, \(cor_val, sil_low, sil_top, sil_high, key){
      row_label <- div(cor_val, style = label_style)

      find_row <- function(sil_val) files_df %>% filter(cor == cor_val, sil == sil_val)

      low_cell <- make_image_cell(
        find_row(sil_low), cor_val, sil_low, label = sil_low,
        is_top = FALSE, is_selected = paste(cor_val, sil_low, sep = "|") %in% keys,
        cluster_annotation = cluster_annotation
      )
      top_cell <- make_image_cell(
        find_row(sil_top), cor_val, sil_top, label = sil_top,
        is_top = TRUE, is_selected = paste(cor_val, sil_top, sep = "|") %in% keys,
        cluster_annotation = cluster_annotation
      )
      high_cell <- make_image_cell(
        find_row(sil_high), cor_val, sil_high, label = sil_high,
        is_top = FALSE, is_selected = paste(cor_val, sil_high, sep = "|") %in% keys,
        cluster_annotation = cluster_annotation
      )

      time_cell <- div(
        style = "min-width: 200px;",
        plotOutput(paste0("top_time_", safe_id(key)), height = "200px")
      )

      list(row_label, low_cell, top_cell, high_cell, time_cell)
    })

    div(
      style = "display: grid; grid-template-columns: auto repeat(3, minmax(150px, 1fr)) minmax(200px, 1fr); gap: 2px; align-items: center;",
      header_cells,
      body_cells
    )
  })

  top_time_data <- reactive({
    rows <- top_rows()
    req(nrow(rows) > 0)
    cluster_annotation <- input$cluster_type
    expr_mat <- expr()

    pmap_dfr(rows, \(cor_val, sil_low, sil_top, sil_high, key){
      res <- scan()[[normalize_key(cor_val)]][[normalize_key(sil_top)]]

      labels <- res[["labels"]][[cluster_annotation]]
      labels <- setNames(labels, row.names(res[["labels"]]))


      clusters <- sort(unique(labels))
      palette <- setNames(color_pal[seq(1:length(clusters))], clusters)

      imap_dfr(split(names(labels), labels), \(genes, cluster_id){
        sub_expr <- expr_mat[intersect(genes, rownames(expr_mat)), , drop = FALSE]
        prepare_time(sub_expr) %>%
          mutate(
            cluster = factor(cluster_id),
            color = palette[[cluster_id]],
            condition = input$condition,
            param_id = key
          )
      })
    })
  })

  observe({
    rows <- top_rows()
    req(nrow(rows) > 0)
    dat_all <- top_time_data()

    for (key in rows$key) {
      local({
        key_local <- key
        id <- safe_id(key_local)

        output[[paste0("top_time_", id)]] <- renderPlot({
          dat <- dat_all %>% filter(param_id == key_local)
          color_res <- setNames(dat$color, dat$cluster)
          suppressWarnings(time_plot(dat, group = cluster, color = color_res, linewidth = 1.8))
        })
      })
    }
  })

  ### Visualize Group tab ----

  selected_results <- reactive({
    req(length(selected_keys()) > 0)
    map(selected_keys(), \(key){
      parts <- strsplit(key, "|", fixed = TRUE)[[1]]
      scan()[[normalize_key(parts[1])]][[normalize_key(parts[2])]]
    }) %>%
      setNames(selected_keys())
  })

  output$n_selected <- renderText({
    req(length(selected_keys()) > 0)
    paste0(length(selected_keys()), " Selected")
  })

  time_data <- reactive({
    cluster_annotation <- input$cluster_type
    expr_mat <- expr()

    imap_dfr(selected_results(), \(res, key){

      labels <- res[["labels"]][[cluster_annotation]]
      labels <- setNames(labels, row.names(res[["labels"]]))

      clusters <- sort(unique(labels))
      palette <- setNames(color_pal[seq(1:length(clusters))], clusters)

      imap_dfr(split(names(labels), labels), \(genes, cluster_id){
        sub_expr <- expr_mat[intersect(genes, rownames(expr_mat)), , drop = FALSE]
        prepare_time(sub_expr) %>%
          mutate(
            cluster = factor(cluster_id),
            color = palette[[cluster_id]],
            condition = input$condition,
            param_id = key
          )
      })
    })
  })

  hist_data <- reactive({
    # Manually set histogram colors by type
    group_colors <- c(
      eigen = "#1B9E77",
      gene  = "#D95F02",
      sil   = "#7570B3",
      var   = "#E7298A"
    )

    # Pull index results and format for ggplot histogram
    index_res <- index() %>%
      mutate(m_gene = (m_gene / 100)) %>%
      select(!"file_id") %>%
      pivot_longer(starts_with("m_"),
                   names_to = "group") %>%
      mutate(group = str_replace(group, "m_", "")) %>%
      mutate(group = factor(group))
    index_plot <- index() %>%
      filter_index() %>%
      mutate(m_gene = (m_gene / 100)) %>%
      select(!"file_id") %>%
      pivot_longer(starts_with("m_"),
                   names_to = "group") %>%
      mutate(group = str_replace(group, "m_", "")) %>%
      mutate(group = factor(group))

    # Map across selected results
    imap(selected_results(), \(res, key){

      # Split results into cor and sil
      parts <- strsplit(key, "|", fixed = TRUE)[[1]]

      # Pull the geom_hline values based off selected result and index
      pull_res <- index_res %>%
        filter(cor == parts[1] & sil == parts[2]) %>%
        mutate(threshold = value)

      # Create 4 plots
       suppressWarnings(map(names(group_colors), \(g_i){
        make_violin_plot(index_plot,
                         pull_res %>% filter(group == g_i),
                         g_i)
      }) %>% set_names(names(group_colors)))
    })
  })

  safe_id <- function(key) gsub("[^A-Za-z0-9]", "_", key)

  output$group_columns <- renderUI({
    keys <- selected_keys()
    req(length(keys) > 0)

    col_width <- 310 # fixed px budget per column: expr_group item 2

    cols <- purrr::map(keys, \(key){
      id <- safe_id(key)
      div(
        style = sprintf("flex: 0 0 %dpx; width: %dpx;", col_width, col_width),
        verbatimTextOutput(paste0("param_text_", id)),
        plotOutput(paste0("expr_plot_", id), height = "250px"),
        plotOutput(paste0("hist_plot_", id), height = "500px")
      )
    })

    div(style = "display: flex; flex-wrap: wrap; gap: 10px; align-items: flex-start;", cols)
  })

  observe({
    keys <- selected_keys()
    req(length(keys) > 0)
    cluster_annotation <- input$cluster_type
    expr_res <- time_data()
    hist_res <- hist_data()
    results <- selected_results()

    for (key in keys) {
      local({
        key_local <- key
        id <- safe_id(key_local)

        output[[paste0("param_text_", id)]] <- renderPrint({
          res <- results[[key_local]]
          parts <- strsplit(key_local, "|", fixed = TRUE)[[1]]

          n_cluster <- res[["labels"]][[cluster_annotation]] %>%
            unique() %>%
            length()

          cat(sprintf(
            "cor = %s, sil = %s\n%d clusters, %d genes\nsilhouette = %.3f",
            parts[1], parts[2], n_cluster, res[["metric"]][["n_gene"]],
            res[["metric"]][["silhouette"]]
          ))
        })

        output[[paste0("expr_plot_", id)]] <- renderPlot({
          dat <- expr_res %>% filter(param_id == key_local)
          color_res <- setNames(dat$color, dat$cluster)
          suppressWarnings(time_plot(dat, group = cluster, color = color_res, linewidth = 1.8))
        })

        output[[paste0("hist_plot_", id)]] <- renderPlot({
          wrap_plots(hist_res[[key_local]], ncol = 2)
        })
      })
    }
  })

  observe({
    keys <- selected_keys()
    current <- isolate(input$download_pick)
    selected <- if (!is.null(current) && current %in% keys) current else if (length(keys) > 0) keys[1] else character(0)
    updateSelectInput(session, "download_pick", choices = keys, selected = selected)
  })

  output$download_result <- downloadHandler(
    filename = function() {
      paste0(setdir(), "_results.tsv")
    },
    content = function(file) {
      req(input$download_pick)
      parts <- strsplit(input$download_pick, "|", fixed = TRUE)[[1]]
      res <- scan()[[normalize_key(parts[1])]][[normalize_key(parts[2])]]

      metric <- res[["metric"]]
      header <- c(
        sprintf("# cor = %s, sil = %s", parts[1], parts[2]),
        sprintf("# silhouette = %s", metric[["silhouette"]]),
        sprintf("# var = %s", metric[["var"]]),
        sprintf("# n_cluster = %s", metric[["n_cluster"]]),
        sprintf("# n_gene = %s", metric[["n_gene"]])
      )

      gene_tbl <- tibble(
        gene = names(res[["labels"]]),
        cluster = res[["labels"]]
      )

      con <- file(file, open = "wt")
      on.exit(close(con))
      writeLines(header, con)
      write.table(gene_tbl, con, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

}

shinyApp(ui, server)
