plot_mediation_gradients <- function(subject_data, treatments, outcomes, covariates, gradient_data, gradients = c(1, 3),
                                     parcels,
                                     gradient_colors = gradient_cols,
                                     plt_title = "",
                                     terms_nice = treatments,
                                     plot_spacing = 0.3,
                                     point_size = 1,
                                     point_alpha = 0.2,
                                     padding_med_str = 0.4,
                                     empty_row_height = -0.1,
                                     base_size_ = 14,
                                     med_lab_size = 11,
                                     l_width = 0.1,
                                     n_lab_size = 14,
                                     add_shade = TRUE,
                                     shade_size = 0.1, 
                                     shade_alpha = 0.02,
                                     rectangle = FALSE) {
  require(ggtext)
  require(patchwork)
  require(scales)
  
  if (length(treatments) != length(outcomes) & length(outcomes) !=1 ) stop("Must have the same number of treatments as outcomes")
  if (length(outcomes) == 1) outcomes <- rep(outcomes, length(treatments))
  
  net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                          col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), #"#DCF8A4"
                          label = c(1:7))
  
  
  
  ests <- tibble()
  for (i in 1:length(treatments)) {
    
    ests_i <- fast_mediation_all(parcel_names = parcels, 
                                 treat = treatments[i], 
                                 outcome = outcomes[i], 
                                 covars = covariates[[i]], 
                                 data = subject_data)
    ests_i <- ests_i |> 
      mutate(
        term = treatments[i], 
        outcome = outcomes[i], 
        covars = covariates[i],
        region = rois
      )
    
    ests <- rbind(ests, ests_i)
  }
  
  
  
  
  atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds")
  perms = readRDS("data/atlas_data/permutations_1000_hungarian.rds")
  
  
  est_plots <- list()
  
  for(i in seq_along(treatments)) {
    
    med_str <- plot_mediation_structure(
      treat = treatments[i],
      mediator = "FCS",
      outcome = outcomes[i],
      lab_size = med_lab_size
    )
    
    plot_ests  <- ests %>% filter(term == treatments[i])
    
    n <- subject_data |> drop_na(all_of(c(treatments[i], outcomes[i], covariates[[i]]))) |> count() |> deframe()
    
    brain_plot <- plot_ests |> 
      right_join(schaef1k$atlas) |> 
      ggplot() +
      geom_sf(aes(geometry = geometry, fill = indirect), linewidth = l_width, show.legend = FALSE) +
      (if (add_shade)
        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
       else NULL) +
      theme_void(base_size = base_size_) +
      scale_fill_gradient2(low = muted("blue"), high = muted("red"))
    
    med_plot <- ggdraw() +
      draw_label(ifelse(n_lab_size != 0, paste0("n = ", n), ""), x = -0.2, y = 0.975, hjust = 0, size = n_lab_size) +
      draw_plot(med_str, x = 0, y = 0.8, width = 1, height = 0.2) +
      draw_plot(brain_plot,  x = 0, y = 0, width = 1, height = 0.8)
    
    est_plots[[treatments[i]]] <- med_plot
    
  }
  
  
  
  grad_char <- paste0("gradient", gradients)
  gradient_plots <- list()
  i = 1
  for (grad in grad_char) {
    
    gradient_plots[[paste0(grad)]] <- gradient_data  %>% 
      right_join(atlas_geometry$atlas, by = "region") %>%
      ggplot() +
      geom_sf(aes(
        fill = .data[[grad]],
        geometry = geometry), linewidth= l_width,
        show.legend = FALSE)+
      (if (add_shade)
        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
       else NULL) +
      theme_void(base_size = base_size_)+
      labs(fill = "", title = str_to_title(str_replace(paste0("_", grad), "_", " "))
      ) + labs(title = case_when(
        grad == "gradient1" ~ "Sens-Assoc",
        grad == "gradient2" ~ "Vis-Mot",
        grad == "gradient3" ~ "Rep-Exec",
        TRUE ~ grad
      )) +
      theme(legend.position = "",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(color = "black", hjust = 0.5, vjust = -5)
      ) +
      scale_fill_gradient2(
        low = gradient_colors[[grad]][1],
        mid = "white",
        high = gradient_colors[[grad]][2] 
      ) 
    
    i = i +1
  }
  
  plots <- list()
  count <- 1
  for (g in grad_char) {
    for (i in seq_along(treatments)) {
      term_ = treatments[i]
      plot_data  <- gradient_data %>% inner_join(ests %>% filter(term == term_), by = "region") 
      
      j <- g
      
      
      lab_grad <- "Gradient score"
      lab_term <- paste0(terms_nice[i], "→", "FCS" , "→", outcomes[i]) 
      
      x_min <- min(plot_data["indirect"])
      x_max <- max(plot_data["indirect"])
      
      y_min <- min(plot_data[j])
      y_max <- max(plot_data[j])
      
      
      
      p <- plot_data %>% 
        ggplot(aes(x = .data[["indirect"]], y =.data[[j]],
                   color = name)) +
        geom_point(alpha =  point_alpha, size = point_size, show.legend = FALSE) +
        stat_poly_line(se = FALSE,
                       color = "#323232") +
        # stat_poly_line(data = plot_data |> filter(.data[[j]] > 0),
        #                se = FALSE,
        #                color = "#323232") +
        # stat_poly_line(data = plot_data |> filter(.data[[j]] < 0),
        #                se = FALSE,
        #                color = "#323232") +
        labs(#title = term_,
          x = lab_term,
          y = "Gradient Score",
          #tag = tag_labs[count],
          color = "Network") +
        xlim(x_min, x_max)+
        #ylim(y_min, y_max)+
        theme_bw(base_size = base_size_) +
        theme(
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          axis.title.x = element_text(size = rel(0.95))
        ) +
        scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
        scale_y_continuous(limits = c(y_min, y_max), position = "left")
      
      r <- cor(plot_data[["indirect"]], plot_data[[j]], method = "pearson")
      p_val <- perm_sphere_p(plot_data[["indirect"]], plot_data[[j]], perm.id = perms, corr.type='pearson')

      #p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
      #label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
      p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("<i>p<sub>spin</sub></i> &lt; ", "<i>p<sub>spin</sub></i> = ", "<i>p<sub>spin</sub></i> &gt; "))(p_val)
      label_corr <- paste0("<i>r</i> = ", r |> round(2), " ,&#8203;&nbsp;" , p_lab)


      p <- p + labs(subtitle = label_corr) +
        theme(plot.subtitle = element_markdown(size = rel(1), hjust = 0.95, margin = margin(b = 3),
                                               color = "black"))

      
      
      require(ggside)
      x_range_df <- data.frame(x = seq(x_min, x_max, length.out = 100))
      colnames(x_range_df)[1] <- "indirect"
      
      y_range_df <- data.frame(y = seq(y_min, y_max, length.out = 100))
      colnames(y_range_df)[1] <- j
      
      
      alpha_switch_xside <- ifelse(g == tail(grad_char, 1), 1, 0)
      alpha_switch_yside <- ifelse(term_ == treatments[1], 1, 0)
      x_axis_colorbar <- c(muted("blue"), muted("red"))
      y_axis_colorbar <- gradient_colors[[g]]
      
      
      
      p <- p +
        geom_xsidetile(data = x_range_df, aes(x = .data[["indirect"]], y = 0, fill = .data[["indirect"]]), 
                       alpha = alpha_switch_xside,
                       show.legend = FALSE, 
                       inherit.aes = FALSE) +
        scale_fill_gradient2(
          low = x_axis_colorbar[1],
          mid = "white",
          high = x_axis_colorbar[2] 
        ) +
        ggnewscale::new_scale_fill() +
        geom_ysidetile(data = y_range_df, aes(y = .data[[j]], x = 0, fill = .data[[j]]), 
                       alpha = alpha_switch_yside,
                       show.legend = FALSE,
                       inherit.aes = FALSE) +
        scale_fill_gradient2(
          low = y_axis_colorbar[1],
          mid = "white",
          high = y_axis_colorbar[2]
        ) +
        theme_ggside_void() +
        ggside(x.pos = "bottom", y.pos = "left") +
        theme(ggside.panel.scale = 0.04)
      
      
      plots[[count]] <- p
      count <- count + 1
    }
  }
  
  n_terms = length(treatments)
  n_gradients = length(gradients)
  n_plot_cols <- n_terms + 1
  col_widths <- c(1, c(rep(1, n_terms), plot_spacing))[-(n_plot_cols+1)]
  empty_cols <- seq(1, n_plot_cols, by = n_terms+1)[-1]
  n_plot_rows <- n_gradients + 1 + (n_gradients)  
  row_heights <- rep(1, n_plot_rows)  
  empty_rows <- seq(2, n_plot_rows, by = 2)  
  row_heights[empty_rows] <- empty_row_height
  
  layout <- c(
    area(3, 1)  # First gradient
  )
  
  if (n_gradients>1) {
    for (g in seq(5, n_plot_rows, by =2)){
      layout <- c(layout, area(g, 1))
    }
  }
  
  for (col in 2:n_plot_cols){
    if (col %in% empty_cols) {
    } else {
      layout <- c(layout, area(1, col))
    }
  }
  
  for(i in seq(3, (n_gradients*2-1)+2, by = 2)){
    for (j in 2:((n_terms) + 1)) {
      if ((j %in% empty_cols)) {
      } else {
        layout <- c(layout, area(i, j))
      }
    }
  }
  
  for(empt in empty_cols){
    layout <- c(layout, area(1, empt , b = n_gradients + 1))
  }
  
  # if (padding != 0) {
  #   for(pad in 1:padding){
  #     layout <- c(layout, area(1, n_plot_cols + pad, b = n_gradients + 1))
  #   }
  # }
  
  plots_to_include <- c(gradient_plots, est_plots, plots)
  
  # if (plt_subtitle) {
  #   f <- mod_formula
  #   rhs <- paste0(f)[2]
  #   rhs <- gsub("\\*", "×", rhs)
  #   expr_str <- paste0("FC[parcel] ~ '~' ~ ", rhs)
  #   #subtit_expr <- parse(text = expr_str)[[1]]
  #   subtit_expr <- bquote(FC[parcel] ~ "~" ~ .(rhs))
  # } else {
  #   subtit_expr <- ""
  # }
  
  # if (!is.null(plt_title)) {
  #   n = ests %>% pull(n) %>% unique()
  #   plt_title <- paste0(plt_title, "(N = ", n, ")")
  # }
  
  p <- Reduce(`+`, plots_to_include) +
    plot_annotation(title = plt_title, 
                    #subtitle = subtit_expr,
                    #tag_levels = c('A', '1'),
                    theme = theme(
                      plot.title = element_text(#size = 28,
                        hjust = 0),
                      plot.subtitle = element_text(size = rel(0.9),
                                                   hjust = 0,
                                                   vjust = -0.05,
                                                   family = "mono",
                                                   face = "italic")
                    )) +
    plot_layout(design = layout, axis_titles = "collect", axes = "collect", guides = "collect",
                widths = col_widths,
                heights = row_heights
    ) 
  if (rectangle) {
    p <- p + plot_annotation(
      theme = theme(plot.background = element_rect(color = "black", fill = NA)))
  }
  p
}

plot_mediation_barplot <- function(subject_data, treatments, outcomes, covariates, gradient_data, gradients = c(1, 3),
                                     parcels,
                                     gradient_colors = gradient_cols,
                                     plt_title = "",
                                   plt_title_size = rel(1),
                                     terms_nice = treatments,
                                     plot_spacing = 0.3,
                                     point_size = 1,
                                     point_alpha = 0.2,
                                     padding_med_str = 0.4,
                                     empty_row_height = -0.1,
                                     base_size_ = 14,
                                     med_lab_size = 11,
                                     l_width = 0.1,
                                     n_lab_size = 14,
                                     add_shade = TRUE,
                                     shade_size = 0.1, 
                                      shade_alpha = 0.02,
                                     rectangle = FALSE,
                                     rasterize = FALSE) {
  require(ggtext)
  require(patchwork)
  require(scales)
  require(png)
  require(ragg)
  
  if (length(treatments) != length(outcomes) & length(outcomes) !=1 ) stop("Must have the same number of treatments as outcomes")
  if (length(outcomes) == 1) outcomes <- rep(outcomes, length(treatments))
  
  net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                          col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), #"#DCF8A4"
                          label = c(1:7))
  
  ests <- tibble()
  for (i in 1:length(treatments)) {
    
    ests_i <- fast_mediation_all(parcel_names = parcels, 
                                 treat = treatments[i], 
                                 outcome = outcomes[i], 
                                 covars = covariates[[i]], 
                                 data = subject_data)
    ests_i <- ests_i |> 
      mutate(
        term = treatments[i], 
        outcome = outcomes[i], 
        covars = covariates[i],
        region = rois
      )
    
    ests <- rbind(ests, ests_i)
  }
  
  
  
  
  atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds")
  perms = readRDS("data/atlas_data/permutations_1000_hungarian.rds")

  make_raster_overlay <- function(plot, bbox,
                                  width_px = 3000,
                                  height_px = 3000,
                                  dpi = 300,
                                  bg = "transparent") {
    tf <- tempfile(fileext = ".png")
    
    ggsave(
      filename = tf,
      plot = plot,
      width = width_px / dpi,
      height = height_px / dpi,
      dpi = dpi,
      bg = bg,
      device = ragg::agg_png
    )
    
    img <- png::readPNG(tf)
    unlink(tf)
    
    annotation_raster(
      raster = as.raster(img),
      xmin = bbox["xmin"],
      xmax = bbox["xmax"],
      ymin = bbox["ymin"],
      ymax = bbox["ymax"]
    )
  }
  
  make_atlas_raster_layer <- function(atlas_sf,
                                      shade_sf = NULL,
                                      shade_size = 0.0005,
                                      shade_alpha = 0.005,
                                      width_px = 3000,
                                      dpi = 300,
                                      fill_limits = NULL,
                                      fill_scale = scale_fill_gradient2(
                                        low = scales::muted("blue"),
                                        mid = "white",
                                        high = scales::muted("red"),
                                        limits = fill_limits
                                      )) {
    bbox <- sf::st_bbox(atlas_sf)
    
    ratio <- as.numeric((bbox["xmax"] - bbox["xmin"]) / (bbox["ymax"] - bbox["ymin"]))
    height_px <- max(1, round(width_px / ratio))
    
    p_raster <- ggplot() +
      geom_sf(
        data = atlas_sf,
        aes(fill = indirect, geometry = geometry),
        linewidth = l_width,
        show.legend = FALSE
      ) +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]),
        expand = FALSE
      ) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)
      ) +
      fill_scale
    
    if (!is.null(shade_sf)) {
      p_raster <- p_raster +
        geom_sf(
          data = shade_sf,
          size = shade_size,
          alpha = shade_alpha,
          show.legend = FALSE
        )
    }
    
    make_raster_overlay(
      plot = p_raster,
      bbox = bbox,
      width_px = width_px,
      height_px = height_px,
      dpi = dpi
    )
  }
  
  
  est_plots <- list()
  
  for(i in seq_along(treatments)) {
    
    med_str <- plot_mediation_structure(
      treat = treatments[i],
      mediator = "FCS",
      outcome = outcomes[i],
      lab_size = med_lab_size
    )
    
    plot_ests  <- ests %>% filter(term == treatments[i])
    
    n <- subject_data |> drop_na(all_of(c(treatments[i], outcomes[i], covariates[[i]]))) |> count() |> deframe()
    
    atlas_plot_df <- plot_ests |>
      right_join(atlas_geometry$atlas, by = "region") |>
      sf::st_as_sf()
    
    fill_limits <- range(atlas_plot_df$indirect, na.rm = TRUE)
    bbox <- sf::st_bbox(atlas_geometry$atlas)
    
    raster_layer <- NULL
    if (rasterize) {
      raster_layer <- make_atlas_raster_layer(
        atlas_sf = atlas_plot_df,
        shade_sf = if (add_shade) atlas_geometry$shade else NULL,
        shade_size = shade_size,
        shade_alpha = shade_alpha,
        width_px = 1500,
        dpi = 300,
        fill_limits = fill_limits
      )
    }
    
    brain_plot <- ggplot() +
      (
        if (rasterize) {
          raster_layer
        } else {
          geom_sf(
            data = atlas_plot_df,
            aes(geometry = geometry, fill = indirect),
            linewidth = l_width,
            show.legend = FALSE
          )
        }
      ) +
      (
        if (!rasterize && add_shade) {
          geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
        } else {
          NULL
        }
      ) +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]),
        expand = FALSE
      ) +
      theme_void(base_size = base_size_) +
      theme(
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)
      ) +
      scale_fill_gradient2(low = muted("blue"), high = muted("red"), limits = fill_limits)
    
    med_plot <- ggdraw() +
      draw_label(ifelse(n_lab_size != 0, paste0("n = ", n), ""), x = 0.0, y = 0.975, hjust = 0, size = n_lab_size) +
      draw_plot(med_str, x = 0, y = 0.8, width = 1, height = 0.2) +
      draw_plot(brain_plot,  x = 0, y = 0, width = 1, height = 0.8)
    
    est_plots[[treatments[i]]] <- med_plot
    
  }
  
  
  
  grad_char <- paste0("gradient", gradients)
  # gradient_plots <- list()
  # i = 1
  # for (grad in grad_char) {
  #   
  #   gradient_plots[[paste0(grad)]] <- gradient_data  %>% 
  #     right_join(atlas_geometry$atlas, by = "region") %>%
  #     ggplot() +
  #     geom_sf(aes(
  #       fill = .data[[grad]],
  #       geometry = geometry), linewidth= l_width,
  #       show.legend = FALSE)+
  #     (if (add_shade)
  #       geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
  #      else NULL) +
  #     theme_void(base_size = base_size_)+
  #     labs(fill = "", title = str_to_title(str_replace(paste0("_", grad), "_", " "))
  #     ) + labs(title = case_when(
  #       grad == "gradient1" ~ "Sens-Assoc",
  #       grad == "gradient2" ~ "Vis-Mot",
  #       grad == "gradient3" ~ "Rep-Exec",
  #       TRUE ~ grad
  #     )) +
  #     theme(legend.position = "",
  #           panel.background = element_rect(fill = "transparent", colour = NA),
  #           plot.background = element_rect(fill = "transparent", colour = NA),
  #           legend.background = element_rect(fill = "transparent", colour = NA),
  #           legend.box.background = element_rect(fill = "transparent", colour = NA),
  #           plot.title = element_text(color = "black", hjust = 0.5, vjust = -5)
  #     ) +
  #     scale_fill_gradient2(
  #       low = gradient_colors[[grad]][1],
  #       mid = "white",
  #       high = gradient_colors[[grad]][2] 
  #     ) 
  #   
  #   i = i +1
  # }
  # 
  plots <- list()
  count <- 1
  #for (g in grad_char) {
    for (i in seq_along(treatments)) {
      term_ = treatments[i]
      if (grepl("age", term_)) g <- "gradient3"
      if (grepl("pathology", term_)) g <- "gradient1"
      
      plot_data  <- gradient_data %>% inner_join(ests %>% filter(term == term_), by = "region") 
      
      j <- g
      
      
      if (g == "gradient3") lab_grad <- "Rep-Exec score"
      if (g == "gradient1") lab_grad <- "Sens-Assoc score"
      lab_term <- "Indirect effect"  #paste0(terms_nice[i], "→", "FCS" , "→", outcomes[i]) 
      
      x_min <- min(plot_data["indirect"])
      x_max <- max(plot_data["indirect"])
      
      y_min <- min(plot_data[j])
      y_max <- max(plot_data[j])
      
      
      
      p <- plot_data %>% 
        ggplot(aes(y = .data[["indirect"]], 
                   x =.data[[j]],
                   fill = name)) +
        geom_col(alpha =  0.6,
                 color = NA,
                 show.legend = FALSE, 
                 width = 0.55) +
        # stat_poly_line(se = FALSE,
        #                color = "#323232") +
        # stat_poly_line(data = plot_data |> filter(.data[[j]] > 0),
        #                se = FALSE,
        #                color = "#323232") +
        # stat_poly_line(data = plot_data |> filter(.data[[j]] < 0),
        #                se = FALSE,
        #                color = "#323232") +
        labs(#title = term_,
          y = lab_term,
          x = lab_grad,
          #tag = tag_labs[count],
          color = "Network") +
        #xlim(x_min, x_max)+
        #ylim(y_min, y_max)+
        theme_bw(base_size = base_size_) +
        theme(
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          axis.title = element_text(size = rel(0.85)),
          axis.text = element_text(size = rel(0.75))
        ) +
        #coord_flip() +
        scale_fill_manual(values = net_names %>% select(name, col) %>% deframe())  +
        ggnewscale::new_scale_fill() 
        #scale_y_continuous(limits = c(y_min, y_max), position = "left")
      
      # # r <- cor(plot_data[["indirect"]], plot_data[[j]], method = "pearson")
      # # p_val <- perm_sphere_p(plot_data[["indirect"]], plot_data[[j]], perm.id = perms, corr.type='pearson')
      # 
      # #p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
      # #label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
      # # p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("<i>p<sub>spin</sub></i> &lt; ", "<i>p<sub>spin</sub></i> = ", "<i>p<sub>spin</sub></i> &gt; "))(p_val)
      # # label_corr <- paste0("<i>r</i> = ", r |> round(2), " ,&#8203;&nbsp;" , p_lab)
      # 
      # 
      # # p <- p + labs(subtitle = label_corr) +
      # #   theme(plot.subtitle = element_markdown(size = rel(1), hjust = 0.95, margin = margin(b = 3),
      # #                                          color = "black"))
      
      
      
      require(ggside)
      y_range_df <- data.frame(y = seq(y_min, y_max, length.out = 100))
      colnames(y_range_df)[1] <- j
      
      x_range_df <- data.frame(x = seq(x_min, x_max, length.out = 100))
      colnames(x_range_df)[1] <- "indirect"
      
      
      alpha_switch_yside <- ifelse(g == tail(grad_char, 1), 1, 0)
      alpha_switch_xside <- ifelse(term_ == treatments[1], 1, 0)
      y_axis_colorbar <- c(muted("blue"), muted("red"))
      x_axis_colorbar <- gradient_colors[[g]]
      
      
      
      p <- p +
        geom_ysidetile(data = x_range_df, aes(y = .data[["indirect"]], 
                                              x = 0, 
                                              fill = .data[["indirect"]]), 
                       #alpha = alpha_switch_xside,
                       show.legend = FALSE, 
                       inherit.aes = FALSE) +
        scale_fill_gradient2(
          low = y_axis_colorbar[1],
          mid = "white",
          high = y_axis_colorbar[2] 
        ) +
        ggnewscale::new_scale_fill() +
        geom_xsidetile(data = y_range_df, aes(x = .data[[j]], y = 0, fill = .data[[j]]), 
                       #alpha = alpha_switch_yside,
                       show.legend = FALSE,
                       inherit.aes = FALSE) +
        scale_fill_gradient2(
          low = x_axis_colorbar[1],
          mid = "white",
          high = x_axis_colorbar[2]
        ) +
        theme_ggside_void() +
        ggside(x.pos = "bottom", y.pos = "left") +
        theme(ggside.panel.scale = 0.04)
      
      
      plots[[count]] <- p
      count <- count + 1
    }
 # }
  
  age_corr_data  <- gradient_data %>% inner_join(ests %>% filter(grepl("age", term)), by = "region") |> 
    select(all_of(colnames(gradient_data)), indirect, outcome)
  wide_data <- age_corr_data |> pivot_wider(names_from = "outcome", values_from = "indirect") 
  
  # r <- cor(wide_data$memory, wide_data$executive, method = "pearson")
  # p_val <- perm_sphere_p(wide_data$memory, wide_data$executive, perm.id = perms, corr.type='pearson')
  # 
  # # p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
  # # label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
  # p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("<i>p<sub>spin</sub></i> &lt; ", "<i>p<sub>spin</sub></i> = ", "<i>p<sub>spin</sub></i> &gt; "))(p_val)
  # label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
  
  r <- cor(wide_data$memory, wide_data$executive, method = "pearson")
  p_val <- perm_sphere_p(wide_data$memory, wide_data$executive, perm.id = perms, corr.type='pearson')
  p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
  label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
  
  
  age_corr <- wide_data |> 
    ggplot(aes(memory, executive, color = name)) +
    geom_point(alpha = 0.2, show.legend = FALSE, size = 0.2) +
    stat_poly_line(se = FALSE, color = "#323232", linewidth = 0.3) +
    ggpp::annotate(geom = "text_npc", label = label_corr,
                   size = 2,
                   npcx = "right", npcy = "top",
                   family = "sans",
                   parse = TRUE) +
    scale_color_manual(values = net_names %>% select(name, col) %>% deframe())  +
    labs(#subtitle = label_corr,
         x = expression(Age %->% FCS %->% Mem),
         y = expression(Age %->% FCS %->% Exec)) +
    theme(axis.title = element_text(size = rel(0.85)))
  
  
  
  pathology_corr_data  <- gradient_data %>% inner_join(ests %>% filter(grepl("pathology", term)), by = "region") |> 
    select(all_of(colnames(gradient_data)), indirect, outcome)
  wide_data <- pathology_corr_data |> pivot_wider(names_from = "outcome", values_from = "indirect") 
  
  # r <- cor(wide_data$memory, wide_data$executive, method = "pearson")
  # p_val <- perm_sphere_p(wide_data$memory, wide_data$executive, perm.id = perms, corr.type='pearson')
  # 
  # p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("<i>p<sub>spin</sub></i> &lt; ", "<i>p<sub>spin</sub></i> = ", "<i>p<sub>spin</sub></i> &gt; "))(p_val)
  # label_corr <- paste0("<i>r</i> = ", r |> round(2), " ,&#8203;&nbsp;" , p_lab)
  
  r <- cor(wide_data$memory, wide_data$executive, method = "pearson")
  p_val <- perm_sphere_p(wide_data$memory, wide_data$executive, perm.id = perms, corr.type='pearson')
  p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
  label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
  
  pathology_corr <- wide_data |> 
    ggplot(aes(memory, executive, color = name)) +
    geom_point(alpha = 0.2, show.legend = FALSE, size = 0.2) +
    stat_poly_line(se = FALSE, color = "#323232", linewidth = 0.3) +
    ggpp::annotate(geom = "text_npc", label = label_corr,
                   size = 2,
                   npcx = "right", npcy = "top",
                   family = "sans",
                   parse = TRUE) +
    scale_color_manual(values = net_names %>% select(name, col) %>% deframe())  +
    labs(#subtitle = label_corr,
         x = expression(ADpath %->% FCS %->% Mem),
         y = expression(ADpath %->% FCS %->% Exec)) +
    theme(axis.title = element_text(size = rel(0.85)))
  
  plots_to_include <- c(est_plots, plots, list(age_corr, pathology_corr))

  lo <- "ABCD
        EFGH
        IIJJ"
  
  big_plot <- Reduce(`+`, plots_to_include) + plot_annotation(title = plt_title, 
                                                            #subtitle = subtit_expr,
                                                            #tag_levels = c('A', '1'),
                                                            theme = theme(
                                                              plot.title = element_text(size = plt_title_size,
                                                                hjust = 0),
                                                              plot.background = element_rect(color = "black")
                                                            )) +
    plot_layout(design =  lo, axis_titles = "collect", 
                axes = "collect",
                guides = "collect",
                # widths = col_widths,
                # heights = row_heights
    ) 
  
  if (rectangle) {
    big_plot  <-  big_plot  + plot_annotation(
      theme = theme(plot.background = element_rect(color = "black", fill = NA, linewidth = 0.5)))
  }
  
  source_data <- gradient_data %>%
    inner_join(ests, by = "region")
  
  list(
    plot = big_plot,
    data = source_data
  )
}

fast_mediation_all <- function(parcel_names, data, treat, outcome, covars = NULL) {
  data <- data |> drop_na(all_of(c(treat, outcome, covars )))
  # --- design pieces ---
  Z   <- model.matrix(reformulate(c(covars)), data)      
  tvec <- data[[treat]]
  y    <- data[[outcome]]
  Y    <- as.matrix(data[parcel_names])
  
  # drop any rows with NA in used columns
  keep <- stats::complete.cases(Z, tvec, y, Y)
  Z    <- Z[keep, , drop = FALSE]
  tvec <- tvec[keep]
  y    <- y[keep]
  Y    <- Y[keep, , drop = FALSE]
  
  n   <- nrow(Z); pZ <- ncol(Z)
  
  # --- (a) paths: FC ~ Z + treat (multivariate) ---
  X_m <- cbind(Z, tvec); colnames(X_m)[ncol(X_m)] <- treat
  fit_m   <- lm.fit(X_m, Y)
  beta_m  <- fit_m$coefficients                         
  rownames(beta_m) <- colnames(X_m)
  a_all   <- beta_m[treat, ]
  
  # SE(a): sigma^2_j * (X'X)^-1[treat,treat]
  resid_m   <- fit_m$residuals
  sigma2_m  <- colSums(resid_m^2) / (n - (pZ + 1))
  XtX_inv_m <- chol2inv(chol(crossprod(X_m)))
  dimnames(XtX_inv_m) <- list(colnames(X_m), colnames(X_m))
  se_a_all  <- sqrt(sigma2_m * XtX_inv_m[treat, treat])
  
  # --- (b) and c' via FWL with Z-only residualization ---
  r_y  <- lm.fit(Z, y)$residuals
  r_t  <- lm.fit(Z, tvec)$residuals
  R_FC <- lm.fit(Z, Y)$residuals                       # n x P
  
  Syy <- sum(r_y^2)
  Syt <- sum(r_y * r_t)
  Stt <- sum(r_t^2)
  Syf <- colSums(r_y * R_FC)
  Stf <- colSums(r_t * R_FC)
  Sff <- colSums(R_FC^2)
  den <- Stt * Sff - Stf^2
  
  # coefficients in r_y ~ r_t + R_FC[,k]
  b_all      <- (Syf * Stt - Syt * Stf) / den          # FC -> cog | treat, Z
  direct_all <- (Syt * Sff - Syf * Stf) / den          # treat -> cog | FC, Z
  
  # SEs for b and c'
  df2   <- n - pZ - 2
  SSE   <- Syy - (direct_all * Syt + b_all * Syf)
  sigma2 <- SSE / df2
  se_b_all      <- sqrt(sigma2 * Stt / den)
  se_direct_all <- sqrt(sigma2 * Sff / den)
  
  # --- total (c): y ~ Z + treat ---
  X_t <- cbind(Z, tvec); colnames(X_t)[ncol(X_t)] <- treat
  total <- lm.fit(X_t, y)$coefficients[treat]
  
  # --- indirect + Sobel SE ---
  indirect <- a_all * b_all
  se_ind   <- sqrt(b_all^2 * se_a_all^2 + a_all^2 * se_b_all^2)
  z        <- indirect / se_ind
  p        <- 2 * pnorm(-abs(z))
  
  tibble::tibble(
    parcel   = parcel_names,
    a        = a_all,
    b        = b_all,
    indirect = indirect,
    p_indirect = p,
    direct   = direct_all,
    total    = total
  )
}


# 
# plot_mediation_structure <- function(treat, mediator, outcome, lab_size = 11, padding = 0.4) {
#   
#   treat <- dplyr::case_when(
#     treat == "pathology_ad" ~ "ADPath",
#     treat == "pathology_ad_" ~ "ADPath",
#     treat == "age" ~ "Age",
#     treat == "age_" ~ "Age",
#     TRUE ~ treat
#   )
#   
#   outcome <- dplyr::case_when(
#     outcome == "memory" ~ "Mem",
#     outcome == "executive" ~ "Exec",
#     TRUE ~ outcome
#   )
#   
#   require(ggforce)
#   require(dplyr)
#   
#   # Node positions
#   nodes <- data.frame(
#     name = c(treat, mediator, outcome),
#     x = c(0, 1, 2),
#     y = c(0, 1, 0)
#   )
#   
#   edges <- data.frame(
#     from = c(treat, mediator, treat),
#     to   = c(mediator, outcome, outcome),
#     label = c("a", "b", "c′ (direct)\n c (total)")
#   )
#   
# 
#   
#   edges2 <- edges %>%
#     rowwise() %>%
#     mutate(
#       x_from = nodes$x[nodes$name == from],
#       y_from = nodes$y[nodes$name == from],
#       x_to   = nodes$x[nodes$name == to],
#       y_to   = nodes$y[nodes$name == to],
#       dx = x_to - x_from,
#       dy = y_to - y_from,
#       dist = sqrt(dx*dx + dy*dy),
#       x_end = x_to - (dx/dist) * padding,
#       y_end = y_to - (dy/dist) * padding
#     )
#   
#   ggplot() +
#     geom_curve(
#       data = edges2,
#       aes(
#         x = x_from, y = y_from,
#         xend = x_end, yend = y_end
#       ),
#       curvature = 0,
#       arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
#       color = "black"
#     ) +
#     geom_label(
#       aes(x = x, y = y, label = name),
#       data = nodes,
#       size = lab_size
#     ) +
#     coord_cartesian(
#       xlim = c(-0.3, 2.3),
#       ylim = c(-0.4, 1.4),
#       clip = "off"
#     ) +
#     theme_void()
# }

plot_mediation_structure <- function(treat, mediator, outcome,
                                            lab_size = 11) {
  
  library(ggraph)
  library(tidygraph)
  library(dplyr)
  
  # --- Rename treat/outcome for readability --------------------------------
  treat <- dplyr::case_when(
    treat %in% c("pathology_ad", "pathology_ad_") ~ "ADpath",
    treat %in% c("age", "age_") ~ "Age",
    TRUE ~ treat
  )
  
  outcome <- dplyr::case_when(
    outcome == "memory" ~ "Mem",
    outcome == "executive" ~ "Exec",
    TRUE ~ outcome
  )
  
  # --- Build node table with manual layout ---------------------------------
  nodes <- tibble(
    name = c(treat, mediator, outcome),
    x = c(0, 1, 2),
    y = c(0, 1, 0)
  )
  
  # --- Edge definitions -----------------------------------------------------
  edges <- tibble(
    from = c(treat, mediator, treat),
    to   = c(mediator, outcome, outcome),
    label = c("a", "b", "c′ (direct)\n c (total)")
  )
  
  # --- Build graph object ---------------------------------------------------
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  # --- Plot with ggraph -----------------------------------------------------
  ggraph(graph, layout = "manual", x = x, y = y) +
    geom_node_text(
      aes(label = name),
      size = lab_size
    ) +
    # Arrows
    geom_edge_link(
      aes(
        end_cap = label_rect(node2.name, fontsize = lab_size*(1/0.5)),
        start_cap = label_rect(node1.name, fontsize = lab_size*(1/0.5))
      ),
      linewidth = 0.1,
      angle_calc = "along",
      label_dodge = unit(0.25, "mm"),
      arrow = arrow(angle = 20, length = unit(0.04, "cm"), type = "closed"),
      color = "black"
    ) +
    
    
    coord_cartesian(
      xlim = c(-0.3, 2.3),
      ylim = c(-0.4, 1.4),
      clip = "off"
    ) +
    theme_void(base_size = lab_size)
}
