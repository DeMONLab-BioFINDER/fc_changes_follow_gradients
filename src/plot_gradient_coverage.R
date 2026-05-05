plot_gradient_coverage <- function(
  grad_df,
  nq_results = readr::read_rds("stuff_for_revisions/schaefer1000_NQ_results.rds"),
  roi_labs = readr::read_delim(
    "stuff_for_revisions/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv",
    show_col_types = FALSE
  ) |>
    dplyr::select(1, 2) |>
    dplyr::rename(region = `ROI Name`, parcel = `ROI Label`),
  gradient_study = "biofinder",
  r_filter = 0.05,
  nq_filter_terms = NULL,
  rasterize = FALSE,
  raster_width_px = 1500,
  raster_dpi = 300,
  base_size = 7,
  wordcloud_max_size = 4.5,
  wordcloud_grid_margin = 4,
  network_legend_text_scale = 0.8,
  network_legend_key_lines = 0.7,
  gradient_cols = data.frame(
    gradient1 = c("#3F596D", "#D38A4E"),
    gradient2 = c("#4682B4", "#781286"),
    gradient3 = c("#8A6081", "#738518")
  ),
  net_names = data.frame(
    name = c("Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default"),
    col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"),
    label = 1:7
  )
) {
  make_gradient_plots_local <- function(
    gradient_data,
    gradient_colors,
    atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds"),
    grad_char = c("gradient1", "gradient2", "gradient3"),
    add_shade = FALSE,
    shade_size = 0.1,
    shade_alpha = 0.3,
    rasterize = FALSE,
    raster_width_px = 1500,
    raster_dpi = 300,
    base_size_ = 10
  ) {
    make_raster_overlay <- function(plot, bbox,
                                    width_px = 3000,
                                    height_px = 3000,
                                    dpi = 300,
                                    bg = "transparent") {
      tf <- tempfile(fileext = ".png")
      
      ggplot2::ggsave(
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
      
      ggplot2::annotation_raster(
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
                                        fill_var = "statistic",
                                        fill_scale = ggplot2::scale_fill_gradient2(
                                          low = scales::muted("blue"),
                                          mid = "white",
                                          high = scales::muted("red")
                                        )) {
      bbox <- sf::st_bbox(atlas_sf)
      ratio <- as.numeric((bbox["xmax"] - bbox["xmin"]) / (bbox["ymax"] - bbox["ymin"]))
      height_px <- max(1, round(width_px / ratio))
      
      p_raster <- ggplot2::ggplot() +
        ggplot2::geom_sf(
          data = atlas_sf,
          ggplot2::aes(fill = .data[[fill_var]], geometry = geometry),
          linewidth = 0.1 * (base_size_ / 10),
          show.legend = FALSE
        ) +
        ggplot2::coord_sf(
          xlim = c(bbox["xmin"], bbox["xmax"]),
          ylim = c(bbox["ymin"], bbox["ymax"]),
          expand = FALSE
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
          plot.background = ggplot2::element_rect(fill = "transparent", colour = NA)
        ) +
        fill_scale
      
      if (!is.null(shade_sf)) {
        p_raster <- p_raster +
          ggplot2::geom_sf(
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
    
    gradient_plots <- list()
    
    for (grad in grad_char) {
      plot_data <- gradient_data |>
        dplyr::right_join(atlas_geometry$atlas, by = "region")
      
      bbox <- sf::st_bbox(atlas_geometry$atlas)
      raster_layer <- NULL
      if (rasterize) {
        raster_layer <- make_atlas_raster_layer(
          atlas_sf = plot_data |> sf::st_as_sf(),
          shade_sf = if (add_shade) atlas_geometry$shade else NULL,
          shade_size = shade_size,
          shade_alpha = shade_alpha,
          width_px = raster_width_px,
          dpi = raster_dpi,
          fill_var = grad,
          fill_scale = ggplot2::scale_fill_gradient2(
            low = gradient_colors[[grad]][1],
            mid = "white",
            high = gradient_colors[[grad]][2]
          )
        )
      }
      
      gradient_plots[[grad]] <- ggplot2::ggplot() +
        (
          if (rasterize) {
            raster_layer
          } else {
            ggplot2::geom_sf(
              data = plot_data,
              ggplot2::aes(fill = .data[[grad]], geometry = geometry),
              linewidth = 0.1 * (base_size_ / 10),
              show.legend = FALSE
            )
          }
        ) +
        (
          if (!rasterize && add_shade) {
            ggplot2::geom_sf(
              data = atlas_geometry$shade,
              size = shade_size * (base_size_ / 10),
              alpha = shade_alpha
            )
          } else {
            NULL
          }
        ) +
        ggplot2::coord_sf(
          xlim = c(bbox["xmin"], bbox["xmax"]),
          ylim = c(bbox["ymin"], bbox["ymax"]),
          expand = FALSE
        ) +
        ggplot2::theme_void(base_size = base_size_) +
        ggplot2::labs(
          fill = "",
          title = dplyr::case_when(
            grad == "gradient1" ~ "Sens-Assoc",
            grad == "gradient2" ~ "Vis-Mot",
            grad == "gradient3" ~ "Rep-Exec",
            TRUE ~ grad
          )
        ) +
        ggplot2::theme(
          legend.position = "",
          panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
          plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
          legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
          legend.box.background = ggplot2::element_rect(fill = "transparent", colour = NA),
          plot.title = ggplot2::element_text(
            color = "black",
            hjust = 0.5,
            size = ggplot2::rel(2.2)
          )
        ) +
        ggplot2::scale_fill_gradient2(
          low = gradient_colors[[grad]][1],
          mid = "white",
          high = gradient_colors[[grad]][2]
        )
    }
    
    gradient_plots
  }
  
  if (is.null(nq_filter_terms)) {
    url <- "http://www.cognitiveatlas.org/api/search"
    resp <- httr2::request(url) |> httr2::req_perform()
    ca <- httr2::resp_body_json(resp)
    ca <- tolower(vapply(ca, function(x) x$name, character(1)))
    
    nq_terms <- nq_results$feature |>
      unique() |>
      stringr::str_remove_all("neuroquery6308_combined_tfidf__")
    
    idx <- stringdist::amatch(nq_terms, ca, method = "jw", maxDist = 0.1)
    nq_match <- data.frame(
      nq_term = nq_terms,
      ca_match = ifelse(is.na(idx), NA, ca[idx])
    ) |>
      dplyr::filter(!is.na(ca_match)) |>
      dplyr::mutate(
        dist = stringdist::stringdist(tolower(nq_term), ca_match, method = "jw")
      )
    
    nq_filter_terms <- nq_match$nq_term
  }
  
  grad_study_df <- grad_df |>
    dplyr::filter(study == gradient_study)
  
  scale_factor <- base_size / 10
  
  nq_g1 <- nq_results |>
    dplyr::arrange(parcel, dplyr::desc(r)) |>
    dplyr::mutate(feature = stringr::str_remove_all(feature, "neuroquery6308_combined_tfidf__")) |>
    dplyr::filter(r > r_filter) |>
    dplyr::filter(feature %in% nq_filter_terms) |>
    dplyr::left_join(roi_labs, by = "parcel") |>
    dplyr::left_join(dplyr::select(grad_study_df, dplyr::starts_with("grad"), region), by = "region") |>
    dplyr::ungroup() |>
    dplyr::group_by(feature) |>
    dplyr::summarise(
      mentions = dplyr::n(),
      mean_g = mean(gradient1),
      mean_r = mean(r),
      .groups = "drop"
    ) |>
    dplyr::filter(mentions > 15)
  
  nq_g3 <- nq_results |>
    dplyr::arrange(parcel, dplyr::desc(r)) |>
    dplyr::mutate(feature = stringr::str_remove_all(feature, "neuroquery6308_combined_tfidf__")) |>
    dplyr::filter(r > r_filter) |>
    dplyr::filter(feature %in% nq_filter_terms) |>
    dplyr::left_join(roi_labs, by = "parcel") |>
    dplyr::left_join(dplyr::select(grad_study_df, dplyr::starts_with("grad"), region), by = "region") |>
    dplyr::ungroup() |>
    dplyr::group_by(feature) |>
    dplyr::summarise(
      mentions = dplyr::n(),
      mean_g = mean(gradient3),
      mean_r = mean(r),
      .groups = "drop"
    ) |>
    dplyr::filter(mentions > 15)
  
  g1_wc_data <- nq_g1 |>
    dplyr::mutate(
      grad_side = ifelse(mean_g > 0, "Association", "Sensory"),
      abs_g = abs(mean_g),
      mean_g = tanh(mean_g / 20)
    ) |>
    dplyr::group_by(grad_side) |>
    dplyr::arrange(mentions, .by_group = TRUE) |>
    dplyr::slice_head(n = 25)
  
  g1_wc <- g1_wc_data |>
    ggplot2::ggplot() +
    ggwordcloud::geom_text_wordcloud(
      ggplot2::aes(label = feature, size = abs(mean_g), color = mean_g, x = mean_g),
      grid_margin = wordcloud_grid_margin * scale_factor,
      seed = 6
    ) +
    ggplot2::scale_color_gradient2(
      low = gradient_cols[1, 1],
      mid = "white",
      high = gradient_cols[2, 1]
    ) +
    ggplot2::xlim(-1, 1.8) +
    ggplot2::scale_size_area(max_size = wordcloud_max_size * scale_factor) +
    ggplot2::theme_void(base_size = 24 * scale_factor)
  
  g3_wc_data <- nq_g3 |>
    dplyr::mutate(
      grad_side = ifelse(mean_g > 0, "Execution", "Representation"),
      mean_g = tanh(mean_g / 20)
    ) |>
    dplyr::group_by(grad_side) |>
    dplyr::arrange(mentions, .by_group = TRUE) |>
    dplyr::slice_head(n = 25)
  
  g3_wc <- g3_wc_data |>
    ggplot2::ggplot() +
    ggwordcloud::geom_text_wordcloud(
      ggplot2::aes(label = feature, size = abs(mean_g), color = mean_g, x = mean_g),
      grid_margin = wordcloud_grid_margin * scale_factor,
      seed = 6
    ) +
    ggplot2::scale_color_gradient2(
      low = gradient_cols[1, 3],
      mid = "white",
      high = gradient_cols[2, 3]
    ) +
    ggplot2::scale_size_area(max_size = wordcloud_max_size * scale_factor) +
    ggplot2::xlim(-1, 1.8) +
    ggplot2::theme_void(base_size = 24 * scale_factor)
  
  grad_plots <- make_gradient_plots_local(
    grad_study_df,
    gradient_colors = gradient_cols,
    add_shade = TRUE,
    shade_alpha = 0.01,
    shade_size = 0.05 * scale_factor,
    rasterize = rasterize,
    raster_width_px = raster_width_px,
    raster_dpi = raster_dpi,
    base_size_ = base_size
  )
  
  net_sep_data <- grad_study_df |>
    tidyr::pivot_longer(dplyr::starts_with("grad"), names_to = "gradient", values_to = "value") |>
    dplyr::filter(gradient != "gradient2") |>
    dplyr::mutate(
      gradient = dplyr::case_when(
        gradient == "gradient1" ~ "Sensory-Association",
        gradient == "gradient3" ~ "Representational-Executive"
      ) |> as_factor()
    )
  
  make_net_sep_plot <- function(plot_df, grad_key, grad_label, show_legend = TRUE) {
    x_span <- range(plot_df$value)#max(abs(plot_df$value), na.rm = TRUE)
    xside_data <- tibble::tibble(
      value = seq(x_span[[1]], x_span[[2]], length.out = 200),
      strip_y = 1
    )
    
    ggplot2::ggplot(plot_df, ggplot2::aes(x = value, y = name, fill = name)) +
      ggplot2::geom_boxplot(linewidth = 0.3 * scale_factor, outliers = FALSE) +
      ggplot2::geom_violin(alpha = 0.3, linewidth = 0.5 * scale_factor) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3 * scale_factor) +
      ggplot2::facet_wrap(~gradient, scales = "free_x") +
      ggrastr::rasterise(
      ggside::geom_xsidetile(
        data = xside_data,
        ggplot2::aes(x = value, y = strip_y, xfill = value),
        inherit.aes = FALSE,
        show.legend = FALSE,
        height = 1
      )) +
      ggplot2::scale_fill_manual(values = net_names |>
        dplyr::select(name, col) |>
        tibble::deframe()) +
      ggside::scale_xfill_gradientn(
        colours = c(
          gradient_cols[1, grad_key],
          "white",
          gradient_cols[2, grad_key]
        ),
        values = scales::rescale(c(x_span[[1]], 0, x_span[[2]])),
        guide = "none"
      ) +
      ggplot2::labs(
        y = "",
        x = "Placement of Yeo Networks on Gradients",
        fill = NULL
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          nrow = 1,
          byrow = TRUE,
          title = NULL
        )
      ) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        legend.position = if (show_legend) "top" else "none",
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.text = ggplot2::element_text(
          hjust = 0.75,
          size = ggplot2::rel(network_legend_text_scale)
        ),
        legend.key.size = grid::unit(network_legend_key_lines * scale_factor, "lines"),
        legend.margin = ggplot2::margin(t = 12 * scale_factor),
        #strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = ggplot2::rel(1)),
        ggside.panel.scale = 0.06
      ) +
      ggside::ggside(x.pos = "bottom")
  }
  
  net_sep_g1 <- make_net_sep_plot(
    plot_df = dplyr::filter(net_sep_data, gradient == "Sensory-Association"),
    grad_key = "gradient1",
    grad_label = "Sensory-Association",
    show_legend = TRUE
  )
  
  net_sep_g3 <- make_net_sep_plot(
    plot_df = dplyr::filter(net_sep_data, gradient == "Representational-Executive"),
    grad_key = "gradient3",
    grad_label = "Representational-Executive",
    show_legend = FALSE
  )
  
  net_sep <- net_sep_g1 + net_sep_g3 +
    patchwork::plot_layout(ncol = 2, guides = "collect", axis_titles = "collect") &
    ggplot2::theme(legend.position = "top")
  
  g1 <- grad_plots[["gradient1"]] +
    ggplot2::ggtitle("Sensory-Association Axis") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 22 * scale_factor))
  
  g3 <- grad_plots[["gradient3"]] +
    ggplot2::ggtitle("Representational-Executive Axis") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 22 * scale_factor))
  
  test <- cowplot::ggdraw() +
    cowplot::draw_plot(g1, width = 0.5, height = 0.38, x = 0, y = 0.62) +
    cowplot::draw_plot(g3, width = 0.5, height = 0.38, x = 0.5, y = 0.62) +
    cowplot::draw_plot(g1_wc, width = 0.5, height = 0.6, x = 0.05, y = 0.17) +
    cowplot::draw_plot(g3_wc, width = 0.5, height = 0.6, x = 0.525, y = 0.19) +
    cowplot::draw_plot(net_sep, width = 1, height = 0.35, x = 0, y = 0) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(colour = "black", linewidth = 1 * scale_factor)
    )
  
  source_data <- list(
    gradients = grad_study_df,
    nq_g1 = nq_g1,
    nq_g3 = nq_g3,
    g1_wordcloud = g1_wc_data,
    g3_wordcloud = g3_wc_data,
    network_separation = net_sep_data
  )
  
  list(
    plot = test,
    data = source_data,
    nq_filter_terms = nq_filter_terms
  )
}
