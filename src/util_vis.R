plot_gradient_relationships <- function(subject_data,
                                        gradient_data,
                                        atlas_geometry = readRDS("data/atlas_data/Schaefer2018_1000Parcels_geometry.rds"),
                                        vect = FALSE,
                                        base_size_ = 11,
                                        gradients = 1:3,
                                        gradient_colors = NULL,
                                        empty_row_height = 0,
                                        padding = 0,
                                        r2_size = rel(4.6),
                                        spintest = TRUE,
                                        perms = readRDS("data/atlas_data/permutations_1000_hungarian.rds"),
                                        list_of_parcel_data, 
                                        id_var = "image_file",
                                        mod_formula = formula(paste0("FC ~ age")),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        filter_criteria = quo(),
                                        plt_title = "",
                                        tag_sep = "",
                                        tag_prefix = "",
                                        layout_construction = "horizontal",
                                        right_term_side = FALSE,
                                        include_gradient_plots = TRUE,
                                        side_color_bar = TRUE,
                                        plot_spacing = 0.3,
                                        show_networks = FALSE,
                                        network_geometry = NULL,
                                        cache_runs = TRUE,
                                        longitudinal = FALSE,
                                        logistic_fit = FALSE,
                                        scale_fc = FALSE,
                                        sub_id = "sid",
                                        longitudinal_formula = formula(paste0("FC ~ age + (1 | ", sub_id, ")"))
                                        ){
  
  source("src/util.R")
  require(tidyverse)
  require(scales)
  require(patchwork)
  require(ggpmisc)
  require(sf)
  
  
  old <- theme_set(theme_bw(base_size = base_size_))
  theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA),
               legend.background = element_rect(fill = "transparent", colour = NA),
               legend.box.background = element_rect(fill = "transparent", colour = NA))

  net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                          col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), #"#DCF8A4"
                          label = c(1:7))
  
  grad_char <- paste0("gradient", gradients)
  analysis_name <- names(list_of_parcel_data)
  
  if (is.null(gradient_colors)) {
    gradient_colors <- matrix(rep(c("#3A3A98", "#832424"), length(grad_char) ), ncol =length(grad_char))
    colnames(gradient_colors) <- grad_char
    gradient_colors <- as.data.frame(gradient_colors)
  } 
  
  
  if (!longitudinal) {
  
    if (cache_runs) {
      prev_mod_formula <- function() {
        tryCatch(
          {
            list(suppressWarnings(readRDS("analysis_cache/model_formula.rds")), 
                 readRDS("analysis_cache/analysis_name.rds"),
                 readRDS("analysis_cache/filter_crit.rds"))
          },
          error = function(cond) {
            message(paste("No cache or file exists"))
            message("Here's the original error message:")
            message(conditionMessage(cond))
            # Choose a return value in case of error
            NA
          })
      }
      prev_values <- prev_mod_formula()
      dir.create(file.path("analysis_cache"), showWarnings = FALSE)
      write_rds(mod_formula, "analysis_cache/model_formula.rds")
      write_rds(analysis_name, "analysis_cache/analysis_name.rds")
      write_rds(as_label(filter_criteria), "analysis_cache/filter_crit.rds")
      
      if (
        identical(mod_formula, prev_values[[1]]) & identical(analysis_name, prev_values[[2]]) & identical(as_label(filter_criteria), prev_values[[3]])
      ){
        list_of_fits <- read_rds("analysis_cache/list_of_fits.rds")
      } else {
        
        list_of_fits <- list()
        for (analysis in analysis_name) {
          print(paste0("Running linear models for ", analysis))
          list_of_fits[[analysis]] <- 
            nodal_regression_fits(
              subject_data %>% filter(!!filter_criteria), 
              list_of_parcel_data[[analysis]], 
              vectorised = vect,
              roi_names = rois,
              id_var = id_var,
              logistic = logistic_fit,
              scale_fc = scale_fc,
              model_formula = mod_formula)
        }
        write_rds(list_of_fits, "analysis_cache/list_of_fits.rds")
      }
      
    } else {
      list_of_fits <- list()
      for (analysis in analysis_name) {
        print(paste0("Running linear models for ", analysis))
        list_of_fits[[analysis]] <- 
          nodal_regression_fits(
            subject_data %>% filter(!!filter_criteria), 
            list_of_parcel_data[[analysis]], 
            vectorised = vect,
            roi_names = rois,
            id_var = id_var,
            model_formula = mod_formula)
      }
    }
    
    list_of_ests <- list()
    print("Getting model estimates")
    for (analysis in names(list_of_fits)){
      list_of_ests[[analysis]] <- get_nodal_ests(list_of_fits[[analysis]], vectorised = vect, mc = TRUE) %>% 
        select(term, region, statistic, n, model_formula)
    }
    
  }
  
  if (longitudinal) {
    list_of_ests <- list()
    print("fitting longitudinal")
    for (analysis in analysis_name){
      print(analysis)
      list_of_ests[[analysis]] <- nodal_lmm_ests(
        subject_data %>% filter(!!filter_criteria),
        list_of_parcel_data[[analysis]], 
        roi_names = rois,
        id_var = id_var,
        subject_id = sub_id,
        model_formula = longitudinal_formula
      )
    }
  }

  
  ests <- list_of_ests[[1]] %>% rename_with(~ names(list_of_ests[1]), statistic)
  if (length(list_of_ests) > 1) {
    for(i in 2:length(list_of_ests)) {
      ests <- inner_join(ests, list_of_ests[[i]] %>% rename_with(~ names(list_of_ests[i]), statistic), by = c("term", "region"))
    }
  }
  
  plot_brain_ests <- function(ests, tag = "a") {
    parcel_line_size = 0.1
    #bsize = 16
    
    terms_of_interest <- unique(ests$term)[!(unique(ests$term) %in% c(covariates, "(Intercept)"))]
    plots_of_terms <- list()
    tag_labs <- paste0(tag_prefix, tag_sep, tag, tag_sep, tolower(as.roman(1:length(terms_of_interest))))
    i = 1
    for (term_of_i in terms_of_interest) {
      
      p <-  ests %>%
        filter(term == term_of_i) %>%
        inner_join(atlas_geometry, by = "region") %>%
        ggplot() +
        geom_sf(aes(
          fill = statistic,
          geometry = geometry), linewidth= 0.1,
          show.legend = FALSE)+
        theme_void(base_size = base_size_)+
        labs(fill = 't', title = str_to_title(str_replace(term_of_i, "_", " ")),
             tag = tag_labs[i]
             ) +
        theme(#legend.position = "",
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(color = "black", hjust = 0.5)
        ) +
        #guides(color = guide_legend(override.aes = list(size = 1))) +
        scale_fill_gradient2(
          low = muted("blue"),
          mid = "white",
          high = muted("red") 
        )
      
      i = i + 1
      
      if (show_networks) {
        p <- p + 
          geom_sf(data = network_geometry %>% drop_na(),
                  aes(
                    #fill = region,
                    color = name,
                    geometry = geometry
                  ), alpha = 0, linewidth = 0.5,
                  show.legend = FALSE) +
          scale_color_manual(
            values = net_names %>% select(name, col) %>% deframe()
          )
      }
      plots_of_terms[[term_of_i]]<- p
    }
    
    return(plots_of_terms)
  }
  
  ## Gradient plots 
  gradient_plots <- list()
  tag_labs <- paste0(tag_prefix, tag_sep, "a",tag_sep, tolower(as.roman(1:length(gradients))))
  i = 1
  for (stud in unique(gradient_data$study)) {
    for (grad in grad_char) {
      
      if (include_gradient_plots) {
        gradient_plots[[paste0(stud,"_", grad)]] <- gradient_data %>% filter(study==stud) %>% 
          #mutate(segregation = ifelse(segregation<0, 0, segregation)) %>% 
          inner_join(atlas_geometry, by = "region") %>%
          ggplot() +
          geom_sf(aes(
            fill = .data[[grad]],
            geometry = geometry), linewidth= 0.1,
            show.legend = FALSE)+
          theme_void(base_size = base_size_)+
          labs(fill = "", title = str_to_title(str_replace(paste0("_", grad), "_", " ")),
               tag = tag_labs[i]
          ) +
          theme(legend.position = "",
                panel.background = element_rect(fill = "transparent", colour = NA),
                plot.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                plot.title = element_text(color = "black", hjust = 0.5)
          ) +
          #guides(color = guide_legend(override.aes = list(size = 1))) +
          scale_fill_gradient2(
            low = gradient_colors[[grad]][1],
            mid = "white",
            high = gradient_colors[[grad]][2] 
          ) 
        i = i +1
      } else {
        gradient_plots[[paste0(stud,"_", grad)]] <- gradient_data %>% filter(study==stud) %>% 
          #mutate(segregation = ifelse(segregation<0, 0, segregation)) %>% 
          inner_join(atlas_geometry, by = "region") %>%
          ggplot() +
          geom_sf(aes(
            fill = .data[[grad]],
            geometry = geometry), alpha = 0, linewidth= 0.1, color = NA,
            show.legend = FALSE)+
          theme_void(base_size = base_size_)+
          labs(fill = "", title = str_to_title(str_replace(paste0("_", grad), "_", " ")),
               tag = tag_labs[i]
          ) +
          theme(legend.position = "",
                panel.background = element_rect(fill = "transparent", colour = NA),
                plot.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                plot.title = element_text(color = NA, hjust = 0.5)
          ) +
          #guides(color = guide_legend(override.aes = list(size = 1))) +
          scale_fill_gradient2(
            low = gradient_colors[[grad]][1],
            mid = "white",
            high = gradient_colors[[grad]][2] 
          ) 
        i = i +1
      }
      
      
    }
  }
  
  
  
  vars_of_interest <- analysis_name
  n_analysis <- length(analysis_name)
  terms_of_interest <- unique(ests$term)[!(unique(ests$term) %in% c(covariates, "(Intercept)"))]
  n_terms <- length(terms_of_interest)
  n_plts_row <- length(terms_of_interest)*length(vars_of_interest)
  tag_labs <- c(length(analysis_name))
  
  n_gradients <- length(grad_char)
  
  labels <- letters[2:(n_analysis+1)]  
  tag_labs  <- c()  
  
  for (g in 0:(n_gradients - 1)) {
    num_seq <- (g * n_terms + 1):((g + 1) * n_terms) 
    for (label in labels) {
      temp <- paste0(tag_prefix, tag_sep, label, tag_sep, num_seq)  
      tag_labs  <- c(tag_labs , temp) 
    }
  }
  
  plots <- list()
  count <- 1
  for (g in grad_char) {
    for (t in vars_of_interest) {
      for (term_ in terms_of_interest){
        plot_data  <- gradient_data %>% inner_join(ests %>% filter(term == term_), by = "region") 
        
        i <- t
        j <- g
        
        # Swap axes if layout is "vertical"
        if (layout_construction == "vertical") {
          i <- g
          j <- t
        }
        
        lab_grad <- "Gradient score"
        lab_term <- paste(str_to_title(str_replace(t, "_", " ")), "t-value")
        
        if (layout_construction == "horizontal") {
          x_lab <- lab_term
          y_lab <- lab_grad
        } else {
          y_lab <- lab_term
          x_lab <- lab_grad
        }

        
          
        x_min <- min(plot_data[i])
        x_max <- max(plot_data[i])
        
        y_min <- min(plot_data[j])
        y_max <- max(plot_data[j])
        
        
        
        p <- plot_data %>% 
          ggplot(aes(x = .data[[i]], y =.data[[j]],
                     color = name)) +
          geom_point(alpha = 0.2) +
          stat_poly_line(se = FALSE, color = "#323232") +
          labs(#title = term_,
            x = x_lab,
            y = y_lab,
            tag = tag_labs[count],
            color = "Network") +
          xlim(x_min, x_max)+
          #ylim(y_min, y_max)+
          theme_bw(base_size = base_size_) +
          theme(
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            #axis.title = element_text(size = rel(1.2))
            ) +
          scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
          scale_y_continuous(limits = c(y_min, y_max), position = ifelse(right_term_side, "right", "left"))
        
        if (spintest) {
          r <- cor(plot_data[[i]], plot_data[[j]], method = "pearson")
          p_val <- perm_sphere_p(plot_data[[i]], plot_data[[j]], perm.id = perms, corr.type='pearson')
          
          p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
          
          label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
          
          
          p <- p + ggpp::annotate(geom = "text_npc", label = label_corr,
                                    size = r2_size,
                                    npcx = "left", npcy = "top",
                                    family = "sans",
                                    parse = TRUE) 
        } else {
          p <- p + stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                                  str_remove(after_stat(rr.confint.label), "95% CI "),
                                                  sep = "*\" \"*")),
                                parse = TRUE, color = "#323232", label.x = "left", label.y = "top", size = r2_size)
        }
        
        
        
        if (side_color_bar) {
          require(ggside)
          x_range_df <- data.frame(x = seq(x_min, x_max, length.out = 100))
          colnames(x_range_df)[1] <- i
          
          y_range_df <- data.frame(y = seq(y_min, y_max, length.out = 100))
          colnames(y_range_df)[1] <- j
          
          if(layout_construction == "horizontal") {
            alpha_switch_xside <- ifelse(g == tail(grad_char, 1), 1, 0)
            alpha_switch_yside <- ifelse(term_ == terms_of_interest[1], 1, 0)
            x_axis_colorbar <- c(muted("blue"), muted("red"))
            y_axis_colorbar <- gradient_colors[[g]]
          }
          if(layout_construction == "vertical") {
            alpha_switch_xside <- ifelse(term_ == tail(terms_of_interest, 1), 1, 0)
            alpha_switch_yside <- ifelse(g == grad_char[1], 1, 0)
            x_axis_colorbar <- gradient_colors[[g]]
            y_axis_colorbar <- c(muted("blue"), muted("red"))
          }

          if(right_term_side & layout_construction == "vertical") alpha_switch_yside <- ifelse(g == tail(grad_char, 1), 1, 0)
          
          
          
          p <- p +
            geom_xsidetile(data = x_range_df, aes(x = .data[[i]], y = 0, fill = .data[[i]]), 
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
            ggside(x.pos = "bottom", y.pos = ifelse(right_term_side, "right", "left")) +
            theme(ggside.panel.scale = 0.04)
        }
        
        plots[[count]] <- p
        count <- count + 1
      }
    }
  }
  
  n_terms = length(unique(ests$term)[!(unique(ests$term) %in% c(covariates, "(Intercept)"))])
  n_analysis <- length(list_of_parcel_data)
  
  list_of_brain_plots_ests <- list()
  letter_tag <- letters[2:(n_analysis+1)]
  i = 1
  for (analysis in names(list_of_ests)){
    list_of_brain_plots_ests[[analysis]] <- plot_brain_ests(list_of_ests[[analysis]], tag = letter_tag[i])
    i = i + 1
  }
  
  
  if (layout_construction == "horizontal") {
    n_plot_cols <- length(1:((n_terms*n_analysis) + (n_analysis - 1) + 1))
    col_widths <- c(1, rep(c(rep(1, n_terms), plot_spacing), n_analysis))[-(n_plot_cols+1)]
    empty_cols <- seq(1, n_plot_cols, by = n_terms+1)[-1]
    
    n_plot_rows <- n_gradients + 1 + (n_gradients)  
    row_heights <- rep(1, n_plot_rows)  
    empty_rows <- seq(2, n_plot_rows, by = 2)  
    
    row_heights[empty_rows] <- empty_row_height
    
    layout <- c(
      area(3, 1)  # First gradient
    )
    
    if (n_gradients>1) {
      for (g in 4:(length(gradients)+2)){
        layout <- c(layout, area(g+1, 1))
      }
    }
    
    for (col in 2:n_plot_cols){
      if (col %in% empty_cols) {
      } else {
        layout <- c(layout, area(1, col))
      }
    }
    
    for(i in seq(3, (n_gradients*2-1)+2, by = 2)){
      for (j in 2:((n_terms*n_analysis) + (n_analysis - 1) + 1)) {
        if ((j %in% empty_cols)) {
        } else {
          layout <- c(layout, area(i, j))
        }
      }
    }
    
    for(empt in empty_cols){
      layout <- c(layout, area(1, empt , b = n_gradients + 1))
    }
    
    if (padding != 0) {
      for(pad in 1:padding){
        layout <- c(layout, area(1, n_plot_cols + pad, b = n_gradients + 1))
      }
    }
    

  }

  
  if (layout_construction == "vertical") {
    n_plot_rows <- (n_terms * n_analysis) + (n_analysis - 1) + 1
    row_heights <- c(1, rep(c(rep(1, n_terms), plot_spacing), n_analysis))[-(n_plot_rows+1)]
    col_widths <- NULL
    empty_rows <- seq(1, n_plot_rows, by = n_terms + 1)[-1]
    
    if (right_term_side) {
      layout <- c(
        area(1, 1)
      )
      
      for (g in 2:(length(gradients))) {
        layout <- c(layout, area(1, g))
      }
      
      for (row in 2:n_plot_rows) {
        if (!(row %in% empty_rows)) {
          layout <- c(layout, area(row, (length(gradients) + 1)))
        }
      }
      
      for (j in 1:(length(gradients))) {
        for (i in 2:n_plot_rows) {
          if (!(i %in% empty_rows)) {
            layout <- c(layout, area(i, j))
          }
        }
      }
      
    } else {
      layout <- c(
        area(1, 2)
      )
      
      for (g in 3:(length(gradients) + 1)) {
        layout <- c(layout, area(1, g))
      }
      
      for (row in 2:n_plot_rows) {
        if (!(row %in% empty_rows)) {
          layout <- c(layout, area(row, 1))
        }
      }
      
      for (j in 2:(length(gradients) + 1)) {
        for (i in 2:n_plot_rows) {
          if (!(i %in% empty_rows)) {
            layout <- c(layout, area(i, j))
          }
        }
      }
    }
    
  }

  
  filt_char = as_label(filter_criteria )

  
  brain_plots <- unlist(list_of_brain_plots_ests, recursive = FALSE)
  
  plots_to_include <- c(gradient_plots, brain_plots, plots)

  
  p <- Reduce(`+`, plots_to_include) +
    plot_annotation(title = plt_title, #subtitle = sub_title, 
                    #tag_levels = c('A', '1'),
                    theme = theme(
                      plot.title = element_text(#size = 28,
                        hjust = 0.5),
                      plot.subtitle = element_text(#size = 28,
                        hjust = 0.5)
                    )) +
    plot_layout(design = layout, axis_titles = "collect", axes = "collect", guides = "collect",
                widths = col_widths,
                heights = row_heights
                ) & 
    theme(legend.position = "",
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.tag.position  = if (right_term_side) c(0.1, 0.96) else c(0.9, 0.96),
          #plot.tag = element_text(size = 8, hjust = 0, vjust = 0)
          )
  list(plot = p, n = ests %>% pull(n) %>% unique(), model_formula = ests %>% pull(model_formula) %>% unique())
  
}



plot_gams_v1 <- function(gam_predictions, grad_df, 
                         atlas_geometry = readRDS("data/atlas_data/Schaefer2018_1000Parcels_geometry.rds"),
                         gradient_cols = data.frame(gradient1 = c("#3F596D", "#D38A4E"), 
                                                    gradient2 =  c("#4682B4", "#781286"),  
                                                    gradient3 =  c("#8A6081", "#738518")),
                         biofinder_data, scale_fac = 3,
                         spintest = TRUE,
                         perms = readRDS("data/atlas_data/permutations_1000_hungarian.rds"),
                         figure_pat = "paper/figures") {
  
  require(ggside)
  require(ggpmisc)
  
  scale_factor <- scale_fac
  old <- theme_set(theme_bw(base_size = 5*scale_factor))
  theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA),
               legend.background = element_rect(fill = "transparent", colour = NA),
               legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  
  legend_labs <-  c("Ab42/40", "Braak12", "Braak34", "Braak56")
  biof_path_plot <- biofinder_df %>% filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) %>% 
    select(pathology_ad, ab_ratio, 
           starts_with("braak"), starts_with("cho")) %>% 
    mutate(across(where(is.numeric) & !pathology_ad, scale)) %>% 
    pivot_longer(-pathology_ad, names_to = "scaled_pat_measures", values_to = "value") %>% 
    ggplot(aes(pathology_ad, value, color = scaled_pat_measures)) +
    geom_smooth() +
    guides(color = guide_legend(
      ncol=1, byrow=TRUE,
      title.position="top", 
      title.hjust = 0)) +
    labs(x = "Pathology score", y = "Scaled value")+
    ggsci::scale_color_nejm(name = "Pathology", labels = legend_labs) +
    theme(legend.position = "right",
          legend.title.position = "top",
          legend.text = element_text(size = rel(0.6)),
          legend.title = element_text(size = rel(0.8)))  
  
  pathology_plot <- biof_path_plot +
    theme(legend.position = "")  
  
  pat_leg <- ggpubr::get_legend(biof_path_plot)
  pat_leg <- ggpubr::as_ggplot(pat_leg)
  
  
  
  gradient_plots <- list()
  i = 1
  grad_char <- c("gradient1", "gradient2", "gradient3")
    for (grad in grad_char) {
        gradient_plots[[paste0(grad)]] <- grad_df %>% filter(study=="biofinder") %>% 
          #mutate(segregation = ifelse(segregation<0, 0, segregation)) %>% 
          inner_join(atlas_geometry, by = "region") %>%
          ggplot() +
          geom_sf(aes(
            fill = .data[[grad]],
            geometry = geometry), linewidth= 0.1,
            show.legend = FALSE)+
          # geom_sf(data = network_geometry %>% drop_na(),
          #         aes(
          #           #fill = region,
          #           color = name,
          #           geometry = geometry
          #         ), alpha = 0, linewidth = 0.5,
          #         show.legend = FALSE) +
          theme_void()+
          labs(fill = "", title = str_to_title(str_replace(paste0("_", grad), "_", " "))
          ) +
          theme(legend.position = "",
                panel.background = element_rect(fill = "transparent", colour = NA),
                plot.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                plot.title = element_text(hjust = 0.5)
          ) +
          #guides(color = guide_legend(override.aes = list(size = 1))) +
          scale_fill_gradient2(
            low = gradient_cols[[grad]][1],
            mid = "white",
            high = gradient_cols[[grad]][2] 
          ) 
        i = i +1
      }
  
  
  brain_pat <- list()
  for (pat_grp in unique(cut(gam_predictions$pat_derivs$pathology_ad, 4))) {
    brain_pat[[pat_grp]] <- gam_predictions$pat_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(pathology_ad = cut(pathology_ad, 4)) %>% 
      filter(pathology_ad == pat_grp) %>% 
      group_by(region, pathology_ad) %>% 
      summarise(value=mean(value)) %>% 
      inner_join(atlas_geometry, by = "region") %>%
      ggplot() +
      geom_sf(aes(
        fill = value,
        geometry = geometry), linewidth= 0.1,
        show.legend = FALSE)+
      theme_void()+
      theme(legend.position = "",
            plot.margin = unit(c(0, 0, 0, 0), "npc"),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(color = "black", hjust = 0.5, size = rel(2)))+
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red"),
        # limits=c(-0.15, 0.15), oob=squish
      ) 
  }
  
  scatter_pat <- list()
  for (pat_grp in unique(cut(gam_predictions$pat_derivs$pathology_ad, 4))) {
    gam_predictions$pat_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(pathology_ad = cut(pathology_ad, 4)) %>% 
      filter(pathology_ad == pat_grp) %>% 
      group_by(region, pathology_ad) %>% 
      summarise(value=mean(value)) %>% 
      inner_join(grad_df) -> plot_df
    
    x_range_df <- data.frame(x = seq(min(plot_df$gradient1), max(plot_df$gradient1), length.out = 100))
    colnames(x_range_df)[1] <- "gradient1"
    y_range_df <- data.frame(y = seq(min(plot_df$value), max(plot_df$value), length.out = 100))
    colnames(y_range_df)[1] <- "value"
    
    x_axis_colorbar <- gradient_cols[[1]]
    
    scatter_pat[[pat_grp]] <- plot_df %>% 
      ggplot(aes(x = value, y = gradient1,
                 color = name)) +
      geom_point(alpha = 0.1, show.legend = FALSE) +
      stat_poly_line(se = FALSE, color = "#323232") +
      #xlim(-0.2, 0.05) +
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE, angle = 30), position = "bottom") +
      labs(
        #title = "Pathology AD {current_frame}",
        color = "Network",
        x = "FC slopes averaged over pathology quartiles", 
        y = "") +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe())+
      scale_y_continuous(position = "right", limits = c(min(grad_df$gradient1), max(grad_df$gradient1))) +
      theme(axis.title.y = element_blank(),
            axis.text = element_text(size = rel(0.6))) +
      geom_xsidetile(data = y_range_df, aes(x = value, y = 0, fill = value), 
                     show.legend = FALSE, 
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red")
      ) +
      ggnewscale::new_scale_fill() +
      geom_ysidetile(data = x_range_df, aes(y = gradient1, x = 0, fill = gradient1), 
                     alpha = ifelse(pat_grp == unique(cut(gam_predictions$pat_derivs$pathology_ad, 4))[4], 1, 0),
                     show.legend = FALSE,
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = x_axis_colorbar[1],
        mid = "white",
        high = x_axis_colorbar[2] 
      ) +
      theme_ggside_void() +
      ggside(x.pos = "bottom", y.pos = "right") +
      theme(ggside.panel.scale = 0.02)
    
    if (spintest) {
      
      r <- cor(plot_df$gradient1, plot_df$value, method = "pearson")
      p_val <- perm_sphere_p(plot_df$gradient1, plot_df$value, perm.id = perms, corr.type='pearson')
      
      p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
      
      label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
      
      
      scatter_pat[[pat_grp]] <- scatter_pat[[pat_grp]] + ggpp::annotate(geom = "text_npc", label = label_corr,
                              size =  3.6,
                              npcx = "left", npcy = "top",
                              family = "sans",
                              parse = TRUE)
    } else {
      scatter_pat[[pat_grp]] <- scatter_pat[[pat_grp]] +       
        stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                       str_remove(after_stat(rr.confint.label), "95% CI "),
                                       sep = "*\" \"*")),
                     parse = TRUE, color = "#323232", label.x = "left", label.y = "top", size = 3.6) 
    }
    
  }
  
  
  brain_age <- list()
  for (age_grp in unique(cut(gam_predictions$age_derivs$age, 4))) {
    brain_age[[age_grp]] <- gam_predictions$age_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(age = cut(age, 4)) %>% 
      filter(age == age_grp) %>% 
      group_by(region, age) %>% 
      summarise(value=mean(value)) %>% 
      inner_join(atlas_geometry, by = "region") %>%
      ggplot() +
      geom_sf(aes(
        #frame = pathology_ad,
        fill = value,
        geometry = geometry), linewidth= 0.1,
        show.legend = FALSE)+
      theme_void()+
      #labs(title = "Pathology AD {current_frame}")+
      theme(legend.position = "",
            plot.margin = unit(c(0, 0, 0, 0), "npc"),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(color = "black", hjust = 0.5, size = rel(2)))+
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red"),
        # limits=c(-0.15, 0.15), oob=squish
      ) 
  }
  
  scatter_age <- list()
  for (age_grp in unique(cut(gam_predictions$age_derivs$age, 4))) {
    gam_predictions$age_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(age = cut(age, 4)) %>% 
      filter(age == age_grp) %>% 
      group_by(region, age) %>% 
      summarise(value=mean(value)) %>% 
      inner_join(grad_df) -> plot_df
    
    x_range_df <- data.frame(x = seq(min(plot_df$gradient3), max(plot_df$gradient3), length.out = 100))
    colnames(x_range_df)[1] <- "gradient3"
    y_range_df <- data.frame(y = seq(min(plot_df$value), max(plot_df$value), length.out = 100))
    colnames(y_range_df)[1] <- "value"
    
    x_axis_colorbar <- gradient_cols[[3]]
    
    scatter_age[[age_grp]] <-  plot_df %>% 
      ggplot(aes(x = value, y = gradient3,
                 color = name)) +
      geom_point(alpha = 0.1, show.legend = FALSE) +
      stat_poly_line(se = FALSE, color = "#323232") +
      ylim(min(grad_df$gradient3), max(grad_df$gradient3)) +
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE, angle = 30)#, position = "top"
                         ) +
      labs(
        color = "Network",
        x = "FC slopes averaged over age quartiles", 
        y = "") +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
      theme(axis.title.y = element_blank(),
            axis.text = element_text(size = rel(0.6))) +
      geom_xsidetile(data = y_range_df, aes(x = value, y = 0, fill = value), 
                     show.legend = FALSE, 
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red")
      ) +
      ggnewscale::new_scale_fill() +
      geom_ysidetile(data = x_range_df, aes(y = gradient3, x = 0, fill = gradient3), 
                     alpha = ifelse(age_grp == unique(cut(gam_predictions$age_derivs$age, 4))[1], 1, 0),
                     show.legend = FALSE,
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = x_axis_colorbar[1],
        mid = "white",
        high = x_axis_colorbar[2] 
      ) +
      theme_ggside_void() +
      ggside(x.pos = "bottom", y.pos = "left") +
      theme(ggside.panel.scale = 0.02)
    
    
    if (spintest) {
      
      r <- cor(plot_df$gradient3, plot_df$value, method = "pearson")
      p_val <- perm_sphere_p(plot_df$gradient3, plot_df$value, perm.id = perms, corr.type='pearson')
      p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
      label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)


      scatter_age[[age_grp]] <- scatter_age[[age_grp]] + ggpp::annotate(geom = "text_npc", label = label_corr,
                              size = 3.6,
                              npcx = "left", npcy = "top",
                              family = "sans",
                              parse = TRUE)
    } else {
      scatter_age[[age_grp]] <- scatter_age[[age_grp]] + 
        stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                       str_remove(after_stat(rr.confint.label), "95% CI "),
                                       sep = "*\" \"*")),
                     parse = TRUE, color = "#323232", label.x = "left", label.y = "top", size = 3.6) 
    }
    
  }
    
  
  quantile_trajectories <- gam_predictions$pat_pred %>% pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "predicted_fc") %>% 
    inner_join(grad_df %>% select(region, gradient1)) %>% 
    mutate(grad_grp = cut(gradient1, breaks = quantile(gradient1, probs = seq(0, 1, length.out = 21)))) %>% 
    group_by(pathology_ad, grad_grp) %>% 
    summarise(mean_pred_fc = mean(predicted_fc),
              mean_grad_value = mean(gradient1)) %>% 
    ggplot(aes(pathology_ad, mean_pred_fc, group = mean_grad_value, color = mean_grad_value)) +
    geom_line(linewidth = 1) +
    labs(x = "Pathology score", y = "Predicted FC") +
    #guides(colour = guide_colorbar(title.position="top", title.hjust = 0.0))+
    scale_color_gradient2(
      high = gradient_cols[[1]][2], mid = "white", low = gradient_cols[[1]][1]) +
    theme(
      legend.key.height = unit(1, "null"),
      legend.key.width = unit(0.02, "npc"),
      legend.box.spacing = unit(0.0025, "npc"),
      legend.title = element_blank(),
      legend.text = element_text(margin = margin(l = 0.8), size = rel(0.7))
    )
  
  r2_pat <- gam_predictions$pat_derivs %>% pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>%
    inner_join(grad_df) %>%
    group_by(pathology_ad) %>%
    summarise(grad_R2 = cor(value, gradient1)) %>%
    ggplot(aes(pathology_ad, grad_R2)) +
    geom_line() +
    labs(y = "Correlation (r)", x = "Pathology score")
  
  
  
  for(x_int in seq_range(gam_predictions$pat_pred$pathology_ad, 5)){
    quantile_trajectories <- quantile_trajectories +
      geom_vline(xintercept = x_int, linetype = "dashed")
    r2_pat <- r2_pat +
      geom_vline(xintercept = x_int, linetype = "dashed")
    pathology_plot <- pathology_plot +
      geom_vline(xintercept = x_int, linetype = "dashed")
  }
  
  
  quantile_trajectories_term2 <- gam_predictions$age_pred %>% 
    pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "predicted_fc") %>% 
    inner_join(grad_df %>% select(region, gradient3)) %>% 
    mutate(grad_grp = cut(gradient3, breaks = quantile(gradient3, probs = seq(0, 1, length.out = 21)))) %>% 
    group_by(age, grad_grp) %>% 
    summarise(mean_pred_fc = mean(predicted_fc),
              mean_grad_value = mean(gradient3)) %>% 
    ggplot(aes(age, mean_pred_fc, group = mean_grad_value, color = mean_grad_value)) +
    geom_line(linewidth = 1) +
    labs(y = "Predicted FC", x = "Age") +
    scale_color_gradient2(
      high = gradient_cols[[3]][2], mid = "white", low = gradient_cols[[3]][1]) +
    labs(color = "") +
    #guides(colour = guide_colorbar(title.position="top", title.hjust = 0.0))+
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    theme(
      legend.key.height = unit(1, "null"),
      legend.key.width = unit(0.02, "npc"),
      legend.box.spacing = unit(0.005, "npc"),
      legend.title = element_blank(),
      legend.text = element_text(margin = margin(l = 0.8), size = rel(0.7))
    )
  
  r2_term2 <- gam_predictions$age_derivs %>% pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>%
    inner_join(grad_df ) %>%
    group_by(age) %>%
    summarise(grad_R2 = cor(value, gradient3)) %>%
    ggplot(aes(age, grad_R2)) +
    geom_line() +
    labs(y = "Correlation (r)", x = "Age")
  
  
  for(x_int in seq_range(gam_predictions$age_pred$age, 5)){
    quantile_trajectories_term2 <- quantile_trajectories_term2 +
      geom_vline(xintercept = x_int, linetype = "dashed")
    r2_term2 <- r2_term2 +
      geom_vline(xintercept = x_int, linetype = "dashed")
  }
  
  
  
  a <- wrap_plots(c(scatter_pat,
                    list(gradient_plots[[1]]),
                    list(plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer()),
                    brain_pat, 
                    list(plot_spacer())
                    ), 
                  nrow = 3) +
    plot_layout(
      guides = "collect",           
      axis_titles = "collect",      
      axes = "collect",
      heights = c(0.5, -0.15, 0.5)          
    )
  
  a <- ggdraw() +
    draw_plot(a) +
    draw_figure_label("A")
  
  
  b <- 
    wrap_plots(list(pathology_plot + labs(tag = "B.1"), plot_spacer(),
                    r2_pat + labs(tag = "B.2"), r2_term2,
                    quantile_trajectories+ labs(tag = "B.3"), quantile_trajectories_term2),
               ncol = 2, byrow = TRUE) +
    plot_layout(
      axis_titles = "collect", axes = "collect"
      #guides = "collect" 
    )  & 
    theme(plot.tag.position  = c(0.175, 1.02),
          plot.tag = element_text(size = rel(0.6), hjust = 0, vjust = 0),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA)
          #legend.box = "horizontal"
    ) 
  
  b <- ggdraw() +
    draw_plot(b) +
    draw_figure_label("B")
  
  c <- wrap_plots(c(list(plot_spacer()),
                    brain_age, 
                    list(plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer()),
                    list(gradient_plots[[3]]),
                    scatter_age), nrow = 3) +
    plot_layout(
      guides = "collect",           
      axis_titles = "collect",      
      axes = "collect",
      heights = c(0.5, -0.2, 0.5)
    )

  c <- ggdraw() +
    draw_plot(c) +
    draw_figure_label("C")
  
  # 
  ggdraw() +
    draw_plot(a, x = 0, y = 0.70, width = 1, height = 0.30) +
    draw_plot(b, x = 0, y = 0.30, width = 1, height = 0.40) +
    draw_plot(pat_leg, y = 0.585, x = 0.4, width = 0.2, height = 0.1) +
    draw_plot(c, x = 0, y = 0.0, width = 1, height = 0.30) +
    draw_line(
      x = c(0.099, 0.01, 0.01),
      y = c(0.685, 0.725, 0.81), 
      linetype = 2
    ) +  
    draw_line(
      x = c(0.099 + 0.082, 0.01 + 0.19, 0.01 + 0.19),
      y = c(0.685, 0.725, 0.81), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.099 + 0.082*2, 0.01 + 0.19*2, 0.01 + 0.19*2),
      y = c(0.685, 0.725, 0.81), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.099 + 0.082*3, 0.01 + 0.19*3, 0.01 + 0.19*3),
      y = c(0.685, 0.725, 0.81), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.099 + 0.082*4, 0.01 + 0.19*4, 0.01 + 0.19*4),
      y = c(0.685, 0.725, 0.81), 
      linetype = 2
    ) + # NEDRE DEL
    draw_line(
      x = c(0.495 + 0.082, 0.04 + 0.19, 0.04 + 0.19),
      y = c(0.35, 0.285, 0.19), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.495 + 0.082*2, 0.04 + 0.19*2, 0.04 + 0.19*2),
      y = c(0.35, 0.285, 0.19), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.495 + 0.082*3, 0.04 + 0.19*3, 0.04 + 0.19*3),
      y = c(0.35, 0.285, 0.19), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.495 + 0.082*4, 0.04 + 0.19*4, 0.04 + 0.19*4),
      y = c(0.35, 0.285, 0.19), 
      linetype = 2
    )+
    draw_line(
      x = c(0.495 + 0.082*5, 0.04 + 0.19*5, 0.04 + 0.19*5),
      y = c(0.35, 0.285, 0.19), 
      linetype = 2
    )

}


figure_one <- function(subject_data,
                       measures_list, measures_list_replication, gradients_df = grad_df %>% filter(study=="biofinder"),
                       gradients_df_replication = grad_df %>% filter(study=="adni"),
                       tag_size = 21,
                       draw_size = 18,
                       b_size = 11,
                       plot_title_size = 1,
                       axes_title_size = 1,
                       r2_sizing1 = 4.6,
                       r2_sizing2 = 4.6,
                       boxed = FALSE,
                       empt_row_height = 0,
                       selected_gradients = c(1, 3), 
                       split = FALSE) {
  library(cowplot)
  # library(showtext)
  # showtext_opts(dpi = 300)
  
  if(!is.list(measures_list)) stop("measures should be in a list")
  if(length(measures_list)>1) stop("list of measures should only contain a single metric")
  if(is.null(names(measures_list))) stop("measures list should be named")
  
  if(!is.list(measures_list_replication)) stop("measures should be in a list")
  if(length(measures_list_replication)>1) stop("list of measures should only contain a single metric")
  if(is.null(names(measures_list_replication))) stop("measures list should be named")
  
  
  l_marg = 0.05
  #text_size = 28
  
  
  
  bf_p <-  plot_gradient_relationships(subject_data %>% filter(fmri_bl), 
                                       gradient_data = gradients_df, 
                                       gradients = selected_gradients,
                                       gradient_colors = gradient_cols,
                                       list_of_parcel_data = measures_list,
                                       empty_row_height = empt_row_height,
                                       base_size_ = b_size,
                                       vect = TRUE,
                                       mod_formula = formula(paste0(" ~ age + pathology_ad + sex + rsqa__MeanFD")),
                                       covariates = c("sex", "rsqa__MeanFD"),
                                       r2_size = rel(r2_sizing1),
                                       filter_criteria = quo(),
                                       show_networks = FALSE,
                                       tag_prefix = "",
                                       tag_sep = "",
                                       layout_construction = "horizontal",
                                       include_gradient_plots = TRUE,
                                       right_term_side = FALSE,
                                       plt_title = "",
                                       cache_runs = FALSE)
  
  n_cs <- bf_p$n
  p_bf <- bf_p$plot &
    theme(#text = element_text(size = text_size),
      plot.tag = element_blank(),
      plot.title = element_text(size = rel(plot_title_size)),
      axis.title = element_text(size = rel(axes_title_size))) 
  
  p_bf <- p_bf + plot_annotation(title = paste0("Cross-sectional", " (N=", n_cs, ")"),
                                 subtitle = expression(italic(FC[parcel] ~ "~" ~ age + pathology + sex + motion)),
                                 theme = theme(plot.subtitle = element_text(size = rel(0.9),
                                                                            hjust = 0,
                                                                            vjust = -0.05,
                                                                            family = "mono",
                                                                            face = "italic",
                                                                            margin = margin(l = l_marg, unit = "npc")),
                                               plot.title.position = "plot",
                                               plot.title = element_text(size = rel(1), hjust =0, margin = margin(l = l_marg, unit = "npc")))
  )
  
  plt_idx <- 4
  if (length(selected_gradients) < 2) plt_idx <- plt_idx-1
  p_bf[[plt_idx[1]]] <- p_bf[[plt_idx[1]]] + labs(title = "AD Pathology")
  
  p_bf <- ggdraw() + 
    draw_plot(p_bf) + 
    draw_plot_label("A", x = l_marg-l_marg, size = tag_size) + 
    draw_label("BioFINDER", x = (1/3.1/2), y = 0.75, hjust = 0.5, size =  draw_size)
  
  if (boxed) p_bf <- p_bf + theme(plot.background = element_rect(color = "black", linewidth = 1))
  

  adni_p <-  plot_gradient_relationships(adni_df %>% filter(fmri_bl), 
                                         gradient_data = gradients_df_replication, 
                                         gradients = selected_gradients,
                                         gradient_colors = gradient_cols,
                                         list_of_parcel_data = measures_list_replication,
                                         mod_formula = formula(paste0(" ~ age + pathology_ad + sex + rsqa__MeanFD")),
                                         empty_row_height = empt_row_height,
                                         base_size_ = b_size,
                                         vect = TRUE,
                                         r2_size = rel(r2_sizing1),
                                         covariates = c("sex", "rsqa__MeanFD"),
                                         id_var = "id_ses",
                                         filter_criteria = quo(),
                                         show_networks = FALSE,
                                         tag_prefix = "",
                                         layout_construction = "horizontal",
                                         include_gradient_plots = TRUE,
                                         right_term_side = FALSE,
                                         plt_title = "",
                                         cache_runs = FALSE)
  
  
  n_adni <- adni_p$n
  p_a <- adni_p$plot &
    theme(#text = element_text(size = text_size),
      plot.tag = element_blank(),
      plot.title = element_text(size = rel(plot_title_size)),
      axis.title = element_text(size = rel(axes_title_size)))
  
  p_a <- p_a + plot_annotation(title = paste0("External Replication", " (N=", n_adni, ")"), 
                               subtitle = expression(italic(FC[parcel] ~ "~" ~ age + pathology + sex + motion)),
                               theme = theme(plot.subtitle = element_text(size = rel(0.9), 
                                                                          hjust = 0, 
                                                                          vjust = -0.05,
                                                                          family = "mono",
                                                                          face = "italic",
                                                                          margin = margin(l = l_marg, unit = "npc")), 
                                             plot.title.position = "plot",
                                             plot.title = element_text(size = rel(1), hjust =0, margin = margin(l = l_marg, unit = "npc")))
  )
  
  
  plt_idx <- 4
  if (length(selected_gradients) < 2) plt_idx <- plt_idx-1
  p_a[[plt_idx[1]]] <- p_a[[plt_idx[1]]] + labs(title = "AD Pathology")
  
  p_a <- ggdraw() + 
    draw_plot(p_a) + 
    draw_plot_label("B", x = l_marg-l_marg, size = tag_size) + 
    draw_label("ADNI", x = (1/3.1/2), y = 0.75, hjust = 0.5, size =  draw_size)
  
  if (boxed) p_a <- p_a + theme(plot.background = element_rect(color = "black", linewidth = 1))
  
  # overlay <- ggdraw() +
  #   draw_plot(p_bf, x = 0, y = 0, width = 0.35, height = 1) +
  #   draw_plot(p_bf_long, x = 0.27, y = 0, width = 0.35, height = 1) +
  #   draw_plot(p_a, x = 0.64, y = 0, width = 0.35, height = 1)  
  
  
  ######
  # Cognition
  #######
  
  
  health_cog <-  plot_gradient_relationships(subject_data %>% filter(fmri_bl, diagnosis=="Normal" | diagnosis=="SCD", abnorm_ab==0, !apoe4),
                                             gradient_data = gradients_df, 
                                             gradients = selected_gradients,
                                             gradient_colors = gradient_cols,
                                             list_of_parcel_data = measures_list,
                                             empty_row_height = empt_row_height,
                                             base_size_ = b_size,
                                             r2_size = rel(r2_sizing2),
                                             mod_formula = formula(paste0("~ scale(age) * scale(-mPACC_v1) + pathology_ad + sex + rsqa__MeanFD")),
                                             logistic_fit = FALSE,
                                             vect = TRUE,
                                             covariates = c("sex", "rsqa__MeanFD"),
                                             filter_criteria = quo(),
                                             show_networks = FALSE,
                                             tag_prefix = "",
                                             tag_sep = "",
                                             layout_construction = "horizontal",
                                             plot_spacing = 0.2,
                                             include_gradient_plots = TRUE,
                                             right_term_side = FALSE,
                                             plt_title = "",
                                             cache_runs = FALSE)
  
  n_health <- health_cog$n
  health_l_marg <- l_marg - 0.016666
  p_health_cog <- health_cog$plot &
    theme(plot.tag = element_blank(),
          plot.title = element_text(size = rel(plot_title_size)),
          axis.title = element_text(size = rel(axes_title_size))
    )
  
  p_health_cog <- p_health_cog + plot_annotation(title = paste0("Cognitively unimpaired, no APOE e4, Ab-", " (N=", n_health, ")"), 
                                                 subtitle = expression(italic(FC[parcel] ~ "~" ~ age * "×" * "-mPACC" + pathology + sex + motion)),
                                                 theme = theme(plot.subtitle = element_text(size = rel(0.9),
                                                                                            hjust = 0, 
                                                                                            vjust = -0.05, 
                                                                                            family = "mono",
                                                                                            face = "italic",
                                                                                            margin = margin(l = health_l_marg, unit = "npc")), 
                                                               plot.title.position = "plot",
                                                               plot.title = element_text(size = rel(1), hjust =0, margin = margin(l = health_l_marg, unit = "npc")))
  )
  
  plt_idx <- 3:6
  if (length(selected_gradients) < 2) plt_idx <- plt_idx-1
  p_health_cog[[plt_idx[1]]] <- p_health_cog[[plt_idx[1]]] + labs(title = "Age") + theme(plot.title = element_text(vjust = -1.5)) 
  p_health_cog[[plt_idx[2]]] <- p_health_cog[[plt_idx[2]]] + labs(title = "-mPACC", subtitle = "(Inverted cognition)") + theme(plot.subtitle = element_text(hjust = 0.5, size = rel(0.6)))
  p_health_cog[[plt_idx[3]]] <- p_health_cog[[plt_idx[3]]] + labs(title = "AD Pathology") + theme(plot.title = element_text(vjust = -1.5)) 
  p_health_cog[[plt_idx[4]]] <- p_health_cog[[plt_idx[4]]] + labs(title = "-mPACC×Age") + theme(plot.title = element_text(vjust = -1.5)) 
  
  p_health_cog <- ggdraw() + draw_plot(p_health_cog) + 
    draw_plot_label(ifelse(split, "A", "C"), x = health_l_marg-health_l_marg, size = tag_size) + 
    draw_label("BioFINDER", x = (1/5/2), y = 0.75, hjust = 0.5, size =  draw_size)
  
  
  
  clinical_cog <-  plot_gradient_relationships(subject_data %>% filter(fmri_bl, diagnosis=="MCI" | diagnosis=="AD", !is.na(mPACC_v1)) %>% 
                                                 mutate(`-mPACC_v1` = -mPACC_v1), 
                                               gradient_data = gradients_df, 
                                               gradients = selected_gradients,
                                               gradient_colors = gradient_cols,
                                               list_of_parcel_data = measures_list,
                                               empty_row_height = empt_row_height,
                                               base_size_ = b_size,
                                               r2_size = rel(r2_sizing2),
                                               vect = TRUE,
                                               mod_formula = formula(paste0(" ~ age + pathology_ad + `-mPACC_v1` +  sex + rsqa__MeanFD")),
                                               logistic_fit = FALSE,
                                               covariates = c("sex", "rsqa__MeanFD"),
                                               filter_criteria = quo(),
                                               show_networks = FALSE,
                                               tag_prefix = "",
                                               tag_sep = "",
                                               layout_construction = "horizontal",
                                               include_gradient_plots = FALSE,
                                               right_term_side = FALSE,
                                               plt_title = "",
                                               cache_runs = FALSE)
  
  n_clin <- clinical_cog$n
  clin_l_marg = 0.3
  p_clinical_cog <-
    clinical_cog$plot  &
    theme(plot.tag = element_blank(),
          plot.title = element_text(size = rel(plot_title_size)),
          axis.title = element_text(size = rel(axes_title_size))
    )
  
  p_clinical_cog <- p_clinical_cog + plot_annotation(title = paste0("Diagnosed MCI/AD (Ab+)", " (N=", n_clin, ")"), 
                                                     subtitle = expression(italic(FC[parcel] ~ "~" ~ age + pathology + "-mPACC" + sex + motion)),
                                                     theme = theme(plot.subtitle = element_text(size = rel(0.9), 
                                                                                                hjust = 0, 
                                                                                                vjust = -0.05, 
                                                                                                family = "mono",
                                                                                                face = "italic",
                                                                                                margin = margin(l = clin_l_marg, unit = "npc")), 
                                                                   plot.title.position = "plot",
                                                                   plot.title = element_text(size = rel(1), hjust =0, margin = margin(l = clin_l_marg, unit = "npc")))
  )
  
  plt_idx <- 3:5
  if (length(selected_gradients) < 2) plt_idx <- plt_idx-1
  p_clinical_cog[[plt_idx[1]]] <- p_clinical_cog[[plt_idx[1]]] + labs(title = "Age") + theme(plot.title = element_text(vjust = -1.5)) 
  p_clinical_cog[[plt_idx[2]]] <- p_clinical_cog[[plt_idx[2]]] + labs(title = "AD Pathology") + theme(plot.title = element_text(vjust = -1.5)) 
  p_clinical_cog[[plt_idx[3]]] <- p_clinical_cog[[plt_idx[3]]] + labs(title = "-mPACC", subtitle = "(Inverted cognition)") + theme(plot.subtitle = element_text(hjust = 0.5, size = rel(0.6)))
  
  p_clinical_cog <- ggdraw() + draw_plot(p_clinical_cog) + draw_plot_label(ifelse(split, "B", "D"), x = clin_l_marg-l_marg, size = tag_size)
  
  
  
  
  # Everyting, everywhere all at once
  get_net_legend <- function(){
    #scale_factor <- 5
    x <- grad_df %>% filter(study == "biofinder") %>% 
      ggplot(aes(gradient1, gradient3, color = name)) +
      geom_point(alpha = 0.5) +
      labs(color = "Yeo Network") +
      theme_bw(base_size = b_size) +
      guides(color = guide_legend(
        label.hjust=0,
        byrow = TRUE,
        nrow = 1, 
        override.aes = list(size = rel(3))
      )) +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe) +
      theme(
        legend.position = "bottom",
        #legend.key.size = unit(0.0, "cm"),
        legend.key.spacing.x = unit(0.75, "cm"),
        legend.direction = "horizontal",
        #legend.title.position = "",
        legend.text.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.8), margin = margin(l = 4, r = 6, unit = "pt")),
        legend.background = element_blank()
      )
    
    leg <- ggpubr::get_legend(x)
    leg <- ggpubr::as_ggplot(leg)
    leg 
  }
  
  net_legend1 <- get_net_legend()
  
  get_net_legend <- function(){
    #scale_factor <- 5
    x <- grad_df %>% filter(study == "biofinder") %>% 
      mutate(name = factor(name, levels = c("DorsAttn", "SomMot", "SalVentAttn", "Default", "Limbic", "Cont", "Vis"))) %>% 
      ggplot(aes(gradient1, gradient3, color = name)) +
      geom_point(alpha = 0.5) +
      theme_bw(base_size = b_size) +
      labs(color = "Yeo Network") +
      guides(color = guide_legend(
        label.hjust=0,
        byrow = TRUE,
        nrow = 2, 
        reverse = FALSE,
        override.aes = list(size = rel(3))
      )) +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe) +
      theme(
        legend.position = "bottom",
        #legend.key.size = unit(0.0, "cm"),
        legend.key.spacing.x = unit(0, "cm"),
        legend.key.spacing.y = unit(-1.35, "cm"),
        legend.direction = "horizontal",
        #legend.title.position = "",
        legend.text.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.8), margin = margin(l = 2, r = 10, unit = "pt")),
        legend.background = element_blank()
      )
    
    leg <- ggpubr::get_legend(x)
    leg <- ggpubr::as_ggplot(leg)
    leg 
  }
  net_legend2 <- get_net_legend()
  
  
  if (split) {
    fig_one_two <- list()
    
    
    if (boxed) {
      fig_one_two[[1]] <- ggdraw() +
        draw_plot(p_bf, x = 0, y = 0.05, width = 0.49, height = 0.95) +
        #draw_plot(p_bf_long, x = 0.26, y = 0.05, width = 0.36, height = 0.95) +
        draw_plot(p_a, x = 0.51, y = 0.05, width = 0.49, height = 0.95) +
        draw_plot(net_legend1, x = 0.3,  y = 0.025, width = 0.4, height = 0.015)
      
    } else {
      fig_one_two[[1]] <- ggdraw() +
        draw_plot(p_bf, x = 0, y = 0.05, width = 0.5, height = 0.95) +
        #draw_plot(p_bf_long, x = 0.26, y = 0.05, width = 0.36, height = 0.95) +
        draw_plot(p_a, x = 0.5, y = 0.05, width = 0.5, height = 0.95) +
        draw_plot(net_legend1, x = 0.4,  y = 0.025, width = 0.3, height = 0.015)
    }
    
    
    if (boxed) {
      fig_one_two[[2]] <- ggdraw() +
        draw_plot(p_health_cog, x = 0.0, y = 0.05, width = 0.6, height = 0.95) +
        draw_plot(p_clinical_cog, x = 0.5, y = 0.05, width = 0.49, height = 0.95) +
        draw_plot(net_legend2, x = 0.4,  y = 0.025, width = 0.4, height = 0.015) +
        theme(plot.background = element_rect(color = "black"))
    } else {
      fig_one_two[[2]] <- ggdraw() +
        draw_plot(p_health_cog, x = 0.0, y = 0.00, width = 0.6, height = 1) +
        draw_plot(p_clinical_cog, x = 0.5, y = 0.00, width = 0.49, height = 1) +
        draw_plot(net_legend2, x = 0.02,  y = 0.035, width = 0.2, height = 0.015) 
    }
    
    
    
    return(fig_one_two)
    
  } else {
    
    if (boxed) {
      
      
      bottom_plots <- ggdraw() +
        draw_plot(p_health_cog, x = 0.0, y = 0.0, width = 0.6, height = 1) +
        draw_plot(p_clinical_cog, x = 0.5, y = 0.0, width = 0.49, height = 1) +
        draw_plot(net_legend2, x = 0.02,  y = 0.035, width = 0.2, height = 0.015) +
        theme(plot.background = element_rect(color = "black", linewidth = 1))
      
      figure1 <- ggdraw() +
        draw_plot(p_bf, x = 0, y = 0.505, width = 0.495, height = 0.495) +
        #draw_plot(p_bf_long, x = 0.26, y = 0.50, width = 0.36, height = 0.5) +
        draw_plot(p_a, x = 0.505, y = 0.505, width = 0.495, height = 0.495) +  
        draw_plot(bottom_plots, x = 0.0, y = 0.0, width = 1, height = 0.49) 

    } else {
      figure1 <- ggdraw() +
        draw_plot(p_bf, x = 0, y = 0.50, width = 0.5, height = 0.5) +
        #draw_plot(p_bf_long, x = 0.26, y = 0.50, width = 0.36, height = 0.5) +
        draw_plot(p_a, x = 0.5, y = 0.50, width = 0.5, height = 0.5) +  
        draw_plot(p_health_cog, x = 0.0, y = 0.01, width = 0.6, height = 0.5) +
        draw_plot(p_clinical_cog, x = 0.5, y = 0.01, width = 0.49, height = 0.5) +
        draw_plot(net_legend, x = 0.4,  y = 0.00, width = 0.4, height = 0.015)
    }
  
    
    
    figure1
    
  }
  
}




net_names_plot <- function(net_names, text_size = 3.8, vertical = TRUE, tile_height = 0.9, alpha_ = 0.8) {
  get_text_color <- function(bg_color) {
    rgb <- col2rgb(bg_color)
    luminance <- (0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]) / 255
    ifelse(luminance > 0.5, "black", "white")
  }
  
  net_names$text_color <- get_text_color(net_names$col)
  if (vertical) {
    net_names$x <- factor(net_names$name, levels = net_names$name)
    net_names$y <- 1
  } else {
    net_names$y <- factor(net_names$name, levels = net_names$name)
    net_names$x <- 1
  }
  
  
  ggplot(net_names, aes(x = x, y = y)) +
    geom_tile(aes(fill = col), alpha = alpha_, width = 0.95, height = tile_height, show.legend = FALSE) +
    geom_text(aes(label = name, color = text_color), size = text_size) +
    scale_fill_identity() +
    scale_color_identity() +
    theme_void() +
    #coord_fixed(ratio = 1 / spacing) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    )
}



overlaid_main_results <- function(subject_data, fc_matrix, 
                                  atlas_geometry = readRDS("data/atlas_data/Schaefer2018_1000Parcels_geometry.rds"),
                                  network_geometry = readRDS("data/atlas_data/Yeo2011_7_geometry.rds"),
                                  base_size_ = 21){
  
  x <- nodal_regression_fits(subject_data %>% filter(fmri_bl),
                             fc_matrix,
                             roi_names = rois,
                             logistic = FALSE,
                             vectorised = TRUE,
                             id_var = "image_file",
                             dep_var = "FC",
                             model_formula = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"))
  
  est <- get_nodal_ests(x)
  
  p_age <-  est %>%
    filter(term == "age") %>%
    inner_join(atlas_geometry, by = "region") %>%
    ggplot() +
    geom_sf(aes(
      fill = statistic,
      geometry = geometry), linewidth= 0.1,
      show.legend = FALSE)+
    theme_void(base_size = base_size_)+
    labs(fill = "T-value", title = "Age"
    ) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.title.position = "top",
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(color = "black", hjust = 0.5),
          legend.key.width = unit(1, "null"),
          plot.margin = margin(1,1,1,1, "cm")
    ) +
    #guides(color = guide_legend(override.aes = list(size = 1))) +
    scale_fill_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red") 
    ) +
    geom_sf(data = network_geometry %>% drop_na(),
            aes(
              #fill = region,
              color = name,
              geometry = geometry
            ), alpha = 0, linewidth = 1,
            show.legend = FALSE) +
    scale_color_manual(
      values = net_names %>% select(name, col) %>% deframe(),
      guide = "none"
    )
  
  
  range_stat <- est %>% filter(term == "age") %>% pull(statistic) %>% range
  y_range_df <- data.frame(statistic = seq(range_stat[1], range_stat[2], length.out = 100))
  
  hist_age <- est %>% filter(term == "age") %>% 
    inner_join(grad_df %>% filter(study=="biofinder")) %>% 
    select(statistic, region, gradient3) %>% 
    inner_join(roi_data) %>% 
    ggplot(aes(x = fct_rev(fct_reorder(region, yeo_label)) #region
               , y = statistic, fill = name)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = net_names %>% select(name, col) %>% deframe()) +
    labs(y = "T-value", x = "Parcels") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      ggside.panel.scale = 0.05
    ) +
    ggnewscale::new_scale_fill() +
    geom_xsidetile(aes(x = region, y = 0, fill = gradient3), show.legend = FALSE) +
    scale_fill_gradient2(
      low = gradient_cols[1, 3],
      mid = "white",
      high = gradient_cols[2, 3]
    ) +
    ggnewscale::new_scale_fill() +
    geom_ysidetile(data = y_range_df,
                   aes(x = 0, y = statistic, fill = statistic), show.legend = FALSE, inherit.aes = FALSE) +
    scale_fill_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red")
    ) +
    theme_ggside_void() +
    ggside(x.pos = "bottom", y.pos = "left")
  
  
  p_path <- est %>%
    filter(term == "pathology_ad") %>%
    inner_join(atlas_geometry, by = "region") %>%
    ggplot() +
    geom_sf(aes(
      fill = statistic,
      geometry = geometry), linewidth= 0.1,
      show.legend = FALSE)+
    theme_void(base_size = base_size_)+
    labs(fill = "T-value", title = "AD pathology"
    ) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.title.position = "top",
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(color = "black", hjust = 0.5),
          legend.key.width = unit(1, "null"),
          plot.margin = margin(1,1,1,1, "cm")
    )  +
    scale_fill_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red")
    )  + geom_sf(data = network_geometry %>% drop_na(),
                 aes(
                   #fill = region,
                   color = name,
                   geometry = geometry
                 ), alpha = 0, linewidth = 1,
                 show.legend = FALSE) +
    scale_color_manual(
      values = net_names %>% select(name, col) %>% deframe(),
      guide = "none"
    )
  
  range_stat <- est %>% filter(term == "pathology_ad") %>% pull(statistic) %>% range
  y_range_df <- data.frame(statistic = seq(range_stat[1], range_stat[2], length.out = 100))
  
  hist_path <- est %>% filter(term == "pathology_ad") %>% 
    inner_join(grad_df %>% filter(study=="biofinder")) %>% 
    select(statistic, region, gradient1) %>% 
    inner_join(roi_data) %>% 
    arrange(yeo_label) %>% 
    ggplot(aes(x = fct_rev(fct_reorder(region, yeo_label))  #region#fct_rev(fct_reorder(region, gradient1))
               , y = statistic, fill = name)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = net_names %>% select(name, col) %>% deframe()) +
    labs(y = "", x = "Parcels") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      ggside.panel.scale = 0.05
    ) +
    ggnewscale::new_scale_fill() +
    geom_xsidetile(aes(x = region, y = 0, fill = gradient1), show.legend = FALSE) +
    scale_fill_gradient2(
      low = gradient_cols[1, 1],
      mid = "white",
      high = gradient_cols[2, 1]
    ) +
    ggnewscale::new_scale_fill() +
    geom_ysidetile(data = y_range_df,
                   aes(x = 0, y = statistic, fill = statistic), show.legend = FALSE, inherit.aes = FALSE) +
    scale_fill_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red")
    ) +
    theme_ggside_void() +
    ggside(x.pos = "bottom", y.pos = "left")
  
  x <- ggplot() +
    geom_sf(data = network_geometry %>% drop_na(),
            aes(
              #fill = region,
              color = name,
              geometry = geometry
            ), alpha = 0, linewidth = 1,
            show.legend = TRUE) +
    scale_color_manual(
      values = net_names %>% select(name, col) %>% deframe()
    ) +  
    guides(color = guide_legend(
      label.hjust=0,
      byrow = TRUE,
      nrow = 1, 
      override.aes = list(size = rel(3))
    )) +
    theme(
      legend.position = "bottom",
      #legend.key.size = unit(0.0, "cm"),
      legend.key.spacing.x = unit(0.75, "cm"),
      legend.direction = "horizontal",
      #legend.title.position = "",
      legend.text.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = rel(0.8), margin = margin(l = 4, r = 6, unit = "pt")),
      legend.background = element_blank()
    )
  
  x <- ggpubr::get_legend(x)
  
  library(cowplot)
  
  ggdraw() +
    draw_plot(p_age, height = 0.6, width = 0.5, x = 0.025, y = 0.425) +
    draw_plot(p_path, height = 0.6, width = 0.5, x = 0.525, y = 0.425) + 
    draw_plot(hist_age, height = 0.35, width = 0.475, x = 0.01, y = 0.1) +
    draw_plot(hist_path, height = 0.35, width = 0.475, x = 0.51, y = 0.1) +
    draw_plot(x, height = 0.1, width = 0.8, x = 0.1, y = 0)
  
}



plot_grads_over_params <- function(connectome_list, 
                                   atlas_geometry = readRDS("data/atlas_data/Schaefer2018_1000Parcels_geometry.rds"),
                                   param_grid = NULL, n_gradients = 1:3, base_size_=18) {

  require(sf)
  require(ggside)
  require(ggpmisc)
  require(patchwork)
  
  old <- theme_set(theme_bw(base_size = base_size_))
  theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA),
               legend.background = element_rect(fill = "transparent", colour = NA),
               legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                          col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), #"#DCF8A4"
                          label = c(1:7))

  
  connectome_est_list <- connectome_list
  
  if (is.null(param_grid)){
    params <- expand_grid(method = c("pca", "diffusion"), affinity = c(FALSE, TRUE), sim_method = "cosine", threshold = c(0.0, 0.25, 0.5, 0.75)) %>% 
      filter(!(method == "diffusion" & !affinity)) %>% 
      filter(!(method == "pca" & affinity)) %>% 
      mutate(sim_method = ifelse(!affinity, NA, sim_method))
  } else {
    params <- param_grid
  }
  
  gradient_data <- c()
  varexp_df <- c()
  for (i in 1:nrow(params)) {
    param_i <- params[i, ]
    grad_list <- get_gradients(connectome_ests = connectome_est_list,
                               n_gradients = n_gradients,
                               threshold = param_i$threshold,
                               similarity_method = param_i$sim_method,
                               on_affinity = param_i$affinity,
                               method = param_i$method,
                               visualize = FALSE,
                               side_density = FALSE)
    
    gradient_data <- rbind(gradient_data, grad_list$gradients)
    varexp_df <- rbind(varexp_df, grad_list$varexp)
  }
  
  gradient_data <- gradient_data %>% distinct()
  grad_char <- paste0("gradient", n_gradients)
  
  std_grad_plots <- list()
  
  for (grad in grad_char) {
    
    low_col <- gradient_cols[1, grad]
    hi_col <- gradient_cols[2, grad]
    
    std_grad_plots[[grad]] <- gradient_data %>% filter(study=="margulies") %>% 
      mutate(region = rois) %>% 
      inner_join(atlas_geometry, by = "region") %>%
      ggplot() +
      geom_sf(aes(
        fill = .data[[grad]],
        geometry = geometry), linewidth= 0.2,
        show.legend = FALSE)+
      theme_void(base_size = base_size_)+
      labs(fill = "", title = "Margulies") +
      theme(legend.position = "",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(color = "black", hjust = 0.5),
            plot.margin = margin(t = 0,  
                                 r = 0, 
                                 b = 0,  
                                 l = 0)
      ) +
      scale_x_continuous(position = "top")+
      scale_fill_gradient2(
        low = low_col,
        mid = "white",
        high = hi_col
      ) 
    
    
  }
  
  grad_text <- list()
  
  for(g_txt in paste0("G", n_gradients)) {
    grad_text[[g_txt]] <- 
      ggplot() +
      theme_void() +  
      annotate("text", x = 0.5, y = 0.5, label = g_txt, size = rel(20), 
               fontface = "bold",
               color = "#323232",
               hjust = 0.5, vjust = 0.5)
  }
  
  
  
  gradient_plots <- list()
  for (stud in names(connectome_est_list)){
    plt_idx <- 0
    for (grad in grad_char) {
      
      low_col <- gradient_cols[1, grad]
      hi_col <- gradient_cols[2, grad]
      
      for (params_row in 1:nrow(params) ){
        
        pars <- params[params_row, ]
        
        plot_pars <- pars %>% mutate(affinity = ifelse(is.na(sim_method), "No", sim_method), 
                                     method = ifelse(method=="pca", "PCA", "Diffusion")) %>% select(-sim_method) %>% 
          rename(thresh = threshold)
        par_char <- paste(names(plot_pars), plot_pars[1, ], sep = ": ") %>% str_to_sentence() 
        
        
        variance_explained <- varexp_df %>% semi_join(pars, by = join_by(method, affinity, sim_method, threshold)) %>% 
          filter(study == stud)
        
        
        plt_idx <- plt_idx + 1
        p <- gradient_data %>% 
          filter(study == stud) %>% 
          semi_join(pars, by = join_by(method, affinity, sim_method, threshold)) %>% 
          inner_join(atlas_geometry, by = "region") %>%
          ggplot() +
          geom_sf(aes(
            fill = .data[[grad]],
            geometry = geometry), linewidth= 0.2,
            show.legend = FALSE)+
          theme_void(base_size = base_size_) +
          labs(fill = "", title = paste0(plot_pars[, "method"],", ", par_char[3]),
               subtitle = ifelse(grad == grad_char[1], pars[1, "threshold"], NA),
               y = ifelse(plot_pars[, "method"] == "PCA" & plot_pars[, "thresh"] == 0, "BioFINDER", ""),
               caption = paste0(round(variance_explained %>% pull(grad)*100), "% var exp")
          ) +
          theme(legend.position = "",
                plot.caption = element_text(vjust = 8),
                panel.background = element_rect(fill = "transparent", colour = NA),
                plot.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                axis.title.y = element_text(angle=90),
                plot.title = element_blank(),
                plot.subtitle = element_text(hjust = 0.5),
                plot.margin = margin(t = 0,  
                                     r = 0, 
                                     b = -10,  
                                     l = 0)
          ) +
          scale_fill_gradient2(
            low = low_col,
            mid = "white",
            high = hi_col 
          ) 
        
        if (grad != grad_char[1]) p <- p + theme(plot.title = element_blank(),
                                                 plot.subtitle = element_blank())
        gradient_plots[[stud]][[plt_idx]] <- p
        
      }
    }
  }
  
  
  plots <- list()
  
  ref_grads <- gradient_data %>% filter(study == "margulies")
  
  for (stud in names(connectome_est_list)){
    plt_idx <- 0
    plot_data_stud <- gradient_data %>% filter(study %in% c(stud))
    for (grad in grad_char) {
      
      low_col <- gradient_cols[1, grad]
      hi_col <- gradient_cols[2, grad]
      
      plot_data_grad <- plot_data_stud %>% pivot_longer(starts_with("gradient"), names_to = "gradient", values_to = "value") %>% 
        filter(gradient == grad) %>% 
        pivot_wider(names_from = "study", values_from = "value")
      
      for (params_row in 1:nrow(params) ){
        pars <- params[params_row, ]
        
        plot_data <- plot_data_grad %>% semi_join(pars) %>% 
          mutate(margulies = ref_grads %>% pull(grad))
        
        x_range_df <- data.frame(x = seq(min(plot_data[, stud]), max(plot_data[, stud]), length.out = 100))
        
        y_range_df <- data.frame(margulies = seq(min(plot_data[, "margulies"]), max(plot_data[, "margulies"]), length.out = 100))
        
        p <- plot_data %>% 
          ggplot(aes(x = .data[[stud]], y = margulies,
                     color = name)) +
          geom_point(alpha = 0.2) +
          stat_poly_eq(color = "#323232", label.x = "left", label.y = "top", size = 6.5) +
          stat_poly_line(se = FALSE, color = "#323232") +
          labs(
            #title = str_to_title(grad),
            y = "",
            x = str_to_title(stud),
            #tag = tag_labs[count],
            color = "Network") +
          theme(
            legend.position = "",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            ggside.panel.scale = 0.05,
            plot.margin = margin(t = 0,  
                                 r = 5, 
                                 b = 0,  
                                 l = 5),
            axis.text.x = element_text(angle =45, hjust = 1)) +
          scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
          scale_y_continuous(limits = c(min(plot_data["margulies"]), max(plot_data["margulies"])))+
          scale_x_continuous(limits = c(min(plot_data[stud]), max(plot_data[stud])), guide = guide_axis(check.overlap = TRUE))
        
        y_side_switch <-  ifelse((pars %>% pull(method) == "pca" & pars %>% pull(threshold) == 0), 1, 0)
        
        p <- p +
          ggnewscale::new_scale_fill() +
          geom_xsidetile(data = x_range_df, 
                         aes(x = x, y = 0, fill = x), 
                         show.legend = FALSE, inherit.aes = FALSE) +
          scale_fill_gradient2(
            low = low_col,
            mid = "white",
            high = hi_col
          ) + ggnewscale::new_scale_fill() +
          geom_ysidetile(data = y_range_df,
                         aes(x = 0, y = margulies, fill = margulies), 
                         alpha = y_side_switch,
                         show.legend = FALSE, 
                         inherit.aes = FALSE) +
          scale_fill_gradient2(
            low = low_col,
            mid = "white",
            high = hi_col
          )
        
        
        p <-  p +
          theme_ggside_void() +
          ggside(x.pos = "bottom", y.pos = "left")
        
        plt_idx <- plt_idx + 1
        plots[[stud]][[plt_idx]] <- p
      }
    }
  }
  
  
  no_of_grads = length(n_gradients)
  n_study <- length(connectome_est_list)
  n_plot_cols <- nrow(params)#((no_of_grads*n_study) + (n_study - 1))
  row_heights <- head(rep(c(1, -0.2, 1, 0.2), no_of_grads), -1)
  
  param_combos <- length(unique(params$threshold))
  n_total_cols <- n_plot_cols + floor(n_plot_cols / param_combos)
  col_widths <- head(c(1, rep(c(rep(1, param_combos), 0.3), floor(n_plot_cols / param_combos))), -1)
  
  layout <- c(
    area(3, 1))
  
  for (g in seq(7, no_of_grads + 8, by = 4)) {
    layout <- c(layout, area(g, 1))
  }
  
  for (g_txt in seq(1, no_of_grads + 8, by = 4) ) {
    layout <- c(layout, area(g_txt, 1))
  }
  
  
  rows <- layout$t - 1
  
  
  for (row in rows) {
    if (row %% 2 == 0) next
    for (col in 2:(n_total_cols + 1)) {
      if (((col - 2) %% (param_combos + 1)) == param_combos) next  
      layout <- c(layout, area(row, col))
    }
  }
  
  for (row in rows + 1) {
    if (row %% 2 == 0) next
    for (col in 2:(n_total_cols + 1)) {
      if (((col - 2) %% (param_combos + 1)) == param_combos) next  
      layout <- c(layout, area(row, col))
    }
  }
  
  
  
  
  patch_plots <- list()
  for (stud in names(connectome_est_list)){
    pp <- 
      std_grad_plots[[1]] + std_grad_plots[2:no_of_grads] +
      grad_text +
      plots[[stud]][1:(nrow(params)*no_of_grads)] +
      gradient_plots[[stud]][1:(nrow(params)*no_of_grads)] +
      plot_layout(design = layout, axis_titles = "collect", axes = "collect_y", guides = "collect",
                  heights = row_heights,
                  widths = col_widths
      ) & 
      theme(legend.position = "",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            axis.title.x =  element_blank(),
            plot.tag.position  = c(.9, .96))  
    
    patch_plots[[stud]] <- pp
  }
  
  
  return(patch_plots)
}



longitudinal_and_window_analysis <- function(long_df, 
                                             b_size = 18,
                                             run_windowing = FALSE,
                                             processed_dir = "data/processed_and_cleaned",
                                             mod_formula = formula(paste0("FC ~ age_bl + time + pathΔ + path_bl + sex + rsqa__MeanFD + (1 | sid)"))) {
  
  
  library(cowplot)
  library(tidyverse)
  source("src/util_vis.R")
  
  long_bf_ <- long_df
  
  old <- theme_set(theme_bw(base_size = b_size))
  theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA),
               legend.background = element_rect(fill = "transparent", colour = NA),
               legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  #######################
  # Linear longitudinal
  #######################
  
  
  bf_longitudinal <-  plot_gradient_relationships(long_bf_, 
                                                  gradient_data = grad_df %>% filter(study=="biofinder"), 
                                                  gradients = c(1, 3),
                                                  empty_row_height = -0.15,
                                                  gradient_colors = gradient_cols,
                                                  r2_size = rel(3.8),
                                                  base_size_ = b_size,
                                                  list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                  covariates = c("time", "sex", "rsqa__MeanFD"),
                                                  filter_criteria = quo(),
                                                  show_networks = FALSE,
                                                  plt_title = "",
                                                  tag_prefix = "",
                                                  tag_sep = "",
                                                  layout_construction = "horizontal",
                                                  include_gradient_plots = TRUE,
                                                  right_term_side = FALSE,
                                                  cache_runs = FALSE,
                                                  longitudinal = TRUE,
                                                  sub_id = "sid",
                                                  longitudinal_formula = formula(paste0("FC ~ time  + age_bl + path_bl + pathΔ + sex + rsqa__MeanFD + (1 | sid)")))
  
  n_long <- bf_longitudinal$n %>% as.numeric()
  long_l_marg <- 0
  p_bf_long <- bf_longitudinal$plot + 
    plot_annotation(title = waiver(), #paste0("Longitudinal", " (N=", n_long, ")"), 
                    subtitle = waiver(), #bquote(FC[parcel] ~ "~" ~ age_bl + path_bl + Δpath + time + sex + motion +"(1|sub)"),
                    theme = theme(plot.background = element_rect(color = "black"),
                                  plot.subtitle = element_text(size = rel(0.8), 
                                                               hjust = 0, 
                                                               vjust = -0.05, 
                                                               margin = margin(l = 0.05, unit = "npc")), 
                                  plot.title = element_text(size = rel(1), 
                                                            hjust =0, margin = margin(l = 0.05, unit = "npc")))
    )  &
    theme(#text = element_text(size = text_size),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 16),
      plot.tag = element_blank())
  
  p_bf_long[[1]] <- p_bf_long[[1]] + theme(plot.title = element_text(vjust = -3))
  p_bf_long[[2]] <- p_bf_long[[2]] + theme(plot.title = element_text(vjust = -3))
  p_bf_long[[3]] <- p_bf_long[[3]] + labs(title = "Age BL") + theme(plot.title = element_text(vjust = -1.5))
  p_bf_long[[4]] <- p_bf_long[[4]] + labs(title = "Pathology BL") + theme(plot.title = element_text(vjust = -1.5))
  p_bf_long[[5]] <- p_bf_long[[5]] + labs(title = bquote(Δ* "Pathology" ~ (t[i]-t[0]) ) ) + theme(plot.title = element_text(vjust = -4))
  
  ###########################
  # Windowing
  ###########################
  
  long_bf_ %>% filter(fmri_bl) %>% 
    select(age_bl, path_bl) %>% 
    mutate(age_bl = cut(age_bl, 30),
           path_bl = cut(path_bl, 30)) %>% 
    group_by(age_bl, path_bl) %>% 
    summarise(n = n(), .groups = "drop")  -> x
  
  num_bins <- 30
  
  long_bf_ %>% 
    filter(fmri_bl) %>% 
    select(age_bl, path_bl) %>% 
    mutate(
      age_bin = cut(age_bl, breaks = num_bins, include.lowest = TRUE, labels = FALSE),
      path_bin = cut(path_bl, breaks = num_bins, include.lowest = TRUE, labels = FALSE)
    ) %>%
    mutate(
      age_mid = min(age_bl) + (age_bin - 0.5) * (max(age_bl) - min(age_bl)) / num_bins,
      path_mid = min(path_bl) + (path_bin - 0.5) * (max(path_bl) - min(path_bl)) / num_bins
    ) %>% 
    group_by(age_mid, path_mid) %>% 
    summarise(n = n(), .groups = "drop") -> x
  
  x %>%
    ggplot() + 
    geom_tile(aes(x=age_mid, y=path_mid, fill=n), color = "black") +  
    annotate("rect", xmin = 18.5, xmax = 88, ymin = 0, ymax = 0.35, fill = NA, color = "black") +
    annotate("rect", xmin = 19, xmax = 88.5, ymin = 0.015, ymax = 0.365, fill = NA, linetype = 2, color = "black") +
    annotate("rect", xmin = 20, xmax = 45, ymin = -0.025, ymax = 1, fill = NA, color = "black") +
    annotate("rect", xmin = 21, xmax = 46, ymin = -0.015, ymax = 1.01, fill = NA, linetype = 2, color = "black") +
    annotate("curve", x = 90, y = 0.20, xend = 90, yend = 0.365,
             arrow = arrow(length = unit(0.1, "inches"))) +
    annotate("curve", x = 35, y = 1.025, xend = 46, yend = 1.025,
             arrow = arrow(length = unit(0.1, "inches")), curvature = -0.3) +
    annotate("text", x = 40, y = 1.1, label = "1 Year", size = 6) +
    annotate("text", x = 95, y = 0.305, label = "0.015 Path", angle = -90, size = 6) +
    theme_gray(base_size = b_size) +
    scale_x_continuous("Age BL") + 
    scale_y_continuous("Pathology BL", breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) + 
    scale_fill_continuous("Count", breaks = c(2, 5, 8)) + 
    theme(#axis.text = element_text(size = 12),
      legend.position = "top",
      #legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-25,-10),
      legend.justification="right",
      legend.title = element_text(vjust = 0),
      legend.text.position = "top",
      legend.text = element_text(vjust = -1),
      plot.background = element_blank(),
      #legend.position.inside = c(0.95, 0.825),
      legend.background = element_blank(),
      #title = element_text(size = 12, face = "bold"),
      panel.border = element_rect(linewidth = 1, color = "black", fill = NA)) -> hist_map
  
  
  ######################
  # Age window analysis
  ######################
  
  long_bf_ %>% filter(fmri_bl) %>% 
    select(age_bl) -> age_df
  
  win_size <- 25
  frames <- data.frame(win_n = seq(min(age_df$age_bl) %>% round(), max(age_df$age_bl) %>% round() - win_size, by = 1)) %>%
    mutate(
      window_min = win_n,
      window_max = win_n + win_size,
    ) %>% 
    rowwise() %>%  
    mutate(sample_size = sum(age_df$age_bl >= window_min & age_df$age_bl <= window_max)) %>%
    mutate(mean_age = mean(age_df$age_bl[age_df$age_bl >= window_min & age_df$age_bl <= window_max])) %>%
    ungroup() %>% 
    mutate(window_range = paste0("[",window_min, ", ", window_max, "]"))
  
  
  if (run_windowing) {
    sliding_window_results <- data.frame()
    for (win_n in frames$win_n) {
      window_df <- long_bf_ %>% filter((age_bl >= win_n) , (age_bl <= win_size + win_n))
      ests_window <- nodal_lmm_ests(window_df, fc_measures_bf$affinity, roi_names = rois, model_formula = mod_formula)
      ests_window$window <- win_n
      sliding_window_results <- rbind(sliding_window_results, ests_window)
    }
    write_rds(sliding_window_results, file.path(processed_dir, "sliding_window_age_res.rds"))
  } else {
    sliding_window_results <- readRDS(file.path(processed_dir, "sliding_window_age_res.rds"))
  }
  
  grad_cor_age <- sliding_window_results %>% 
    filter(!(term %in% c("(Intercept)", "time", "rsqa__MeanFD", "sex"))) %>% 
    group_by(window, term, n) %>% 
    summarise(Gradient1 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient1), statistic), 
              Gradient3 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient3), statistic)) %>% 
    ungroup() %>% 
    inner_join(frames, join_by(window == win_n))
  
  ###########################
  # Pathology window analysis
  ###########################
  
  long_bf_ %>% filter(fmri_bl) %>% 
    select(path_bl) -> path_df
  #seq(min(path_df$path_bl) %>% round(), max(path_df$path_bl) %>% round() - win_size, by = 0.01)
  win_size <- 0.35
  steps <- seq(min(path_df$path_bl) %>% round(), (max(path_df$path_bl) %>% round() - win_size) + 0.01, by = 0.015)
  steps %>% length()
  frames_path <- data.frame(win_n = steps) %>%
    mutate(
      window_min = win_n,
      window_max = win_n + win_size,
    ) %>% 
    rowwise() %>%  # Ensures calculations are done row-by-row
    mutate(sample_size = sum(path_df$path_bl >= window_min & path_df$path_bl <= window_max)) %>%
    mutate(mean_path = mean(path_df$path_bl[path_df$path_bl >= window_min & path_df$path_bl <= window_max])) %>%
    ungroup() %>% 
    mutate(window_range = paste0("[",window_min, ", ", window_max, "]"))
  
  if (run_windowing) {
    sliding_window_results_path <- data.frame()
    for (win_n in frames_path$win_n) {
      window_df <- long_bf_ %>% filter((path_bl >= win_n) , (path_bl <= win_size + win_n)) 
      ests_window <- nodal_lmm_ests(window_df, fc_measures_bf$affinity, roi_names = rois, model_formula = mod_formula)
      ests_window$window <- win_n
      sliding_window_results_path <- rbind(sliding_window_results_path, ests_window)
    }
    write_rds(sliding_window_results_path, file.path(processed_dir, "sliding_window_path_res.rds"))
  } else {
    sliding_window_results_path <- readRDS(file.path(processed_dir, "sliding_window_path_res.rds"))
  }
  
  grad_cor_path <- sliding_window_results_path %>% 
    filter(!(term %in% c("(Intercept)", "time", "rsqa__MeanFD", "sex"))) %>% 
    group_by(window, term, n) %>% 
    summarise(Gradient1 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient1), statistic), 
              Gradient3 = cor(grad_df %>% filter(study == "biofinder") %>% pull(gradient3), statistic)) %>% 
    inner_join(frames_path, join_by(window == win_n))
  
  ###########################
  # Plotting window analyses
  ###########################
  
  
  k = 5
  grad_cor_p_age <- grad_cor_age %>% 
    pivot_longer(starts_with("Grad"), names_to = "gradient", values_to = "grad_corr") %>% 
    filter((gradient == "Gradient3" & term == "age_bl") | (gradient == "Gradient1" & term == "path_bl") | (gradient == "Gradient1" & term == "pathΔ") ) %>% 
    mutate(term = case_when(
      term == "age_bl" ~ "Age~at~baseline",
      term == "path_bl" ~ "Pathology~at~baseline",
      term == "pathΔ" ~ "Δ*Pathology~(t[i]-t[0])",
    )) %>%
    ggplot(aes(
      window,
      grad_corr,
      color = gradient,
      fill = gradient
    )) +
    geom_point(aes(alpha = n), show.legend = FALSE) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = k, bs = "cs")) +
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5 )+
    facet_wrap(~term, nrow = 1, labeller = label_parsed) +
    #theme_bw(base_size = 14) +
    scale_x_continuous(labels = frames$window_range[seq(1, length(frames$win_n), by = 5)], 
                       breaks = frames$win_n[seq(1, length(frames$win_n), by = 5)],
                       guide = guide_axis(check.overlap = TRUE)) +
    labs(y = "Gradient corr (r)",
         x = "Baseline Age Window") +
    theme(
      legend.position = "",
      legend.title = element_blank()
    ) +
    geom_xsidecol(aes(y = sample_size, fill = NULL, color = NULL)) +
    scale_xsidey_continuous(breaks = c(0, 250))+
    #theme_ggside_void() +
    ggside(x.pos = "bottom") +
    theme(ggside.panel.scale = 0.2,
          axis.text.x = element_text(angle =-30, hjust =0, size = rel(0.8))) +
    scale_color_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_fill_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_alpha_continuous(range = c(0.3, 1))
  
  ####### PATHOLOGY
  
  
  grad_cor_p_path <- grad_cor_path %>% 
    pivot_longer(starts_with("Grad"), names_to = "gradient", values_to = "grad_corr") %>% 
    filter((gradient == "Gradient3" & term == "age_bl") | (gradient == "Gradient1" & term == "path_bl") | (gradient == "Gradient1" & term == "pathΔ") ) %>%
    mutate(term = case_when(
      term == "age_bl" ~ "Age~at~baseline",
      term == "path_bl" ~ "Pathology~at~baseline",
      term == "pathΔ" ~ "Δ*Pathology~(t[i]-t[0])",
    )) %>%
    ggplot(aes(
      window,
      grad_corr,
      color = gradient,
      fill = gradient
    )) +
    geom_point(aes(alpha = n), show.legend = FALSE) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs")) +
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5 )+
    facet_wrap(~term, nrow = 1, labeller = label_parsed) +
    #theme_bw(base_size = 14) +
    scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75) )+
    scale_x_continuous(labels = frames_path$window_range[seq(1, length(frames_path$win_n), by = 10)], 
                       breaks = frames_path$win_n[seq(1, length(frames_path$win_n), by = 10)],
                       guide = guide_axis(check.overlap = TRUE)) +
    labs(y = "Gradient corr (r)",
         x = "Baseline Pathology Window") +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.12, -0.45),
      legend.direction = "horizontal",
      legend.text = element_text(
        margin = margin(r = 20, l = 10,  unit = "pt")),
      legend.title = element_blank()
    ) +
    geom_xsidecol(aes(y = sample_size, fill = NULL, color = NULL)) +
    #theme_ggside_void() +
    ggside(x.pos = "bottom") +
    scale_xsidey_continuous(breaks = c(0, 250))+
    theme(ggside.panel.scale = 0.2,
          axis.text.x = element_text(angle =-30, hjust =0, size = rel(0.8))) +
    scale_color_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_fill_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_alpha_continuous(range = c(0.3, 1))
  
  window_plots <-  wrap_plots(list(grad_cor_p_age, grad_cor_p_path), nrow = 2)+ 
    plot_annotation(
      theme = theme(
        plot.background = element_rect(color = "black") 
      )
    ) 
  
  #############
  # Fonky way to get legend
  #############
  
  get_net_legend <- function(){
    #scale_factor <- 5
    x <- grad_df %>% filter(study == "biofinder") %>% 
      ggplot(aes(gradient1, gradient3, color = name)) +
      geom_point(alpha = 0.5) +
      labs(color = "Yeo Network") +
      guides(color = guide_legend(
        label.hjust=0,
        byrow = TRUE,
        nrow = 2, 
        reverse = TRUE,
        override.aes = list(size = rel(4))
      )) +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe) +
      theme(
        legend.position = "bottom",
        #legend.key.size = unit(0.0, "cm"),
        legend.key.spacing.x = unit(0, "cm"),
        legend.key.spacing.y = unit(-1.25, "cm"),
        legend.direction = "horizontal",
        #legend.title.position = "",
        legend.text.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.65), margin = margin(l = 5, r = 10, unit = "pt")),
        legend.background = element_blank()
      )
    
    leg <- ggpubr::get_legend(x)
    leg <- ggpubr::as_ggplot(leg)
    leg 
  }
  net_legend <- get_net_legend()
  
  ##########################
  # Put together all
  #########################
  
  p_bf_long_cp <- ggdraw() + 
    draw_plot(p_bf_long) + 
    draw_plot(net_legend, x = 0.15, y = 0.025, width = 0.2, height = 0.05) +
    draw_plot_label("A", size = 20) +
    draw_label(paste0("BioFINDER Longitudinal", " (N=", n_long, ")"), x = 0.065, y = 0.965, hjust = 0, size =  18) +
    draw_label("FC ~ age_bl + path_bl + Δpath + \ntime + sex + motion + (1|sub)", x = 0.015, y = 0.89, hjust = 0, size =  12)
  
  window_plots_cp <- ggdraw() + 
    draw_plot(window_plots) + 
    draw_plot_label("C", size = 20)
  
  lmm_tab <- image_read("paper/figures/conceptual_plot/LMM_table.png")
  
  ggdraw() +
    draw_plot(hist_map, x = 0.01, y = 0.3, width = 0.98, height = 0.7) +
    draw_plot(p_bf_long[[5]] + labs(title = "ΔPath") + theme(plot.title = element_text(vjust = 0)),
              x = 0.39, y = 0.02,  width = 0.3, height = 0.3) +
    draw_plot(p_bf_long[[1]] + theme(plot.title = element_text(vjust = 0)),
              x = 0.69, y = 0.025, width = 0.29, height = 0.29) +
    draw_image(lmm_tab, x = 0.02, y = 0.00, width = 0.35, height = 0.25) +
    draw_plot_label("B", size = 20) + 
    annotate("curve",  x = 0.54, y = 0.09, xend = 0.83, yend = 0.09, linewidth = 1,
             arrow = arrow(length = unit(0.13, "inches"), ends = "both")) +
    annotate("segment",  x = 0.33, y = 0.425, xend = 0.225, yend = 0.22, linewidth = 1,
             arrow = arrow(length = unit(0.13, "inches"))) +
    annotate("label", label = "Parcel-wise LMM \n in age window", x = 0.25, y = 0.3, size = 6) +
    annotate("text", label = "Pearson Correlation", x = 0.7, y = 0.02, size = 6) +
    annotate("segment",  x = 0.38, y = 0.145, xend = 0.43, yend = 0.145, 
             arrow = arrow(length = unit(0.13, "inches")),
             linewidth = 1) +
    theme(plot.background = element_rect(color = "black", linewidth = 1)) -> slide_meth
  
  final_fig <- ggdraw() +
    draw_plot(p_bf_long_cp, x = 0, y = 0.51, width = 0.645, height = 0.49) +
    draw_plot(slide_meth, x = 0.655, y = 0.51, width = 0.345, height = 0.49) +
    draw_plot(window_plots_cp, x = 0, y = 0.0, width = 1, height = 0.5) +
    annotate("curve",  x = 0.89, y = 0.5105, xend = 0.706, yend = 0.39,
             linewidth = 3,
             color = "white",
             curvature = 0.29) +
    annotate("curve",  x = 0.89, y = 0.5105, xend = 0.706, yend = 0.39, linewidth = 1,
             arrow = arrow(length = unit(0.13, "inches")), 
             curvature = 0.29) 
  
  ###################
  # Supplementary fig
  ###################
  
  grad_cor_p_path_supp <- grad_cor_path %>% 
    pivot_longer(starts_with("Grad"), names_to = "gradient", values_to = "grad_corr") %>% 
    mutate(term = case_when(
      term == "age_bl" ~ "Age~at~baseline",
      term == "path_bl" ~ "Pathology~at~baseline",
      term == "pathΔ" ~ "Δ*Pathology~(t[i]-t[0])",
    )) %>%
    ggplot(aes(
      window,
      grad_corr,
      color = gradient,
      fill = gradient
    )) +
    geom_point(aes(alpha = n), show.legend = FALSE) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = k, bs = "cs")) +
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5 )+
    facet_wrap(~term, nrow = 1, labeller = label_parsed) +
    #theme_bw(base_size = 14) +
    scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75) )+
    scale_x_continuous(labels = frames_path$window_range[seq(1, length(frames_path$win_n), by = 10)], 
                       breaks = frames_path$win_n[seq(1, length(frames_path$win_n), by = 10)],
                       guide = guide_axis(check.overlap = TRUE)) +
    labs(y = "Gradient corr (r)",
         x = "Baseline Pathology Window") +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.12, -0.45),
      legend.direction = "horizontal",
      legend.text = element_text(
        margin = margin(r = 20, l = 10,  unit = "pt")),
      legend.title = element_blank()
    ) +
    geom_xsidecol(aes(y = sample_size/2, fill = NULL, color = NULL)) +
    #theme_ggside_void() +
    ggside(x.pos = "bottom") +
    scale_xsidey_continuous(breaks = c(0, 250))+
    theme(ggside.panel.scale = 0.2,
          axis.text.x = element_text(angle =-30, hjust =0, size = rel(0.8))) +
    scale_color_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_fill_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_alpha_continuous(range = c(0.3, 1))
  
  
  grad_cor_p_age_supp <- grad_cor_age %>% 
    pivot_longer(starts_with("Grad"), names_to = "gradient", values_to = "grad_corr") %>% 
    mutate(term = case_when(
      term == "age_bl" ~ "Age~at~baseline",
      term == "path_bl" ~ "Pathology~at~baseline",
      term == "pathΔ" ~ "Δ*Pathology~(t[i]-t[0])",
    )) %>%
    ggplot(aes(
      window,
      grad_corr,
      color = gradient,
      fill = gradient
    )) +
    geom_point(aes(alpha = n), show.legend = FALSE) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = k, bs = "cs")) +
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.5 )+
    facet_wrap(~term, nrow = 1, labeller = label_parsed) +
    #theme_bw(base_size = 14) +
    scale_x_continuous(labels = frames$window_range[seq(1, length(frames$win_n), by = 5)], 
                       breaks = frames$win_n[seq(1, length(frames$win_n), by = 5)],
                       guide = guide_axis(check.overlap = TRUE)) +
    labs(y = "Gradient corr (r)",
         x = "Baseline Age Window") +
    theme(
      legend.position = "",
      legend.title = element_blank()
    ) +
    geom_xsidecol(aes(y = sample_size, fill = NULL, color = NULL)) +
    scale_xsidey_continuous(breaks = c(0, 250))+
    #theme_ggside_void() +
    ggside(x.pos = "bottom") +
    theme(ggside.panel.scale = 0.2,
          axis.text.x = element_text(angle =-30, hjust =0, size = rel(0.8))) +
    scale_color_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_fill_manual(values = c(Gradient3 = "#8A6081", Gradient1 = "#D38A4E")) +
    scale_alpha_continuous(range = c(0.3, 1))
  
  window_plots_both_G <-  wrap_plots(list(grad_cor_p_age_supp, grad_cor_p_path_supp), nrow = 2)+ 
    plot_annotation(
      theme = theme(
        plot.background = element_rect(color = "black") 
      )
    ) 
  
  # Return figures:
  list(main_fig = final_fig, supp_fig = window_plots_both_G)
  
}


