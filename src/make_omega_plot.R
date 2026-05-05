library(psych)
library(ggraph)
library(tidygraph)

plot_cognition_factor_model <- function() {
  
  scale_fac <- 7 / 20
  
  sl <- cog_fac$schmid$sl
  sl <- sl[, 1:5]
  
  latent_factors <- c("global_cog",   "executive",    "memory",       "language",     "visuospat")
  #colnames(domain_scores)
  
  latent_factors <- str_to_title(str_replace(latent_factors, "_", " "))
  domains <- latent_factors[-1]
  colnames(sl) <- latent_factors
  
  sl <- as.data.frame(sl) |>
    rownames_to_column("test") |>
    mutate(
      max_domain = domains[ max.col(across(all_of(domains)), ties.method = "first") ],
      max_loading = pmax(!!!rlang::syms(domains))
    ) |>
    filter(max_loading > 0.2) |>
    mutate(
      max_domain = factor(max_domain, levels = latent_factors)
    ) |>
    arrange(max_domain, desc(max_loading)) |> 
    select(-max_domain, -max_loading) |> 
    column_to_rownames("test")
  
  item_names <- c(adas_delayed_word_recall = "ADAS Delayed Recall",
                  adas_immediate_word_recall_average = "ADAS Immediate Recall",
                  bnt_15_2 = "Boston Naming",
                  vosp_cube = "VOSP Cube",
                  animal_fluency = "Animal Fluency",
                  letter_s = "Letter Fluency",
                  trailmaking_a = "Trail Making A",
                  trailmaking_b = "Trail Making B",
                  symbol_digit = "Symbol Digit",
                  aqt_color_form = "AQT Color-Form")
  
  item_names <- item_names[rownames(sl)]
  
  sl <- sl[enframe(item_names)$name, ]
  rownames(sl) <- item_names
  
  df <- as_tibble(sl, rownames = NA) |> rownames_to_column("item") 
  
  edges <- df %>%
    pivot_longer(cols = -item, names_to = "factor", values_to = "loading") %>%
    filter(!is.na(loading), loading > 0.2) %>%
    mutate(loading = round(loading, 1)) |> 
    rename(from = factor, to = item)
  
  source_data <- edges |>
    transmute(
      factor = from,
      item = to,
      loading = loading
    )
  
  nodes <- tibble(
    name = c(latent_factors, df$item),
    type = c(
      rep("latent", length(latent_factors)),
      rep("item", nrow(df))
    ),
    x = c(
      -0.15, rep(1.15, 4),
      rep(0.5, nrow(df))
    ),
    y = c(0.5,
          seq(0.9, 0.1, length.out = 4),
          seq(1, 0, length.out = nrow(df))
    )
  )
  
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  make_push_vector <- function(edges, push_df, default = 0, unit = "mm") {
    merged <- edges %>%
      left_join(push_df, by = c("to", "from")) %>%
      mutate(push = ifelse(is.na(push), default, push))
    unit_list <- lapply(merged$push, function(x) unit(x, unit))
    push_vec <- do.call(grid::unit.c, unit_list)
    return(push_vec)
  }
  
  push_vec <- make_push_vector(edges = edges, push_df = tibble(
    to = c("VOSP Cube", "ADAS Delayed Recall", "ADAS Immediate Recall"),
    from = c("Executive", "Memory", "Memory"),
    push = c(15, 5, 5) * scale_fac
  ))
  
  
  omega_plot <- ggraph(graph, layout = "manual", x = x, y = y) +
    geom_edge_link(
      aes(label = loading, 
          # Shrink the item-label collision box so long labels do not truncate arrows early.
          end_cap = label_rect(
            node2.name,
            padding = margin(
              0,
              -20 * scale_fac,
              0,
              -20 * scale_fac,
              "mm"
            )
          )
      ),
      angle_calc = "along",
      arrow = arrow(angle = 15, length = unit(3 * scale_fac, "mm"), type = "closed"),
      label_dodge = unit(3 * scale_fac, "mm"),
      label_push = push_vec,
      color = "grey30",
      linewidth = 0.5 * scale_fac,
      label_size = 5 * scale_fac
    ) +
    geom_node_label(
      aes(label = name, filter = type == "item"), 
      label.padding = unit(0.7, "mm"),
      label.r = unit(1.2, "mm"),
      label.size = 0.15,
      size = 5 * scale_fac
    )+
    geom_node_circle(aes(filter = type == "latent", r = 0.11), fill = "white", linewidth = 0.3) +
    geom_node_text(aes(label = name, filter = type == "latent"), size = 5 * scale_fac) +
    theme_void(base_size = 7) +
    ggtitle(" Bifactor model of cognition (Schmid-Leiman)") +
    theme(plot.background = element_rect(color = "black", linewidth = 0.5), plot.title = element_text(margin = margin(l = 10, t = 3), hjust = 0))
  
  list(
    plot = omega_plot,
    data = source_data
  )
}
