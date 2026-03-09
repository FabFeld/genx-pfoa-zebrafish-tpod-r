#' Multi-group density plot for BMD distributions
#'
#' Creates a combined density plot for multiple substances/groups with optional
#' threshold markers (percentiles, 20th gene, first peak, LCRC).
#'
#' @param df_list list Named list of data.frames containing BMD columns.
#' @param BMD_type character Column name for BMD values (e.g., "BMD.zSD").
#' @param group_var character Column name used to identify groups.
#' @param dens_adj numeric Density adjustment parameter.
#' @param unit character Unit label for the x-axis.
#' @param colors character Optional vector of colors.
#' @param expand_factor numeric Axis expansion factor.
#' @param proportional logical Whether to normalize densities by group size.
#' @param show_elements character Which threshold lines to include.
#' @param vline_size numeric Line width for thresholds.
#' @param dens_alpha numeric Density fill alpha.
#' @param dens_line_width numeric Density line width.
#' @param z_value numeric Optional z value used in the subtitle.
#'
#' @return ggplot2 object.
#' @export
dens_plot_multi <- function(df_list,
                            BMD_type,
                            group_var = "Substance",
                            dens_adj = 2,
                            unit = "µg/L",
                            colors = NULL,
                            expand_factor = 0.1,
                            proportional = FALSE,
                            show_elements = c("percentiles", "gene_20", "first_peak", "LCRC"),
                            vline_size = 2,
                            dens_alpha = 0.5,
                            dens_line_width = 1,
                            z_value = NULL) {
  library(ggplot2)
  library(dplyr)
  library(purrr)
  
  data <- bind_rows(df_list, .id = group_var)
  data <- data %>% filter(!is.na(.data[[BMD_type]]))
  
  group_counts <- table(data[[group_var]])
  total_count <- sum(group_counts)
  
  group_counts <- data %>%
    group_by(.data[[group_var]]) %>%
    summarise(n = sum(!is.na(.data[[BMD_type]]))) %>%
    mutate(label = paste0(.data[[group_var]], " (", n, ")"))
  subtitle_counts <- paste(group_counts$label, collapse = ", ")
  subtitle_counts <- paste0("Number of BMDs: ", subtitle_counts)
  
  if (is.null(colors)) {
    colors <- scales::viridis_pal()(length(unique(data[[group_var]])))
  }
  
  color_names <- setNames(colors[seq_along(unique(data[[group_var]]))],
                          unique(data[[group_var]]))
  
  x_values <- log10(data[[BMD_type]])
  x_range <- range(x_values, na.rm = TRUE)
  expansion <- diff(x_range) * expand_factor
  x_limits <- c(x_range[1] - expansion, x_range[2] + expansion)
  
  # ---- Density estimation ----
  density_data <- map_dfr(unique(data[[group_var]]), function(g) {
    group_data <- data %>% filter(.data[[group_var]] == g)
    log_x <- log10(group_data[[BMD_type]])
    
    if (proportional) {
      weights <- rep(1, length(log_x))
    } else {
      weights <- rep(1 / length(log_x), length(log_x))
    }
    
    dens <- density(log_x,
                    weights = weights / sum(weights),
                    adjust = dens_adj,
                    na.rm = TRUE,
                    n = 512)
    
    tibble(
      x = dens$x,
      y = dens$y * if (proportional) length(log_x) / total_count else 1,
      group = g
    )
  })
  
  # ---- Base plot ----
  if (is.null(z_value)) {
    z_value <- if (exists("z", inherits = TRUE)) z else "NA"
  }

  gg_dens <- ggplot(density_data, aes(x = x, y = y, fill = group, color = group)) +
    geom_area(alpha = dens_alpha, position = "identity", size = dens_line_width) +
    scale_fill_manual(values = color_names, name = group_var) +
    scale_color_manual(values = color_names, name = group_var ) +
    labs(
      x = glue::glue("log10(Gene-level BMD) ({unit})"),
      y = "Density",
      title = paste("Density plot by", group_var),
      subtitle = paste(
        paste0("BMDtype = ",  z_value, "SD"),
        paste0("Dens_adj = ", dens_adj),
        paste0("Proportional =", proportional),
        subtitle_counts,
        sep = "; "
      )
    )+
    theme_bw() +
    theme(legend.position = "right") +
    xlim(x_limits)+
    scale_y_continuous(
      limits =   c(0, NA),
      expand = expansion(mult = c(0,0.05))
    )
  
  # ---- Compute and plot statistical lines ----
  stats_list <- map(unique(data[[group_var]]), function(g) {
    grp_data <- data %>% filter(.data[[group_var]] == g)
    log_x <- log10(grp_data[[BMD_type]])
    
    list(
      group = g,
      percentiles = quantile(log_x, 0.1, na.rm = TRUE),
      gene_20 = if (length(log_x) >= 20) sort(log_x)[20] else NA,
      first_peak = {
        dens <- density(log_x, adjust = dens_adj)
        d_y <- dens$y
        peaks <- which(diff(sign(diff(d_y))) == -2) + 1
        if (length(peaks)) dens$x[peaks[1]] else median(log_x)
      },
      LCRC = {
        df <- grp_data %>% arrange(.data[[BMD_type]])
        if (nrow(df) <= 1) return(log10(df[[BMD_type]][1]))
        df <- df %>%
          mutate(bmc_ratio = lead(.data[[BMD_type]]) / .data[[BMD_type]]) %>%
          filter(!is.na(bmc_ratio))
        above_thresh <- df %>% filter(bmc_ratio > 1.66)
        if (nrow(above_thresh) > 0) {
          log10(above_thresh %>% tail(1) %>% pull(!!sym(BMD_type)))
        } else {
          log10(df[[BMD_type]][1])
        }
      }
    )
  }) %>% setNames(unique(data[[group_var]]))
  
  # Filter requested elements
  all_elements <- c("percentiles", "gene_20", "first_peak", "LCRC")
  valid_elements <- intersect(show_elements, all_elements)
  line_types <- c("solid", "dashed", "dotted", "dotdash")[match(valid_elements, all_elements)]
  labels_base <- c("10th percentile", "20th gene", "1st peak", "LCRC")[match(valid_elements, all_elements)]
  
  # Build vertical line dataframe
  tPOD <- map_dfr(names(stats_list), function(g) {
    stats <- stats_list[[g]]
    vals <- stats[valid_elements]
    tibble(
      xintercept = unlist(vals),
      group = g,
      label = paste0(g, ": ", labels_base[seq_along(vals)]),
      linetype = line_types[seq_along(vals)],
      color = color_names[g],
      rounded = round(unlist(vals), 2),
      linear = round(10^unlist(vals), 2)
    )
  })
  
  if (nrow(tPOD) > 0) {
    gg_dens <- gg_dens +
      geom_vline(
        data = tPOD,
        aes(xintercept = xintercept, color = group, linetype = label),
        size = vline_size,
        key_glyph = "path"
      ) +
      scale_linetype_manual(
        values = setNames(tPOD$linetype, tPOD$label),
        labels = setNames(
          paste0(tPOD$label, ": ", tPOD$rounded, " (", tPOD$linear, " µg/L)"),
          tPOD$label
        ),
        name = "Thresholds"
      ) +
      theme(legend.key.width = unit(3, "cm"))
  }
  
  # Assign thresholds to global env (optional)
  for (grp in names(stats_list)) {
    for (elem in valid_elements) {
      assign(paste0(elem, "_", grp), stats_list[[grp]][[elem]], envir = .GlobalEnv)
    }
  }
  
  message("Density plot created.")
  return(gg_dens)
}
