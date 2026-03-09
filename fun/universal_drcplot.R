universal_drcplot <- function(
    res,
    gene_no = 36,
    sorted_by = "maxychange",
    gene_filter = NULL,         # z.B. function(x) x$maxychange > 1
    facet_by = c("id", "SYMBOL"),
    plot.type = c("dose_fitted", "dose_residuals", "fitted_residuals"),
    dose_log_transfo = TRUE,
    npts = 500,
    nr = NULL,
    nc = NULL,
    showCI = TRUE,
    CIribbon_alpha = 0.1
) {
  plot.type <- match.arg(plot.type)
  facet_by <- match.arg(facet_by)
  
  fitres <- res$drcfit$fitres
  omicdata <- res$omicdata
  
  
 # Frühe Filter: erst genefilter, dann (optional) CI-Filter wenn showCI = TRUE
  if (!is.null(gene_filter)) {
    if (is.function(gene_filter)) {
      # erwartet: gene_filter liefert einen logischen Vektor der Länge nrow(fitres)
      keep <- gene_filter(fitres)
      if (!is.logical(keep) || length(keep) != nrow(fitres)) {
        stop("gene_filter function must return a logical vector of length nrow(fitres)")
      }
      fitres <- fitres[keep, , drop = FALSE]
    } else {
      fitres <- fitres[fitres$id %in% gene_filter, , drop = FALSE]
    }
  }

  if (showCI && !is.null(res$bmdboot)) {
    bmdres_all <- res$bmdboot$res
    if (!is.null(bmdres_all$BMD.zSD.lower) && !is.null(bmdres_all$BMD.zSD.upper)) {
      valid_ci <- !is.na(bmdres_all$BMD.zSD.lower) & !is.na(bmdres_all$BMD.zSD.upper) &
                  is.finite(bmdres_all$BMD.zSD.lower) & is.finite(bmdres_all$BMD.zSD.upper) &
                  (bmdres_all$BMD.zSD.lower > 0) & (bmdres_all$BMD.zSD.upper > 0) &
                  (bmdres_all$BMD.zSD.lower <= bmdres_all$BMD.zSD.upper)
      valid_ci_ids <- unique(bmdres_all$id[valid_ci])
      fitres <- fitres[fitres$id %in% valid_ci_ids, , drop = FALSE]
    }
  }
  
  # Auswahl & Sortierung nach sorted_by
    if (!is.null(gene_no) && !is.null(sorted_by)) {
    fitres <- fitres[!is.na(fitres[[sorted_by]]), ] # NA entfernen
    fitres <- fitres[order(-as.numeric(fitres[[sorted_by]])), ]
    fitres <- head(fitres, gene_no)
  }
  # IDs, SYMBOLs in sortierter Reihenfolge
  ids <- fitres$id
  symbols <- fitres$SYMBOL
  
  # Bringe omicdata in die gleiche Reihenfolge wie fitres
  idx <- match(ids, omicdata$item)
  data <- omicdata$data[idx, , drop = FALSE]
  data.mean <- omicdata$data.mean[idx, , drop = FALSE]
  dose <- omicdata$dose
  
  # Facet-Levels in sortierter Reihenfolge
  lev <- if (!is.null(nr) && !is.null(nc) && (nrow(fitres) < (nr * nc))) {
    c(if (facet_by == "id") ids else symbols,
      strrep(" ", 1:(nr * nc - nrow(fitres))))
  } else {
    if (facet_by == "id") ids else symbols
  }
  
  nitems <- nrow(fitres)
  nobs <- length(dose)
  doseu <- as.numeric(colnames(data.mean))
  ndose <- length(doseu)
  
  # Plotdaten vorbereiten
  dataobs <- data.frame(dose = numeric(), signal = numeric(), facet = character())
  dataobsmean <- data.frame(dose = numeric(), signal = numeric(), facet = character())
  datatheo <- data.frame(dose = numeric(), signal = numeric(), facet = character())
  
  # x-Achse für Fitted
  if (dose_log_transfo) {
    dose_pos <- sort(unique(dose[dose > 0]))
    spacing_factor <- dose_pos[2] / dose_pos[1]
    minx <- if (any(dose == 0)) dose_pos[1] / (3 * spacing_factor) else dose_pos[1]
    maxx <- max(dose)
    xplot <- 10^seq(log10(minx), log10(maxx), length.out = npts)
  } else {
    xplot <- seq(0, max(dose), length.out = npts)
  }
  
  for (i in seq_len(nitems)) {
    facet_val <- if (facet_by == "id") fitres$id[i] else fitres$SYMBOL[i]
    datai <- data[i, ]
    datameani <- data.mean[i, ]
    # Fitted curves
    model <- as.character(fitres$model[i])
    if (model == "exponential") datapred <- fExpo(x = xplot, d = fitres$d[i], b = fitres$b[i], e = fitres$e[i])
    if (model == "Hill") datapred <- fHill(x = xplot, c = fitres$c[i], d = fitres$d[i], b = fitres$b[i], e = fitres$e[i])
    if (model == "log-Gauss-probit") datapred <- fLGauss5p(x = xplot, c = fitres$c[i], d = fitres$d[i], b = fitres$b[i], e = fitres$e[i], f = fitres$f[i])
    if (model == "Gauss-probit") datapred <- fGauss5p(x = xplot, c = fitres$c[i], d = fitres$d[i], b = fitres$b[i], e = fitres$e[i], f = fitres$f[i])
    if (model == "linear") datapred <- xplot * fitres$b[i] + fitres$d[i]
    if (model == "const") datapred <- rep(mean(datai), length(xplot))
    
  dataobs <- rbind(dataobs, 
           data.frame(dose = dose, signal = unname(datai), facet = rep(facet_val, nobs)))
  dataobsmean <- rbind(dataobsmean, 
             data.frame(dose = doseu, signal = unname(datameani), facet = rep(facet_val, ndose)))
    datatheo <- rbind(datatheo,
                      data.frame(dose = xplot, signal = datapred, facet = rep(facet_val, npts)))
  }
  
  dataobs$facet <- factor(dataobs$facet, levels = lev)
  dataobsmean$facet <- factor(dataobsmean$facet, levels = lev)
  datatheo$facet <- factor(datatheo$facet, levels = lev)
  
  # Plot
  if (plot.type == "dose_fitted") {
    g <- ggplot(dataobs, aes(x = dose, y = signal)) +
      geom_point(shape = 1) +
      facet_wrap(~ facet, scales = "free_y", nrow = nr, ncol = nc, drop = FALSE) +
      geom_point(data = dataobsmean, shape = 19) +
      geom_line(data = datatheo, colour = "blue", linewidth =1)
    
    if (dose_log_transfo) g <- g + scale_x_log10()
  } 
  
  # BMD und Konfidenzintervall
  if (!is.null(res$bmdboot) && plot.type == "dose_fitted") {
    bmdres <- res$bmdboot$res
    bmdres <- bmdres[bmdres$id %in% ids, ]
    # Angenommen, g ist dein ggplot-Objekt
    facet_order <- as.character(levels(g$data$facet))
    bmdres$facet <- as.character(bmdres[[facet_by]])
    bmdres$facet <- factor(bmdres$facet, levels = facet_order)
    bmdres <- bmdres[match(facet_order, bmdres$facet), ]
    # BMD-Linie
    g <- g + geom_vline(
      data = bmdres,
      aes(xintercept = BMD.zSD),
      colour = "red", linetype = 1
    )
    # Konfidenzintervall als Schattierung
    if (showCI &&
        !is.null(bmdres$BMD.zSD.lower) &&
        !is.null(bmdres$BMD.zSD.upper)) {
      g <- g + geom_rect(
        data = bmdres,
        aes(xmin = BMD.zSD.lower,
            xmax = BMD.zSD.upper,
            ymin = -Inf, ymax = Inf),
        fill = "red", alpha = CIribbon_alpha,
        inherit.aes = FALSE
      )
    }
  }
   # Beispiel (auskommentiert): zusätzliche Titel / Achsen-Settings innerhalb der Funktion
  g <- g +
    ggtitle(substance) +
    labs(subtitle = paste0("BMD: no of SD = ", z, "; ", method, " models (n = ", sum(!is.na(res$bmdcalc$res$BMD.zSD)), "); niter = ", niter)) +
    xlab(paste0("concentration [", unit, "]")) +
    theme_bw() +
    scale_x_log10(labels = scales::trans_format("log10", math_format(10^.x)),
                  limits = c(min(xplot), NA)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  return(g)
}
