## Flower plots ---------------------------------------------------------------

#' Flower plot of a statistical unit
#'
#' A flower plot shows all the distributions observed on a unit as the petals
#' of a flower. Each variable is a stem going from the centre outwards; along
#' the stem, the distance from the centre is the value (its position in the
#' range of the variable, as in the other plots of the package), and at each
#' value a petal opens sideways, with a width proportional to its probability
#' and the colour of the value. For discrete distributions each support value
#' is a petal (an ellipse); for histograms the petals of the classes are
#' joined into a violin-shaped petal whose width follows the density.
#' Compared with a Green Eyes plot, where the radius is the cumulative
#' probability, a flower plot shows the \emph{shape} of the distributions
#' (modes, gaps, skewness) at a glance.
#'
#' @param x a \code{dvdata} object, a \code{MatH}, a \code{MatD} or a
#'   \code{data.frame} (passed to \code{\link{dv_data}} with the arguments in
#'   \code{...}).
#' @param row the unit to draw (index or name).
#' @param petal \code{"auto"} (ellipses for discrete distributions, violins
#'   for histograms), \code{"ellipse"} or \code{"violin"}.
#' @param petal_width the width of a petal is proportional to its
#'   probability; \code{petal_width} is the half-width (as a fraction of the
#'   half-sector) of the largest petal: with \code{width_scale = "common"} the
#'   largest probability among the drawn units, with \code{"unit"} the largest
#'   probability of each unit.
#' @param width_scale which petal has the width \code{petal_width}:
#'   \code{"common"} (default) the largest petal among all the drawn units
#'   and variables, so that all the widths are comparable; \code{"unit"} the
#'   largest petal of each unit; \code{"variable"} the largest petal of each
#'   variable among the drawn units, so that widths are comparable across
#'   units for the same variable (useful for histograms of variables whose
#'   distributions are very different: with a common scale the most
#'   concentrated variable would make all the other petals thin).
#' @param width_type how the width of the petals is measured:
#'   \code{"angular"} (default, as in the original design: the half-width is
#'   an angle, so that a petal far from the centre is wider than a petal of
#'   the same probability near the centre) or \code{"linear"} (the half-width
#'   is a distance from the stem, the same at every radius: the width is
#'   exactly proportional to the probability, and the area of the petals is
#'   comparable along the stem).
#' @param center \code{"flower"} (green receptacle with orange sepals, as in
#'   the original design) or \code{"plain"} (a grey disk).
#' @param center_size radius of the centre (the radius of the flower is 1).
#' @param show_median logical; mark the median of each distribution on the
#'   stem?
#' @param n_levels for histograms: number of classes (see \code{\link{dv_data}});
#'   more classes give smoother violins.
#' @param smooth for violin petals: standard deviation (in classes) of the
#'   Gaussian kernel used to smooth the profile of the petals (0 = no
#'   smoothing). The colours are not smoothed.
#' @param style \code{"dark"} (default, black background) or \code{"light"}.
#' @param alpha opacity of the petals.
#' @inheritParams gei_plot
#' @return a \code{ggplot} object.
#' @examples
#' dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
#' flower_plot(dv, row = "Hotel_01")
#' flower_plot(dv, row = 2, var_order = "tsp", style = "light", show_median = TRUE)
#' flower_plot(HistDAWass::BLOOD, row = 1, palette = "Zissou 1")
#' @export
flower_plot <- function(x, row = 1, palette = "Paquin", reverse = FALSE,
                        petal = c("auto", "ellipse", "violin"), petal_width = 0.9,
                        width_scale = c("common", "unit", "variable"), width_type = c("angular", "linear"),
                        center = c("flower", "plain"),
                        center_size = 0.22, show_median = FALSE, var_order = NULL,
                        labels = TRUE, title = NULL, style = c("dark", "light"),
                        alpha = 0.95, n_levels = 30, smooth = 1.5, ...) {
  x <- .flower_data(x, n_levels, ...)
  if (is.character(row)) row <- match(row, x$rows)
  if (length(row) != 1 || is.na(row)) stop("row must identify a single unit")
  g <- flower_grid(x, rows = row, palette = palette, reverse = reverse, petal = petal,
                   petal_width = petal_width, width_scale = width_scale, width_type = width_type,
                   center = center,
                   center_size = center_size, show_median = show_median, var_order = var_order,
                   labels = labels, style = style, alpha = alpha, strip = FALSE, smooth = smooth)
  style <- match.arg(style)
  g + ggplot2::labs(title = if (is.null(title)) x$rows[row] else title) +
    ggplot2::theme(plot.title = ggplot2::element_text(
      hjust = 0.5, colour = if (style == "dark") "grey90" else "grey15"))
}

#' Flower plots of several units
#'
#' Draws the flower plots (see \code{\link{flower_plot}}) of several units
#' in a grid, optionally ordered or grouped as in \code{\link{gei_grid}}.
#' @inheritParams flower_plot
#' @inheritParams gei_grid
#' @return a \code{ggplot} object.
#' @examples
#' dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
#' flower_grid(dv, rows = 1:6, ncol = 3, var_order = "tsp")
#' @export
flower_grid <- function(x, rows = NULL, order = c("data", "cluster"), groups = NULL,
                        ncol = NULL, palette = "Paquin", reverse = FALSE,
                        petal = c("auto", "ellipse", "violin"), petal_width = 0.9,
                        width_scale = c("common", "unit", "variable"), width_type = c("angular", "linear"),
                        center = c("flower", "plain"),
                        center_size = 0.22, show_median = FALSE, var_order = NULL,
                        labels = FALSE, style = c("dark", "light"), alpha = 0.95,
                        strip = TRUE, n_levels = 30, smooth = 1.5, ...) {
  .flower_core(x = x, rows = rows, order = order, groups = groups, ncol = ncol,
               palette = palette, reverse = reverse, petal = petal, petal_width = petal_width,
               width_scale = width_scale, width_type = width_type, center = center,
               center_size = center_size, show_median = show_median, var_order = var_order,
               labels = labels, style = style, alpha = alpha, strip = strip,
               n_levels = n_levels, smooth = smooth, ...)
}

.flower_core <- function(x, rows = NULL, order = c("data", "cluster"), groups = NULL,
                        ncol = NULL, palette = "Paquin", reverse = FALSE,
                        petal = c("auto", "ellipse", "violin"), petal_width = 0.9,
                        width_scale = c("common", "unit", "variable"), width_type = c("angular", "linear"),
                        center = c("flower", "plain"),
                        center_size = 0.22, show_median = FALSE, var_order = NULL,
                        labels = FALSE, style = c("dark", "light"), alpha = 0.95,
                        strip = TRUE, n_levels = 30, smooth = 1.5, place = NULL, ...) {
  x <- .flower_data(x, n_levels, ...)
  order <- match.arg(order); petal <- match.arg(petal); width_scale <- match.arg(width_scale)
  width_type <- match.arg(width_type)
  center <- match.arg(center); style <- match.arg(style)
  if (petal == "auto") petal <- if (is(x$source, "MatH")) "violin" else "ellipse"
  nr <- length(x$rows); p <- length(x$vars)
  if (is.null(rows)) rows <- seq_len(nr)
  if (is.character(rows)) rows <- match(rows, x$rows)
  if (anyNA(rows)) stop("unknown rows")
  if (order == "cluster") {
    ro <- dv_cluster(x, rows = TRUE, cols = FALSE)$row_order
    rows <- ro[ro %in% rows]
  }
  grp_lab <- NULL
  if (!is.null(groups)) {
    if (length(groups) != nr) stop("groups must have one value per row of the data")
    rows <- rows[order(groups[rows], match(rows, rows))]
    grp_lab <- as.character(groups)
  }
  vo <- .var_order(x, var_order)
  facet_lab <- if (is.null(grp_lab)) x$rows else paste0("[", grp_lab, "] ", x$rows)
  fg <- if (style == "dark") "grey85" else "grey25"
  bgc <- if (style == "dark") "black" else "white"
  stem_col <- if (style == "dark") "grey55" else "grey65"

  df <- x$df[x$df$IDr %in% rows & x$df$freq > 0, ]
  df$pos <- match(df$IDc, vo)
  r0 <- center_size
  rr <- function(u) r0 + (1 - r0) * u                 # radius of a value
  half <- pi / p                                      # half-sector (angle)
  if (petal == "violin") {
    prof <- lapply(split(x$df[x$df$IDr %in% rows, ], list(x$df$IDr[x$df$IDr %in% rows],
                                                          x$df$IDc[x$df$IDr %in% rows]), drop = TRUE),
                   function(f) {
                     f <- f[order(f$u), ]
                     f$fs <- .smooth_profile(f$freq, smooth)
                     f
                   })
    allf <- do.call(rbind, prof)
  }
  # largest petal used to scale the widths: overall, of each unit or of each
  # variable (fmax(i, j) for unit i and variable j)
  ref <- if (petal == "violin") data.frame(IDr = allf$IDr, IDc = allf$IDc, f = allf$fs) else
    data.frame(IDr = df$IDr, IDc = df$IDc, f = df$freq)
  fmax <- switch(width_scale,
                 common = function(i, j) rep(max(ref$f), length(i)),
                 unit = { m <- tapply(ref$f, ref$IDr, max); function(i, j) m[as.character(i)] },
                 variable = { m <- tapply(ref$f, ref$IDc, max); function(i, j) m[as.character(j)] })
  df$a <- petal_width * half * df$freq / fmax(df$IDr, df$IDc)   # half-width (as an angle)
  df$theta <- (df$pos - 0.5) / p * 2 * pi            # direction of the stem

  # polar (radius, angle) -> cartesian, angle clockwise from the top
  pc <- function(r, th) data.frame(x = r * sin(th), y = r * cos(th))
  # a point of a petal: radius r on the stem of direction th, half-width a
  # (angular: a is an angle; linear: a is converted into a constant distance
  # from the stem, perpendicular to it)
  wmax <- 0.6 * sin(min(half, pi / 2))   # linear half-width of the widest petal
  pw <- if (width_type == "angular") {
    function(r, th, a) pc(r, th + a)
  } else {
    function(r, th, a) {
      w <- a / half * wmax
      data.frame(x = r * sin(th) + w * cos(th), y = r * cos(th) - w * sin(th))
    }
  }
  polys <- list()
  if (petal == "ellipse") {
    # radial semi-axis: 45% of the smallest gap between consecutive values
    allu <- tapply(x$df$u, x$df$IDc, function(u) { d <- diff(sort(unique(u))); if (length(d)) min(d) else 0.5 })
    b <- 0.45 * (1 - r0) * allu[as.character(df$IDc)]
    phi <- seq(0, 2 * pi, length.out = 41)
    k <- 0
    for (i in seq_len(nrow(df))) {
      e <- pw(rr(df$u[i]) + b[i] * cos(phi), df$theta[i], df$a[i] * sin(phi))
      k <- k + 1
      polys[[k]] <- data.frame(id = i, e, u = df$u[i], IDr = df$IDr[i])
    }
  } else {
    # violin: slices between consecutive class centres, width interpolated
    k <- 0
    for (f in prof) {
      if (!any(f$freq > 0)) next
      a <- petal_width * half * f$fs / fmax(f$IDr[1], f$IDc[1])
      vb <- .violin_ends(f$u, a)
      th <- (match(f$IDc[1], vo) - 0.5) / p * 2 * pi
      for (s in seq_len(length(vb$u) - 1)) {
        if (vb$a[s] < 1e-9 && vb$a[s + 1] < 1e-9) next
        r1 <- rr(vb$u[s]); r2 <- rr(vb$u[s + 1])
        sl <- pw(c(r1, r2, r2, r1), th, c(-vb$a[s], -vb$a[s + 1], vb$a[s + 1], vb$a[s]))
        k <- k + 1
        polys[[k]] <- data.frame(id = paste(f$IDr[1], f$IDc[1], s), sl,
                                 u = (vb$uc[s] + vb$uc[s + 1]) / 2, IDr = f$IDr[1])
      }
    }
  }
  pet <- do.call(rbind, polys)
  pet$facet <- factor(facet_lab[pet$IDr], levels = facet_lab[rows])

  # stems: from the centre to the largest value with positive probability
  st <- do.call(rbind, lapply(split(df, list(df$IDr, df$IDc), drop = TRUE), function(d) {
    e <- pc(c(r0, rr(max(d$u))), d$theta[1])
    data.frame(x = e$x[1], y = e$y[1], xend = e$x[2], yend = e$y[2], IDr = d$IDr[1])
  }))
  st$facet <- factor(facet_lab[st$IDr], levels = facet_lab[rows])

  # outline of the violins
  ol <- NULL
  if (petal == "violin") {
    ol <- do.call(rbind, lapply(prof, function(f) {
      a <- petal_width * half * f$fs / fmax(f$IDr[1], f$IDc[1])
      vb <- .violin_ends(f$u, a)
      th <- (match(f$IDc[1], vo) - 0.5) / p * 2 * pi
      e <- pw(c(rr(vb$u), rev(rr(vb$u))), th, c(-vb$a, rev(vb$a)))
      data.frame(id = paste(f$IDr[1], f$IDc[1]), e, IDr = f$IDr[1])
    }))
    ol$facet <- factor(facet_lab[ol$IDr], levels = facet_lab[rows])
  }
  # centre of the flower (and sepals)
  thc <- seq(0, 2 * pi, length.out = 121)
  cen <- data.frame(x = rep(r0 * sin(thc), length(rows)), y = rep(r0 * cos(thc), length(rows)),
                    IDr = rep(rows, each = length(thc)),
                    facet = factor(rep(facet_lab[rows], each = length(thc)), levels = facet_lab[rows]))
  sep <- NULL
  if (center == "flower") {
    phi <- seq(0, 2 * pi, length.out = 41)
    sep1 <- do.call(rbind, lapply(seq_len(p), function(j) {
      th <- (j - 0.5) / p * 2 * pi
      e <- pc(r0 / 2 + r0 * 0.42 * cos(phi), th + min(half, 0.6) * 0.75 * sin(phi))
      data.frame(id = j, e)
    }))
    sep <- do.call(rbind, lapply(rows, function(i) { s1 <- sep1; s1$id <- paste(i, s1$id); s1$IDr <- i
      s1$facet <- factor(facet_lab[i], levels = facet_lab[rows]); s1 }))
  }
  md <- NULL
  if (show_median) {
    md <- do.call(rbind, lapply(split(df, list(df$IDr, df$IDc), drop = TRUE), function(d) {
      d <- d[order(d$u), ]
      um <- d$u[which(d$cdf >= 0.5)[1]]
      e <- pc(rr(um), d$theta[1])
      data.frame(x = e$x, y = e$y, IDr = d$IDr[1])
    }))
    md$facet <- factor(facet_lab[md$IDr], levels = facet_lab[rows])
  }
  cols <- dv_palette(palette, 64, reverse)
  # flowers placed at given positions (flower_pca): move and scale the
  # geometry and return the layers
  if (!is.null(place)) {
    return(list(st = place(st), pet = place(pet), ol = place(ol), cen = place(cen),
                sep = place(sep), md = place(md), cols = cols, petal = petal, alpha = alpha,
                center = center, fg = fg, stem_col = stem_col, vo = vo))
  }
  g <- .flower_layers(ggplot2::ggplot(), st, pet, ol, cen, sep, md, cols, petal, alpha, center,
                      fg, stem_col, group_var = "facet")
  lim <- if (labels) {
    thm <- (seq_len(p) - 0.5) / p * 2 * pi
    lb <- data.frame(x = 1.07 * sin(thm), y = 1.07 * cos(thm), lab = x$vars[vo],
                     hjust = ifelse(sin(thm) > 0.15, 0, ifelse(sin(thm) < -0.15, 1, 0.5)),
                     vjust = ifelse(cos(thm) > 0.15, 0, ifelse(cos(thm) < -0.15, 1, 0.5)))
    g <- g + ggplot2::geom_text(data = lb, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab,
                                                        hjust = .data$hjust, vjust = .data$vjust),
                                size = 2.8, colour = fg)
    w <- 1.1 + 0.075 * max(nchar(x$vars))
    list(x = c(-w, w), y = c(-1.2, 1.2))
  } else list(x = c(-1.05, 1.05), y = c(-1.05, 1.05))
  g <- g + ggplot2::coord_equal(xlim = lim$x, ylim = lim$y, expand = FALSE, clip = "off")
  if (length(rows) > 1 || strip) {
    g <- g + ggplot2::facet_wrap(ggplot2::vars(.data$facet), ncol = ncol, strip.position = "bottom")
  }
  g + ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bgc, colour = bgc),
      panel.background = ggplot2::element_rect(fill = bgc, colour = NA),
      strip.text = if (strip) ggplot2::element_text(size = 9, colour = fg, margin = ggplot2::margin(2, 0, 4, 0))
                   else ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8, colour = fg),
      panel.spacing = ggplot2::unit(0.6, "lines"),
      plot.margin = ggplot2::margin(6, 6, 6, 6)) +
    ggplot2::labs(caption = paste("Petals, clockwise from the top:", paste(x$vars[vo], collapse = ", ")))
}

# dv_data with more classes by default for histograms
.flower_data <- function(x, n_levels, ...) {
  if (inherits(x, "dvdata")) return(x)
  if (is(x, "MatH")) return(dv_data(x, n_levels = n_levels, ...))
  dv_data(x, ...)
}

# Gaussian smoothing of a profile (sd in number of classes), keeping the total
.smooth_profile <- function(f, sd) {
  n <- length(f)
  if (sd <= 0 || n < 3) return(f)
  h <- ceiling(3 * sd)
  w <- stats::dnorm(-h:h, sd = sd)
  out <- vapply(seq_len(n), function(i) {
    idx <- i + (-h:h); ok <- idx >= 1 & idx <= n
    sum(f[idx[ok]] * w[ok]) / sum(w[ok])
  }, 0)
  out * sum(f) / sum(out)
}

# close a violin at both ends (half a class beyond the first/last class)
.violin_ends <- function(u, a) {
  n <- length(u)
  du <- if (n > 1) diff(u) else 0.05
  u2 <- c(max(u[1] - 0.5 * du[1], 0), u, min(u[n] + 0.5 * du[length(du)], 1))
  list(u = u2, a = c(0, a, 0), uc = c(u[1], u, u[n]))
}

# the layers of the flowers (shared by flower_grid and flower_pca)
.flower_layers <- function(g, st, pet, ol, cen, sep, md, cols, petal, alpha, center, fg, stem_col,
                           group_var = "facet", pt_size = 1.6, lw = 1, stroke_map = TRUE) {
  # stroke_map: the border of the petals takes their colour (hides the seams
  # between the slices of the violins); FALSE leaves the colour scale free
  pmap <- if (stroke_map) {
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$id, fill = .data$u, colour = .data$u)
  } else ggplot2::aes(x = .data$x, y = .data$y, group = .data$id, fill = .data$u)
  g <- g +
    ggplot2::geom_segment(data = st, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                                                  yend = .data$yend), colour = stem_col,
                          linewidth = 0.4 * lw) +
    (if (stroke_map) {
      ggplot2::geom_polygon(data = pet, pmap, linewidth = if (petal == "violin") 0.15 * lw else 0,
                            alpha = alpha)
    } else ggplot2::geom_polygon(data = pet, pmap, colour = NA, alpha = alpha)) +
    ggplot2::scale_fill_gradientn(colours = cols, limits = c(0, 1), guide = "none",
                                  aesthetics = if (stroke_map) c("fill", "colour") else "fill")
  if (!is.null(ol)) {
    g <- g + ggplot2::geom_path(data = ol, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
                                colour = fg, linewidth = 0.25 * lw, alpha = 0.7)
  }
  g <- g + ggplot2::geom_polygon(data = cen, ggplot2::aes(x = .data$x, y = .data$y,
                                                          group = .data[[group_var]]),
                                 fill = if (center == "flower") "#2e7d32" else "grey40", colour = NA)
  if (!is.null(sep)) {
    g <- g + ggplot2::geom_polygon(data = sep, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
                                   fill = "darkorange", colour = "#8d3b12", linewidth = 0.25 * lw)
  }
  if (!is.null(md)) {
    g <- g + ggplot2::geom_point(data = md, ggplot2::aes(x = .data$x, y = .data$y),
                                 shape = 21, size = pt_size, fill = "white", colour = "black",
                                 stroke = 0.4)
  }
  g
}

