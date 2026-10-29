## Difference eyes ----------------------------------------------------------

#' Difference eyes: departures of the units from a reference distribution
#'
#' A variant of the Green Eyes plot (see \code{\link{gei_grid}}) that shows
#' \emph{where} the distributions of a unit depart from a reference. Each
#' variable is a sector and, as in a Green Eyes plot, the radius is the
#' probability level \eqn{t} (from \eqn{t = 0} at the pupil to \eqn{t = 1} at
#' the border); the colour is the difference between the quantile function
#' of the unit and that of the reference,
#' \deqn{D_{ij}(t) = Q_{ij}(t) - Q^{ref}_j(t),}
#' red where the unit is below the reference, green where it is above, white
#' where they coincide. By default the reference is the Wasserstein
#' barycenter of all the units, whose quantile function is the mean of the
#' quantile functions: \eqn{\int_0^1 D_{ij}(t)^2 dt} is then the contribution
#' of the cell to the Wasserstein variance of variable \eqn{j}, so a
#' difference eye is a picture of how each unit contributes to the
#' variability of the data.
#'
#' Quantile functions are computed exactly from the original data (step
#' functions for discrete distributions, piecewise linear functions for
#' histograms) on a grid of \code{levels} probability levels.
#'
#' @inheritParams gei_grid
#' @param reference the reference distribution: \code{"mean"} (default, the
#'   Wasserstein barycenter of all the units), the name or index of a unit,
#'   or a vector of names/indices of units (their barycenter, e.g. a group).
#' @param scale how the differences are expressed: \code{"raw"} (in the
#'   units of the variables), \code{"range"} (divided by the range of each
#'   variable), \code{"sd"} (divided by the Wasserstein standard deviation of
#'   each variable). By default \code{"raw"} if all the variables have the
#'   same range, \code{"sd"} otherwise.
#' @param limit (optional) the colour scale spans
#'   \eqn{[-}\code{limit}\eqn{, }\code{limit}\eqn{]}; by default the largest
#'   absolute difference in the drawn units.
#' @param palette a diverging palette: a vector of colours (low to high) or a
#'   palette name (see \code{\link{dv_palette}}); by default red - white -
#'   green, consistent with the Paquin palette.
#' @param levels number of probability levels of the grid.
#' @param legend logical; draw the colour legend?
#' @return a \code{ggplot} object, with the array of the differences
#'   (units x variables x levels) in the attribute \code{"diff"}.
#' @examples
#' dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
#' gei_diff(dv, rows = 1:8, ncol = 4)
#' gei_diff(dv, rows = 1:4, reference = "Hotel_01", ncol = 4)
#' gei_diff(HistDAWass::BLOOD, rows = 1:6, ncol = 3)
#' @export
gei_diff <- function(x, rows = NULL, reference = "mean", scale = NULL, order = c("data", "cluster"),
                     groups = NULL, ncol = NULL, var_order = NULL, limit = NULL,
                     palette = NULL, reverse = FALSE, pupil = 0.3, pupil_color = "grey40",
                     median_line = TRUE, levels = 200, labels = FALSE, legend = TRUE,
                     strip = TRUE, bg = "transparent", ...) {
  x <- dv_data(x, ...)
  order <- match.arg(order)
  nr <- length(x$rows); p <- length(x$vars)
  if (is.null(rows)) rows <- seq_len(nr)
  if (is.character(rows)) rows <- match(rows, x$rows)
  if (anyNA(rows)) stop("unknown rows")

  # quantile functions on a grid of levels (mid points of the bands)
  tg <- (seq_len(levels) - 0.5) / levels
  Q <- .quantile_array(x$source, tg)                       # n x p x levels
  # reference
  if (identical(reference, "mean")) {
    ref_rows <- seq_len(nr); ref_lab <- "the barycenter of all the units"
  } else {
    ref_rows <- if (is.character(reference)) match(reference, x$rows) else as.integer(reference)
    if (anyNA(ref_rows) || any(ref_rows < 1 | ref_rows > nr)) stop("unknown reference units")
    ref_lab <- if (length(ref_rows) == 1) x$rows[ref_rows] else
      sprintf("the barycenter of %d units", length(ref_rows))
  }
  Qref <- apply(Q[ref_rows, , , drop = FALSE], c(2, 3), mean)   # barycenter: mean of quantiles
  D <- sweep(Q, c(2, 3), Qref)

  # scale
  rg <- .var_ranges(x)
  if (is.null(scale)) scale <- if (nrow(unique(round(rg, 10))) == 1) "raw" else "sd"
  scale <- match.arg(scale, c("raw", "range", "sd"))
  den <- switch(scale,
                raw = rep(1, p),
                range = pmax(rg[, 2] - rg[, 1], 1e-12),
                sd = {
                  sdv <- sqrt(apply(sweep(Q, c(2, 3), apply(Q, c(2, 3), mean))^2, 2, mean))
                  pmax(sdv, 1e-12)
                })
  D <- sweep(D, 2, den, "/")
  dimnames(D) <- list(x$rows, x$vars, NULL)
  if (is.null(limit)) limit <- max(abs(D[rows, , ]), 1e-12)

  # order of the units and of the sectors (as in gei_grid)
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

  # bands: consecutive levels with the same difference are merged
  lev <- c(0, seq_len(levels) / levels)
  step <- limit / 50
  bands <- vector("list", length(rows) * p); k <- 0
  for (i in rows) for (j in seq_len(p)) {
    # differences quantized to the resolution of the colour scale (101 levels)
    v <- round(pmin(pmax(D[i, j, ], -limit), limit) / step) * step
    r <- rle(v); hi <- cumsum(r$lengths); lo <- c(0, hi[-length(hi)])
    k <- k + 1
    bands[[k]] <- data.frame(IDr = i, xpos = match(j, vo), ymin = lev[lo + 1], ymax = lev[hi + 1],
                             d = r$values)
  }
  df <- do.call(rbind, bands)
  df$facet <- factor(facet_lab[df$IDr], levels = facet_lab[rows])

  cols <- if (is.null(palette)) {
    c("#831818", "#c62320", "#f78462", "#f7f2ec", "#98ab76", "#47632a", "#275024")
  } else if (length(palette) > 1) palette else dv_palette(palette, 64)
  if (reverse) cols <- rev(cols)
  leg_name <- switch(scale, raw = "Q(t) - Qref(t)", range = "(Q(t) - Qref(t)) / range",
                     sd = "(Q(t) - Qref(t)) / sd")

  # eyes drawn as polygons in cartesian coordinates (much faster than
  # coord_polar with many bands); radius 1, pupil as in gei_plot
  rad <- function(y) (y + pupil) / (1 + pupil)
  npt <- max(2, ceiling(120 / p))
  poly <- .band_polygons(cx = rep(0, nrow(df)), cy = rep(0, nrow(df)),
                         th1 = (df$xpos - 1) / p * 2 * pi, th2 = df$xpos / p * 2 * pi,
                         r1 = rad(df$ymin), r2 = rad(df$ymax), npt = npt)
  poly$d <- df$d[poly$id]; poly$facet <- df$facet[poly$id]
  th <- seq(0, 2 * pi, length.out = 121)
  circ <- function(r) {
    data.frame(x = rep(r * sin(th), length(rows)), y = rep(r * cos(th), length(rows)),
               facet = factor(rep(facet_lab[rows], each = length(th)), levels = facet_lab[rows]))
  }
  ths <- (seq_len(p) - 1) / p * 2 * pi
  seps <- data.frame(x = rad(0) * sin(ths), y = rad(0) * cos(ths), xend = sin(ths), yend = cos(ths))
  g <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = poly, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                                    fill = .data$d, colour = .data$d), linewidth = 0.1) +
    ggplot2::geom_polygon(data = circ(rad(0)), ggplot2::aes(x = .data$x, y = .data$y, group = .data$facet),
                          fill = pupil_color, colour = NA) +
    ggplot2::geom_segment(data = seps, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                                                    yend = .data$yend), colour = "white", linewidth = 0.6) +
    ggplot2::scale_fill_gradientn(colours = cols, limits = c(-limit, limit),
                                  oob = function(v, ...) pmin(pmax(v, -limit), limit),
                                  aesthetics = c("fill", "colour"), name = leg_name,
                                  guide = if (legend) ggplot2::guide_colourbar(
                                    barwidth = ggplot2::unit(10, "lines"), title.position = "top",
                                    title.hjust = 0.5) else "none")
  if (median_line) {
    g <- g + ggplot2::geom_path(data = circ(rad(0.5)), ggplot2::aes(x = .data$x, y = .data$y,
                                                                    group = .data$facet),
                                linetype = "dashed", alpha = 0.6, linewidth = 0.3)
  }
  if (labels) {
    thm <- (seq_len(p) - 0.5) / p * 2 * pi
    lb <- data.frame(x = 1.08 * sin(thm), y = 1.08 * cos(thm), lab = x$vars[vo],
                     hjust = ifelse(sin(thm) > 0.15, 0, ifelse(sin(thm) < -0.15, 1, 0.5)),
                     vjust = ifelse(cos(thm) > 0.15, 0, ifelse(cos(thm) < -0.15, 1, 0.5)))
    g <- g + ggplot2::geom_text(data = lb, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab,
                                                        hjust = .data$hjust, vjust = .data$vjust),
                                size = 2.8, colour = "grey20")
  }
  lim <- if (labels) {
    # room for the names of the variables
    w <- 1.1 + 0.075 * max(nchar(x$vars))
    list(x = c(-w, w), y = c(-1.22, 1.22))
  } else list(x = c(-1.02, 1.02), y = c(-1.02, 1.02))
  g <- g + ggplot2::coord_equal(xlim = lim$x, ylim = lim$y, expand = FALSE, clip = "off")
  if (length(rows) > 1 || strip) {
    g <- g + ggplot2::facet_wrap(ggplot2::vars(.data$facet), ncol = ncol, strip.position = "bottom")
  }
  g <- g + .gei_theme(bg, FALSE, strip) +
    ggplot2::theme(legend.position = if (legend) "bottom" else "none",
                   plot.caption = ggplot2::element_text(hjust = 0, size = 8, colour = "grey30")) +
    ggplot2::labs(caption = paste0("Reference: ", ref_lab, ". Sectors, clockwise from the top: ",
                                   paste(x$vars[vo], collapse = ", "), "."))
  if (length(rows) == 1 && !strip) g <- g + ggplot2::labs(title = x$rows[rows]) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  attr(g, "diff") <- D
  g
}

# array (units x variables x levels) of quantile values
.quantile_array <- function(src, tg) {
  M <- src@M
  n <- nrow(M); p <- ncol(M)
  Q <- array(NA_real_, c(n, p, length(tg)))
  hist <- is(src, "MatH")
  for (i in seq_len(n)) for (j in seq_len(p)) {
    d <- M[i, j][[1]]
    Q[i, j, ] <- if (hist) {
      stats::approx(d@p, d@x, xout = tg, ties = "ordered", rule = 2)$y
    } else {
      cdf <- cumsum(d@p)
      d@x[pmin(findInterval(tg - 1e-12, cdf) + 1, length(d@x))]
    }
  }
  Q
}

# ranges of the variables (for the scale of the differences)
.var_ranges <- function(x) {
  if (!is.null(x$ranges) && x$scheme != "standardized") return(as.matrix(x$ranges))
  src <- x$source@M
  t(vapply(seq_len(ncol(src)), function(j) {
    r <- range(unlist(lapply(src[, j], function(d) d@x)))
    r
  }, numeric(2)))
}
