## Green Eyes (GEI) plots ---------------------------------------------------

#' Green Eyes plot of a statistical unit
#'
#' A Green Eyes plot (Green Eye Iris plot, Crocetta et al., 2025) shows all the distributions observed on
#' a statistical unit as a circular plot: each variable is a sector, and the
#' sector is filled from the inner circle (the "pupil") outwards with the
#' cumulative distribution of the variable, each value being coloured
#' according to its position in the range of the variable (low values in
#' red, high values in green with the default palette). The dashed circle
#' marks the median level (cdf = 0.5): if the colour crossing it is green the
#' median of the variable is high. A unit with "green eyes" has high values
#' in all the variables.
#'
#' @param x a \code{dvdata} object, a \code{MatH}, a \code{MatD} or a
#'   \code{data.frame} (passed to \code{\link{dv_data}} with the arguments in
#'   \code{...}).
#' @param row the unit to draw (index or name).
#' @param palette palette name or vector of colours (see
#'   \code{\link{dv_palette}}).
#' @param reverse logical; reverse the palette?
#' @param pupil size of the pupil (relative to the radius of the cdf
#'   ring, default 0.3).
#' @param pupil_color colour of the pupil.
#' @param median_line logical; draw the dashed circle at cdf = 0.5?
#' @param show_sd logical; draw, for each variable, a bar whose length is
#'   proportional to the standard deviation of the distribution (relative to
#'   the largest standard deviation observed for that variable).
#' @param show_skewness logical; draw a tick at a radial position related to
#'   the skewness of the distribution (\eqn{(\sqrt[3]{\gamma} + 3)/6}).
#' @param var_order order of the sectors (clockwise from the top):
#'   \code{NULL} (the order of the data), a vector of indices/names,
#'   \code{"cluster"} (order of a hierarchical clustering of the variables,
#'   see \code{\link{dv_cluster}}) or \code{"tsp"} (alias
#'   \code{"hamiltonian"}: the shortest Hamiltonian cycle through the
#'   variables with the correlation distance \eqn{\sqrt{2(1-r)}}, so that
#'   correlated variables are adjacent, see \code{\link{dv_var_order}}).
#' @param labels logical; write the variable names around the eye?
#' @param title (optional) title; by default the name of the unit.
#' @param polar logical; if \code{FALSE} the plot is drawn as stacked bars
#'   (without the pupil).
#' @param bg background colour.
#' @param ... arguments passed to \code{\link{dv_data}}.
#' @return a \code{ggplot} object.
#' @references Crocetta C., Irpino A., Antonucci L., Marin C. (2025). A
#'   performance indicator and its decomposition according to the impacts of
#'   different aspects based on distributional data. \emph{Quality & Quantity},
#'   59(3), 2091--2110. \doi{10.1007/s11135-024-01956-0}
#' @examples
#' dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
#' gei_plot(dv, row = "Hotel_01")
#' gei_plot(dv, row = 2, show_sd = TRUE, show_skewness = TRUE)
#' gei_plot(HistDAWass::BLOOD, row = 3, n_levels = 15, palette = "Zissou 1")
#' @export
gei_plot <- function(x, row = 1, palette = "Paquin", reverse = FALSE, pupil = 0.3,
                     pupil_color = "grey40", median_line = TRUE, show_sd = FALSE,
                     show_skewness = FALSE, var_order = NULL, labels = TRUE,
                     title = NULL, polar = TRUE, bg = "transparent", ...) {
  x <- dv_data(x, ...)
  if (is.character(row)) row <- match(row, x$rows)
  if (length(row) != 1 || is.na(row)) stop("row must identify a single unit")
  g <- gei_grid(x, rows = row, palette = palette, reverse = reverse, pupil = pupil,
                pupil_color = pupil_color, median_line = median_line, show_sd = show_sd,
                show_skewness = show_skewness, var_order = var_order, labels = labels,
                polar = polar, bg = bg, strip = FALSE)
  g + ggplot2::labs(title = if (is.null(title)) x$rows[row] else title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Green Eyes plots of several units
#'
#' Draws the Green Eyes plots (see \code{\link{gei_plot}}) of several units
#' in a grid (one facet per unit), optionally ordered or grouped by a
#' clustering of the units.
#' @inheritParams gei_plot
#' @param rows the units to draw (default all).
#' @param order \code{"data"} (default) or \code{"cluster"}: order the
#'   units according to a hierarchical clustering (Wasserstein distance).
#' @param groups (optional) a vector of group labels for the units (e.g. the
#'   output of \code{cutree} or of \code{DDAWASS::WD_kmeans}); units are
#'   then sorted by group and the group is reported in the facet label.
#' @param ncol number of columns of the grid.
#' @param strip logical; write the unit names under each eye?
#' @return a \code{ggplot} object.
#' @examples
#' dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
#' gei_grid(dv, rows = 1:8, ncol = 4)
#' hc <- DDAWASS::WD_hclust(DDAWASS::Satisfaction, method = "ward.D2")
#' gei_grid(dv, rows = 1:12, groups = cutree(hc, 3), ncol = 6)
#' @export
gei_grid <- function(x, rows = NULL, order = c("data", "cluster"), groups = NULL,
                     ncol = NULL, palette = "Paquin", reverse = FALSE, pupil = 0.3,
                     pupil_color = "grey40", median_line = TRUE, show_sd = FALSE,
                     show_skewness = FALSE, var_order = NULL, labels = FALSE,
                     polar = TRUE, bg = "transparent", strip = TRUE, ...) {
  x <- dv_data(x, ...)
  order <- match.arg(order)
  if (!polar) pupil <- 0
  nr <- length(x$rows); p <- length(x$vars)
  if (is.null(rows)) rows <- seq_len(nr)
  if (is.character(rows)) rows <- match(rows, x$rows)
  # order of units
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
  # order of variables
  vo <- .var_order(x, var_order)
  df <- x$df[x$df$IDr %in% rows, ]
  df$xpos <- match(df$IDc, vo)
  facet_lab <- if (is.null(grp_lab)) x$rows else paste0("[", grp_lab, "] ", x$rows)
  df$facet <- factor(facet_lab[df$IDr], levels = facet_lab[rows])
  cols <- dv_palette(palette, 64, reverse)

  g <- ggplot2::ggplot(df) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$xpos - 0.5, xmax = .data$xpos + 0.5,
                                    ymin = .data$cdf_lo, ymax = .data$cdf, fill = .data$u),
                       colour = NA) +
    ggplot2::annotate("rect", xmin = 0.5, xmax = p + 0.5, ymin = -pupil, ymax = 0,
                      fill = pupil_color, colour = NA) +
    ggplot2::scale_fill_gradientn(colours = cols, limits = c(0, 1), guide = "none")
  if (median_line) {
    g <- g + ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.6)
  }
  if (show_sd || show_skewness) {
    st <- .cell_stats(x, rows, vo, facet_lab)
    if (show_sd) {
      g <- g + ggplot2::geom_segment(data = st, ggplot2::aes(x = .data$xpos, xend = .data$xpos,
                                                             y = 0.5 - .data$sdprop / 2,
                                                             yend = 0.5 + .data$sdprop / 2),
                                     colour = "grey15", linewidth = 1.2)
    }
    if (show_skewness) {
      g <- g + ggplot2::geom_segment(data = st, ggplot2::aes(x = .data$xpos - 0.3, xend = .data$xpos + 0.3,
                                                             y = .data$skpos, yend = .data$skpos),
                                     colour = "black", linewidth = 0.9)
    }
  }
  g <- g + ggplot2::scale_x_continuous(breaks = seq_len(p), labels = x$vars[vo],
                                       limits = c(0.5, p + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-pupil, 1), expand = c(0, 0))
  if (polar) g <- g + ggplot2::coord_polar(theta = "x", clip = "off")
  if (length(rows) > 1 || strip) {
    g <- g + ggplot2::facet_wrap(ggplot2::vars(.data$facet), ncol = ncol, strip.position = "bottom")
  }
  g + .gei_theme(bg, labels, strip)
}

.var_order <- function(x, var_order) {
  p <- length(x$vars)
  if (is.null(var_order)) return(seq_len(p))
  if (identical(var_order, "cluster")) return(dv_cluster(x, rows = FALSE, cols = TRUE)$col_order)
  if (length(var_order) == 1 && var_order %in% c("tsp", "hamiltonian")) return(dv_var_order(x)$order)
  vo <- if (is.character(var_order)) match(var_order, x$vars) else var_order
  if (anyNA(vo) || length(vo) != p || anyDuplicated(vo)) stop("var_order must be a permutation of the variables")
  vo
}

.cell_stats <- function(x, rows, vo, facet_lab) {
  sd <- x$stats$std; sk <- x$stats$skewness
  maxsd <- apply(sd, 2, max, na.rm = TRUE)
  maxsd[!is.finite(maxsd) | maxsd == 0] <- 1
  out <- expand.grid(IDr = rows, IDc = seq_along(x$vars))
  out$xpos <- match(out$IDc, vo)
  out$sdprop <- sd[cbind(out$IDr, out$IDc)] / maxsd[out$IDc] * 0.5
  s <- sk[cbind(out$IDr, out$IDc)]
  s[!is.finite(s)] <- 0
  out$skpos <- pmin(pmax((sign(s) * abs(s)^(1 / 3) + 3) / 6, 0), 1)
  out$facet <- factor(facet_lab[out$IDr], levels = facet_lab[rows])
  out
}

.gei_theme <- function(bg, labels, strip) {
  ggplot2::theme_void() +
    ggplot2::theme(
      axis.text.x = if (labels) ggplot2::element_text(size = 8, colour = "grey20") else ggplot2::element_blank(),
      strip.text = if (strip) ggplot2::element_text(size = 9, margin = ggplot2::margin(2, 0, 4, 0)) else ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0.4, "lines"),
      plot.background = ggplot2::element_rect(fill = bg, colour = NA),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}
