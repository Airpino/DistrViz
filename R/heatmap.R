## Distributional heatmaps ---------------------------------------------------

#' Distributional heatmap
#'
#' Draws a heatmap of a distributional data table: each cell (unit x
#' variable) is a horizontal strip split into segments whose widths are the
#' probabilities of the support values (or of the histogram classes) and
#' whose colours are given by the position of the values in the range of the
#' variable. Rows and columns can be reordered by hierarchical clustering
#' (Wasserstein distance between rows, Wasserstein correlation between
#' columns) and the dendrograms are drawn on the left and on the top.
#'
#' @param x a \code{dvdata} object, a \code{MatH}, a \code{MatD} or a
#'   \code{data.frame} (passed to \code{\link{dv_data}} with the arguments in
#'   \code{...}).
#' @param cluster_rows,cluster_cols logical; reorder rows/columns by
#'   hierarchical clustering and draw the dendrograms? \code{cluster_cols}
#'   can also be \code{"tsp"}: the columns are ordered by a short
#'   Hamiltonian path through the variables (see \code{\link{dv_var_order}},
#'   with \code{col_dist} as distance) and no column dendrogram is drawn.
#' @param row_method,col_method agglomeration methods (see
#'   \code{\link{dv_cluster}}).
#' @param col_dist distance between variables: \code{"sqrt(1-r^2)"}
#'   (default), \code{"1-r"} or \code{"sqrt(2(1-r))"}.
#' @param cluster_on \code{"source"} or \code{"discretized"} (see
#'   \code{\link{dv_cluster}}).
#' @param palette,reverse colour palette (see \code{\link{dv_palette}}).
#' @param gap_x,gap_y fraction of each cell left empty between cells.
#' @param dendro_size relative size of the row and column dendrograms.
#' @param style \code{"light"} (default) or \code{"dark"}.
#' @param label_size font size of the labels.
#' @param angle_x angle of the column labels.
#' @param legend logical; draw the colour legend?
#' @param legend_labels (optional) labels of the legend at positions 0, 0.5,
#'   1.
#' @param interactive logical; return an interactive \pkg{plotly} plot (with
#'   tooltips)? Requires the \pkg{plotly} package.
#' @param ... arguments passed to \code{\link{dv_data}} (e.g.
#'   \code{n_levels}, \code{scheme}, \code{ranges}).
#' @references Irpino A., Verde R. (2015). Basic statistics for distributional
#'   symbolic variables: a new metric-based approach. \emph{Advances in Data
#'   Analysis and Classification}, 9(2), 143--175.
#'   \doi{10.1007/s11634-014-0176-4}
#' @return a \code{ggplot} object (or a \code{plotly} object), with the
#'   clustering results in the attribute \code{"clustering"}.
#' @examples
#' dist_heatmap(DDAWASS::Satisfaction, ranges = c(0, 10))
#' dist_heatmap(HistDAWass::BLOOD, n_levels = 20, scheme = "standardized",
#'              palette = "Zissou 1", style = "dark")
#' @export
dist_heatmap <- function(x, cluster_rows = TRUE, cluster_cols = TRUE, row_method = "ward.D2",
                         col_method = "complete", col_dist = c("sqrt(1-r^2)", "1-r", "sqrt(2(1-r))"),
                         cluster_on = c("source", "discretized"), palette = "Paquin",
                         reverse = FALSE, gap_x = 0.04, gap_y = 0.12,
                         dendro_size = c(0.18, 0.15), style = c("light", "dark"),
                         label_size = 8, angle_x = 90, legend = TRUE, legend_labels = NULL,
                         interactive = FALSE, ...) {
  x <- dv_data(x, ...)
  style <- match.arg(style)
  col_dist <- match.arg(col_dist)
  cluster_on <- match.arg(cluster_on)
  n <- length(x$rows); p <- length(x$vars)
  tsp_cols <- identical(cluster_cols, "tsp")
  cl <- dv_cluster(x, rows = cluster_rows, cols = isTRUE(cluster_cols), row_method = row_method,
                   col_method = col_method, col_dist = col_dist, on = cluster_on)
  if (tsp_cols) {
    cl$col_order <- dv_var_order(x, dist = col_dist, type = "path", on = cluster_on)$order
  }
  rpos <- match(seq_len(n), cl$row_order)   # position of each row in the order
  cpos <- match(seq_len(p), cl$col_order)
  df <- x$df
  df$xc <- cpos[df$IDc]
  df$yc <- n - rpos[df$IDr] + 1             # first row of the order on top
  df$xmin <- df$xc - 0.5 + gap_x / 2 + df$cdf_lo * (1 - gap_x)
  df$xmax <- df$xc - 0.5 + gap_x / 2 + df$cdf * (1 - gap_x)
  df$ymin <- df$yc - 0.5 + gap_y / 2
  df$ymax <- df$yc + 0.5 - gap_y / 2
  df$text <- paste0(df$row_label, "<br>", df$var_label, "<br>value: ", df$label,
                    "<br>p = ", signif(df$freq, 3), "  F = ", signif(df$cdf, 3))

  fg <- if (style == "dark") "white" else "grey15"
  bgc <- if (style == "dark") "black" else "white"
  cols <- dv_palette(palette, 64, reverse)
  if (is.null(legend_labels)) legend_labels <- .legend_labels(x)

  map <- ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin,
                      ymax = .data$ymax, fill = .data$u)
  if (interactive) map[["text"]] <- ggplot2::aes(text = .data$text)[["text"]]
  g <- ggplot2::ggplot(df) +
    suppressWarnings(ggplot2::geom_rect(map, colour = NA)) +
    ggplot2::scale_fill_gradientn(colours = cols, limits = c(0, 1), breaks = c(0, 0.5, 1),
                                  labels = legend_labels, name = NULL,
                                  guide = if (legend) ggplot2::guide_colourbar(barheight = ggplot2::unit(0.4, "npc"),
                                                                                barwidth = ggplot2::unit(0.6, "lines")) else "none")
  xlim_lo <- 0.5; ylim_hi <- n + 0.5
  if (!is.null(cl$hc_rows)) {
    s <- .dendro_segments(cl$hc_rows)
    w <- dendro_size[1] * p
    s <- data.frame(x = 0.5 - s$y * w - 0.1, xend = 0.5 - s$yend * w - 0.1,
                    y = n - s$x + 1, yend = n - s$xend + 1)
    g <- g + ggplot2::geom_segment(data = s, ggplot2::aes(x = .data$x, xend = .data$xend,
                                                          y = .data$y, yend = .data$yend),
                                   colour = fg, linewidth = 0.3)
    xlim_lo <- 0.5 - w - 0.2
  }
  if (!is.null(cl$hc_cols)) {
    s <- .dendro_segments(cl$hc_cols)
    h <- dendro_size[2] * n
    s <- data.frame(x = s$x, xend = s$xend, y = n + 0.6 + s$y * h, yend = n + 0.6 + s$yend * h)
    g <- g + ggplot2::geom_segment(data = s, ggplot2::aes(x = .data$x, xend = .data$xend,
                                                          y = .data$y, yend = .data$yend),
                                   colour = fg, linewidth = 0.3)
    ylim_hi <- n + 0.7 + h
  }
  row_lab <- unique(df[, c("IDr", "row_label")])
  row_lab <- row_lab$row_label[order(row_lab$IDr)]
  var_lab <- unique(df[, c("IDc", "var_label")])
  var_lab <- var_lab$var_label[order(var_lab$IDc)]
  g <- g +
    ggplot2::scale_x_continuous(breaks = seq_len(p), labels = var_lab[cl$col_order],
                                limits = c(xlim_lo, p + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = n:1, labels = row_lab[cl$row_order],
                                limits = c(0.5, ylim_hi), expand = c(0, 0), position = "right") +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = label_size, angle = angle_x, colour = fg,
                                          hjust = if (angle_x == 0) 0.5 else 1,
                                          vjust = if (angle_x == 90) 0.5 else 1),
      axis.text.y = ggplot2::element_text(size = label_size, colour = fg, hjust = 0),
      plot.background = ggplot2::element_rect(fill = bgc, colour = bgc),
      legend.text = ggplot2::element_text(colour = fg, size = label_size),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
  attr(g, "clustering") <- cl
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) stop("install the 'plotly' package")
    pl <- plotly::ggplotly(g, tooltip = "text")
    attr(pl, "clustering") <- cl
    return(pl)
  }
  g
}

.legend_labels <- function(x) {
  if (x$scheme == "standardized") return(c("-3 sd", "mean", "+3 sd"))
  if (is(x$source, "MatH")) {
    b <- x$breaks
    same <- length(b) == 1 || all(vapply(b, function(z) isTRUE(all.equal(z, b[[1]])), logical(1)))
    if (!same) return(c("low", "mid", "high"))
    # centres of the first, middle and last classes (u = 0, 0.5, 1)
    v <- (b[[1]][-1] + b[[1]][-length(b[[1]])]) / 2
    v <- c(v[1], mean(range(v)), v[length(v)])
    return(format(signif(v, 3), trim = TRUE))
  }
  rg <- x$ranges
  if (nrow(unique(rg)) == 1) {
    v <- c(rg[1, 1], mean(rg[1, ]), rg[1, 2])
    return(format(signif(v, 3), trim = TRUE))
  }
  c("min", "mid", "max")
}
