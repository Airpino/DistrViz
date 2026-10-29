## Flowers on a factorial plane or on a map ---------------------------------

#' Flower plots on a factorial plane or on a map
#'
#' Draws the flower plot (see \code{\link{flower_plot}}) of each unit at the
#' position of the unit on a plane: the first factorial plane of a
#' quantile-based PCA/MFA (as in \code{\link{gei_pca}}) or any pair of
#' coordinates, e.g. longitude and latitude, possibly over a map
#' (\code{basemap}). The diameter of each flower is \code{size} times the
#' smallest of the width and the height of the cloud of points. A key on the
#' right shows the flower of the barycenter of the units, with the names of
#' the variables around the petals.
#'
#' @inheritParams gei_pca
#' @inheritParams flower_plot
#' @param key \code{"mean"} (default: the flower of the barycenter of the
#'   drawn units, see \code{\link{gei_pca}}) or \code{"none"}.
#' @param key_position (optional) \code{c(x, y)}: centre of the key flower
#'   inside the plot (see \code{\link{gei_pca}}).
#' @param style \code{"light"} (default, white background) or \code{"dark"}.
#' @return a \code{ggplot} object, with the coordinates used in the attribute
#'   \code{"coord"}.
#' @examples
#' dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
#' flower_pca(dv, rows = 1:12, var_order = "tsp")
#' \donttest{
#' # on a map: monthly temperatures of the Chinese stations
#' if (requireNamespace("maps", quietly = TRUE)) {
#'   library(HistDAWass)
#'   Temp <- China_Month[, 13:24]
#'   lonlat <- with(stations_coordinates, cbind(LON.deg + LON.min / 60,
#'                                              LAT.deg + LAT.min / 60))
#'   china <- ggplot2::map_data("world", region = c("China", "Taiwan"))
#'   flower_pca(Temp, lonlat, basemap = china, scheme = "global", n_levels = 60,
#'              palette = "Blue-Red 3", axis_titles = c("longitude", "latitude"),
#'              labels = FALSE)
#' }
#' }
#' @export
flower_pca <- function(x, pca = NULL, axes = c(1, 2), rows = NULL, size = 0.1, quantiles = 10,
                       groups = NULL, labels = TRUE, label_size = 2.3, basemap = NULL,
                       axis_titles = NULL, key = c("mean", "none"), key_size = 0.3,
                       key_position = NULL,
                       palette = "Paquin", reverse = FALSE,
                       petal = c("auto", "ellipse", "violin"), petal_width = 0.9,
                       width_scale = c("common", "unit", "variable"),
                       width_type = c("angular", "linear"), center = c("flower", "plain"),
                       center_size = 0.22, show_median = FALSE, var_order = NULL,
                       style = c("light", "dark"), alpha = 0.95, n_levels = 30, smooth = 1.5,
                       ...) {
  x <- .flower_data(x, n_levels, ...)
  key <- match.arg(key); style <- match.arg(style)
  nr <- length(x$rows); p <- length(x$vars)
  if (length(axes) != 2) stop("axes must give two axes")
  if (is.null(pca)) pca <- .dv_pca(x, quantiles, max(axes))
  coord <- .pca_coord(pca, x$rows)
  if (max(axes) > ncol(coord)) stop("the factorial analysis has only ", ncol(coord), " axes")
  xy <- coord[, axes, drop = FALSE]
  if (is.null(rows)) rows <- seq_len(nr)
  if (is.character(rows)) rows <- match(rows, x$rows)
  if (anyNA(rows)) stop("unknown rows")
  if (!is.null(groups) && length(groups) != nr) stop("groups must have one value per row of the data")
  ww <- c(diff(range(xy[rows, 1])), diff(range(xy[rows, 2])))
  w <- if (any(ww > 0)) min(ww[ww > 0]) else 1
  R <- size * w / 2
  pad <- R * 1.1 + if (labels) R * 0.5 else 0
  xlim <- range(xy[rows, 1]) + c(-1, 1) * R * 1.1
  ylim <- range(xy[rows, 2]) + c(-pad, R * 1.1)
  bm <- .basemap_df(basemap)
  if (!is.null(bm)) { xlim <- range(c(xlim, bm$x)); ylim <- range(c(ylim, bm$y)) }

  # key: the barycenter is added to the data as an extra unit
  kid <- NULL
  if (key == "mean") {
    bc <- .bary_cells(x, rows)
    kid <- nr + 1
    add <- data.frame(row = "Mean", variable = x$vars[bc$IDc], IDr = kid, IDc = bc$IDc,
                      value = bc$value, label = bc$label, u = bc$u, freq = bc$freq,
                      cdf_lo = bc$cdf_lo, cdf = bc$cdf, row_label = "Mean",
                      var_label = x$vars[bc$IDc])
    x$df$row <- as.character(x$df$row); x$df$variable <- as.character(x$df$variable)
    x$df <- rbind(x$df[, names(add)], add)
    x$rows <- c(x$rows, "Mean")
    # position and size of the key flower, on the right of the plane
    h <- diff(ylim)
    Rk <- key_size * h / 2
    L <- max(nchar(x$vars))
    kc <- 0.55 * 2.6 / 200
    if (is.null(key_position)) {
      a <- min(2 * L * kc, 0.5)
      b0 <- xlim[2] + 0.3 * Rk
      x1 <- (b0 + 2.56 * Rk - a * xlim[1]) / (1 - a)
      cw <- kc * (x1 - xlim[1])
      xk <- b0 + 0.2 * Rk + Rk * 1.08 + L * cw
      yk <- mean(ylim)
      by <- ylim
    } else {
      cw <- kc * diff(xlim)
      xk <- key_position[1]; yk <- key_position[2]
      b0 <- xk - Rk * 1.1 - L * cw - 0.2 * Rk
      x1 <- xk + Rk * 1.1 + L * cw + 0.2 * Rk
      by <- c(yk - Rk * 1.1 - 2.4 * cw - 0.2 * Rk, yk + Rk * 1.1 + 6 * cw + 0.2 * Rk)
    }
  }
  pos <- rbind(xy, if (!is.null(kid)) c(xk, yk))
  rad <- c(rep(R, nr), if (!is.null(kid)) Rk)
  place <- function(d) {
    if (is.null(d)) return(NULL)
    cx <- pos[d$IDr, 1]; cy <- pos[d$IDr, 2]; s <- rad[d$IDr]
    d$x <- cx + s * d$x; d$y <- cy + s * d$y
    if (!is.null(d$xend)) { d$xend <- cx + s * d$xend; d$yend <- cy + s * d$yend }
    d
  }
  fl <- .flower_core(x, rows = c(rows, kid), palette = palette, reverse = reverse, petal = petal,
                     petal_width = petal_width, width_scale = width_scale, width_type = width_type,
                     center = center, center_size = center_size, show_median = show_median,
                     var_order = var_order, labels = FALSE, style = style, alpha = alpha,
                     strip = FALSE, smooth = smooth, place = place)
  fl$cen$grp <- fl$cen$IDr
  fg <- if (style == "dark") "grey85" else "grey20"
  bgc <- if (style == "dark") "black" else "white"

  g <- ggplot2::ggplot()
  if (is.null(bm)) {
    g <- g +
      ggplot2::annotate("segment", x = min(xy[rows, 1]) - R * 1.1, xend = max(xy[rows, 1]) + R * 1.1,
                        y = 0, yend = 0, colour = "grey75", linetype = "dashed") +
      ggplot2::annotate("segment", y = min(xy[rows, 2]) - R * 1.6, yend = max(xy[rows, 2]) + R * 1.1,
                        x = 0, xend = 0, colour = "grey75", linetype = "dashed")
  } else {
    g <- g + .basemap_layer(bm, fill = if (style == "dark") "grey15" else "grey95",
                            colour = if (style == "dark") "grey45" else "grey60")
  }
  if (!is.null(kid)) {
    g <- g + ggplot2::annotate("rect", xmin = b0, xmax = x1, ymin = by[1], ymax = by[2],
                               fill = bgc, colour = "grey80", linewidth = 0.3)
  }
  g <- .flower_layers(g, fl$st, fl$pet, fl$ol, fl$cen, fl$sep, fl$md, fl$cols, fl$petal,
                      fl$alpha, fl$center, fl$fg, fl$stem_col, group_var = "grp",
                      pt_size = 0.9, lw = 0.7, stroke_map = is.null(groups))
  if (!is.null(groups)) {
    th <- seq(0, 2 * pi, length.out = 121)
    oc <- data.frame(id = rep(rows, each = length(th)),
                     x = rep(xy[rows, 1], each = length(th)) + R * 1.02 * sin(th),
                     y = rep(xy[rows, 2], each = length(th)) + R * 1.02 * cos(th))
    oc$group <- factor(groups[oc$id])
    g <- g + ggplot2::geom_path(data = oc, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                                        colour = .data$group), linewidth = 0.4) +
      ggplot2::labs(colour = "group")
  }
  if (labels) {
    lb <- data.frame(x = xy[rows, 1], y = xy[rows, 2] - R * 1.05, lab = x$rows[rows])
    g <- g + ggplot2::geom_text(data = lb, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
                                size = label_size, vjust = 1, colour = fg)
  }
  if (!is.null(kid)) {
    thm <- (seq_len(p) - 0.5) / p * 2 * pi
    kl <- data.frame(x = xk + 1.06 * Rk * sin(thm), y = yk + 1.06 * Rk * cos(thm), lab = x$vars[fl$vo],
                     hjust = ifelse(sin(thm) > 0.15, 0, ifelse(sin(thm) < -0.15, 1, 0.5)),
                     vjust = ifelse(cos(thm) > 0.15, 0, ifelse(cos(thm) < -0.15, 1, 0.5)))
    g <- g + ggplot2::geom_text(data = kl, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab,
                                                        hjust = .data$hjust, vjust = .data$vjust),
                                size = 2.6, colour = fg) +
      ggplot2::annotate("text", x = xk, y = yk + Rk * 1.06 + 3.2 * cw,
                        label = sprintf("Mean (n = %d)", length(rows)), fontface = "bold",
                        size = 2.9, colour = fg, vjust = 0)
    xlim <- range(c(xlim, b0, x1)); ylim <- range(c(ylim, by))
  }
  pct <- .pca_eig(pca)
  lab_ax <- if (!is.null(axis_titles)) rep(axis_titles, length.out = 2) else
    paste0("Dim ", axes, if (!is.null(pct)) sprintf(" (%.1f%%)", pct[axes]) else "")
  g <- g + ggplot2::coord_equal(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
    ggplot2::labs(x = lab_ax[1], y = lab_ax[2],
                  caption = paste("Petals, clockwise from the top:",
                                  paste(x$vars[fl$vo], collapse = ", "))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = bgc, colour = NA),
                   axis.text = ggplot2::element_text(colour = fg),
                   axis.title = ggplot2::element_text(colour = fg),
                   legend.text = ggplot2::element_text(colour = fg),
                   legend.title = ggplot2::element_text(colour = fg),
                   panel.grid.major = ggplot2::element_line(
                     colour = if (style == "dark") "grey25" else "grey90"),
                   plot.caption = ggplot2::element_text(hjust = 0, size = 8, colour = fg))
  attr(g, "coord") <- xy
  g
}
