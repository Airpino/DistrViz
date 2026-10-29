#' @import methods
#' @importFrom rlang .data
NULL

#' Prepare distributional data for visualization
#'
#' Transforms a matrix of histograms (\code{HistDAWass::MatH}), a matrix of
#' discrete distributions (\code{DDAWASS::MatD}) or a table of discrete data
#' into a \code{dvdata} object, the common input of \code{\link{gei_plot}},
#' \code{\link{gei_grid}} and \code{\link{dist_heatmap}}.
#'
#' Every support value (or histogram class) is mapped to a colour position
#' \eqn{u \in [0, 1]}:
#' \itemize{
#'   \item \code{MatH}: each histogram is discretized into \code{n_levels}
#'     equal-width classes (see \code{DDAWASS::MatH2MatD}) spanning the range
#'     of the variable, or \code{ranges} if given; class \eqn{k}
#'     gets \eqn{u = (k - 1)/(n_{levels} - 1)}.
#'   \item \code{MatD}: \eqn{u = (x - a_j)/(b_j - a_j)} where
#'     \eqn{[a_j, b_j]} is the range of variable \eqn{j}: taken from
#'     \code{ranges} if given (e.g. the admissible range of a rating scale,
#'     \code{c(0, 10)}), otherwise from the data. Variables measured on
#'     different scales (e.g. 1-5 and 0-10) are thus coloured consistently.
#' }
#' With \code{scheme = "global"} the same range is used for all variables;
#' with \code{scheme = "standardized"} each variable is standardised with its
#' Wasserstein mean and standard deviation and the colours span \eqn{[-3, 3]}
#' standard deviations.
#'
#' @param x a \code{MatH}, a \code{MatD} or a \code{data.frame}.
#' @param n_levels number of classes used to discretize histograms (only for
#'   \code{MatH}).
#' @param scheme \code{"variable"} (default), \code{"global"} or
#'   \code{"standardized"}.
#' @param ranges (optional) the ranges of the variables: a vector
#'   \code{c(min, max)} used for all variables, a two-column matrix (one row
#'   per variable) or a list of \code{c(min, max)}. For a \code{MatD} the
#'   colour position of a value is its position in the range; for a
#'   \code{MatH} the range is split into \code{n_levels} classes of equal
#'   width (the mass outside the range goes to the first or last class).
#'   With an odd \code{n_levels} and a range symmetric around a value
#'   \eqn{v}, the central class is centred on \eqn{v} (useful with diverging
#'   palettes). Not allowed with \code{scheme = "standardized"}.
#' @param unit,vars when \code{x} is a \code{data.frame} of raw data (one row
#'   per observation): the grouping column and the variables (see
#'   \code{DDAWASS::data2MatD}).
#' @param support,min.n when \code{x} is a \code{data.frame} of raw data:
#'   the admissible values of the variables (e.g. \code{1:7} for a 1-7 rating
#'   scale, so that values never observed get probability zero) and the
#'   minimum number of observations required to keep a unit (see
#'   \code{DDAWASS::data2MatD}).
#' @param value,freq,variable when \code{x} is a long frequency table:
#'   the columns with the values, the frequencies and the variable names
#'   (see \code{DDAWASS::long2MatD}). The long format is used when
#'   \code{freq} is given.
#' @param row_labels,var_labels (optional) longer labels for rows and
#'   variables (used in tooltips and axes).
#' @return an object of class \code{dvdata}: a list with the long table
#'   \code{df} (columns \code{row}, \code{variable}, \code{IDr}, \code{IDc},
#'   \code{value}, \code{label}, \code{u}, \code{freq}, \code{cdf_lo},
#'   \code{cdf}), the summary statistics \code{stats}, and the original data
#'   \code{source}.
#' @examples
#' dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
#' dv
#' dvh <- dv_data(HistDAWass::BLOOD, n_levels = 12)
#' @export
dv_data <- function(x, n_levels = 10, scheme = c("variable", "global", "standardized"),
                    ranges = NULL, unit = NULL, vars = NULL, support = NULL,
                    min.n = 1, value = "value",
                    freq = NULL, variable = "variable", row_labels = NULL,
                    var_labels = NULL) {
  scheme <- match.arg(scheme)
  if (inherits(x, "dvdata")) return(x)
  if (is.data.frame(x)) {
    if (is.null(unit)) stop("for a data.frame you must give the 'unit' column")
    x <- if (!is.null(freq)) {
      DDAWASS::long2MatD(x, unit = unit, variable = variable, value = value, freq = freq)
    } else {
      if (is.null(vars)) vars <- setdiff(names(x)[vapply(x, is.numeric, logical(1))], unit)
      if (is.null(ranges) && !is.null(support)) ranges <- range(support)
      DDAWASS::data2MatD(x, unit = unit, vars = vars, support = support, min.n = min.n)
    }
  }
  if (is(x, "MatH") || is(x, "MatD")) x <- .fix_names(x)
  if (is(x, "MatH")) {
    res <- .dv_from_MatH(x, n_levels, scheme, ranges)
  } else if (is(x, "MatD")) {
    res <- .dv_from_MatD(x, scheme, ranges)
  } else {
    stop("x must be a MatH, a MatD or a data.frame")
  }
  rn <- rownames(x@M); cn <- colnames(x@M)
  res$df$row <- factor(res$df$row, levels = rn)
  res$df$variable <- factor(res$df$variable, levels = cn)
  res$df$IDr <- as.integer(res$df$row)
  res$df$IDc <- as.integer(res$df$variable)
  res$df$row_label <- if (is.null(row_labels)) as.character(res$df$row) else row_labels[res$df$IDr]
  res$df$var_label <- if (is.null(var_labels)) as.character(res$df$variable) else var_labels[res$df$IDc]
  res$rows <- rn
  res$vars <- cn
  res$source <- x
  res$scheme <- scheme
  class(res) <- "dvdata"
  res
}

.ranges_matrix <- function(ranges, p) {
  rg <- if (is.list(ranges) && !is.data.frame(ranges)) do.call(rbind, ranges) else
    if (is.matrix(ranges) || is.data.frame(ranges)) as.matrix(ranges) else
      matrix(ranges, p, 2, byrow = TRUE)
  if (nrow(rg) != p || ncol(rg) != 2) stop("ranges must have one row (min, max) per variable")
  if (any(rg[, 2] <= rg[, 1])) stop("in ranges, max must be greater than min")
  rg
}

# missing or NA row/column names (e.g. after subsetting a MatH by names in
# some versions of HistDAWass) are replaced by I1, I2, ... and V1, V2, ...
.fix_names <- function(x) {
  rn <- rownames(x@M); cn <- colnames(x@M)
  if (is.null(rn) || anyNA(rn)) rn <- paste0("I", seq_len(nrow(x@M)))
  if (is.null(cn) || anyNA(cn)) {
    cn <- paste0("V", seq_len(ncol(x@M)))
    warning("the variables have no names: V1, V2, ... are used")
  }
  dimnames(x@M) <- list(make.unique(rn), make.unique(cn))
  x
}

.dv_from_MatH <- function(x, n_levels, scheme, ranges = NULL) {
  brk <- NULL
  if (!is.null(ranges)) {
    if (scheme == "standardized") stop("'ranges' cannot be used with scheme = \"standardized\"")
    rg <- .ranges_matrix(ranges, ncol(x@M))
    brk <- lapply(seq_len(nrow(rg)), function(j) seq(rg[j, 1], rg[j, 2], length.out = n_levels + 1))
  }
  dd <- DDAWASS::MatH2MatD(x, n_levels = n_levels, scheme = scheme, breaks = brk,
                           values = "index", return.breaks = TRUE)
  M <- dd$MatD
  df <- DDAWASS::MatD2long(M)
  j <- as.integer(df$variable)
  br <- dd$breaks
  names(br) <- colnames(x@M)
  lab <- vapply(seq_len(nrow(df)), function(r) {
    b <- br[[j[r]]]; k <- df$x[r]
    sprintf("[%s, %s]", format(signif(b[k], 4)), format(signif(b[k + 1], 4)))
  }, character(1))
  out <- data.frame(row = df$row, variable = df$variable, value = df$x, label = lab,
                    u = (df$x - 1) / (n_levels - 1), freq = df$p,
                    cdf_lo = df$cdf - df$p, cdf = df$cdf)
  st <- .dv_stats_MatH(x)
  list(df = out, stats = st, breaks = br, n_levels = n_levels, matD = M)
}

.dv_from_MatD <- function(x, scheme, ranges) {
  df <- DDAWASS::MatD2long(x)
  p <- ncol(x@M)
  j <- as.integer(df$variable)
  if (scheme == "standardized") {
    mu <- vapply(seq_len(p), function(k) DDAWASS::meanD(DDAWASS::WD.vec.mean(x[, k])), 0)
    sg <- sqrt(diag(DDAWASS::WD.var.covar(x)))
    sg[sg == 0] <- 1
    u <- pmin(pmax(((df$x - mu[j]) / sg[j] + 3) / 6, 0), 1)
    rg <- cbind(mu - 3 * sg, mu + 3 * sg)
  } else {
    if (is.null(ranges)) {
      rg <- t(vapply(seq_len(p), function(k) range(df$x[j == k]), numeric(2)))
      if (scheme == "global") rg <- matrix(range(df$x), p, 2, byrow = TRUE)
    } else {
      rg <- .ranges_matrix(ranges, p)
    }
    w <- rg[, 2] - rg[, 1]
    w[w == 0] <- 1
    u <- pmin(pmax((df$x - rg[j, 1]) / w[j], 0), 1)
  }
  out <- data.frame(row = df$row, variable = df$variable, value = df$x,
                    label = format(signif(df$x, 4), trim = TRUE), u = u, freq = df$p,
                    cdf_lo = df$cdf - df$p, cdf = df$cdf)
  st <- .dv_stats_MatD(x)
  list(df = out, stats = st, ranges = rg, matD = x)
}

.dv_stats_MatD <- function(x) {
  g <- function(s) DDAWASS::get.MatD.stats(x, s)$mat
  list(mean = g("mean"), std = g("std"), skewness = g("skewness"), median = g("median"))
}

.dv_stats_MatH <- function(x) {
  g <- function(s) HistDAWass::get.MatH.stats(x, stat = s)$mat
  list(mean = g("mean"), std = g("std"), skewness = g("skewness"), median = g("median"))
}

#' @export
print.dvdata <- function(x, ...) {
  cat("DistrViz data:", length(x$rows), "rows x", length(x$vars), "variables",
      sprintf("(from a %s, scheme '%s')\n", class(x$source)[1], x$scheme))
  print(utils::head(x$df[, c("row", "variable", "value", "label", "u", "freq", "cdf")]))
  invisible(x)
}

## clustering helpers -------------------------------------------------------

#' Orders of rows and columns from hierarchical clustering
#'
#' Rows are clustered with the L2 Wasserstein distance (\code{WH_hclust} for
#' histograms, \code{WD_hclust} for discrete distributions), columns with the
#' distance \eqn{\sqrt{1 - r^2}} (or \eqn{1 - r}, or \eqn{\sqrt{2(1-r)}}) where \eqn{r} is the
#' Wasserstein correlation between variables.
#' @param x a \code{dvdata} object (or anything accepted by \code{dv_data}).
#' @param rows,cols logical; cluster rows / columns?
#' @param row_method,col_method agglomeration methods.
#' @param col_dist \code{"sqrt(1-r^2)"} (default), \code{"1-r"} or
#'   \code{"sqrt(2(1-r))"}.
#' @param on \code{"source"} (default: cluster the original data) or
#'   \code{"discretized"} (cluster the discretized distributions that are
#'   drawn).
#' @return a list with \code{hc_rows}, \code{hc_cols} (\code{hclust} objects
#'   or \code{NULL}), \code{row_order} and \code{col_order}.
#' @examples
#' cl <- dv_cluster(dv_data(DDAWASS::Satisfaction))
#' cl$row_order
#' @export
dv_cluster <- function(x, rows = TRUE, cols = TRUE, row_method = "ward.D2",
                       col_method = "complete", col_dist = c("sqrt(1-r^2)", "1-r", "sqrt(2(1-r))"),
                       on = c("source", "discretized")) {
  x <- dv_data(x)
  col_dist <- match.arg(col_dist)
  on <- match.arg(on)
  obj <- if (on == "source") x$source else x$matD
  hr <- hcc <- NULL
  if (rows && length(x$rows) > 2) {
    hr <- if (is(obj, "MatH")) HistDAWass::WH_hclust(obj, method = row_method) else
      DDAWASS::WD_hclust(obj, method = row_method)
    hr$labels <- x$rows
  }
  if (cols && length(x$vars) > 2) {
    D <- .var_dist(x, col_dist, on)
    hcc <- stats::hclust(stats::as.dist(D), method = col_method)
    hcc$labels <- x$vars
  }
  list(hc_rows = hr, hc_cols = hcc,
       row_order = if (is.null(hr)) seq_along(x$rows) else hr$order,
       col_order = if (is.null(hcc)) seq_along(x$vars) else hcc$order)
}

# rectangular dendrogram segments from an hclust object
# leaves at positions 1..n (in hc$order), heights in [0, 1]
.dendro_segments <- function(hc) {
  n <- length(hc$order)
  pos_leaf <- numeric(n)
  pos_leaf[hc$order] <- seq_len(n)
  h <- hc$height / max(hc$height)
  xn <- numeric(n - 1); yn <- h
  getxy <- function(k) if (k < 0) c(pos_leaf[-k], 0) else c(xn[k], yn[k])
  segs <- vector("list", n - 1)
  for (i in seq_len(n - 1)) {
    a <- getxy(hc$merge[i, 1]); b <- getxy(hc$merge[i, 2])
    xn[i] <- (a[1] + b[1]) / 2
    segs[[i]] <- data.frame(x = c(a[1], b[1], a[1]), xend = c(a[1], b[1], b[1]),
                            y = c(a[2], b[2], yn[i]), yend = c(yn[i], yn[i], yn[i]))
  }
  do.call(rbind, segs)
}

## palettes -----------------------------------------------------------------

#' Colour palettes for DistrViz
#'
#' Returns a vector of colours. \code{palette} can be a vector of colours,
#' one of the built-in palettes (\code{"Paquin"}, the default red-to-green
#' palette of Green Eyes plots, \code{"Hiroshige"}, \code{"Grey"}) or any
#' name in \code{\link[grDevices]{hcl.pals}()} (e.g. \code{"RdYlGn"},
#' \code{"Zissou 1"}, \code{"Blue-Red 3"}, \code{"Viridis"}).
#' @param palette a palette name or a vector of colours.
#' @param n number of colours.
#' @param reverse logical; reverse the order?
#' @return a character vector of colours.
#' @examples
#' dv_palette("Paquin", 11)
#' dv_palette("Zissou 1", 5, reverse = TRUE)
#' @export
dv_palette <- function(palette = "Paquin", n = 11, reverse = FALSE) {
  builtin <- list(
    Paquin = c("#831818", "#c62320", "#f05b43", "#f78462", "#feac81", "#f7dea3",
               "#ced1af", "#98ab76", "#748f46", "#47632a", "#275024"),
    Hiroshige = c("#e76254", "#ef8a47", "#f7aa58", "#ffd06f", "#ffe6b7", "#aadce0",
                  "#72bcd5", "#528fad", "#376795", "#1e466e"),
    Grey = c("#f7f7f7", "#000000")
  )
  cols <- if (length(palette) > 1) palette else if (palette %in% names(builtin)) {
    builtin[[palette]]
  } else if (palette %in% grDevices::hcl.pals()) {
    grDevices::hcl.colors(max(n, 2), palette, rev = FALSE)
  } else stop("unknown palette: see ?dv_palette")
  cols <- grDevices::colorRampPalette(cols)(n)
  if (reverse) rev(cols) else cols
}
