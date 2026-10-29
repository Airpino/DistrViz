## Ordering of the variables: Hamiltonian cycle (TSP) ---------------------

#' Order of the variables by a shortest Hamiltonian cycle
#'
#' Finds an order of the variables such that variables with similar
#' (strongly correlated) distributions are adjacent. Variables are the nodes
#' of a complete graph whose edges are weighted by a correlation distance
#' between variables; the order is the shortest Hamiltonian cycle on this
#' graph, i.e. the solution of a Travelling Salesman Problem (TSP). Since the
#' sectors of a Green Eyes plot are arranged on a circle, the last variable is
#' adjacent to the first one, and the cycle is the natural choice; for a
#' linear arrangement (e.g. the columns of a heatmap) the cycle is cut at its
#' longest edge, giving a short Hamiltonian path.
#'
#' The correlations are the Wasserstein correlations between the
#' distributional variables (\code{DDAWASS::WD.correlation} or
#' \code{HistDAWass::WH.correlation}). With \code{dist = "sqrt(2(1-r))"}
#' (default) negatively correlated variables are far apart; with
#' \code{dist = "sqrt(1-r^2)"} strongly correlated variables are close
#' whatever the sign of the correlation.
#'
#' The cycle is computed exactly (by enumeration) for up to 9 variables;
#' otherwise it is found heuristically: nearest-neighbour tours from every
#' variable and from \code{repetitions} random tours, each improved by the
#' 2-opt local search, keeping the shortest one.
#'
#' @param x a \code{dvdata} object, or anything accepted by
#'   \code{\link{dv_data}}.
#' @param dist the correlation distance: \code{"sqrt(2(1-r))"} (default),
#'   \code{"sqrt(1-r^2)"} or \code{"1-r"}.
#' @param type \code{"cycle"} (default) or \code{"path"}.
#' @param on \code{"source"} (default) or \code{"discretized"} (see
#'   \code{\link{dv_cluster}}).
#' @param repetitions number of random starting tours of the heuristic (more
#'   than 9 variables only).
#' @param seed (optional) seed for the random starting tours.
#' @return a list with \code{order} (the indices of the variables, starting
#'   from the first variable for a cycle), \code{labels}, \code{length} (the
#'   length of the cycle or path), \code{dist} (the distance matrix) and
#'   \code{exact} (logical).
#' @references Crocetta C., Irpino A., Antonucci L., Marin C. (2025). A
#'   performance indicator and its decomposition according to the impacts of
#'   different aspects based on distributional data. \emph{Quality & Quantity},
#'   59(3), 2091--2110. \doi{10.1007/s11135-024-01956-0}
#' @examples
#' tsp <- dv_var_order(DDAWASS::Satisfaction)
#' tsp$labels
#' gei_plot(DDAWASS::Satisfaction, row = 1, var_order = "tsp", ranges = c(0, 10))
#' @export
dv_var_order <- function(x, dist = c("sqrt(2(1-r))", "sqrt(1-r^2)", "1-r"),
                         type = c("cycle", "path"), on = c("source", "discretized"),
                         repetitions = 100, seed = NULL) {
  x <- dv_data(x)
  dist <- match.arg(dist)
  type <- match.arg(type)
  on <- match.arg(on)
  D <- .var_dist(x, dist, on)
  p <- nrow(D)
  if (p <= 3) {
    ord <- seq_len(p)
    exact <- TRUE
  } else if (p <= 9) {
    ord <- .tsp_exact(D)
    exact <- TRUE
  } else {
    ord <- .tsp_heuristic(D, repetitions, seed)
    exact <- FALSE
  }
  len_cycle <- .tour_length(D, ord)
  if (type == "path" && p > 2) {
    # cut the cycle at its longest edge
    nxt <- c(ord[-1], ord[1])
    e <- which.max(D[cbind(ord, nxt)])
    ord <- if (e == p) ord else c(ord[(e + 1):p], ord[seq_len(e)])
    len <- len_cycle - D[ord[p], ord[1]]
  } else {
    ord <- .canonical_cycle(ord, D)
    len <- len_cycle
  }
  list(order = ord, labels = x$vars[ord], length = len, dist = D, exact = exact)
}

.var_dist <- function(x, dist, on = "source") {
  obj <- if (on == "source") x$source else x$matD
  R <- if (is(obj, "MatH")) HistDAWass::WH.correlation(obj) else DDAWASS::WD.correlation(obj)
  R <- pmin(pmax(as.matrix(R), -1), 1)
  D <- switch(dist,
              "sqrt(2(1-r))" = sqrt(2 * (1 - R)),
              "sqrt(1-r^2)" = sqrt(pmax(1 - R^2, 0)),
              "1-r" = 1 - R)
  diag(D) <- 0
  dimnames(D) <- list(x$vars, x$vars)
  D
}

.tour_length <- function(D, ord) sum(D[cbind(ord, c(ord[-1], ord[1]))])

# start the cycle at variable 1 and follow the direction in which the
# neighbour of variable 1 has the smaller index
.canonical_cycle <- function(ord, D) {
  p <- length(ord)
  k <- which(ord == 1)
  ord <- c(ord[k:p], ord[seq_len(k - 1)])
  if (p > 2 && ord[p] < ord[2]) ord <- c(ord[1], rev(ord[-1]))
  ord
}

# all the permutations of a vector (rows of a matrix)
.perms <- function(v) {
  if (length(v) <= 1) return(matrix(v, 1))
  do.call(rbind, lapply(seq_along(v), function(i) cbind(v[i], .perms(v[-i]))))
}

.tsp_exact <- function(D) {
  p <- nrow(D)
  P <- cbind(1, .perms(2:p))            # cycles starting from node 1
  P <- P[P[, 2] < P[, p], , drop = FALSE] # each cycle once (one direction)
  len <- numeric(nrow(P))
  for (k in seq_len(p)) {
    nx <- if (k < p) P[, k + 1] else P[, 1]
    len <- len + D[cbind(P[, k], nx)]
  }
  P[which.min(len), ]
}

.two_opt <- function(D, ord) {
  p <- length(ord)
  improved <- TRUE
  while (improved) {
    improved <- FALSE
    for (i in 1:(p - 2)) {
      for (j in (i + 2):p) {
        if (i == 1 && j == p) next
        a <- ord[i]; b <- ord[i + 1]; c <- ord[j]; d <- ord[if (j == p) 1 else j + 1]
        delta <- D[a, c] + D[b, d] - D[a, b] - D[c, d]
        if (delta < -1e-12) {
          ord[(i + 1):j] <- rev(ord[(i + 1):j])
          improved <- TRUE
        }
      }
    }
  }
  ord
}

.nearest_neighbour <- function(D, start) {
  p <- nrow(D)
  ord <- start
  left <- setdiff(seq_len(p), start)
  while (length(left)) {
    nx <- left[which.min(D[ord[length(ord)], left])]
    ord <- c(ord, nx)
    left <- setdiff(left, nx)
  }
  ord
}

.tsp_heuristic <- function(D, repetitions = 100, seed = NULL) {
  p <- nrow(D)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  starts <- c(lapply(seq_len(p), function(s) .nearest_neighbour(D, s)),
              lapply(seq_len(repetitions), function(r) sample.int(p)))
  best <- NULL; best_len <- Inf
  for (s in starts) {
    t <- .two_opt(D, s)
    l <- .tour_length(D, t)
    if (l < best_len - 1e-12) { best <- t; best_len <- l }
  }
  best
}
