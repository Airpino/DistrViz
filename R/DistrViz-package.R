#' DistrViz: Green Eyes, flower plots and heatmaps for distributional data
#'
#' Visualization tools for distributional data tables. The input can be a
#' matrix of histograms (\code{HistDAWass::MatH}), a matrix of discrete
#' distributions (\code{DDAWASS::MatD}), a data.frame of raw discrete data
#' (one row per observation, with grouping columns) or a long frequency
#' table.
#' \itemize{
#'   \item \code{\link{dv_data}} prepares the data and maps values to colours;
#'     \code{\link{dv_palette}} gives the palettes.
#'   \item \code{\link{gei_plot}} and \code{\link{gei_grid}} draw Green Eyes
#'     plots; \code{\link{gei_pca}} places them on a factorial plane (also
#'     interactively); \code{\link{gei_diff}} draws difference eyes.
#'   \item \code{\link{flower_plot}} and \code{\link{flower_grid}} draw
#'     flower plots.
#'   \item \code{\link{dist_heatmap}} draws distributional heatmaps with
#'     Wasserstein-based dendrograms.
#'   \item \code{\link{dv_cluster}} and \code{\link{dv_var_order}} order units
#'     and variables (hierarchical clustering, shortest Hamiltonian cycle).
#' }
#' See \code{browseVignettes("DistrViz")} for six worked examples.
#' @keywords internal
"_PACKAGE"
