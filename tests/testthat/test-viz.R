test_that("dv_data accepts MatD, MatH, raw and long data", {
  d1 <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
  expect_s3_class(d1, "dvdata")
  expect_true(all(d1$df$u >= 0 & d1$df$u <= 1))
  sums <- tapply(d1$df$freq, list(d1$df$IDr, d1$df$IDc), sum)
  expect_equal(as.numeric(sums), rep(1, length(sums)), tolerance = 1e-10)

  d2 <- dv_data(HistDAWass::BLOOD, n_levels = 7)
  expect_equal(sort(unique(d2$df$value)), 1:7)

  d3 <- dv_data(DDAWASS::Satisfaction_raw, unit = "hotel", vars = c("Staff", "Price"))
  expect_equal(length(d3$rows), 30)

  lg <- DDAWASS::MatD2long(DDAWASS::Satisfaction[1:3, 1:2])
  d4 <- dv_data(lg, unit = "row", variable = "variable", value = "x", freq = "p")
  expect_equal(d4$vars, c("Staff", "Cleanliness"))

  d5 <- dv_data(DDAWASS::Satisfaction, scheme = "standardized")
  expect_true(all(d5$df$u >= 0 & d5$df$u <= 1))
})

test_that("plots are built", {
  dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
  expect_s3_class(gei_plot(dv, 1, show_sd = TRUE, show_skewness = TRUE), "ggplot")
  expect_s3_class(gei_grid(dv, rows = 1:6, order = "cluster"), "ggplot")
  expect_s3_class(gei_grid(dv, groups = rep(1:3, 10), var_order = "cluster"), "ggplot")
  h <- dist_heatmap(dv)
  expect_s3_class(h, "ggplot")
  expect_length(attr(h, "clustering")$row_order, 30)
  expect_s3_class(dist_heatmap(HistDAWass::BLOOD, cluster_cols = FALSE, style = "dark"), "ggplot")
  expect_silent(ggplot2::ggplot_build(gei_plot(HistDAWass::BLOOD, 2)))
})

test_that("dendrogram segments are well formed", {
  hc <- stats::hclust(stats::dist(1:6))
  s <- DistrViz:::.dendro_segments(hc)
  expect_equal(nrow(s), 3 * 5)
  expect_true(all(s$y >= 0 & s$yend <= 1))
})

test_that("palettes", {
  expect_length(dv_palette("Paquin", 5), 5)
  expect_length(dv_palette("RdYlGn", 7), 7)
  expect_error(dv_palette("nonexistent"))
})

test_that("dv_data handles support/min.n for raw data", {
  d <- data.frame(g = rep(c("a", "b", "c"), c(5, 5, 2)),
                  y = c(1, 2, 2, 3, 3, 4, 5, 5, 5, 4, 1, 1))
  dv <- dv_data(d, unit = "g", vars = "y", support = 1:5, min.n = 3)
  expect_equal(dv$rows, c("a", "b"))
  expect_equal(range(dv$df$u), c(0, 1))
  expect_equal(sum(dv$df$freq[dv$df$row == "a"]), 1)
})

test_that("ranges work for histograms", {
  dv <- dv_data(HistDAWass::BLOOD[, 1], n_levels = 5, ranges = c(0, 500))
  expect_equal(dv$breaks[[1]], seq(0, 500, length.out = 6))
  expect_equal(names(dv$breaks), colnames(HistDAWass::BLOOD@M)[1])
  tot <- tapply(dv$df$freq, dv$df$row, sum)
  expect_true(all(abs(tot - 1) < 1e-8))
  expect_error(dv_data(HistDAWass::BLOOD, ranges = c(0, 1), scheme = "standardized"))
  lab <- DistrViz:::.legend_labels(dv)
  expect_equal(lab, c("50", "250", "450"))
})

test_that("missing column names are repaired", {
  B <- HistDAWass::BLOOD
  colnames(B@M) <- NULL
  expect_warning(dv <- dv_data(B, n_levels = 5))
  expect_equal(dv$vars, c("V1", "V2", "V3"))
})

test_that("polar = FALSE draws no pupil", {
  g <- gei_plot(DDAWASS::Satisfaction, row = 1, polar = FALSE)
  expect_s3_class(g, "ggplot")
})

test_that("Hamiltonian cycle ordering", {
  set.seed(3)
  P <- matrix(runif(16), 8, 2)
  D <- as.matrix(dist(P))
  ex <- DistrViz:::.tsp_exact(D)
  he <- DistrViz:::.tsp_heuristic(D, 30, 1)
  expect_equal(sort(ex), 1:8)
  expect_equal(DistrViz:::.tour_length(D, he), DistrViz:::.tour_length(D, ex))
  # brute-force check of the exact solver on 6 nodes
  D6 <- D[1:6, 1:6]
  P6 <- cbind(1, DistrViz:::.perms(2:6))
  best <- min(apply(P6, 1, function(o) DistrViz:::.tour_length(D6, o)))
  expect_equal(DistrViz:::.tour_length(D6, DistrViz:::.tsp_exact(D6)), best)

  tsp <- dv_var_order(DDAWASS::Satisfaction)
  expect_equal(sort(tsp$order), seq_along(colnames(DDAWASS::Satisfaction@M)))
  expect_equal(tsp$order[1], 1L)
  pth <- dv_var_order(DDAWASS::Satisfaction, type = "path")
  expect_lte(pth$length, tsp$length)
  g <- gei_plot(DDAWASS::Satisfaction, row = 1, var_order = "tsp")
  expect_s3_class(g, "ggplot")
  h <- dist_heatmap(DDAWASS::Satisfaction, cluster_cols = "tsp")
  expect_null(attr(h, "clustering")$hc_cols)
})

test_that("gei_pca places the eyes on the plane", {
  dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
  g <- gei_pca(dv)
  expect_s3_class(g, "ggplot")
  xy <- attr(g, "coord")
  expect_equal(nrow(xy), length(dv$rows))
  # eye diameter <= 1/10 of width and height
  co <- matrix(rnorm(2 * length(dv$rows)), ncol = 2, dimnames = list(dv$rows, NULL))
  g2 <- gei_pca(dv, co, groups = rep(1:2, length.out = length(dv$rows)))
  pts <- ggplot2::layer_data(g2, 3)
  d <- tapply(pts$x, pts$group, function(v) diff(range(v)))
  expect_true(all(d <= 0.1 * min(diff(range(co[, 1])), diff(range(co[, 2]))) + 1e-8))
  expect_error(gei_pca(dv, co[1:3, ]))
  gh <- gei_pca(HistDAWass::BLOOD, n_levels = 8)
  expect_s3_class(gh, "ggplot")
})

test_that("interactive gei_pca", {
  skip_if_not_installed("ggiraph")
  dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
  w <- gei_pca(dv, interactive = TRUE, hidden = c("Hotel_01", "Hotel_03"),
               groups = rep(1:3, length.out = length(dv$rows)))
  expect_s3_class(w, "girafe")
  expect_equal(w$x$settings$select$selected, c("Hotel_01", "Hotel_03"))
  expect_true(grepl("data-id", w$x$html))
  expect_true(any(grepl("gei-controls", unlist(w$jsHooks))))
  w2 <- gei_pca(dv, interactive = TRUE, controls = FALSE)
  expect_length(w2$jsHooks, 0)
})

test_that("key eyes of gei_pca", {
  dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
  bc <- DistrViz:::.bary_cells(dv, seq_along(dv$rows))
  expect_equal(as.numeric(tapply(bc$freq, bc$IDc, sum)), rep(1, length(dv$vars)))
  expect_true(all(bc$u >= 0 & bc$u <= 1))
  g0 <- gei_pca(dv, key = "none"); g1 <- gei_pca(dv)
  expect_gt(length(g1$layers), length(g0$layers))
  gr <- rep(1:3, length.out = length(dv$rows))
  g2 <- gei_pca(dv, groups = gr, key = "groups")
  expect_s3_class(g2, "ggplot")
  dvh <- dv_data(HistDAWass::BLOOD, n_levels = 8)
  bh <- DistrViz:::.bary_cells(dvh, 1:5)
  expect_equal(as.numeric(tapply(bh$freq, bh$IDc, sum)), rep(1, 3))
})

test_that("difference eyes", {
  g <- gei_diff(DDAWASS::Satisfaction, rows = 1:6)
  expect_s3_class(g, "ggplot")
  D <- attr(g, "diff")
  M <- DDAWASS::Satisfaction
  expect_equal(dim(D), c(nrow(M@M), ncol(M@M), 200))
  # mean of the differences from the barycenter is zero at each level
  expect_lt(max(abs(apply(D, c(2, 3), mean))), 1e-10)
  # mean squared difference = Wasserstein variance
  v <- colMeans(apply(attr(gei_diff(M, levels = 2000, scale = "raw"), "diff")^2, c(1, 2), mean))
  expect_equal(unname(v), unname(diag(DDAWASS::WD.var.covar(M))), tolerance = 0.01)
  # reference unit: its own differences are zero
  g2 <- gei_diff(M, reference = 1)
  expect_equal(max(abs(attr(g2, "diff")[1, , ])), 0)
  gh <- gei_diff(HistDAWass::BLOOD, rows = 1:3)
  expect_s3_class(gh, "ggplot")
  expect_error(gei_diff(M, reference = "nope"))
})

test_that("flower plots", {
  dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
  g <- flower_plot(dv, row = 1, show_median = TRUE)
  expect_s3_class(g, "ggplot")
  g2 <- flower_grid(dv, rows = 1:4, var_order = "tsp", style = "light", center = "plain",
                    width_scale = "unit")
  expect_s3_class(g2, "ggplot")
  gh <- flower_plot(HistDAWass::BLOOD, row = 2, n_levels = 20)
  expect_s3_class(gh, "ggplot")
  expect_error(flower_plot(dv, row = "nope"))
  # smoothing keeps the total mass
  f <- c(0, 0.1, 0.5, 0.3, 0.1, 0)
  expect_equal(sum(DistrViz:::.smooth_profile(f, 1.5)), 1)
  expect_equal(DistrViz:::.smooth_profile(f, 0), f)
})

test_that("flower widths: linear and per variable", {
  A <- HistDAWass::Age_Pyramids_2014[1:4, 2:3]
  g <- flower_grid(A, n_levels = 21, ranges = c(0, 105), width_type = "linear",
                   width_scale = "variable")
  expect_s3_class(g, "ggplot")
  # linear widths: points of a petal at the same distance from the stem
  pd <- ggplot2::layer_data(g, 2)
  expect_true(all(is.finite(pd$x)))
  g2 <- flower_plot(DDAWASS::Satisfaction, row = 1, width_type = "linear", width_scale = "unit")
  expect_s3_class(g2, "ggplot")
})

test_that("eyes and flowers on given coordinates and on a map", {
  dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
  co <- matrix(stats::rnorm(2 * length(dv$rows)), ncol = 2, dimnames = list(dv$rows, NULL))
  bm <- data.frame(long = c(-3, 3, 3, -3), lat = c(-3, -3, 3, 3), group = 1)
  g <- gei_pca(dv, co, basemap = bm, axis_titles = c("lon", "lat"), key_position = c(-2.5, -2))
  expect_s3_class(g, "ggplot")
  expect_equal(g$labels$x, "lon")
  f <- flower_pca(dv, co, basemap = bm, key_position = c(-2.5, -2),
                  groups = rep(1:2, length.out = length(dv$rows)))
  expect_s3_class(f, "ggplot")
  f2 <- flower_pca(HistDAWass::BLOOD, n_levels = 15, key = "none")
  expect_s3_class(f2, "ggplot")
  expect_error(gei_pca(dv, co, basemap = data.frame(a = 1, b = 2)))
})
