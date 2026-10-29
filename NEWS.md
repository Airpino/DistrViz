# DistrViz 1.0.0

* First release.
* Green Eyes plots (`gei_plot()`, `gei_grid()`) and distributional heatmaps
  (`dist_heatmap()`, also interactive with plotly) for `MatH`, `MatD`, raw
  data and long frequency tables.
* New `gei_pca()`: Green Eyes plots placed on a factorial plane (quantile
  PCA/MFA of DDAWASS or HistDAWass, any FactoMineR result or a matrix of
  coordinates); each eye takes at most `size` (default 1/10) of the width and
  of the height of the plane. A key on the right shows the eye of the
  barycenter of all the units (or of each group) with the names of the
  variables around the sectors. With `interactive = TRUE` (ggiraph) the eyes
  have tooltips, can be hidden by clicking or with the lasso (the point of
  application stays) and the plane can be zoomed; switches for the median
  circles and the names and a slider for the size of all the eyes.
* New `flower_plot()` and `flower_grid()`: flower plots, with one petal for
  each value (discrete distributions) or violin petals (histograms), widths
  proportional to the probabilities and colours of the values (port of
  `flower_plot`, `flower_plot2` and `flower_plot_discr`). `width_type =
  "linear"` makes the width of the petals a distance from the stem (exactly
  proportional to the probability); `width_scale = "variable"` scales the
  widths variable by variable.
* Six vignettes: the sixth is on flower plots of histogram data (population
  pyramids, summer climate of Chinese stations).
* New `flower_pca()`: flower plots on a factorial plane or at any coordinates.
* `gei_pca()` and `flower_pca()` draw a map under the eyes/flowers
  (`basemap`, e.g. from `ggplot2::map_data()`), with `axis_titles` and
  `key_position` (the key inside the plot, e.g. in an empty corner of a map).
* New `gei_diff()`: difference eyes, coloured by the difference between the
  quantile functions of each unit and of a reference (by default the
  Wasserstein barycenter; or a unit, or the barycenter of a group). The mean
  squared difference is the contribution of each cell to the Wasserstein
  variance.
* New `dv_var_order()` and `var_order = "tsp"`: the sectors are ordered by the
  shortest Hamiltonian cycle through the variables with the correlation
  distance sqrt(2(1-r)) (exact up to 9 variables, nearest neighbour + 2-opt
  otherwise). `dist_heatmap(cluster_cols = "tsp")` orders the columns by a
  short Hamiltonian path.
* `col_dist = "sqrt(2(1-r))"` in `dv_cluster()` and `dist_heatmap()`.
* `dv_data()` accepts `support` and `min.n` for raw data (passed to
  `DDAWASS::data2MatD()`); the support range is used as colour range.
* `ranges` in `dv_data()` also works for histograms (`MatH`): the classes
  split the given range, e.g. to centre a diverging palette on a value.
* Legends of histogram heatmaps report the values when the classes are common
  to all the variables.
* Missing or `NA` names of rows/columns are replaced by `I1, ...`/`V1, ...`.
* `gei_plot(polar = FALSE)` draws stacked bars without the pupil.
