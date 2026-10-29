## Green Eyes on a factorial plane -----------------------------------------

#' Green Eyes plots on a factorial plane
#'
#' Draws the Green Eyes plot (see \code{\link{gei_plot}}) of each unit at the
#' position of the unit on a factorial plane, e.g. the first two components
#' of a quantile-based PCA or MFA of the distributional data. The eyes are
#' drawn directly in the coordinates of the plane (as polygons), so the plot
#' is a single \code{ggplot} with an isometric scale. The diameter of each eye
#' is \code{size} times the smallest of the width and the height of the cloud
#' of points: with the default \code{size = 0.1} an eye never takes more than
#' one tenth of the width or of the height of the factorial plane.
#'
#' @inheritParams gei_plot
#' @param pca (optional) the factorial analysis: an object with the
#'   coordinates of the units in \code{$ind$coord} (as returned by
#'   \code{DDAWASS::WD.MultiplePCA}, \code{DDAWASS::WD.1d.PCA},
#'   \code{HistDAWass::WH.MultiplePCA}, \code{HistDAWass::WH.1d.PCA} or
#'   by \pkg{FactoMineR} functions), or a matrix/data frame of coordinates
#'   with one row per unit. Rows are matched by name when possible. If
#'   \code{NULL}, a quantile-based MFA of all the variables (PCA if there is a
#'   single variable) is computed on the data.
#' @param axes the two axes of the plane.
#' @param rows (optional) the units to draw (default all).
#' @param size diameter of the eyes as a fraction of the smallest of the
#'   width and the height of the cloud of points (default 0.1).
#' @param quantiles number of quantiles used when \code{pca} is computed.
#' @param groups (optional) a vector of group labels for the units (e.g. a
#'   partition from \code{cutree}); the outline of each eye is coloured by
#'   group.
#' @param labels logical; write the unit names under the eyes?
#' @param label_size font size of the unit names.
#' @param arc_points number of points used to draw a full circle.
#' @param show_points logical; draw the points of application (the
#'   coordinates of the units) under the eyes? Always \code{TRUE} in the
#'   interactive version.
#' @param point_size size of the points of application.
#' @param interactive logical; return an interactive plot (an
#'   \pkg{ggiraph} \code{girafe} widget)? Requires the \pkg{ggiraph}
#'   package. See Details.
#' @param hidden (interactive version) units whose eyes are hidden at start
#'   (names or indices).
#' @param hover_fade (interactive version) opacity of the other eyes while
#'   the mouse is over an eye (1 = no fading).
#' @param key the key drawn on the right of the plane, telling which
#'   variable is in each sector: \code{"mean"} (default) the eye of the
#'   barycenter of all the drawn units, \code{"groups"} one eye for the
#'   barycenter (prototype) of each group (requires \code{groups}), or
#'   \code{"none"}. Variable names are written around the key eyes. For
#'   discrete distributions the barycenter is the Wasserstein barycenter
#'   projected on the observed support (\code{DDAWASS::WD.vec.mean.simplified});
#'   for histograms it is the Wasserstein barycenter
#'   (\code{HistDAWass::WH.vec.mean}), discretized with the same classes as
#'   the eyes.
#' @param key_size diameter of the key eye as a fraction of the height of the
#'   plane (with \code{key = "groups"} the key eyes are reduced if needed
#'   to fit in the height).
#' @param basemap (optional) a map drawn under the eyes, when the coordinates
#'   in \code{pca} are geographical (e.g. longitude and latitude): a data
#'   frame of polygons with columns \code{long}, \code{lat} (or \code{x},
#'   \code{y}) and \code{group}, as returned by
#'   \code{ggplot2::map_data()}. When a map is given, the dashed axes of the
#'   factorial plane are not drawn and the limits of the plot include the map.
#' @param key_position (optional) \code{c(x, y)}: the key is drawn inside the
#'   plot, centred at this point (in the coordinates of the plane), instead of
#'   in a box on the right; useful on maps, where some corners are empty.
#' @param axis_titles (optional) the titles of the two axes (e.g.
#'   \code{c("longitude", "latitude")}), also used in the tooltips.
#' @param controls (interactive version) logical; add the controls above the
#'   plot (switches for the median circles and the names, slider for the size
#'   of the eyes)?
#' @param size_range (interactive version) range of the size slider (same
#'   unit as \code{size}).
#' @param width_svg,height_svg (interactive version) size of the drawing, in
#'   inches; by default the height follows the aspect ratio of the plane.
#' @details In the interactive version (\code{interactive = TRUE}):
#' \itemize{
#'   \item hovering an eye (or its point, label or outline) shows a tooltip
#'     with the name of the unit, its group and its coordinates; over a band
#'     of the eye the tooltip also reports the variable, the value (or the
#'     class) and its probability. The other eyes are faded;
#'   \item clicking an eye hides it, leaving its point of application
#'     (drawn in red); clicking the point shows the eye again. The lasso
#'     buttons of the toolbar hide (or show) all the eyes in an area, and the
#'     reset button shows all the eyes;
#'   \item the zoom buttons of the toolbar activate zooming and panning with
#'     the mouse wheel and drag;
#'   \item above the plot, a switch shows or hides the dashed median circles
#'     of all the eyes (initial state: \code{median_line}), a switch shows or
#'     hides the names of the units, and a slider changes the size of all the
#'     eyes (each eye is rescaled around its point of application; the
#'     initial value is \code{size}).
#' }
#' In a Shiny app the hidden units are available in
#' \code{input$<outputId>_selected}.
#' @return a \code{ggplot} object (a \code{girafe} widget if
#'   \code{interactive = TRUE}), with the coordinates used in the attribute
#'   \code{"coord"} and the factorial analysis in the attribute \code{"pca"}.
#' @examples
#' dv <- dv_data(DDAWASS::Satisfaction, ranges = c(0, 10))
#' gei_pca(dv, var_order = "tsp")
#' \donttest{
#' res <- DDAWASS::WD.MultiplePCA(DDAWASS::Satisfaction, 1:4, plots = FALSE)
#' gei_pca(dv, res, axes = c(1, 3), size = 0.08)
#' if (requireNamespace("ggiraph", quietly = TRUE)) {
#'   gei_pca(dv, res, interactive = TRUE, hidden = c("Hotel_01", "Hotel_02"))
#' }
#' }
#' @export
gei_pca <- function(x, pca = NULL, axes = c(1, 2), rows = NULL, size = 0.1,
                    quantiles = 10, palette = "Paquin", reverse = FALSE, pupil = 0.3,
                    pupil_color = "grey40", median_line = TRUE, var_order = NULL,
                    groups = NULL, labels = TRUE, label_size = 2.5, arc_points = 120,
                    show_points = FALSE, point_size = 1.2, interactive = FALSE,
                    hidden = NULL, hover_fade = 0.5, width_svg = 8, height_svg = NULL,
                    controls = TRUE, size_range = c(0.02, 0.4),
                    key = c("mean", "groups", "none"), key_size = 0.3,
                    basemap = NULL, axis_titles = NULL, key_position = NULL, ...) {
  x <- dv_data(x, ...)
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

  # size of the eyes
  wx <- diff(range(xy[rows, 1])); wy <- diff(range(xy[rows, 2]))
  ww <- c(wx, wy)
  w <- if (any(ww > 0)) min(ww[ww > 0]) else 1
  R <- size * w / 2
  vo <- .var_order(x, var_order)

  # polygons of the cdf bands
  df <- x$df[x$df$IDr %in% rows & x$df$freq > 0, ]
  df$xpos <- match(df$IDc, vo)
  rad <- function(y) R * (y + pupil) / (1 + pupil)
  npt <- max(2, ceiling(arc_points / p))
  poly <- .band_polygons(cx = xy[df$IDr, 1], cy = xy[df$IDr, 2],
                         th1 = (df$xpos - 1) / p * 2 * pi, th2 = df$xpos / p * 2 * pi,
                         r1 = rad(df$cdf_lo), r2 = rad(df$cdf), npt = npt)
  poly$u <- df$u[poly$id]

  circ <- function(r, id_rows) {
    th <- seq(0, 2 * pi, length.out = arc_points + 1)
    data.frame(id = rep(id_rows, each = length(th)),
               x = rep(xy[id_rows, 1], each = length(th)) + r * sin(th),
               y = rep(xy[id_rows, 2], each = length(th)) + r * cos(th))
  }
  pup <- circ(rad(0), rows)
  cols <- dv_palette(palette, 64, reverse)

  pct <- .pca_eig(pca)
  lab_ax <- if (!is.null(axis_titles)) rep(axis_titles, length.out = 2) else
    paste0("Dim ", axes, if (!is.null(pct)) sprintf(" (%.1f%%)", pct[axes]) else "")
  ax_tip <- if (!is.null(axis_titles)) lab_ax else paste("Dim", axes)
  bm <- .basemap_df(basemap)

  # interactive elements: every element of an eye has the unit as data_id
  if (interactive && !requireNamespace("ggiraph", quietly = TRUE)) {
    stop("install the 'ggiraph' package for interactive plots")
  }
  geom <- function(name) {
    if (interactive) getExportedValue("ggiraph", paste0(name, "_interactive"))
    else getExportedValue("ggplot2", name)
  }
  imap <- function(map) {
    if (interactive) {
      map[["tooltip"]] <- rlang::quo(.data$tooltip)
      map[["data_id"]] <- rlang::quo(.data$data_id)
    }
    map
  }
  unit_tip <- vapply(seq_len(nr), function(i) {
    paste0("<b>", .html_escape(x$rows[i]), "</b>",
           if (!is.null(groups)) paste0("<br>group: ", .html_escape(groups[i])) else "",
           sprintf("<br>%s: %.3g<br>%s: %.3g", ax_tip[1], xy[i, 1], ax_tip[2], xy[i, 2]))
  }, character(1))
  did <- function(id) x$rows[id]
  poly$data_id <- did(df$IDr[poly$id])
  poly$tooltip <- paste0(unit_tip[df$IDr], "<br><i>", .html_escape(x$vars[df$IDc]), "</i>: ",
                         .html_escape(df$label), sprintf(" (p = %.3g)", df$freq))[poly$id]
  pup$data_id <- did(pup$id); pup$tooltip <- unit_tip[pup$id]
  centres <- data.frame(x = xy[rows, 1], y = xy[rows, 2], data_id = did(rows),
                        tooltip = unit_tip[rows])

  g <- ggplot2::ggplot()
  if (is.null(bm)) {
    g <- g +
      ggplot2::annotate("segment", x = min(xy[rows, 1]) - R * 1.1, xend = max(xy[rows, 1]) + R * 1.1,
                        y = 0, yend = 0, colour = "grey75", linetype = "dashed") +
      ggplot2::annotate("segment", y = min(xy[rows, 2]) - R * 1.6, yend = max(xy[rows, 2]) + R * 1.1,
                        x = 0, xend = 0, colour = "grey75", linetype = "dashed")
  } else {
    g <- g + .basemap_layer(bm)
  }
  if (interactive || show_points) {
    # points of application (under the eyes: visible when an eye is hidden)
    g <- g + geom("geom_point")(data = centres, imap(ggplot2::aes(x = .data$x, y = .data$y)),
                                shape = 21, size = point_size, fill = "grey30", colour = "white",
                                stroke = 0.3)
  }
  g <- g +
    geom("geom_polygon")(data = poly, imap(ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                                         fill = .data$u)), colour = NA) +
    geom("geom_polygon")(data = pup, imap(ggplot2::aes(x = .data$x, y = .data$y, group = .data$id)),
                         fill = pupil_color, colour = NA) +
    ggplot2::scale_fill_gradientn(colours = cols, limits = c(0, 1), guide = "none")
  if (median_line || interactive) {   # interactive: drawn anyway, shown/hidden by a switch
    mc <- circ(rad(0.5), rows)
    mc$data_id <- did(mc$id); mc$tooltip <- unit_tip[mc$id]
    g <- g + geom("geom_path")(data = mc, imap(ggplot2::aes(x = .data$x, y = .data$y, group = .data$id)),
                               linetype = "dashed", linewidth = 0.2, alpha = 0.6)
  }
  if (!is.null(groups)) {
    oc <- circ(rad(1), rows)
    oc$group <- factor(groups[oc$id])
    oc$data_id <- did(oc$id); oc$tooltip <- unit_tip[oc$id]
    g <- g + geom("geom_path")(data = oc, imap(ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                                             colour = .data$group)), linewidth = 0.5) +
      ggplot2::labs(colour = "group")
  }
  if (labels || interactive) {   # interactive: drawn anyway, shown/hidden by a switch
    lb <- data.frame(x = xy[rows, 1], y = xy[rows, 2] - R * 1.08, lab = x$rows[rows],
                     data_id = did(rows), tooltip = unit_tip[rows])
    g <- g + geom("geom_text")(data = lb, imap(ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab)),
                               size = label_size, vjust = 1, colour = "grey20")
  }
  pad <- R * 1.1 + if (labels || interactive) R * 0.5 else 0
  xlim <- range(xy[rows, 1]) + c(-1, 1) * R * 1.1
  ylim <- range(xy[rows, 2]) + c(-pad, R * 1.1)
  if (!is.null(bm)) {
    xlim <- range(c(xlim, bm$x)); ylim <- range(c(ylim, bm$y))
  }
  key <- match.arg(key)
  if (key == "groups" && is.null(groups)) key <- "mean"
  if (key != "none") {
    kk <- .gei_key(x, g, key, groups, rows, vo, xlim, ylim, key_size, pupil, pupil_color,
                   median_line || interactive, arc_points, label_size, interactive,
                   key_position = key_position)
    g <- kk$g; xlim <- kk$xlim; ylim <- kk$ylim
  }
  g <- g + ggplot2::coord_equal(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
    ggplot2::labs(x = lab_ax[1], y = lab_ax[2],
                  caption = paste("Sectors, clockwise from the top:",
                                  paste(x$vars[vo], collapse = ", "))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.caption = ggplot2::element_text(hjust = 0, size = 8, colour = "grey30"))
  if (interactive) {
    if (is.character(hidden)) hidden <- match(hidden, x$rows)
    sel <- if (length(hidden)) x$rows[hidden[hidden %in% rows]] else character(0)
    if (is.null(height_svg)) {
      # follow the aspect ratio of the plane (plus room for axes, legend, caption)
      xr <- diff(range(xy[rows, 1])) + 2.2 * R
      yr <- diff(range(xy[rows, 2])) + 1.1 * R + pad
      leg <- if (!is.null(groups)) 1 else 0
      height_svg <- min(max((width_svg - 0.8 - leg) * yr / xr + 1.3, 3), 3 * width_svg)
    }
    g <- ggiraph::girafe(
      ggobj = g, width_svg = width_svg, height_svg = height_svg,
      options = list(
        ggiraph::opts_tooltip(css = paste0("background:white;color:#222;border:1px solid #999;",
                                           "border-radius:4px;padding:4px 6px;font-size:11px;",
                                           "font-family:sans-serif;"), opacity = 0.95),
        ggiraph::opts_hover(css = ggiraph::girafe_css(
          css = "", area = "stroke:black;stroke-width:0.6px;",
          point = "fill:#d7301f;", text = "font-weight:bold;")),
        ggiraph::opts_hover_inv(css = paste0("opacity:", hover_fade, ";")),
        ggiraph::opts_selection(type = "multiple", only_shiny = FALSE, selected = sel,
                                css = ggiraph::girafe_css(
                                  css = "", area = "opacity:0;pointer-events:none;",
                                  line = "opacity:0;pointer-events:none;",
                                  text = "opacity:0;pointer-events:none;",
                                  point = "fill:#d7301f;stroke:#d7301f;")),
        ggiraph::opts_zoom(min = 1, max = 20),
        ggiraph::opts_toolbar(saveaspng = TRUE, position = "topright")
      ))
    if (controls) {
      g <- htmlwidgets::onRender(g, .gei_pca_js, data = list(
        median = isTRUE(median_line), labels = isTRUE(labels), size = size,
        smin = size_range[1], smax = size_range[2]))
    }
  }
  attr(g, "coord") <- xy
  attr(g, "pca") <- pca
  g
}

.html_escape <- function(s) {
  s <- gsub("&", "&amp;", as.character(s), fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}

# annular sectors as polygons (one polygon per band)
.band_polygons <- function(cx, cy, th1, th2, r1, r2, npt) {
  n <- length(cx)
  t <- seq(0, 1, length.out = npt)
  m <- 2 * npt
  id <- rep(seq_len(n), each = m)
  th_out <- outer(th2 - th1, t) + th1                # n x npt
  th <- cbind(th_out, th_out[, npt:1, drop = FALSE])
  r <- cbind(matrix(r2, n, npt), matrix(r1, n, npt))
  data.frame(id = id,
             x = as.vector(t(cx + r * sin(th))),
             y = as.vector(t(cy + r * cos(th))))
}

.dv_pca <- function(x, quantiles, ncp) {
  src <- x$source
  p <- length(x$vars)
  ncp <- max(ncp, 2)
  if (is(src, "MatH")) {
    if (p == 1) {
      HistDAWass::WH.1d.PCA(src, var = 1, quantiles = quantiles, plots = FALSE)
    } else {
      grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
      HistDAWass::WH.MultiplePCA(src, seq_len(p), quantiles = quantiles, graph = FALSE)
    }
  } else {
    if (p == 1) DDAWASS::WD.1d.PCA(src, var = 1, quantiles = quantiles, plots = FALSE, ncp = ncp)
    else DDAWASS::WD.MultiplePCA(src, seq_len(p), quantiles = quantiles, ncp = ncp, plots = FALSE)
  }
}

.pca_coord <- function(pca, rn) {
  co <- if (is.list(pca) && !is.data.frame(pca) && !is.null(pca$ind$coord)) pca$ind$coord else pca
  co <- as.matrix(co)
  if (!is.numeric(co) || ncol(co) < 2) stop("pca must contain the coordinates of the units on at least two axes")
  if (!is.null(rownames(co)) && all(rn %in% rownames(co))) {
    co <- co[rn, , drop = FALSE]
  } else if (nrow(co) != length(rn)) {
    stop("the coordinates must have one row per unit")
  }
  co
}

.pca_eig <- function(pca) {
  if (is.list(pca) && !is.data.frame(pca) && !is.null(pca$eig)) return(pca$eig[, 2])
  NULL
}

# JavaScript for the controls of the interactive gei_pca (htmlwidgets::onRender)
.gei_pca_js <- "
function(el, x, opt) {
  var setup = function() {
    var svg = el.querySelector('svg');
    if (!svg) { setTimeout(setup, 50); return; }
    if (el.querySelector('.gei-controls')) return;
    // centres of the eyes: the points of application
    var ctr = {};
    svg.querySelectorAll('circle[data-id]').forEach(function(c) {
      ctr[c.getAttribute('data-id')] = [parseFloat(c.getAttribute('cx')), parseFloat(c.getAttribute('cy'))];
    });
    var shapes = [], texts = [];
    svg.querySelectorAll('polygon[data-id], polyline[data-id], path[data-id]').forEach(function(e) {
      var c = ctr[e.getAttribute('data-id')];
      if (c) { e.style.vectorEffect = 'non-scaling-stroke'; shapes.push([e, c]); }
    });
    svg.querySelectorAll('text[data-id]').forEach(function(e) {
      var c = ctr[e.getAttribute('data-id')];
      if (c) { var b = e.getBBox(); texts.push([e, c, b.y - c[1], e.getAttribute('transform') || '']); }
    });
    var medians = svg.querySelectorAll('polyline[data-id][stroke-dasharray]');
    svg.querySelectorAll('[data-id=__key__]').forEach(function(e) { e.style.pointerEvents = 'none'; });
    // stack the elements eye by eye (ggplot draws them layer by layer), so
    // that overlapping eyes do not mix; the eye under the mouse comes on top
    var byUnit = {}, order = [];
    shapes.forEach(function(a) {
      var id = a[0].getAttribute('data-id');
      if (!byUnit[id]) { byUnit[id] = []; order.push(id); }
      byUnit[id].push(a[0]);
    });
    var firstText = texts.length ? texts[0][0] : null;
    var raise = function(id) {
      byUnit[id].forEach(function(e) {
        var par = e.parentNode;
        if (firstText && firstText.parentNode === par) par.insertBefore(e, firstText);
        else par.appendChild(e);
      });
    };
    order.forEach(raise);
    var current = null;
    shapes.forEach(function(a) {
      a[0].addEventListener('mouseenter', function() {
        var id = a[0].getAttribute('data-id');
        if (id !== current) { current = id; raise(id); }
      });
    });

    var box = document.createElement('div');
    box.className = 'gei-controls';
    box.style.cssText = 'font:12px sans-serif;color:#333;display:flex;flex-wrap:wrap;gap:16px;' +
                        'align-items:center;padding:4px 8px;';
    box.innerHTML =
      '<label><input type=checkbox class=gei-med> median circles</label>' +
      '<label><input type=checkbox class=gei-lab> names</label>' +
      '<label>eye size <input type=range class=gei-size style=vertical-align:middle;width:180px></label>' +
      '<span class=gei-val></span>';
    el.insertBefore(box, el.firstChild);
    var med = box.querySelector('.gei-med'), lab = box.querySelector('.gei-lab'),
        sz = box.querySelector('.gei-size'), val = box.querySelector('.gei-val');
    sz.min = opt.smin; sz.max = opt.smax; sz.step = 'any'; sz.value = opt.size;
    med.checked = opt.median; lab.checked = opt.labels;
    var showMed = function() {
      medians.forEach(function(m) { m.style.display = med.checked ? '' : 'none'; });
    };
    var showLab = function() {
      texts.forEach(function(t) { t[0].style.display = lab.checked ? '' : 'none'; });
    };
    var resize = function() {
      var s = parseFloat(sz.value) / opt.size;
      shapes.forEach(function(a) {
        var c = a[1];
        a[0].setAttribute('transform', 'translate(' + c[0] + ',' + c[1] + ') scale(' + s +
                          ') translate(' + (-c[0]) + ',' + (-c[1]) + ')');
      });
      texts.forEach(function(t) {
        t[0].setAttribute('transform', 'translate(0,' + ((s - 1) * t[2]) + ') ' + t[3]);
      });
      val.textContent = 'size = ' + parseFloat(sz.value).toFixed(3) + ' (x' + s.toFixed(2) + ')';
    };
    med.addEventListener('change', showMed);
    lab.addEventListener('change', showLab);
    sz.addEventListener('input', resize);
    showMed(); showLab(); resize();
  };
  setup();
}
"

# key eyes: barycenter of all the units, or of each group ----------------------
.gei_key <- function(x, g, key, groups, rows, vo, xlim, ylim, key_size, pupil, pupil_color,
                     median_line, arc_points, label_size, interactive = FALSE,
                     key_position = NULL) {
  p <- length(x$vars)
  sets <- if (key == "groups") {
    gr <- groups[rows]
    lv <- if (is.factor(gr)) levels(droplevels(gr)) else sort(unique(gr))
    stats::setNames(lapply(lv, function(l) rows[gr == l]), paste("Group", lv))
  } else list("Mean" = rows)
  nk <- length(sets)
  h <- diff(ylim)
  # room for: key eye, variable names around it, title above it
  # approximate width of a character as a fraction of the total width
  # (assuming a plot about 165 mm wide), to leave room for the names
  L <- max(nchar(x$vars))
  kc <- label_size * 1.1 * 0.55 / 200
  a <- min(2 * L * kc, 0.5)
  with_lab <- 1   # names around the first key eye only (same order in all the eyes)
  geom_key <- function(Rk) {
    b0 <- xlim[2] + 0.3 * Rk
    x1 <- (b0 + 2.56 * Rk - a * xlim[1]) / (1 - a)
    cw <- kc * (x1 - xlim[1])
    list(b0 = b0, x1 = x1, cw = cw, ch = cw * 1.7,
         xk = b0 + 0.2 * Rk + Rk * 1.08 + L * cw)
  }
  Rk <- min(key_size * h / 2, h / nk / 2 / 1.7)
  for (it in 1:3) {   # the height of the text depends on the total width
    gk <- geom_key(Rk)
    need <- nk * 2.25 * Rk + with_lab * 2.4 * gk$ch + nk * 1.6 * gk$ch
    if (need > h) Rk <- Rk * h / need
  }
  gk <- geom_key(Rk)
  b0 <- gk$b0; x1 <- gk$x1; xk <- gk$xk; ch <- gk$ch
  # heights of the slots: eye + names (if any) + title
  hs <- 2.25 * Rk + 1.6 * ch + ifelse(seq_len(nk) <= with_lab, 2.4 * ch, 0)
  gap <- (h - sum(hs)) / nk
  top <- ylim[2] - gap / 2 - c(0, cumsum(hs + gap))[seq_len(nk)]
  yk <- top - 1.6 * ch - ifelse(seq_len(nk) <= with_lab, 1.2 * ch, 0) - 1.12 * Rk
  box_y <- ylim
  if (!is.null(key_position)) {
    # key inside the plot: key_position is the centre of the first key eye
    Rk <- key_size * h / 2 / max(1, nk / 2)
    cw <- kc * diff(xlim); ch <- 1.7 * cw
    xk <- key_position[1]
    step <- 2.3 * Rk + 1.8 * ch
    yk <- key_position[2] - c(0, cumsum(rep(step, nk - 1)) + 2.4 * ch * (seq_len(nk - 1) >= 1))[seq_len(nk)]
    b0 <- xk - Rk * 1.1 - L * cw - 0.2 * Rk
    x1 <- xk + Rk * 1.1 + L * cw + 0.2 * Rk
    box_y <- c(yk[nk] - Rk * 1.1 - 1.4 * ch - 0.2 * Rk, yk[1] + Rk * 1.1 + 3 * ch + 0.2 * Rk)
  }
  rad <- function(y) Rk * (y + pupil) / (1 + pupil)
  th_mid <- (seq_len(p) - 0.5) / p * 2 * pi
  npt <- max(2, ceiling(arc_points / p))
  polys <- list(); pups <- list(); meds <- list(); outl <- list(); labs <- list(); tit <- list()
  for (k in seq_len(nk)) {
    cell <- .bary_cells(x, sets[[k]])
    cell <- cell[cell$freq > 0, ]
    cell$xpos <- match(cell$IDc, vo)
    pg <- .band_polygons(cx = rep(xk, nrow(cell)), cy = rep(yk[k], nrow(cell)),
                         th1 = (cell$xpos - 1) / p * 2 * pi, th2 = cell$xpos / p * 2 * pi,
                         r1 = rad(cell$cdf_lo), r2 = rad(cell$cdf), npt = npt)
    pg$u <- cell$u[pg$id]; pg$id <- paste(k, pg$id)
    polys[[k]] <- pg
    th <- seq(0, 2 * pi, length.out = arc_points + 1)
    pups[[k]] <- data.frame(id = k, x = xk + rad(0) * sin(th), y = yk[k] + rad(0) * cos(th))
    meds[[k]] <- data.frame(id = k, x = xk + rad(0.5) * sin(th), y = yk[k] + rad(0.5) * cos(th))
    outl[[k]] <- data.frame(id = k, group = sub("^Group ", "", names(sets)[k]),
                            x = xk + rad(1) * sin(th), y = yk[k] + rad(1) * cos(th))
    if (k <= with_lab) {
      r_lab <- Rk * 1.08
      sn <- sin(th_mid); cs <- cos(th_mid)
      labs[[k]] <- data.frame(x = xk + r_lab * sn, y = yk[k] + r_lab * cs, lab = x$vars[vo],
                              hjust = ifelse(sn > 0.15, 0, ifelse(sn < -0.15, 1, 0.5)),
                              vjust = ifelse(cs > 0.15, 0, ifelse(cs < -0.15, 1, 0.5)))
    }
    tit[[k]] <- data.frame(x = xk, y = yk[k] + Rk * 1.1 + if (k <= with_lab) 1.3 * ch else 0,
                           lab = sprintf("%s (n = %d)", names(sets)[k], length(sets[[k]])))
  }
  poly <- do.call(rbind, polys); pup <- do.call(rbind, pups); lab <- do.call(rbind, labs)
  tt <- do.call(rbind, tit)
  g <- g +
    ggplot2::annotate("rect", xmin = b0, xmax = x1, ymin = box_y[1], ymax = box_y[2],
                      fill = "white", colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_polygon(data = poly, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                                    fill = .data$u), colour = NA) +
    ggplot2::geom_polygon(data = pup, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
                          fill = pupil_color, colour = NA)
  if (median_line) {
    # interactive: data_id "__key__" so that the median switch reaches it too
    md <- do.call(rbind, meds)
    g <- g + if (interactive) {
      ggiraph::geom_path_interactive(data = md, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
                                     data_id = "__key__", linetype = "dashed", linewidth = 0.25, alpha = 0.7)
    } else {
      ggplot2::geom_path(data = md, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
                         linetype = "dashed", linewidth = 0.25, alpha = 0.7)
    }
  }
  if (key == "groups") {
    oc <- do.call(rbind, outl)
    oc$group <- factor(oc$group, levels = levels(factor(groups)))
    g <- g + ggplot2::geom_path(data = oc, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                                        colour = .data$group), linewidth = 0.8)
  }
  g <- g +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab,
                                                hjust = .data$hjust, vjust = .data$vjust),
                       size = label_size * 1.1, colour = "grey15") +
    ggplot2::geom_text(data = tt, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
                       size = label_size * 1.15, fontface = "bold", vjust = 0, colour = "grey15")
  if (!is.null(key_position)) {
    return(list(g = g, xlim = range(c(xlim, b0, x1)), ylim = range(c(ylim, box_y))))
  }
  list(g = g, xlim = c(xlim[1], x1), ylim = ylim)
}

# discretized barycenter of a set of units, in the format of x$df
.bary_cells <- function(x, set) {
  src <- x$source
  p <- length(x$vars)
  out <- vector("list", p)
  if (is(src, "MatH")) {
    for (j in seq_len(p)) {
      m <- HistDAWass::WH.vec.mean(src[set, j])
      br <- x$breaks[[j]]
      cdf <- vapply(br, function(q) if (q == -Inf) 0 else if (q == Inf) 1 else HistDAWass::compP(m, q), 0)
      cdf[1] <- 0; cdf[length(cdf)] <- 1
      pr <- pmax(diff(cdf), 0); L <- length(pr)
      fin <- signif(br, 4)
      out[[j]] <- data.frame(IDc = j, value = seq_len(L),
                             label = sprintf("[%s, %s]", fin[-length(fin)], fin[-1]),
                             u = (seq_len(L) - 1) / (L - 1), freq = pr,
                             cdf_lo = cumsum(pr) - pr, cdf = cumsum(pr))
    }
  } else {
    sub <- src[set, , drop = FALSE]
    rg <- x$ranges
    for (j in seq_len(p)) {
      sup <- sort(unique(unlist(lapply(src@M[, j], function(d) d@x))))
      m <- DDAWASS::WD.vec.mean.simplified(sub[, j], support = sup)
      w <- rg[j, 2] - rg[j, 1]; if (w == 0) w <- 1
      out[[j]] <- data.frame(IDc = j, value = m@x, label = format(signif(m@x, 4), trim = TRUE),
                             u = pmin(pmax((m@x - rg[j, 1]) / w, 0), 1), freq = m@p,
                             cdf_lo = cumsum(m@p) - m@p, cdf = cumsum(m@p))
    }
  }
  do.call(rbind, out)
}

# map polygons: a data frame with long/lat (or x/y) and group
.basemap_df <- function(basemap) {
  if (is.null(basemap)) return(NULL)
  b <- as.data.frame(basemap)
  if (all(c("long", "lat") %in% names(b))) { b$x <- b$long; b$y <- b$lat }
  if (!all(c("x", "y") %in% names(b))) stop("basemap must have columns long, lat (or x, y)")
  if (!"group" %in% names(b)) b$group <- 1
  b[, c("x", "y", "group")]
}

.basemap_layer <- function(bm, fill = "grey95", colour = "grey60") {
  ggplot2::geom_polygon(data = bm, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
                        fill = fill, colour = colour, linewidth = 0.25, inherit.aes = FALSE)
}

