# Node x (horizontal position) and y (vertical position) for a rectangular tree
# layout. Tips are placed at integer y positions in cladewise order; internal
# nodes sit at the mean of their children.
# inputs: tree (phylo); use_edge_length (if FALSE, a cladogram layout with tips
#   aligned and internal nodes spread by topology, avoiding a long root edge)
# output: list(x, y) numeric vectors indexed by node id
.tree_xy <- function(tree, use_edge_length = FALSE) {
  n_tip <- length(tree$tip.label)
  n_all <- n_tip + tree$Nnode

  if (use_edge_length) {
    if (is.null(tree$edge.length)) {
      tree$edge.length <- rep(1, nrow(tree$edge))
    }
    x <- ape::node.depth.edgelength(tree)
  } else {
    # cladogram: each node's height is the max edges to a descendant tip; place
    # x so tips align on the right and internal nodes spread evenly by topology
    he <- numeric(n_all)
    po <- stats::reorder(tree, "postorder")
    for (k in seq_len(nrow(po$edge))) {
      pk <- po$edge[k, 1]
      ck <- po$edge[k, 2]
      he[pk] <- max(he[pk], he[ck] + 1)
    }
    x <- max(he) - he
  }

  cw <- stats::reorder(tree, "cladewise")
  tip_order <- cw$edge[cw$edge[, 2] <= n_tip, 2]
  y <- numeric(n_all)
  y[tip_order] <- seq_along(tip_order)

  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  internal <- (n_tip + 1L):n_all
  # assign deepest internal nodes first so children are always set beforehand
  for (nd in internal[order(x[internal], decreasing = TRUE)]) {
    y[nd] <- mean(y[child[parent == nd]])
  }
  list(x = x, y = y)
}

# Elbow segments (one horizontal and one vertical per edge) for a tree, with x
# scaled to [0, 1] and y scaled to [0, 1].
# inputs: tree; xy from .tree_xy()
# output: data frame with x, y, xend, yend
.tree_segments <- function(tree, xy) {
  n_tip <- length(tree$tip.label)
  xs <- if (max(xy$x) > 0) xy$x / max(xy$x) else xy$x
  ys <- if (n_tip > 1) (xy$y - 1) / (n_tip - 1) else rep(0.5, length(xy$y))
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  rbind(
    data.frame(x = xs[parent], y = ys[child],
               xend = xs[child], yend = ys[child]),   # horizontal
    data.frame(x = xs[parent], y = ys[parent],
               xend = xs[parent], yend = ys[child])    # vertical
  )
}

# Tip label / position table with x and y scaled to [0, 1].
# inputs: tree; xy
# output: data frame with label, x (actual tip depth), y
.tree_tips <- function(tree, xy) {
  n_tip <- length(tree$tip.label)
  xs <- if (max(xy$x) > 0) xy$x[seq_len(n_tip)] / max(xy$x) else xy$x[seq_len(n_tip)]
  ys <- if (n_tip > 1) (xy$y[seq_len(n_tip)] - 1) / (n_tip - 1) else 0.5
  data.frame(label = tree$tip.label, x = xs, y = ys, stringsAsFactors = FALSE)
}

# Symbionts belonging to significant nodes, read from the per-node Newick stored
# in a codiv result. inputs: codiv_results; threshold. output: character vector.
.significant_symbionts <- function(codiv_results, threshold) {
  params <- attr(codiv_results, "codiv_params")
  methods <- if (!is.null(params)) params$methods else
    intersect(c("hommola", "paco", "parafit"),
              sub("_.*$", "", colnames(codiv_results)))
  # significance from the first method's per-node permutation p-value
  sig_col <- NULL
  for (m in methods) {
    mc <- .codiv_method_cols(m, colnames(codiv_results))
    sig_col <- mc[["pvalue"]]
    if (!is.null(sig_col) && !is.na(sig_col)) break
  }
  if (is.null(sig_col) || !("Symbiont_Tree" %in% colnames(codiv_results))) {
    return(character(0))
  }
  sig_rows <- which(codiv_results[[sig_col]] < threshold)
  tips <- unlist(lapply(sig_rows, function(i) {
    nwk <- codiv_results[["Symbiont_Tree"]][i]
    if (is.na(nwk) || !nzchar(nwk)) return(character(0))
    tryCatch(ape::read.tree(text = nwk)$tip.label, error = function(e) character(0))
  }))
  unique(tips)
}

#' Visualize Host-Symbiont Codiversification
#'
#' Draws a tanglegram: the host tree on the left, the symbiont tree mirrored on
#' the right, and a line for each host-symbiont association. Optionally rotates
#' the symbiont tree to align it with the host tree and colours the association
#' lines.
#'
#' @param Host_tree A phylo object; host phylogenetic tree
#' @param Symbiont_tree A phylo object; symbiont phylogenetic tree
#' @param Host_to_Symbiont_df A data frame with columns "Host" and "Symbiont"
#'   linking symbionts to hosts
#' @param codiv_results Optional [codiv()] result; required for
#'   `color_by = "significance"`
#' @param color_by How to colour the association lines: "host" (default, by
#'   host), "none", or "significance" (symbionts in significant nodes)
#' @param significance_threshold p-value threshold for `color_by =
#'   "significance"`; default 0.05
#' @param line_width Width of the association lines; default 0.5
#' @param line_alpha Transparency of the association lines (0-1); default 0.6
#' @param title Plot title; default NULL
#' @param use_procrustes If TRUE (default), rotate the symbiont tree to align it
#'   with the host tree and reduce line crossings
#' @param show_host_labels If TRUE (default), label host tips
#' @param show_symbiont_labels If TRUE, label symbiont tips; default FALSE
#'   (symbiont trees are often large)
#' @param label_size Tip label text size; default 3.2
#' @param use_branch_lengths If FALSE (default), use a cladogram layout (tips
#'   aligned, nodes spread by topology) for an uncluttered view; if TRUE, scale
#'   horizontal position by branch length
#'
#' @return A ggplot object, customizable with further ggplot2 layers and
#'   saveable with [ggplot2::ggsave()].
#'
#' @details
#' Requires the ggplot2 package. The tree layout is computed with ape and the
#' alignment uses [ape::rotateConstr()]. By default a cladogram layout is used
#' (tips aligned, nodes spread by topology); set `use_branch_lengths = TRUE` to
#' position nodes by branch length instead.
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   sim <- simulate_codiv_data(n_hosts = 8, n_clades = 2, seed = 1)
#'   plot_codiv_trees(sim$host_tree, sim$symbiont_tree, sim$links,
#'                    color_by = "host")
#' }
#' }
#'
plot_codiv_trees <- function(Host_tree, Symbiont_tree, Host_to_Symbiont_df,
                             codiv_results = NULL,
                             color_by = c("host", "none", "significance"),
                             significance_threshold = 0.05, line_width = 0.5,
                             line_alpha = 0.6, title = NULL,
                             use_procrustes = TRUE, show_host_labels = TRUE,
                             show_symbiont_labels = FALSE, label_size = 3.2,
                             use_branch_lengths = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_codiv_trees() requires the 'ggplot2' package. ",
         "Install it with install.packages('ggplot2').")
  }
  color_by <- match.arg(color_by)

  Host_to_Symbiont_df$Host <- as.character(Host_to_Symbiont_df$Host)
  Host_to_Symbiont_df$Symbiont <- as.character(Host_to_Symbiont_df$Symbiont)

  if (color_by == "significance" && is.null(codiv_results)) {
    stop("color_by = 'significance' requires `codiv_results`")
  }

  # keep only tips that participate in associations
  Host_tree <- ape::keep.tip(Host_tree,
    intersect(Host_tree$tip.label, Host_to_Symbiont_df$Host))
  Symbiont_tree <- ape::keep.tip(Symbiont_tree,
    intersect(Symbiont_tree$tip.label, Host_to_Symbiont_df$Symbiont))

  # host layout gives the target tip order used to align the symbiont tree
  host_xy <- .tree_xy(Host_tree, use_branch_lengths)
  host_tips <- .tree_tips(Host_tree, host_xy)

  if (use_procrustes) {
    host_y <- stats::setNames(host_tips$y, host_tips$label)
    # order symbionts by their host's vertical position, then rotate to match
    sym_host_y <- host_y[Host_to_Symbiont_df$Host[
      match(Symbiont_tree$tip.label, Host_to_Symbiont_df$Symbiont)]]
    target <- Symbiont_tree$tip.label[order(sym_host_y)]
    Symbiont_tree <- tryCatch(ape::rotateConstr(Symbiont_tree, target),
                              error = function(e) Symbiont_tree)
  }

  sym_xy <- .tree_xy(Symbiont_tree, use_branch_lengths)
  sym_tips <- .tree_tips(Symbiont_tree, sym_xy)

  # Layout across x: host tree [0, tw], label gutter [tw, 1], association band
  # [1, 2], label gutter [2, 3 - tw], symbiont tree [3 - tw, 3]. Pulling the
  # trees back to width `tw` leaves room for tip labels that do not overlap them.
  tw <- 0.78
  host_seg <- .tree_segments(Host_tree, host_xy)
  host_seg$x <- host_seg$x * tw
  host_seg$xend <- host_seg$xend * tw
  sym_seg <- .tree_segments(Symbiont_tree, sym_xy)
  sym_seg$x <- 3 - sym_seg$x * tw
  sym_seg$xend <- 3 - sym_seg$xend * tw

  # tip alignment: extend every tip out to a common line (host at x = 1,
  # symbiont at x = 2) so tips line up and the association lines meet them
  host_ext <- data.frame(x = host_tips$x * tw, y = host_tips$y,
                         xend = 1, yend = host_tips$y)
  sym_ext <- data.frame(x = 3 - sym_tips$x * tw, y = sym_tips$y,
                        xend = 2, yend = sym_tips$y)
  host_tips$x <- 1
  sym_tips$x <- 2

  # association lines from host tip line to symbiont tip line
  links <- Host_to_Symbiont_df
  links$y <- host_tips$y[match(links$Host, host_tips$label)]
  links$yend <- sym_tips$y[match(links$Symbiont, sym_tips$label)]
  links$x <- 1
  links$xend <- 2
  links <- links[!is.na(links$y) & !is.na(links$yend), ]

  if (color_by == "host") {
    links$link_color <- links$Host
  } else if (color_by == "significance") {
    sig_symb <- .significant_symbionts(codiv_results, significance_threshold)
    links$link_color <- ifelse(links$Symbiont %in% sig_symb,
                               "significant", "not significant")
  } else {
    links$link_color <- "association"
  }

  p <- ggplot2::ggplot() +
    # tip-alignment extensions (light dotted guides)
    ggplot2::geom_segment(
      data = host_ext,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.25, linetype = "dotted", color = "grey60") +
    ggplot2::geom_segment(
      data = sym_ext,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.25, linetype = "dotted", color = "grey60") +
    ggplot2::geom_segment(
      data = host_seg,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.4) +
    ggplot2::geom_segment(
      data = sym_seg,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.4)

  if (color_by == "none") {
    p <- p + ggplot2::geom_segment(
      data = links,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = line_width, alpha = line_alpha, color = "#4477AA")
  } else {
    p <- p + ggplot2::geom_segment(
      data = links,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = link_color),
      linewidth = line_width, alpha = line_alpha) +
      ggplot2::labs(color = if (color_by == "host") "Host" else "Significance")
  }

  if (show_host_labels) {
    # labels sit in the gutter just left of the host tip line, clear of the tree
    p <- p + ggplot2::geom_text(
      data = host_tips,
      ggplot2::aes(x = x, y = y, label = label),
      hjust = 1, nudge_x = -0.03, nudge_y = 0.012, size = label_size)
  }
  if (show_symbiont_labels) {
    p <- p + ggplot2::geom_text(
      data = sym_tips,
      ggplot2::aes(x = x, y = y, label = label),
      hjust = 0, nudge_x = 0.03, nudge_y = 0.012, size = label_size)
  }

  p +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
}

# columns referenced inside ggplot2::aes(); declared to avoid R CMD check NOTES
utils::globalVariables(c("x", "y", "xend", "yend", "label", "link_color"))
