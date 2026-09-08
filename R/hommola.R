
# Melt a distance object/matrix into flat off-diagonal ordered pairs.
# inputs: d (a dist, distTips, or matrix)
# output: list with row (labels), col (labels), dist (numeric vector)
# Auto-labels rows/cols with integers when dimnames are absent, matching the
# behaviour of the previous tibble-based implementation.
.melt_dist <- function(d) {
  m <- as.matrix(d)
  if (is.null(rownames(m))) rownames(m) <- as.character(seq_len(nrow(m)))
  if (is.null(colnames(m))) colnames(m) <- as.character(seq_len(ncol(m)))
  off <- which(row(m) != col(m))
  list(row = rownames(m)[row(m)[off]],
       col = colnames(m)[col(m)[off]],
       dist = m[off])
}

# Coerce a distance object to a labelled matrix, adding integer dimnames when
# they are missing so downstream label-based indexing is well defined.
# inputs: d (a dist, distTips, or matrix)
# output: numeric matrix with non-null row/col names
.as_labelled_matrix <- function(d) {
  m <- as.matrix(d)
  if (is.null(rownames(m))) rownames(m) <- as.character(seq_len(nrow(m)))
  if (is.null(colnames(m))) colnames(m) <- as.character(seq_len(ncol(m)))
  m
}

# Fast Hommola correlation core shared by hommola() and hommola_wf().
# inputs:
#   sym: precomputed melt from .melt_dist() of the symbiont distances
#   host_mat: labelled host distance matrix from .as_labelled_matrix()
#   symbiont_labels, host_labels: aligned vectors defining the symbiont -> host
#     mapping (host_labels[i] is the host of symbiont_labels[i])
# output: Pearson correlation, or 0 when host distances are degenerate
.hommola_core <- function(sym, host_mat, symbiont_labels, host_labels) {
  s2h <- stats::setNames(host_labels, symbiont_labels)
  host1 <- s2h[sym$row]
  host2 <- s2h[sym$col]
  # index by matched integer position so labels absent from the host matrix
  # yield NA rather than a subscript error
  ri <- match(host1, rownames(host_mat))
  ci <- match(host2, colnames(host_mat))
  host_dist <- host_mat[cbind(ri, ci)]
  # guard against a degenerate mapping (all host distances equal or unresolved),
  # which would make the correlation 0 or undefined
  if (length(unique(host_dist)) == 1) {
    return(0)
  }
  cor(sym$dist, host_dist, method = "pearson")
}

# Vectorized permutation correlations for hommola_wf.
# Computes the Hommola correlation for every column of permutations_mat at once
# using integer indexing (no per-permutation string lookups) and a single
# column-wise cor() call.
# inputs:
#   sym: precomputed melt from .melt_dist() of the symbiont distances
#   host_mat: labelled host distance matrix from .as_labelled_matrix()
#   symbiont_labels, host_labels: the (unpermuted) Symbiont/Host columns
#   permutations_mat: N x P integer matrix, each column a permutation of 1:N
# output: numeric vector of P correlations (degenerate host mappings give 0,
#   matching .hommola_core)
.hommola_perm_cors <- function(sym, host_mat, symbiont_labels, host_labels,
                               permutations_mat) {
  n <- length(symbiont_labels)
  p <- ncol(permutations_mat)
  m <- length(sym$dist)
  nrow_h <- nrow(host_mat)

  # fixed per-pair symbiont positions and per-row host-matrix indices
  dr <- match(sym$row, symbiont_labels)
  dc <- match(sym$col, symbiont_labels)
  hidx <- match(host_labels, rownames(host_mat))

  # inverse of each permutation: inv[perm[j], p] = j
  inv_mat <- matrix(0L, n, p)
  inv_mat[cbind(as.vector(permutations_mat), rep(seq_len(p), each = n))] <-
    rep(seq_len(n), p)

  # host-matrix index assigned to each symbiont position, per permutation
  g_mat <- matrix(hidx[inv_mat], n, p)

  # PERFORMANCE: the permuted host distances form an m x p matrix, and m grows
  # with the square of the clade size. A 3,135-tip clade at 999 permutations
  # asks for 4.9M x 999 doubles, or 39 TB, and the allocation fails. Chunking
  # the permutation columns bounds peak memory at the budget below; the
  # arithmetic is unchanged, so results match an unchunked run exactly.
  chunk_elements <- 8e6                            # about 64 MB per matrix
  chunk_p <- max(1L, min(p, as.integer(chunk_elements %/% max(1, m))))

  cors <- numeric(p)
  for (start in seq(1L, p, by = chunk_p)) {
    idx <- start:min(start + chunk_p - 1L, p)
    ri <- g_mat[dr, idx, drop = FALSE]             # m x length(idx) row indices
    ci <- g_mat[dc, idx, drop = FALSE]             # m x length(idx) col indices
    host_dist <- matrix(host_mat[(ci - 1L) * nrow_h + ri], m, length(idx))
    cors[idx] <- suppressWarnings(as.numeric(cor(sym$dist, host_dist)))
  }
  # a constant host-distance column gives NA; match .hommola_core by mapping
  # those degenerate permutations to 0 (only when the symbiont distances vary)
  if (stats::sd(sym$dist) > 0) {
    cors[is.na(cors)] <- 0
  }
  cors
}

#' Calculate Hommola correlation coefficient
#'
#' Computes Pearson correlation between pairwise phylogenetic distances in
#' host and symbiont trees. Tests the null hypothesis of no association
#' between host and symbiont topology (Hommola et al. 2011).
#'
#' @param i_host_subtree_dist A distance matrix (e.g., from `dist()` or
#'   `adephylo::distTips()`) for the host subtree tips
#' @param i_symbiont_subtree_dist A distance matrix for the symbiont subtree tips
#' @param i_host_to_symbiont_df A data frame with "Host" and "Symbiont" columns
#'   linking each symbiont to its host
#'
#' @return A numeric correlation coefficient (Pearson's r) ranging from -1 to 1.
#'   Higher positive values indicate symbiont topology more closely matches
#'   host topology. Returns 0 if all host distances are identical (degenerate case).
#'
#' @importFrom stats cor setNames
#' @export
#'
#' @examples
#' \donttest{
#' library(adephylo)
#' h_tree <- ape::rtree(5)
#' s_tree <- ape::rtree(15)
#' hs_df <- data.frame(
#'   Host = rep(h_tree$tip.label, 3),
#'   Symbiont = s_tree$tip.label
#' )
#' h_dist <- distTips(h_tree, method = "patristic")
#' s_dist <- distTips(s_tree, method = "patristic")
#' r <- hommola(h_dist, s_dist, hs_df)
#' }
#'
hommola <- function(i_host_subtree_dist, i_symbiont_subtree_dist, i_host_to_symbiont_df) {
  sym <- .melt_dist(i_symbiont_subtree_dist)
  host_mat <- .as_labelled_matrix(i_host_subtree_dist)
  .hommola_core(sym, host_mat,
                i_host_to_symbiont_df$Symbiont,
                i_host_to_symbiont_df$Host)
}
