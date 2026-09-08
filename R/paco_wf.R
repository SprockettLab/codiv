
#' Phylogenetic Association with Co-diversification (PACo) analysis
#'
#' Performs PACo (Phylogenetic Association with Co-diversification) analysis
#' to test host-symbiont co-evolution. Uses principal coordinates to reduce
#' phylogenetic distance matrices and computes a goodness-of-fit statistic.
#'
#' @param i_host_subtree A phylo object for the host subtree
#' @param i_symbiont_subtree A phylo object for the symbiont subtree
#' @param i_host_to_symbiont_mat A binary matrix with hosts as rows and
#'   symbionts as columns, where 1 indicates an association
#' @param permutations Number of randomizations for significance testing;
#'   at least 99 recommended
#' @param seed Random seed for reproducibility; if NA, no seed is set
#'
#' @return A PACo object containing:
#'   - gof: Goodness-of-fit statistics including ss (squared sum) and
#'     p-value from permutation test
#'   - Additional fields from paco::PACo()
#'
#' @references
#' Balbuena JA, Míguez-Argüello R, Blasco-Costa I, Llopis-Beltrán A (2013).
#' PACo: A novel procedure to estimate the packedness of a given association
#' matrix. Journal of Biogeography 40: 948-960.
#'
#' @importFrom paco prepare_paco_data add_pcoord PACo
#' @importFrom stats cophenetic
#' @export
#'
#' @examples
#' \donttest{
#' h_tree <- ape::rtree(5)
#' s_tree <- ape::rtree(15)
#' links <- data.frame(Host = rep(h_tree$tip.label, 3),
#'                     Symbiont = s_tree$tip.label)
#' hs_mat <- host_symbiont_links(links)   # binary host-by-symbiont matrix
#' result <- paco_wf(h_tree, s_tree, hs_mat, permutations = 99, seed = 123)
#' }
#'
paco_wf <- function(i_host_subtree, i_symbiont_subtree, i_host_to_symbiont_mat, permutations, seed) {
  gdist <- cophenetic(i_host_subtree)
  ldist <- cophenetic(i_symbiont_subtree)

  # DECISION: rescale the symbiont distances to unit maximum before the
  # principal coordinate step. add_pcoord() discards every eigenvalue below a
  # fixed tolerance, so a clade of near-identical symbionts (branch lengths
  # around 1e-4, common for MAGs of the same strain) collapses to zero
  # dimensions and paco fails with "'dims' cannot be of length 0". add_pcoord()
  # normalizes the symbiont side internally, so scaling only that matrix leaves
  # gof$ss and gof$p unchanged; scaling the host side would not.
  max_ldist <- max(ldist)

  # No symbiont divergence at all means there is nothing for the ordination to
  # place: every tip sits at the same point. codiv() filters these clades out
  # before scanning, but paco_wf() is exported, so guard here too rather than
  # letting add_pcoord() fail with "'dims' cannot be of length 0".
  if (!is.finite(max_ldist) || max_ldist <= 0) {
    warning("Symbiont distances carry no variation; PACo is undefined for ",
            "this clade", call. = FALSE)
    return(list(gof = list(ss = NA_real_, p = NA_real_, n = nrow(ldist))))
  }

  ldist <- ldist / max_ldist

  # add_pcoord() can still fail on a distance matrix that is rank deficient
  # rather than merely small: it drops eigenvalues below a tolerance, and if
  # too few survive the ordination has no dimensions and sweep() errors with
  # "'dims' cannot be of length 0". Rescaling above fixes the case where every
  # distance is tiny, but not the case where the matrix genuinely has too few
  # independent dimensions.
  #
  # Failing here used to abort the whole node, discarding the Hommola, ParaFit
  # and TreeDist results computed alongside it. PACo alone is returned as NA
  # instead, so one undefined ordination costs one statistic rather than four.
  D <- tryCatch({
    D <- prepare_paco_data(gdist, ldist, i_host_to_symbiont_mat)
    D <- add_pcoord(D)
    suppressWarnings(PACo(D, nperm = permutations, seed = seed, method = "swap"))
  }, error = function(e) {
    warning("PACo could not be computed for this clade (", conditionMessage(e),
            "); returning NA for PACo only", call. = FALSE)
    list(gof = list(ss = NA_real_, p = NA_real_, n = nrow(ldist)))
  })
  return(D)
}

