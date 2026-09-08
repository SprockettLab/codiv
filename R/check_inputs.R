#' Validate inputs to codiv()
#'
#' Checks that phylogenetic trees are valid, parameters are in acceptable
#' ranges, and the host-to-symbiont linking data frame has the correct structure.
#'
#' @param Host_tree A phylo object; binary phylogenetic tree of hosts
#' @param Symbiont_tree A phylo object; rooted, binary phylogenetic tree of
#'   symbionts. Rooting is required: see Details.
#' @param Host_to_Symbiont_df A data frame with columns "Host" and "Symbiont"
#'   that links each symbiont tip to the host it was isolated from
#' @param min_hosts Minimum number of hosts required for a node to be included
#'   in the scan; must be >= 3
#' @param min_symbiont_tips Minimum number of symbiont tips in a subtree;
#'   default 7 is recommended
#' @param span_fraction Fraction of total tree span; filters nodes with span
#'   <= this fraction of maximum span (0-1)
#' @param permutations Number of permutations for significance testing;
#'   default 99 recommended, higher values more accurate
#'
#' @return Invisibly returns TRUE if all checks pass; stops with error message
#'   if validation fails
#'
#' @details
#' `Symbiont_tree` must be rooted. [codiv()] tests one clade per internal node,
#' so the rooting decides which clades exist to be tested. An unrooted tree read
#' by [ape::read.tree()] still has a node sitting in the root position, picked
#' by the order of the Newick string, so a scan would run and report an
#' arbitrary set of clades rather than failing. Root the tree before resolving
#' polytomies: [ape::multi2di()] applied to an unrooted tree splits the basal
#' polytomy with a zero-length branch, after which [ape::is.rooted()] returns
#' TRUE for a tree nobody rooted.
#' @export
#'
#' @examples
#' Host_tree <- ape::rtree(10)
#' Symbiont_tree <- ape::rtree(100)
#' Host_to_Symbiont_df <- data.frame(
#'   "Host" = rep(Host_tree$tip.label, 10),
#'   "Symbiont" = Symbiont_tree$tip.label
#' )
#' check_inputs(Host_tree, Symbiont_tree, Host_to_Symbiont_df, 4, 7, 0.25, 99)
#'
check_inputs <- function(Host_tree, Symbiont_tree, Host_to_Symbiont_df,
                         min_hosts, min_symbiont_tips, span_fraction,
                         permutations) {
  if (!inherits(Host_tree, "phylo") || !inherits(Symbiont_tree, "phylo")) {
    stop("`Host_tree` and `Symbiont_tree` must be phylo objects")
  }

  # count parameters must be single, positive whole numbers
  .is_count <- function(x) {
    is.numeric(x) && length(x) == 1 && !is.na(x) && x > 0 && x == round(x)
  }
  if (!.is_count(min_hosts)) {
    stop("`min_hosts` must be a single positive integer")
  }
  if (!.is_count(min_symbiont_tips)) {
    stop("`min_symbiont_tips` must be a single positive integer")
  }
  if (!.is_count(permutations)) {
    stop("`permutations` must be a single positive integer")
  }
  if (!is.numeric(span_fraction) || length(span_fraction) != 1 ||
      is.na(span_fraction)) {
    stop("`span_fraction` must be a single number between 0 and 1")
  }

  if (min_hosts < 3) {
    stop("`min_hosts` must be >= 3")
  }

  if (span_fraction <= 0 || span_fraction > 1) {
    stop("`span_fraction` must be between 0 (exclusive) and 1")
  }

  if (min_symbiont_tips < 7) {
    warning("Recommend setting `min_symbiont_tips` to 7 or more")
  }

  if (permutations < 99) {
    warning("Recommend setting `permutations` to 99 or more for statistical validity")
  }

  # both methods rely on branch lengths for patristic distances
  if (is.null(Host_tree$edge.length) || is.null(Symbiont_tree$edge.length)) {
    stop("Both trees must have branch lengths (`edge.length`)")
  }

  # tip labels must be unique so symbionts and hosts map unambiguously
  if (anyDuplicated(Host_tree$tip.label)) {
    stop("`Host_tree` has duplicated tip labels")
  }
  if (anyDuplicated(Symbiont_tree$tip.label)) {
    stop("`Symbiont_tree` has duplicated tip labels")
  }

  # The scan walks the symbiont tree node by node, so the root decides which
  # clades exist to be tested at all. An unrooted tree read by ape still has a
  # node in the root position, chosen by the order of the Newick string rather
  # than by anything biological, so the scan would run and quietly report an
  # arbitrary set of clades. Rooting is a judgement call, so it is left to the
  # caller instead of being applied here.
  if (!ape::is.rooted(Symbiont_tree)) {
    stop("`Symbiont_tree` is unrooted. The scan tests one clade per internal ",
         "node, so the rooting determines which clades exist and cannot be ",
         "chosen automatically. Root it first, for example with ",
         "castor::root_at_midpoint(tree), ape::root(tree, outgroup = ...), or ",
         "phangorn::midpoint(tree). Root before resolving polytomies, not ",
         "after.")
  }

  # ape::multi2di() on an unrooted tree splits the basal polytomy with a
  # zero-length branch, after which is.rooted() reports TRUE for a tree nobody
  # rooted. A real root normally has positive length on both sides, so a
  # zero-length root edge is worth flagging.
  root_node <- length(Symbiont_tree$tip.label) + 1L
  root_edges <- Symbiont_tree$edge.length[Symbiont_tree$edge[, 1] == root_node]
  if (length(root_edges) > 0 && any(root_edges == 0)) {
    warning("`Symbiont_tree` has a zero-length branch at the root, which is ",
            "what ape::multi2di() leaves behind when it is used on an ",
            "unrooted tree. If that is what happened, the root is arbitrary: ",
            "root the tree first, then resolve polytomies.")
  }

  if (!ape::is.rooted(Host_tree)) {
    warning("`Host_tree` is unrooted. Patristic distances are unaffected, but ",
            "the Host_Colless and Host_Sackin shape metrics assume a root.")
  }

  # tree-shape and co-phylogeny methods assume fully bifurcating trees
  if (!ape::is.binary(Host_tree) || !ape::is.binary(Symbiont_tree)) {
    warning("One or both trees are not fully bifurcating; resolve polytomies ",
            "with ape::multi2di() for correct tree-shape metrics")
  }

  if (!all(c("Host", "Symbiont") %in% colnames(Host_to_Symbiont_df))) {
    stop("`Host_to_Symbiont_df` must contain columns named 'Host' and 'Symbiont'")
  }

  if (anyNA(Host_to_Symbiont_df$Host) || anyNA(Host_to_Symbiont_df$Symbiont)) {
    stop("`Host_to_Symbiont_df` must not contain NA in 'Host' or 'Symbiont'")
  }

  if (!all(Host_to_Symbiont_df$Symbiont %in% Symbiont_tree$tip.label)) {
    stop("Some symbionts in `Host_to_Symbiont_df` are not tips in `Symbiont_tree`")
  }

  if (!all(Host_to_Symbiont_df$Host %in% Host_tree$tip.label)) {
    stop("Some hosts in `Host_to_Symbiont_df` are not tips in `Host_tree`")
  }

  invisible(TRUE)
}
