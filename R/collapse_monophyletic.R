#' Collapse monophyletic groups in a phylogenetic tree
#'
#' Identifies groups of symbiont tips that are monophyletic (all descended
#' from the same host) and reduces each such group to a single representative
#' tip. Used during codiversification analysis to remove within-host variation
#' and focus on between-host branching patterns.
#'
#' Collapsing does not depend on branch lengths. A monophyletic group of tips
#' from one host is reduced to a single tip whether its members sit on long
#' branches, on branches of equal length, on zero-length branches (identical
#' sequences), or on a tree with no branch lengths at all.
#'
#' @param symbiont_subtree A phylo object representing a subtree of the
#'   symbiont phylogeny
#' @param Host_to_Symbiont_df A data frame with "Host" and "Symbiont" columns
#'   linking symbionts to their source hosts
#'
#' @return A phylo object in which every monophyletic group of tips from a
#'   single host has been reduced to one representative tip. Tips from a host
#'   that do not form a single monophyletic group (polyphyly) are kept, one
#'   per monophyletic group. When a group carries branch lengths the tip
#'   furthest from the group's root is kept; otherwise the first tip is kept.
#'
#' @export
#'
#' @examples
#' host_tree <- ape::rtree(5)
#' symbiont_tree <- ape::rtree(30)
#' hs_df <- data.frame(
#'   Host = rep(host_tree$tip.label, 6),
#'   Symbiont = symbiont_tree$tip.label
#' )
#' collapsed <- collapse_monophyletic(symbiont_tree, hs_df)
#'
collapse_monophyletic <- function(symbiont_subtree, Host_to_Symbiont_df) {
  host_of <- Host_to_Symbiont_df$Host[match(symbiont_subtree$tip.label,
                                            Host_to_Symbiont_df$Symbiont)]
  names(host_of) <- symbiont_subtree$tip.label

  # every internal node's tip set, plus the subtree it roots, widest first so
  # each tip is handled by the largest monophyletic single-host group that
  # contains it and nested groups are skipped
  node_subtrees <- lapply(seq_len(ape::Nnode(symbiont_subtree)), function(i)
    castor::get_subtree_at_node(symbiont_subtree, i)$subtree)
  node_subtrees <- node_subtrees[order(-vapply(node_subtrees,
    function(s) length(s$tip.label), integer(1L)))]

  tips_to_drop <- character(0)
  for (s in node_subtrees) {
    grp <- s$tip.label
    if (length(grp) < 2L) next
    if (any(grp %in% tips_to_drop)) next          # already inside a wider group
    grp_hosts <- host_of[grp]
    if (anyNA(grp_hosts) || length(unique(grp_hosts)) != 1L) next

    # keep exactly one representative, drop the rest, regardless of whether the
    # group carries branch lengths: the tip furthest from the group's root, or
    # the first tip when branch lengths are absent or all equal
    if (is.null(s$edge.length) || all(s$edge.length == 0)) {
      keep <- grp[1L]
    } else {
      root_dist <- castor::get_all_distances_to_root(s)[seq_along(grp)]
      keep <- grp[which.max(root_dist)]
    }
    tips_to_drop <- c(tips_to_drop, setdiff(grp, keep))
  }

  if (length(tips_to_drop) == 0L) {
    return(symbiont_subtree)
  }
  ape::drop.tip(symbiont_subtree, tips_to_drop)
}
