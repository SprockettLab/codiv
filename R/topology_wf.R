# Choose one representative symbiont tip per host for the topology test.
# For each host, its symbionts may form several maximal single-host clades. We
# keep the tip from the LARGEST such clade (tie-broken by smallest total clade
# branch length). A host whose symbionts are entirely scattered (every clade is
# a single tip) has no defined position and is excluded.
# inputs:
#   tree: symbiont subtree (phylo)
#   host_of: named character vector, tip label -> host, covering all tips
# output: list(keep, hosts_used, fraction_kept)
#   keep = tip labels to retain (one per included host)
.host_representatives <- function(tree, host_of) {
  n_tip <- length(tree$tip.label)
  n_all <- n_tip + tree$Nnode
  if (is.null(tree$edge.length)) {
    tree$edge.length <- rep(1, nrow(tree$edge))
  }
  edge_len <- numeric(n_all)
  edge_len[tree$edge[, 2]] <- tree$edge.length
  children <- split(tree$edge[, 2], tree$edge[, 1])

  pure <- logical(n_all)
  pure_host <- rep(NA_character_, n_all)
  n_tips <- integer(n_all)
  brlen <- numeric(n_all)

  # tips
  pure[seq_len(n_tip)] <- TRUE
  pure_host[seq_len(n_tip)] <- host_of[tree$tip.label]
  n_tips[seq_len(n_tip)] <- 1L

  # internal nodes, deepest first (postorder by node depth)
  depths <- ape::node.depth.edgelength(tree)
  internal <- (n_tip + 1L):n_all
  for (nd in internal[order(depths[internal], decreasing = TRUE)]) {
    ch <- children[[as.character(nd)]]
    n_tips[nd] <- sum(n_tips[ch])
    brlen[nd] <- sum(brlen[ch] + edge_len[ch])
    ch_hosts <- pure_host[ch]
    if (all(pure[ch]) && length(unique(ch_hosts)) == 1L) {
      pure[nd] <- TRUE
      pure_host[nd] <- ch_hosts[1]
    }
  }

  # each tip's maximal single-host clade root = highest pure ancestor
  parent_of <- integer(n_all)
  parent_of[tree$edge[, 2]] <- tree$edge[, 1]
  clade_root <- integer(n_tip)
  for (tp in seq_len(n_tip)) {
    cur <- tp
    repeat {
      p <- parent_of[cur]
      if (p == 0L || !pure[p]) break
      cur <- p
    }
    clade_root[tp] <- cur
  }

  # one entry per clade: host, size, branch length, a representative tip
  clade_ids <- unique(clade_root)
  clades <- lapply(clade_ids, function(cr) {
    tips_in <- which(clade_root == cr)
    list(host = pure_host[cr], size = n_tips[cr], brlen = brlen[cr],
         rep = tree$tip.label[tips_in[1]])
  })
  clade_host <- vapply(clades, function(x) x$host, character(1))

  hosts <- unique(host_of[tree$tip.label])
  hosts_total <- length(hosts)
  keep <- character(0)
  kept_tip_count <- 0L
  for (h in hosts) {
    idx <- which(clade_host == h)
    sizes <- vapply(clades[idx], function(x) x$size, numeric(1))
    if (length(idx) >= 2 && max(sizes) == 1) {
      next  # fully scattered host: no defined position, exclude
    }
    best <- idx[order(-sizes, vapply(clades[idx], function(x) x$brlen, numeric(1)))][1]
    keep <- c(keep, clades[[best]]$rep)
    kept_tip_count <- kept_tip_count + clades[[best]]$size
  }

  list(keep = keep,
       hosts_used = length(keep),
       fraction_kept = if (n_tip > 0) kept_tip_count / n_tip else NA_real_)
}

#' Topology-congruence test with permutation significance
#'
#' Tests whether a symbiont tree's branching order matches a host tree's,
#' ignoring branch lengths. Each host is reduced to one representative tip (its
#' largest single-host clade, tie-broken by smallest total clade branch length);
#' hosts whose symbionts are entirely scattered have no defined position and are
#' excluded. Both trees are then pruned to the shared hosts and compared with
#' normalized mutual clustering information (1 = identical branching order). The
#' null distribution is built by shuffling the host labels.
#'
#' This complements the distance-based methods (Hommola, PACo, ParaFit): it is
#' insensitive to branch-length differences, so it can detect congruent topology
#' even when evolutionary rates differ between the trees.
#'
#' @param host_tree A phylo object; the host (sub)tree
#' @param symbiont_subtree A phylo object; the symbiont (sub)tree
#' @param host_to_symbiont_df A data frame with "Host" and "Symbiont" columns
#' @param permutations Number of permutations for the null; at least 99
#'   recommended
#' @param seed Random seed for reproducibility; if NA, no seed is set
#'
#' @return A list with:
#'   - `stat`: normalized mutual clustering information (0-1; 1 = identical
#'     branching order), or NA when fewer than 3 hosts have a defined position
#'   - `pvalue`: one-tailed permutation p-value
#'   - `hosts_used`: number of hosts with a defined position
#'   - `fraction_kept`: fraction of symbiont tips represented after reducing each
#'     host to its dominant clade (low values flag heavy pruning)
#'
#' @export
#'
#' @examples
#' \donttest{
#' h <- ape::rcoal(8)
#' s <- h; s$tip.label <- paste0("s", seq_along(s$tip.label))
#' hs <- data.frame(Host = h$tip.label, Symbiont = s$tip.label)
#' topology_wf(h, s, hs, permutations = 99, seed = 1)
#' }
#'
topology_wf <- function(host_tree, symbiont_subtree, host_to_symbiont_df,
                        permutations, seed) {
  host_of <- stats::setNames(host_to_symbiont_df$Host,
                             host_to_symbiont_df$Symbiont)[symbiont_subtree$tip.label]

  reps <- .host_representatives(symbiont_subtree, host_of)
  if (reps$hosts_used < 3) {
    return(list(stat = NA_real_, pvalue = NA_real_,
                hosts_used = reps$hosts_used, fraction_kept = reps$fraction_kept))
  }

  # symbiont tree reduced to one tip per host, relabeled with the host
  g <- ape::keep.tip(symbiont_subtree, reps$keep)
  g$tip.label <- host_of[g$tip.label]
  shared <- intersect(host_tree$tip.label, g$tip.label)

  # reorder so both trees carry the "order" attribute TreeDist expects
  h <- stats::reorder(ape::keep.tip(host_tree, shared), "cladewise")
  g <- stats::reorder(ape::keep.tip(g, shared), "cladewise")

  obs <- as.numeric(TreeDist::MutualClusteringInfo(h, g, normalize = TRUE))

  # permutation null: shuffle host labels on the symbiont tree; seed locally and
  # restore the global RNG so there are no side effects
  if (!is.na(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
    }
    set.seed(seed)
  }
  labels <- g$tip.label
  perm_vals <- vapply(seq_len(permutations), function(k) {
    gp <- g
    gp$tip.label <- sample(labels)
    as.numeric(TreeDist::MutualClusteringInfo(h, gp, normalize = TRUE))
  }, numeric(1))

  list(stat = obs,
       pvalue = (sum(perm_vals >= obs) + 1) / (permutations + 1),
       hosts_used = reps$hosts_used,
       fraction_kept = reps$fraction_kept)
}
