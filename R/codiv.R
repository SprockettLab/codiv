# Maximum root-to-tip distance within the clade rooted at each internal node.
# Computed in a single post-order pass over the edges, which matches
# castor::get_tree_span()'s max_distance for every node without extracting a
# subtree per node.
# input: tree (phylo with edge lengths)
# output: numeric vector of length Nnode, indexed by internal node 1..Nnode
.node_heights <- function(tree) {
  po <- reorder(tree, "postorder")
  n_tip <- length(tree$tip.label)
  heights <- numeric(n_tip + tree$Nnode)  # tips start at 0
  parent <- po$edge[, 1]
  child <- po$edge[, 2]
  el <- po$edge.length
  for (k in seq_along(parent)) {
    cand <- el[k] + heights[child[k]]
    if (cand > heights[parent[k]]) heights[parent[k]] <- cand
  }
  heights[(n_tip + 1L):(n_tip + tree$Nnode)]
}

# Calculate Colless index: measures tree imbalance.
# For each bifurcating internal node, add |left tips - right tips|; non-binary
# nodes are skipped (as in the reference apTreeshape definition).
# Vectorized via castor::count_tips_per_node instead of per-node traversal.
# Verified against apTreeshape 1.5-0 across 100+ random trees.
.calc_colless <- function(tree) {
  if (tree$Nnode <= 1) return(0)
  n_tips <- length(tree$tip.label)

  # tip count beneath every node: 1 for a tip, clade size for an internal node
  tipcount <- integer(n_tips + tree$Nnode)
  tipcount[seq_len(n_tips)] <- 1L
  tipcount[(n_tips + 1L):(n_tips + tree$Nnode)] <- castor::count_tips_per_node(tree)

  # group child tip counts by parent; a bifurcating node contributes |c1 - c2|
  child_tipcounts <- split(tipcount[tree$edge[, 2]], tree$edge[, 1])
  sum(vapply(child_tipcounts,
             function(cs) if (length(cs) == 2) abs(cs[1] - cs[2]) else 0,
             numeric(1)))
}

# Calculate Sackin index: sum of depths (edges from root) of all leaves.
# Depth is filled in a single pre-order pass (reverse post-order edges).
# Verified against apTreeshape 1.5-0 across 100+ random trees.
.calc_sackin <- function(tree) {
  n_tips <- length(tree$tip.label)
  if (n_tips == 0) return(0)

  depth <- integer(n_tips + tree$Nnode)
  po <- reorder(tree, "postorder")
  # reverse post-order visits parents before children (pre-order)
  ord <- rev(seq_len(nrow(po$edge)))
  parent <- po$edge[ord, 1]
  child <- po$edge[ord, 2]
  for (k in seq_along(parent)) {
    depth[child[k]] <- depth[parent[k]] + 1L
  }
  sum(depth[seq_len(n_tips)])
}

#' Prepare symbiont tree with unique node labels
#'
#' Ensures symbiont tree has unique, non-null node labels for tracking nodes
#' through the codiversification scan.
#'
#' @param tree A phylo object to prepare
#' @param verbose If TRUE, print status messages
#'
#' @return A phylo object with guaranteed unique node labels
#'
#' @keywords internal
.prepare_symbiont_tree <- function(tree, verbose = TRUE) {
  if (is.null(tree$node.label)) {
    # no labels at all: create simple unique ones
    if (verbose) message("Creating unique node labels")
    tree$node.label <- paste0("Node_", seq_len(tree$Nnode))
  } else if (anyDuplicated(tree$node.label) ||
             anyNA(tree$node.label) ||
             any(tree$node.label == "")) {
    # some labels are missing, blank, or duplicated: prepend the node index,
    # which guarantees uniqueness while preserving any original label text
    if (verbose) message("Making node labels unique")
    tree$node.label <- paste0(seq_len(tree$Nnode), "_", tree$node.label)
  }
  tree
}

# Default ceiling on clade size, scaled to the tree.
#
# A fixed number does not travel between trees: 1500 tips is 2% of a 77k-tip
# bacterial tree and three times the whole 473-tip archaeal one. A pure
# fraction fails the other way, capping a small tree at a handful of tips. The
# floor keeps the cap inert on small trees, and the fraction keeps it
# meaningful as a catalog grows.
#
# 2% was chosen against a 76,697-tip tree, where it excludes 5 of 8,593
# scanned nodes whose best correlation was 0.559 -- clades that are
# statistically significant only because thousands of tips give the
# permutation test enormous power at trivial effect sizes.
# inputs: n_tips (integer, tips in the whole symbiont tree)
# output: single integer cap
.default_max_symbiont_tips <- function(n_tips) {
  max(500L, as.integer(ceiling(0.02 * n_tips)))
}


#' Identify nodes to scan and apply filtering criteria
#'
#' Builds list of symbiont tree subtrees and filters by span, tips, and hosts.
#'
#' @param Symbiont_tree A phylo object with unique node labels
#' @param Host_to_Symbiont_df Host-symbiont linkage data frame
#' @param min_hosts Minimum hosts threshold
#' @param min_symbiont_tips Minimum symbiont tips threshold
#' @param span_fraction Span cutoff as fraction of max span
#' @param verbose If TRUE, print filtering statistics
#'
#' @return A list with elements:
#'   - nodes_to_scan: Character vector of node labels to analyze
#'   - subtree_list: Named list of subtrees
#'   - Symbiont_df: Data frame with filtering metadata
#'
#' @keywords internal
.identify_nodes_to_scan <- function(Symbiont_tree, Host_to_Symbiont_df,
                                     min_hosts, min_symbiont_tips,
                                     span_fraction, max_symbiont_tips = NULL,
                                     verbose = TRUE) {
  n_nodes <- Nnode(Symbiont_tree)

  # PERFORMANCE: compute the cheap per-node metrics (span and symbiont tip
  # counts) for all nodes without extracting any subtrees. Subtrees are only
  # built later for the nodes that survive these first two filters.
  span_list <- .node_heights(Symbiont_tree)
  symbiont_tip_count_list <- castor::count_tips_per_node(Symbiont_tree)

  n1 <- n_nodes
  span_cutoff <- span_fraction * max(span_list)

  # A clade whose tips all sit at zero distance from one another carries no
  # phylogenetic signal: every symbiont pairwise distance is 0, so the Hommola
  # correlation is undefined and the PACo ordination has no dimensions to work
  # with. FastTree produces these whenever a group of MAGs has an identical
  # marker alignment. They are dropped rather than scanned and reported as
  # failures.
  span_ok <- span_list > 0 & span_list <= span_cutoff
  n_zero_span <- sum(span_list == 0)

  # span_fraction limits how DEEP a clade is; it does not limit how many tips
  # it holds, and the two are only loosely related (Spearman 0.66 on a
  # 76,697-tip bacterial tree, where clades of 1000+ tips ranged from 5% to
  # 100% of the tree's maximum span). A densely sampled recent radiation is
  # enormous but shallow and passes the span filter untouched. Capping tips is
  # therefore a separate control, not a substitute.
  if (is.null(max_symbiont_tips)) {
    max_symbiont_tips <- .default_max_symbiont_tips(length(Symbiont_tree$tip.label))
  }
  tips_ok <- symbiont_tip_count_list >= min_symbiont_tips &
             symbiont_tip_count_list <= max_symbiont_tips
  n_too_big <- sum(symbiont_tip_count_list > max_symbiont_tips)
  n2 <- sum(span_ok)
  n3 <- sum(span_ok & tips_ok)
  candidate_nodes <- which(span_ok & tips_ok)

  # PERFORMANCE: extract subtrees and compute the (more expensive) distinct-host
  # counts only for candidate nodes that already passed the span/tip filters
  subtree_list <- vector("list", length(candidate_nodes))
  host_counts <- integer(length(candidate_nodes))
  for (j in seq_along(candidate_nodes)) {
    st <- castor::get_subtree_at_node(Symbiont_tree, candidate_nodes[j])$subtree
    subtree_list[[j]] <- st
    host_counts[j] <- length(unique(Host_to_Symbiont_df$Host[
      match(st$tip.label, Host_to_Symbiont_df$Symbiont)]))
  }
  names(subtree_list) <- Symbiont_tree$node.label[candidate_nodes]

  hosts_ok <- host_counts >= min_hosts
  n4 <- sum(hosts_ok)

  # keep only surviving nodes and their already-extracted subtrees
  subtree_list <- subtree_list[hosts_ok]
  Symbiont_df <- data.frame(
    Node_Label = Symbiont_tree$node.label[candidate_nodes][hosts_ok],
    node = candidate_nodes[hosts_ok],
    Subtree_Span = span_list[candidate_nodes][hosts_ok],
    Subtree_Symbiont_Tips = symbiont_tip_count_list[candidate_nodes][hosts_ok],
    Subtree_Host_Tips = host_counts[hosts_ok],
    stringsAsFactors = FALSE)

  if (verbose) {
    message(paste0("Of the ", n1, " internal nodes in the symbiont tree,"))
    message(paste0("   ", n2, " (", scales::label_percent(big.mark = ",")(n2/n1), ") have span > 0 and <= ",
                   scales::label_percent(big.mark = ",")(span_fraction), " of max (",
                   n_zero_span, " dropped for zero span)"))
    message(paste0("   ", n3, " (", scales::label_percent(big.mark = ",")(n3/n1), ") have ",
                   min_symbiont_tips, "-", max_symbiont_tips, " symbiont tips (",
                   n_too_big, " dropped as too large)"))
    message(paste0("   ", n4, " (", scales::label_percent(big.mark = ",")(n4/n1), ") have >= ",
                   min_hosts, " hosts"))
    message(paste0("Scanning ", n4, " nodes for codiversification."))
  }

  # Sort by size (process smaller trees first for early termination benefit)
  Symbiont_df <- Symbiont_df[order(Symbiont_df$Subtree_Symbiont_Tips,
                                    Symbiont_df$Subtree_Host_Tips), ]
  nodes_to_scan <- Symbiont_df$Node_Label
  subtree_list <- subtree_list[nodes_to_scan]

  list(nodes_to_scan = nodes_to_scan, subtree_list = subtree_list,
       Symbiont_df = Symbiont_df)
}

# Scan a single symbiont node for codiversification signal.
# inputs:
#   i: node label to scan
#   subtree_list: named list of symbiont subtrees
#   Host_tree, Host_to_Symbiont_df: full host tree and linkage data
#   methods: character vector of methods to run
#   permutations, seed: permutation settings
#   subtree_features: whether to compute TreeDist metrics
#   focus_hosts, focus_hosts_cols: optional focus-host tracking
# output: named list keyed by output column name (a subset of all_cols), or a
#   list with Node_ID plus a `.error` string when the node fails. Each call is
#   wrapped so one bad node cannot abort a whole scan.
.scan_one_node <- function(i, subtree_list, Host_tree, Host_to_Symbiont_df,
                           methods, permutations, seed, subtree_features,
                           focus_hosts, focus_hosts_cols,
                           uncollapsed_hommola = FALSE) {
  out <- list(Node_ID = i)
  result <- tryCatch({
    # isolate each node as a subtree
    i_symbiont_subtree <- subtree_list[[i]]
    # reduce linkage data and host tree to hosts relevant for this subtree
    i_host_to_symbiont_df <- Host_to_Symbiont_df[
      Host_to_Symbiont_df$Symbiont %in% i_symbiont_subtree$tip.label, ]
    i_host_subtree <- keep.tip(Host_tree, unique(i_host_to_symbiont_df$Host))

    # binary host-by-symbiont association matrix
    i_host_to_symbiont_mat <- host_symbiont_links(i_host_to_symbiont_df)

    # patristic distance matrices
    # PERFORMANCE: the full tip-by-tip matrix is only needed for the
    # uncollapsed Hommola test. Everything else works from the collapsed tree,
    # so when that test is off the distances are computed on the collapsed
    # tree directly. On a clade of a few thousand near-identical MAGs the
    # collapsed tree can be an order of magnitude smaller.
    i_symbiont_subtree_dist <- if (uncollapsed_hommola) {
      as.matrix(adephylo::distTips(i_symbiont_subtree, tips = "all",
                                   method = "patristic", useC = TRUE))
    } else NULL
    i_host_subtree_dist <- adephylo::distTips(i_host_subtree, tips = "all",
                                              method = "patristic", useC = TRUE)

    # collapse symbionts from the same host that form a monophyletic group
    i_collapsed_symbiont_subtree <- collapse_monophyletic(i_symbiont_subtree,
                                                          i_host_to_symbiont_df)
    # PERFORMANCE: patristic distances between retained tips are unchanged when
    # other tips are dropped, so the collapsed distance matrix is just a subset
    # of the full one (avoids a second distTips call per node)
    collapsed_tips <- i_collapsed_symbiont_subtree$tip.label
    i_collapsed_symbiont_subtree_dist <- if (is.null(i_symbiont_subtree_dist)) {
      as.matrix(adephylo::distTips(i_collapsed_symbiont_subtree, tips = "all",
                                   method = "patristic", useC = TRUE))[
        collapsed_tips, collapsed_tips, drop = FALSE]
    } else {
      i_symbiont_subtree_dist[collapsed_tips, collapsed_tips, drop = FALSE]
    }
    i_collapsed_host_to_symbiont_df <- Host_to_Symbiont_df[
      Host_to_Symbiont_df$Symbiont %in% i_collapsed_symbiont_subtree$tip.label, ]
    i_collapsed_host_to_symbiont_mat <- host_symbiont_links(
      i_collapsed_host_to_symbiont_df)

    if ("hommola" %in% methods) {
      # DECISION: the uncollapsed test is off by default. It keeps every MAG
      # from a host that formed a monophyletic group, so a host sampled deeply
      # contributes many near-identical tips and the permutation test treats
      # them as independent. On a 76,697-tip scan it called 2,017 nodes
      # significant that the collapsed test did not, against 57 the other way,
      # while the two correlations agreed at Spearman 0.83. That asymmetry is
      # pseudoreplication, not sensitivity. It is also the expensive test: its
      # permuted distances scale with the square of the uncollapsed clade
      # size, which is what made the largest clades take hours each.
      if (uncollapsed_hommola) {
        hommola_results <- hommola_wf(i_host_subtree_dist, i_symbiont_subtree_dist,
                                      i_host_to_symbiont_df, permutations, seed)
        out[["Uncollapsed_Hommola_r"]] <- hommola_results$Uncollapsed_Hommola_r
        out[["Uncollapsed_Hommola_pvalue"]] <- hommola_results$Uncollapsed_Hommola_pvalue
      }

      collapsed_hommola <- hommola_wf(i_host_subtree_dist,
                                      i_collapsed_symbiont_subtree_dist,
                                      i_host_to_symbiont_df, permutations, seed)
      out[["Hommola_r"]] <- collapsed_hommola$Uncollapsed_Hommola_r
      out[["Hommola_pvalue"]] <- collapsed_hommola$Uncollapsed_Hommola_pvalue
    }

    if ("paco" %in% methods) {
      paco_results <- paco_wf(i_host_subtree, i_collapsed_symbiont_subtree,
                              i_collapsed_host_to_symbiont_mat, permutations, seed)
      out[["PACo_ss"]] <- paco_results$gof$ss
      out[["PACo_pvalue"]] <- paco_results$gof$p
    }

    if ("parafit" %in% methods) {
      # ape::parafit matches host.D, para.D and the association matrix by
      # position, not by name, so the association matrix must be reordered to
      # the distance-matrix row orders or the statistic is computed on
      # mismatched host-symbiont pairs
      parafit_host_D <- as.matrix(i_host_subtree_dist)
      parafit_para_D <- i_collapsed_symbiont_subtree_dist
      parafit_hp <- i_collapsed_host_to_symbiont_mat[
        rownames(parafit_host_D), rownames(parafit_para_D), drop = FALSE]
      parafit_results <- ape::parafit(parafit_host_D,
                                      parafit_para_D,
                                      parafit_hp,
                                      nperm = permutations, test.links = FALSE,
                                      seed = seed, correction = "cailliez",
                                      silent = TRUE)
      out[["ParaFitGlobal"]] <- parafit_results$ParaFitGlobal
      out[["ParaFit_pvalue"]] <- parafit_results$p.global
    }

    if ("topology" %in% methods) {
      topo <- topology_wf(i_host_subtree, i_symbiont_subtree,
                          i_host_to_symbiont_df, permutations, seed)
      out[["Topology_congruence"]] <- topo$stat
      out[["Topology_pvalue"]] <- topo$pvalue
      out[["Topology_hosts_used"]] <- topo$hosts_used
      out[["Topology_fraction_kept"]] <- topo$fraction_kept
    }

    # tree shape and size features
    out[["Symbiont_Tree"]] <- write.tree(i_symbiont_subtree)
    out[["N_Symbionts"]] <- length(i_symbiont_subtree$tip.label)
    out[["Symbiont_Colless"]] <- .calc_colless(i_symbiont_subtree)
    out[["Symbiont_Sackin"]] <- .calc_sackin(i_symbiont_subtree)

    out[["Host_Tree"]] <- write.tree(i_host_subtree)
    out[["N_Hosts"]] <- length(i_host_subtree$tip.label)
    out[["Host_Colless"]] <- .calc_colless(i_host_subtree)
    out[["Host_Sackin"]] <- .calc_sackin(i_host_subtree)

    if (length(focus_hosts_cols) > 0) {
      present <- vapply(focus_hosts,
                        function(h) h %in% i_host_to_symbiont_df$Host,
                        logical(1))
      out[focus_hosts_cols] <- present
    }

    if (subtree_features) {
      # rename collapsed symbiont tips with their host so the two trees share a
      # tip-label set, then compute TreeDist metrics (see github.com/ms609/TreeDist)
      f_symbiont_subtree <- reorder(i_collapsed_symbiont_subtree,
                                    order = "cladewise")
      f_symbiont_subtree$tip.label <- i_collapsed_host_to_symbiont_df$Host[
        match(f_symbiont_subtree$tip.label,
              i_collapsed_host_to_symbiont_df$Symbiont)]

      out[["TreeDistance"]] <- TreeDistance(i_host_subtree, f_symbiont_subtree)
      out[["SharedPhylogeneticInfo"]] <- SharedPhylogeneticInfo(i_host_subtree, f_symbiont_subtree)
      out[["DifferentPhylogeneticInfo"]] <- DifferentPhylogeneticInfo(i_host_subtree, f_symbiont_subtree)
      out[["NyeSimilarity"]] <- NyeSimilarity(i_host_subtree, f_symbiont_subtree)
      out[["JaccardRobinsonFoulds"]] <- JaccardRobinsonFoulds(i_host_subtree, f_symbiont_subtree)
      out[["MatchingSplitDistance"]] <- MatchingSplitDistance(i_host_subtree, f_symbiont_subtree)
      out[["MatchingSplitInfoDistance"]] <- MatchingSplitInfoDistance(i_host_subtree, f_symbiont_subtree)
      out[["MutualClusteringInfo"]] <- MutualClusteringInfo(i_host_subtree, f_symbiont_subtree)
    }
    out
  }, error = function(e) {
    warning("Node '", i, "' failed and was skipped: ", conditionMessage(e),
            call. = FALSE)
    list(Node_ID = i, .error = conditionMessage(e))
  })
  result
}

#' Perform codiversification scan between host and symbiont trees
#'
#' Scans a symbiont phylogenetic tree for evidence of codiversification with
#' a host tree. For each internal node of the symbiont tree meeting size
#' thresholds, calculates correlation coefficients and p-values using one or
#' more methods (Hommola, PACo, ParaFit) to test if the symbiont topology
#' matches the host topology at that node.
#'
#' @param Host_tree A phylo object; binary phylogenetic tree of hosts
#' @param Symbiont_tree A phylo object; binary phylogenetic tree of symbionts
#' @param Host_to_Symbiont_df A data frame with columns "Host" and "Symbiont"
#'   linking each symbiont tip to its source host
#' @param min_hosts Minimum number of hosts required at a node for inclusion
#'   in the scan; default 3
#' @param min_symbiont_tips Minimum number of symbiont tips in a subtree for
#'   inclusion; default 7 (higher values recommended for statistical power)
#' @param max_symbiont_tips Maximum number of symbiont tips in a subtree. If
#'   NULL (default), `max(500, 2% of the tree's tips)`. This is not the same
#'   control as `span_fraction`: that limits how deep a clade is, this limits
#'   how many genomes it holds, and the two are only loosely related. A
#'   densely sampled recent radiation is enormous but shallow, so it passes
#'   any reasonable span cutoff. Very large clades are also where the
#'   permutation test reaches significance at trivial effect sizes simply
#'   because n is large.
#' @param uncollapsed_hommola If TRUE, also run the Hommola test on the
#'   uncollapsed subtree and report `Uncollapsed_Hommola_r`/`Uncollapsed_Hommola_pvalue`. Default
#'   FALSE. The uncollapsed test counts every MAG from a host that formed a
#'   monophyletic group as an independent observation, so a deeply sampled
#'   host is pseudoreplicated; PACo and ParaFit already use the collapsed
#'   tree. It is also by far the most expensive statistic on large clades.
#' @param span_fraction Filters nodes to only those with phylogenetic span
#'   <= this fraction of total tree span (0-1); default 0.1 focuses on distal nodes
#' @param permutations Number of randomizations for significance testing;
#'   default 99; higher values (999+) recommended for publication
#' @param seed Random seed for reproducible permutations; default 8675309
#' @param verbose If TRUE (default), prints progress and filtering statistics
#' @param Save_fp File path where results are saved after each node (enables
#'   resumable runs via `continue`); if NULL (default), a unique temporary file
#'   is used, so separate calls never share a checkpoint. Supply an explicit
#'   path to enable resuming a specific run
#' @param methods Character vector of methods to use: any combination of
#'   "hommola", "paco", "parafit" (default uses these three, all distance-based),
#'   and the opt-in "topology" (a branching-order congruence test that ignores
#'   branch lengths; each host is reduced to its largest single-host clade and
#'   hosts with no coherent position are excluded, reported via the
#'   Topology_hosts_used and Topology_fraction_kept columns)
#' @param focus_hosts Optional character vector of host names. For each name a
#'   logical `<host>_PRESENT` column is added to the results, TRUE at nodes whose
#'   subtree contains that host. Useful for sorting/filtering the scan when you
#'   care about co-diversifying symbionts in particular hosts. Names not found in
#'   `Host_to_Symbiont_df` trigger a warning and are FALSE everywhere.
#' @param subtree_features If TRUE (default), calculates multiple tree distance
#'   metrics (TreeDistance, SharedPhylogeneticInfo, etc.) for each node
#' @param cores Number of cores for parallelizing the node scan. If NULL
#'   (default), uses all detected cores minus one. Forking is unavailable on
#'   Windows, where the scan always runs on a single core.
#' @param batch_size Number of nodes to scan per checkpoint batch. Results are
#'   written to `Save_fp` after each batch, so smaller values checkpoint more
#'   often (more restart safety, slightly more disk I/O). If NULL (default), a
#'   size of `max(200, cores * 50)` is used. Mutually exclusive with `n_batches`.
#' @param n_batches Alternative to `batch_size`: split the scanned nodes into
#'   this many roughly equal checkpoint batches. Mutually exclusive with
#'   `batch_size`.
#' @param continue If TRUE (default), resumes an interrupted scan by skipping
#'   completed nodes; requires Save_fp file to exist
#' @return A `codiv` object: a data frame (with one row per symbiont tree node
#'   scanned) carrying the run parameters as a `codiv_params` attribute, with
#'   `print` and `summary` methods. It behaves like a data frame for subsetting
#'   and column access. Columns include:
#'   - Node_ID: Node identifier
#'   - N_Symbionts, N_Hosts: Tip counts in each subtree
#'   - Symbiont_Colless, Symbiont_Sackin: Tree shape metrics
#'   - Hommola_r, Hommola_pvalue: Hommola Pearson correlation on the collapsed
#'     subtree and its permutation p-value (if the hommola method is used)
#'   - Uncollapsed_Hommola_r, Uncollapsed_Hommola_pvalue: the same test on the
#'     uncollapsed subtree (only when `uncollapsed_hommola = TRUE`)
#'   - PACo_ss, PACo_pvalue: PACo results (if method used)
#'   - ParaFitGlobal, ParaFit_pvalue: ParaFit results (if method used)
#'   - Topology_congruence, Topology_pvalue, Topology_hosts_used,
#'     Topology_fraction_kept: topology-congruence results (if method used)
#'   - Tree distance metrics if subtree_features=TRUE
#'
#'   P-values are the per-node permutation p-values. They are not corrected for
#'   multiple testing across nodes: nested clades are non-independent, so a
#'   per-node correction (e.g. Benjamini-Hochberg) would be invalid. Use
#'   [codiv_null_scans()] for a scan-level false-discovery assessment.
#'
#' @importFrom ape write.tree extract.clade node.depth.edgelength keep.tip Nnode
#' @importFrom stats cor reorder
#' @importFrom utils read.table write.table
#' @importFrom pbmcapply pbmclapply
#' @importFrom TreeDist TreeDistance SharedPhylogeneticInfo DifferentPhylogeneticInfo NyeSimilarity JaccardRobinsonFoulds MatchingSplitDistance MatchingSplitInfoDistance MutualClusteringInfo
#' @export
#'
#' @examples
#' \donttest{
#' # simulate a small host/symbiont dataset with known co-diversification
#' sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
#'
#' results <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
#'                  methods = "hommola", permutations = 99, verbose = FALSE)
#' results            # compact overview
#' head(as.data.frame(results))
#'
#' # Bundled example trees (real data) load with system.file():
#' # host <- ape::read.tree(system.file("extdata", "Host_Tree.nwk",
#' #                                    package = "codiv"))
#' }
#'
codiv <- function(Host_tree, Symbiont_tree, Host_to_Symbiont_df,
                  min_hosts = 3, min_symbiont_tips = 7,
                  max_symbiont_tips = NULL, span_fraction = 0.1,
                  uncollapsed_hommola = FALSE,
                  permutations = 99, seed = 8675309,
                  verbose = TRUE, Save_fp = NULL, methods = c("hommola", "paco", "parafit"),
                  focus_hosts = NULL, subtree_features = TRUE, cores = NULL,
                  batch_size = NULL, n_batches = NULL,
                  continue = TRUE) {

  # default: use all available cores but one, leaving the machine responsive
  if (is.null(cores)) {
    nc <- parallel::detectCores()
    cores <- if (is.na(nc)) 1L else max(1L, nc - 1L)
  }

  # preserve the caller's global RNG state: codiv seeds each node internally
  # but must not leave the user's random stream changed
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    saved_rng <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", saved_rng, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }

  # validate requested methods
  valid_methods <- c("hommola", "paco", "parafit", "topology")
  bad_methods <- setdiff(methods, valid_methods)
  if (length(bad_methods) > 0) {
    stop("Unknown method(s): ", paste(bad_methods, collapse = ", "),
         ". Choose from: ", paste(valid_methods, collapse = ", "))
  }

  # cores must be a single positive whole number
  if (!is.numeric(cores) || length(cores) != 1 || is.na(cores) ||
      cores < 1 || cores != round(cores)) {
    stop("`cores` must be a single positive integer")
  }

  # batch_size and n_batches control checkpoint frequency; they are mutually
  # exclusive and each, if given, must be a single positive whole number
  .check_pos_int <- function(x, nm) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 1 ||
        x != round(x)) {
      stop("`", nm, "` must be a single positive integer")
    }
  }
  if (!is.null(batch_size) && !is.null(n_batches)) {
    stop("Supply only one of `batch_size` or `n_batches`, not both")
  }
  if (!is.null(batch_size)) .check_pos_int(batch_size, "batch_size")
  if (!is.null(n_batches)) .check_pos_int(n_batches, "n_batches")

  # focus_hosts, if given, must be host names present in the linking table
  if (!is.null(focus_hosts)) {
    if (!is.character(focus_hosts) || length(focus_hosts) < 1) {
      stop("`focus_hosts` must be a character vector of host names")
    }
    focus_hosts <- unique(focus_hosts)
    missing_focus <- setdiff(focus_hosts, as.character(Host_to_Symbiont_df$Host))
    if (length(missing_focus) > 0) {
      warning("These focus_hosts are not in Host_to_Symbiont_df and will be ",
              "FALSE at every node: ", paste(missing_focus, collapse = ", "))
    }
  }

  check_inputs(Host_tree, Symbiont_tree, Host_to_Symbiont_df, min_hosts,
               min_symbiont_tips, span_fraction, permutations)

  # coerce link columns to character so factor levels cannot corrupt matching
  Host_to_Symbiont_df$Host <- as.character(Host_to_Symbiont_df$Host)
  Host_to_Symbiont_df$Symbiont <- as.character(Host_to_Symbiont_df$Symbiont)

  # every symbiont tip must map to a host; drop tips with no association so the
  # methods (which need matching symbiont and association sets) do not fail
  orphan_symbionts <- setdiff(Symbiont_tree$tip.label,
                              Host_to_Symbiont_df$Symbiont)
  if (length(orphan_symbionts) > 0) {
    if (verbose) {
      message("Dropping ", length(orphan_symbionts),
              " symbiont tip(s) with no host association")
    }
    Symbiont_tree <- ape::drop.tip(Symbiont_tree, orphan_symbionts)
  }

  # Set Save_fp default and create directory if needed. Use a unique temp file
  # so two codiv() calls in one session cannot cross-contaminate through a
  # shared checkpoint when `continue = TRUE`.
  if (is.null(Save_fp)) {
    Save_fp <- tempfile(pattern = "codiv_results_", fileext = ".tsv")
  } else {
    save_dir <- dirname(Save_fp)
    if (save_dir != "." && !dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }

  # Prepare symbiont tree and identify nodes to scan
  Symbiont_tree <- .prepare_symbiont_tree(Symbiont_tree, verbose = verbose)
  if (is.null(max_symbiont_tips)) {
    max_symbiont_tips <- .default_max_symbiont_tips(length(Symbiont_tree$tip.label))
  }
  scan_info <- .identify_nodes_to_scan(Symbiont_tree, Host_to_Symbiont_df,
                                        min_hosts, min_symbiont_tips,
                                        span_fraction, max_symbiont_tips,
                                        verbose = verbose)
  nodes_to_scan <- scan_info$nodes_to_scan
  subtree_list <- scan_info$subtree_list

  # Per-node seeds: derive a distinct, reproducible seed for each node from the
  # base seed and the node's (stable) tree index. This makes the permutation
  # nulls independent across nodes yet fully deterministic and independent of
  # core count and scan order. A base seed of NA keeps permutations random.
  if (is.na(seed)) {
    node_seeds <- stats::setNames(rep(NA_real_, nrow(scan_info$Symbiont_df)),
                                  scan_info$Symbiont_df$Node_Label)
  } else {
    node_seeds <- stats::setNames(seed + scan_info$Symbiont_df$node,
                                  scan_info$Symbiont_df$Node_Label)
  }

  # Build output column structure
  basic_cols <- c("Node_ID", "Symbiont_Tree", "N_Symbionts", "Symbiont_Colless",
                  "Symbiont_Sackin", "Host_Tree", "N_Hosts", "Host_Colless", "Host_Sackin")
  method_cols <- c()
  if ("hommola" %in% methods) {
    if (uncollapsed_hommola) {
      method_cols <- c(method_cols, "Uncollapsed_Hommola_r", "Uncollapsed_Hommola_pvalue")
    }
    method_cols <- c(method_cols, "Hommola_r",
                     "Hommola_pvalue")
  }
  if ("paco" %in% methods) {
    method_cols <- c(method_cols, "PACo_ss", "PACo_pvalue")
  }
  if ("parafit" %in% methods) {
    method_cols <- c(method_cols, "ParaFitGlobal", "ParaFit_pvalue")
  }
  if ("topology" %in% methods) {
    method_cols <- c(method_cols, "Topology_congruence", "Topology_pvalue",
                     "Topology_hosts_used", "Topology_fraction_kept")
  }
  tree_metric_cols <- c()
  if (subtree_features) {
    tree_metric_cols <- c("TreeDistance", "SharedPhylogeneticInfo",
                          "DifferentPhylogeneticInfo", "NyeSimilarity",
                          "JaccardRobinsonFoulds", "MatchingSplitDistance",
                          "MatchingSplitInfoDistance", "MutualClusteringInfo")
  }
  focus_hosts_cols <- if (!is.null(focus_hosts)) paste0(focus_hosts, "_PRESENT") else c()
  all_cols <- c(basic_cols, method_cols, tree_metric_cols, focus_hosts_cols)

  # Initialize results dataframe
  Results_df <- data.frame(matrix(nrow = length(nodes_to_scan), ncol = length(all_cols)))
  colnames(Results_df) <- all_cols
  rownames(Results_df) <- nodes_to_scan
  Results_df$Node_ID <- nodes_to_scan

  # Handle resuming from a previous run (checkpoint file). A node counts as
  # complete only when every completion column is filled, so partially written
  # rows are re-scanned rather than trusted.
  if (continue && file.exists(Save_fp)) {
    prev_results <- tryCatch(
      suppressWarnings(read.table(Save_fp, sep = "\t", header = TRUE,
                                  fill = TRUE, stringsAsFactors = FALSE)),
      error = function(e) NULL)

    if (is.null(prev_results) || !("Node_ID" %in% colnames(prev_results)) ||
        !all(all_cols %in% colnames(prev_results))) {
      if (verbose) {
        message("Checkpoint file present but unreadable or incompatible; ",
                "starting a fresh scan.")
      }
    } else {
      # columns whose presence marks a node as scanned
      completion_cols <- if (length(method_cols) > 0) method_cols else "N_Symbionts"
      completion_cols <- intersect(completion_cols, colnames(prev_results))
      complete_mask <- !is.na(prev_results$Node_ID) &
        Reduce(`&`, lapply(completion_cols,
                           function(cc) !is.na(prev_results[[cc]])))
      completed <- prev_results$Node_ID[complete_mask]

      # carry completed rows forward into the results table
      carry <- complete_mask & (prev_results$Node_ID %in% Results_df$Node_ID)
      if (any(carry)) {
        src <- prev_results[carry, all_cols, drop = FALSE]
        Results_df[match(src$Node_ID, Results_df$Node_ID), all_cols] <- src
      }
      nodes_to_scan <- nodes_to_scan[!(nodes_to_scan %in% completed)]
      if (verbose) {
        message("Resuming from checkpoint: ", length(completed),
                " nodes already complete, ", length(nodes_to_scan),
                " remaining.")
      }
    }
  }

  write.table(Results_df, Save_fp, sep = "\t", quote = FALSE, row.names = FALSE)
  # Fill one scanned node's results (a named list from .scan_one_node) into the
  # results table, column by column, matching on Node_ID.
  fill_row <- function(df, row_list) {
    node <- row_list$Node_ID
    cols <- setdiff(intersect(names(row_list), all_cols), "Node_ID")
    for (cc in cols) {
      df[node, cc] <- row_list[[cc]]
    }
    df
  }

  # Determine the parallel backend. Forking is unavailable on Windows, where
  # the scan falls back to a single core.
  eff_cores <- if (.Platform$OS.type == "windows") 1L else max(1L, as.integer(cores))
  if (isTRUE(cores > 1) && .Platform$OS.type == "windows" && verbose) {
    message("Parallel scanning is not supported on Windows; running serially.")
  }
  # respect the 2-core limit enforced during R CMD check / on CRAN
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
    eff_cores <- min(eff_cores, 2L)
  }
  if (verbose && length(nodes_to_scan) > 0) {
    message("Scanning ", length(nodes_to_scan), " nodes",
            if (eff_cores > 1L) paste0(" across ", eff_cores, " cores") else "",
            " ...")
  }

  scan_fun <- function(i) {
    .scan_one_node(i, subtree_list, Host_tree, Host_to_Symbiont_df, methods,
                   permutations, node_seeds[[i]], subtree_features, focus_hosts,
                   focus_hosts_cols, uncollapsed_hommola)
  }

  # Process nodes in chunks so results are checkpointed to disk periodically.
  # Each chunk runs in parallel via pbmcapply, which draws its own live progress
  # bar; the message before each bar makes clear these are sequential
  # checkpointed batches of nodes, not permutations.
  if (length(nodes_to_scan) > 0) {
    # checkpoint granularity: n_batches splits the nodes into that many batches
    # (a fractional divisor makes the count exact); batch_size fixes the nodes
    # per batch; otherwise a sensible default.
    n_nodes <- length(nodes_to_scan)
    divisor <- if (!is.null(n_batches)) {
      n_nodes / n_batches
    } else if (!is.null(batch_size)) {
      as.numeric(batch_size)
    } else {
      max(200, eff_cores * 50)
    }
    chunks <- split(nodes_to_scan, ceiling(seq_along(nodes_to_scan) / divisor))
    n_chunks <- length(chunks)
    done <- 0L
    for (k in seq_len(n_chunks)) {
      ch <- chunks[[k]]
      if (verbose && n_chunks > 1L) {
        message(sprintf("Batch %d/%d: nodes %d-%d of %d (checkpointed after each)",
                        k, n_chunks, done + 1L, done + length(ch),
                        length(nodes_to_scan)))
      }
      done <- done + length(ch)
      rows <- pbmcapply::pbmclapply(ch, scan_fun, mc.cores = eff_cores)
      for (r in rows) {
        Results_df <- fill_row(Results_df, r)
      }
      # checkpoint after each chunk
      write.table(Results_df, Save_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  # P-values are left as per-node permutation p-values. They are not corrected
  # across nodes: nested clades are non-independent, so a per-node multiple-
  # testing correction would be invalid. Use codiv_null_scans() for a scan-level
  # false-discovery assessment.
  write.table(Results_df, Save_fp, sep = "\t", quote = FALSE, row.names = FALSE)

  # return a self-describing S3 object (still a data frame underneath)
  attr(Results_df, "codiv_params") <- list(
    methods = intersect(valid_methods, methods),
    permutations = permutations,
    min_hosts = min_hosts,
    min_symbiont_tips = min_symbiont_tips,
    span_fraction = span_fraction,
    seed = seed,
    subtree_features = subtree_features
  )
  class(Results_df) <- c("codiv", "data.frame")
  Results_df
}
