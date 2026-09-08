# Event-based cophylogeny simulation via treeducken (robustness check)
# ------------------------------------------------------------------------
# treeducken simulates host-symbiont cophylogenies under a birth-death process
# with cospeciation, host-switching, within-host duplication, and loss -- i.e.
# congruence EMERGES from the rates rather than being planted. Use this as a
# realism check alongside the constructive simulate_codiv_data(): if the scan
# still localizes signal when data come from a mechanistic process, the result
# is not an artifact of how we built the trees.
#
# treeducken was archived from CRAN; install it from GitHub:
#   remotes::install_github("wadedismukes/treeducken")
#
# This wrapper returns the SAME shape as codiv::simulate_codiv_data():
#   list(host_tree, symbiont_tree, links, truth)
# so it is a drop-in for the analyses in validation_plan.R. Note: treeducken's
# ground truth is the event history (cospeciation vs host-switch), not a clean
# per-clade co-diversifying label -- so `truth` here marks each symbiont by the
# host it was sampled on, and the "co-diversifying" question is addressed at the
# whole-scan level (higher cospeciation rate -> more detected signal) rather
# than per node. Verify argument names against your installed treeducken.

simulate_codiv_treeducken <- function(host_birth = 1.0, host_death = 0.2,
                                      cosp_rate = 2.0, host_switch_rate = 0.5,
                                      symb_birth = 1.0, symb_death = 0.2,
                                      time = 1.0, seed = NULL) {
  if (!requireNamespace("treeducken", quietly = TRUE)) {
    stop("simulate_codiv_treeducken() requires the 'treeducken' package.\n",
         "It is archived on CRAN; install from GitHub:\n",
         "  remotes::install_github('wadedismukes/treeducken')")
  }
  if (!is.null(seed)) set.seed(seed)

  # sim_cophyBD: cophylogenetic birth-death. Argument names vary by version;
  # check ?treeducken::sim_cophyBD and adjust if needed.
  sim <- treeducken::sim_cophyBD(
    hbr = host_birth, hdr = host_death,
    cosp_rate = cosp_rate, host_exp_rate = host_switch_rate,
    sbr = symb_birth, sdr = symb_death,
    numbsim = 1, time_to_sim = time)[[1]]

  host_tree <- sim$host_tree
  symbiont_tree <- sim$symb_tree
  assoc <- sim$association_mat            # hosts x symbionts (0/1)

  # convert the association matrix to a Host/Symbiont data frame
  # (rows = hosts, cols = symbionts, as in treeducken output)
  idx <- which(assoc > 0, arr.ind = TRUE)
  links <- data.frame(
    Host = rownames(assoc)[idx[, 1]],
    Symbiont = colnames(assoc)[idx[, 2]],
    stringsAsFactors = FALSE)

  # a symbiont associated with exactly one host is what codiv expects; if
  # treeducken produced multi-host symbionts, keep the first association
  links <- links[!duplicated(links$Symbiont), ]

  truth <- data.frame(
    Symbiont = links$Symbiont,
    host = links$Host,
    stringsAsFactors = FALSE)

  list(host_tree = host_tree, symbiont_tree = symbiont_tree,
       links = links, truth = truth)
}

# Example workflow (requires treeducken + codiv):
#
#   library(codiv)
#   sim <- simulate_codiv_treeducken(host_switch_rate = 0.2, cosp_rate = 3,
#                                    time = 1, seed = 1)
#   res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
#                min_hosts = 3, min_symbiont_tips = 7, methods = "hommola")
#   null <- codiv_null_scans(sim$host_tree, sim$symbiont_tree, sim$links,
#                            n_permutations = 100)
#   # Robustness sweep: vary host_switch_rate (low = high congruence) and show
#   # the number of detected co-diversifying clades / global permutation p-value
#   # tracks the true amount of cospeciation.
