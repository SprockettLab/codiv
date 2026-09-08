# Simulation-based validation plan for the codiv methods paper
# ------------------------------------------------------------
# Uses the package's simulate_codiv_data() to generate host/symbiont data with
# KNOWN ground truth, then runs four analyses that address the gaps not covered
# by the published empirical validation (Sanders et al. 2023; Sprockett 2025):
#   1. Type I error - null data: are per-node p-values uniform?
#   2. Power        - detection rate vs. congruence
#   3. Localization - node-level AUC (does the scan flag the right clades?)
#   4. FDR calibration - naive per-node counting vs. the second-order
#                        permutation test (codiv_null_scans) under the null
# Plus the headline congruence sweep: node-level AUC vs. congruence, for the
# uncollapsed vs. collapsed statistic and with/without within-host copies.
#
# Run: Rscript simulations/validation_plan.R   (scale up reps/permutations for
# the paper). A realistic event-based generator (treeducken) is provided in
# simulations/simulate_treeducken.R as a robustness check.

suppressMessages({ library(ape); library(codiv) })

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# Annotate EVERY scanned node with the known ground truth, so you can compare
# what codiv calculated (Hommola_r, p-value, ...) against the truth for each
# node. Adds: majority_clade, clade_purity (fraction of the node's tips in that
# clade), truth_codiversifying, truth_congruence, truth_n_per_host.
annotate_nodes <- function(res_df, sim) {
  clade_of <- stats::setNames(sim$truth$clade, sim$truth$Symbiont)
  clade_summary <- sim$truth[!duplicated(sim$truth$clade),
                             c("clade", "codiversifying", "congruence",
                               "n_per_host")]
  ann <- lapply(res_df$Symbiont_Tree, function(nw) {
    tips <- tryCatch(read.tree(text = nw)$tip.label, error = function(e) NA)
    if (all(is.na(tips))) {
      return(data.frame(majority_clade = NA, clade_purity = NA,
                        truth_codiversifying = NA, truth_congruence = NA,
                        truth_n_per_host = NA))
    }
    tab <- sort(table(clade_of[tips]), decreasing = TRUE)
    cl <- names(tab)[1]
    row <- clade_summary[clade_summary$clade == cl, ]
    data.frame(majority_clade = cl, clade_purity = tab[[1]] / length(tips),
               truth_codiversifying = row$codiversifying,
               truth_congruence = row$congruence,
               truth_n_per_host = row$n_per_host)
  })
  cbind(res_df, do.call(rbind, ann))
}

# convenience: the ground-truth co-diversifying label per scanned node
node_truth <- function(res_df, sim) annotate_nodes(res_df, sim)$truth_codiversifying

# AUC via the Mann-Whitney relationship (no extra dependency)
auc <- function(score, truth) {
  pos <- score[truth]; neg <- score[!truth]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
}

run_scan <- function(sim, span_fraction = 0.6, permutations = 199,
                     methods = "hommola") {
  as.data.frame(codiv(sim$host_tree, sim$symbiont_tree, sim$links,
                      min_hosts = 3, min_symbiont_tips = 7,
                      span_fraction = span_fraction, permutations = permutations,
                      methods = methods, cores = 2, seed = 1,
                      subtree_features = FALSE, verbose = FALSE))
}

# ---------------------------------------------------------------------------
# Analysis 1: Type I error under the null (all clades non-co-diversifying)
# ---------------------------------------------------------------------------
analysis_type1 <- function(reps = 20, n_hosts = 10, alpha = 0.01) {
  pvals <- c()
  for (r in seq_len(reps)) {
    sim <- simulate_codiv_data(n_hosts = n_hosts, clades = replicate(4,
      list(codiversifying = FALSE, n_per_host = 1), simplify = FALSE),
      seed = r)
    pvals <- c(pvals, run_scan(sim, span_fraction = 1)$Hommola_pvalue)
  }
  cat(sprintf("Type I: FPR at alpha=%.3g is %.3g (target %.3g)\n",
              alpha, mean(pvals < alpha, na.rm = TRUE), alpha))
  invisible(pvals)  # for the paper: QQ-plot vs uniform, KS test
}

# ---------------------------------------------------------------------------
# Analysis 2: Power vs. congruence
# ---------------------------------------------------------------------------
analysis_power <- function(reps = 10, congruence_grid = c(0, 0.25, 0.5, 0.75, 1),
                           n_hosts = 12, alpha = 0.01, r_threshold = 0.75) {
  out <- data.frame()
  for (cg in congruence_grid) {
    detected <- logical(reps)
    for (r in seq_len(reps)) {
      sim <- simulate_codiv_data(n_hosts = n_hosts, clades = list(
        list(codiversifying = TRUE, congruence = cg, n_per_host = 1)),
        seed = r)
      res <- run_scan(sim, span_fraction = 1)
      detected[r] <- any(res$Hommola_r > r_threshold &
                           res$Hommola_pvalue < alpha, na.rm = TRUE)
    }
    out <- rbind(out, data.frame(congruence = cg, power = mean(detected)))
  }
  print(out); invisible(out)
}

# ---------------------------------------------------------------------------
# Analysis 3: Localization AUC at a single setting
# ---------------------------------------------------------------------------
analysis_localization <- function(n_hosts = 10, congruence = 1, n_per_host = 4) {
  sim <- simulate_codiv_data(n_hosts = n_hosts, clades = list(
    list(codiversifying = TRUE,  congruence = congruence, n_per_host = 1),
    list(codiversifying = TRUE,  congruence = congruence, n_per_host = n_per_host),
    list(codiversifying = FALSE, n_per_host = 1),
    list(codiversifying = FALSE, n_per_host = n_per_host)), seed = 1)
  res <- run_scan(sim)
  nt <- node_truth(res, sim); keep <- !is.na(nt)
  cat(sprintf("localization: %d nodes (%d co-diversifying)\n",
              sum(keep), sum(nt[keep])))
  cat(sprintf("  AUC uncollapsed: %.3f | collapsed: %.3f\n",
              auc(res$Hommola_r[keep], nt[keep]),
              auc(res$Collapsed_Hommola_r[keep], nt[keep])))
  invisible(res)
}

# ---------------------------------------------------------------------------
# Headline: congruence sweep -> node-level AUC curves
# ---------------------------------------------------------------------------
# For each congruence value, build a tree with co-diversifying clades at that
# congruence plus random controls, and measure how well the scan separates
# co-diversifying from random nodes. Two statistics (uncollapsed vs collapsed)
# and two within-host settings (n_per_host 1 vs 4). The collapsed-vs-uncollapsed
# gap at n_per_host > 1 is the within-host pseudoreplication effect.
analysis_congruence_sweep <- function(congruence_grid = seq(0, 1, 0.2),
                                      n_per_host_vals = c(1, 4), reps = 5,
                                      n_hosts = 10, permutations = 199,
                                      out_png = "simulations/auc_sweep.png") {
  rows <- list()
  for (nph in n_per_host_vals) {
    for (cg in congruence_grid) {
      unc <- col <- numeric(reps)
      for (r in seq_len(reps)) {
        sim <- simulate_codiv_data(n_hosts = n_hosts, clades = list(
          list(codiversifying = TRUE,  congruence = cg, n_per_host = nph),
          list(codiversifying = TRUE,  congruence = cg, n_per_host = nph),
          list(codiversifying = FALSE, n_per_host = nph),
          list(codiversifying = FALSE, n_per_host = nph)),
          seed = r * 1000 + round(cg * 100) + nph)
        res <- run_scan(sim, permutations = permutations)
        nt <- node_truth(res, sim); keep <- !is.na(nt)
        unc[r] <- auc(res$Hommola_r[keep], nt[keep])
        col[r] <- auc(res$Collapsed_Hommola_r[keep], nt[keep])
      }
      rows[[length(rows) + 1]] <- data.frame(
        congruence = cg, n_per_host = nph,
        statistic = c("uncollapsed", "collapsed"),
        auc = c(mean(unc, na.rm = TRUE), mean(col, na.rm = TRUE)))
    }
  }
  out <- do.call(rbind, rows)
  print(out)

  if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(out_png)) {
    out$facet <- paste0("n_per_host = ", out$n_per_host)
    p <- ggplot2::ggplot(out, ggplot2::aes(congruence, auc,
                                           color = statistic)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::facet_wrap(~ facet) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted") +
      ggplot2::labs(x = "Congruence of co-diversifying clades",
                    y = "Node-level AUC (localization)",
                    title = "Localization accuracy vs. congruence") +
      ggplot2::ylim(0.4, 1) + ggplot2::theme_bw()
    ggplot2::ggsave(out_png, p, width = 8, height = 4, dpi = 110)
    cat("saved", out_png, "\n")
  }
  invisible(out)
}

# ---------------------------------------------------------------------------
# Analysis 4: FDR calibration - naive per-node vs. second-order permutation
# ---------------------------------------------------------------------------
analysis_fdr_calibration <- function(reps = 10, n_hosts = 12, n_per_host = 4) {
  for (r in seq_len(reps)) {
    sim <- simulate_codiv_data(n_hosts = n_hosts, clades = replicate(4,
      list(codiversifying = FALSE, n_per_host = n_per_host), simplify = FALSE),
      seed = r)
    res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, min_hosts = 3,
                 min_symbiont_tips = 7, span_fraction = 1, permutations = 999,
                 methods = "hommola", cores = 2, seed = r,
                 subtree_features = FALSE, verbose = FALSE)
    naive <- sum(as.data.frame(res)$Hommola_pvalue < 0.05, na.rm = TRUE)
    null <- codiv_null_scans(sim$host_tree, sim$symbiont_tree, sim$links,
                             n_permutations = 50, observed = res,
                             stat_threshold = 0.75, p_threshold = 0.01,
                             min_hosts = 3, min_symbiont_tips = 7,
                             span_fraction = 1, permutations = 999,
                             methods = "hommola", cores = 2, seed = r,
                             verbose = FALSE)
    cat(sprintf("rep %d | naive p<0.05: %d | observed: %d | null mean: %.1f | global p: %.3f\n",
                r, naive, null$observed, mean(null$null_counts),
                null$global_pvalue))
  }
}

# ---------------------------------------------------------------------------
# Run small versions (scale up for the paper)
# ---------------------------------------------------------------------------
if (identical(environment(), globalenv())) {
  set.seed(1)
  analysis_type1(reps = 10)
  analysis_power(reps = 5)
  analysis_localization()
  analysis_congruence_sweep(reps = 3)
  analysis_fdr_calibration(reps = 3)
}
