#!/usr/bin/env Rscript
# Scaled-up congruence sweep for the codiv methods paper
# ------------------------------------------------------
# Node-level localization AUC (uncollapsed vs collapsed statistic) as a function
# of the REALIZED congruence between co-diversifying clades and the host tree,
# faceted by within-host copy number (n_per_host). Realized congruence is
# measured per dataset (normalized mutual clustering information of the clade's
# host-relabeled topology vs the host tree), which is a cleaner x-axis than the
# input SPR-based congruence knob.
#
# Cluster-friendly: results are appended to a CSV after every dataset, so a
# preempted job resumes where it left off. Configure via environment variables
# or edit the defaults below.
#
#   REPS=50 NHOSTS=20 CORES=8 PERMS=999 Rscript simulations/run_congruence_sweep.R
#
# Then plot:  Rscript simulations/run_congruence_sweep.R plot

suppressMessages({ library(ape); library(codiv); library(TreeDist) })

cfg <- list(
  reps        = as.integer(Sys.getenv("REPS", "50")),
  n_hosts     = as.integer(Sys.getenv("NHOSTS", "20")),
  cores       = as.integer(Sys.getenv("CORES", "4")),
  permutations= as.integer(Sys.getenv("PERMS", "999")),
  n_per_host  = as.integer(strsplit(Sys.getenv("NPERHOST", "1,4"), ",")[[1]]),
  # input congruence values spread the realized congruence across [0,1]
  congruence  = seq(0, 1, by = 0.1),
  out_csv     = Sys.getenv("OUT_CSV", "simulations/congruence_sweep.csv"),
  out_png     = Sys.getenv("OUT_PNG", "simulations/congruence_sweep.png")
)

# --- helpers ---------------------------------------------------------------

# realized congruence: mean normalized mutual clustering info between each
# co-diversifying clade (one tip per host, relabeled with host) and the host
# tree (1 = identical branching order)
realized_congruence <- function(sim) {
  codiv_clades <- unique(sim$truth$clade[sim$truth$codiversifying])
  vals <- vapply(codiv_clades, function(cl) {
    symb <- sim$truth$Symbiont[sim$truth$clade == cl]
    lk <- sim$links[sim$links$Symbiont %in% symb, ]
    reps <- lk$Symbiont[!duplicated(lk$Host)]
    g <- keep.tip(sim$symbiont_tree, reps)
    g$tip.label <- lk$Host[match(g$tip.label, lk$Symbiont)]
    shared <- intersect(sim$host_tree$tip.label, g$tip.label)
    if (length(shared) < 3) return(NA_real_)
    h <- reorder(keep.tip(sim$host_tree, shared), "cladewise")
    g <- reorder(keep.tip(g, shared), "cladewise")
    as.numeric(MutualClusteringInfo(h, g, normalize = TRUE))
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

node_truth <- function(res_df, sim) {
  clade_of <- stats::setNames(sim$truth$clade, sim$truth$Symbiont)
  truth_of <- stats::setNames(
    tapply(sim$truth$codiversifying, sim$truth$clade, `[`, 1),
    unique(sim$truth$clade))
  vapply(res_df$Symbiont_Tree, function(nw) {
    tips <- tryCatch(read.tree(text = nw)$tip.label, error = function(e) NA)
    if (all(is.na(tips))) return(NA)
    unname(truth_of[names(sort(table(clade_of[tips]), decreasing = TRUE))[1]])
  }, logical(1))
}

auc <- function(score, truth) {
  pos <- score[truth]; neg <- score[!truth]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
}

# --- run one dataset -------------------------------------------------------
run_one <- function(nph, cg, rep) {
  sim <- simulate_codiv_data(n_hosts = cfg$n_hosts, clades = list(
    list(codiversifying = TRUE,  congruence = cg, n_per_host = nph),
    list(codiversifying = TRUE,  congruence = cg, n_per_host = nph),
    list(codiversifying = FALSE, n_per_host = nph),
    list(codiversifying = FALSE, n_per_host = nph)),
    seed = rep * 10000 + round(cg * 100) + nph)
  res <- as.data.frame(codiv(sim$host_tree, sim$symbiont_tree, sim$links,
    min_hosts = 3, min_symbiont_tips = 7, span_fraction = 0.6,
    permutations = cfg$permutations, methods = "hommola", cores = cfg$cores,
    seed = 1, subtree_features = FALSE, verbose = FALSE))
  nt <- node_truth(res, sim); keep <- !is.na(nt)
  data.frame(
    n_per_host = nph, input_congruence = cg, rep = rep,
    realized_congruence = realized_congruence(sim),
    n_nodes = sum(keep),
    auc_uncollapsed = auc(res$Hommola_r[keep], nt[keep]),
    auc_collapsed   = auc(res$Collapsed_Hommola_r[keep], nt[keep]))
}

# --- sweep with incremental CSV (resumable) --------------------------------
run_sweep <- function() {
  done <- if (file.exists(cfg$out_csv)) {
    d <- utils::read.csv(cfg$out_csv)
    paste(d$n_per_host, d$input_congruence, d$rep)
  } else character(0)

  grid <- expand.grid(nph = cfg$n_per_host, cg = cfg$congruence,
                      rep = seq_len(cfg$reps))
  for (i in seq_len(nrow(grid))) {
    key <- paste(grid$nph[i], grid$cg[i], grid$rep[i])
    if (key %in% done) next
    row <- run_one(grid$nph[i], grid$cg[i], grid$rep[i])
    if (file.exists(cfg$out_csv)) {
      utils::write.table(row, cfg$out_csv, sep = ",", row.names = FALSE,
                         col.names = FALSE, append = TRUE)
    } else {
      utils::write.csv(row, cfg$out_csv, row.names = FALSE)
    }
    cat(sprintf("[%d/%d] nph=%d cg=%.1f rep=%d | realized=%.2f auc: unc=%.2f col=%.2f\n",
                i, nrow(grid), grid$nph[i], grid$cg[i], grid$rep[i],
                row$realized_congruence, row$auc_uncollapsed, row$auc_collapsed))
  }
  cat("sweep complete:", cfg$out_csv, "\n")
}

# --- plot AUC vs realized congruence ---------------------------------------
make_plot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plotting")
  }
  d <- utils::read.csv(cfg$out_csv)
  long <- rbind(
    data.frame(d[c("n_per_host", "realized_congruence")],
               statistic = "uncollapsed", auc = d$auc_uncollapsed),
    data.frame(d[c("n_per_host", "realized_congruence")],
               statistic = "collapsed", auc = d$auc_collapsed))
  long$facet <- paste0("n_per_host = ", long$n_per_host)
  p <- ggplot2::ggplot(long, ggplot2::aes(realized_congruence, auc,
                                          color = statistic)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggplot2::geom_smooth(se = TRUE, method = "loess", formula = y ~ x) +
    ggplot2::facet_wrap(~ facet) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted") +
    ggplot2::labs(x = "Realized congruence (normalized mutual clustering info)",
                  y = "Node-level localization AUC",
                  title = "Localization accuracy vs. realized congruence") +
    ggplot2::coord_cartesian(ylim = c(0.4, 1)) + ggplot2::theme_bw()
  ggplot2::ggsave(cfg$out_png, p, width = 8, height = 4, dpi = 120)
  cat("saved", cfg$out_png, "\n")
}

# --- entry point -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] == "plot") {
  make_plot()
} else {
  run_sweep()
  if (requireNamespace("ggplot2", quietly = TRUE)) make_plot()
}
