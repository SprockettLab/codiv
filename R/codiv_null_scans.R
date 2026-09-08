# Format a duration in seconds as a compact string (s / m / h).
.fmt_duration <- function(secs) {
  if (!is.finite(secs)) return("?")
  if (secs < 90) return(sprintf("%.0fs", secs))
  if (secs < 5400) return(sprintf("%.1fm", secs / 60))
  sprintf("%.1fh", secs / 3600)
}

#' Second-order permutation test for scan-wide co-diversification
#'
#' Assesses whether a co-diversification scan detects more significant clades
#' than expected by chance, accounting for the phylogenetic non-independence and
#' pseudoreplication of nodes that per-node FDR corrections (e.g. Benjamini-
#' Hochberg) ignore. The host tree's tip labels are permuted and the entire
#' [codiv()] scan is rerun; repeating this builds a null distribution for the
#' *number* of significant clades. Because each permutation reruns the whole
#' pipeline, the null captures exactly the dependence structure that a per-node
#' correction does not.
#'
#' This is the recommended way to gauge false discoveries for a whole scan
#' (Sanders et al. 2023; Sprockett et al. 2025). [codiv()] deliberately does not
#' apply a per-node multiple-testing correction, which would be invalid for
#' non-independent nested clades.
#'
#' @param Host_tree A phylo object; binary host tree
#' @param Symbiont_tree A phylo object; binary symbiont tree
#' @param Host_to_Symbiont_df A data frame with "Host" and "Symbiont" columns
#' @param n_permutations Number of host-label permutations (outer null); the
#'   analyses in the source papers used 100
#' @param statistic Name of the test-statistic column used to call a node
#'   significant. NULL (default) uses `Hommola_r`, the collapsed Hommola
#'   correlation that [codiv()] always writes when the Hommola method runs
#' @param stat_threshold A node counts as significant when its statistic exceeds
#'   this value; default 0.75
#' @param pvalue Name of the p-value column used to call a node significant;
#'   NULL (default) resolves alongside `statistic`
#' @param p_threshold A node counts as significant when its p-value is below this
#'   value; default 0.01
#' @param observed Optional precomputed [codiv()] result for the real (unpermuted)
#'   data; if NULL it is computed
#' @param seed Random seed for reproducibility
#' @param verbose If TRUE (default), report progress across permutations
#' @param ... Additional arguments passed to [codiv()] (e.g. min_hosts,
#'   min_symbiont_tips, span_fraction, permutations, methods, cores)
#'
#' @return A list with:
#'   - `observed`: number of significant clades in the real data
#'   - `null_counts`: number of significant clades in each permuted scan
#'   - `global_pvalue`: (permutations with count >= observed + 1) / (n + 1)
#'   - `empirical_fdr`: expected false-discovery proportion at the chosen
#'     thresholds, `mean(null_counts) / observed` (capped at 1)
#'   - `thresholds`: the statistic/p-value columns and cutoffs used
#'
#' @details
#' Runtime is roughly `(n_permutations + 1) x` a single scan, so use `cores` and
#' consider fewer inner `permutations` for a first pass. Each rerun uses its own
#' checkpoint file and does not resume.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sim <- simulate_codiv_data(n_hosts = 8, n_clades = 2, seed = 1)
#' null <- codiv_null_scans(sim$host_tree, sim$symbiont_tree, sim$links,
#'                          n_permutations = 10, min_hosts = 3,
#'                          min_symbiont_tips = 7, span_fraction = 0.6,
#'                          permutations = 49, methods = "hommola",
#'                          verbose = FALSE)
#' null$observed
#' null$empirical_fdr
#' }
#'
codiv_null_scans <- function(Host_tree, Symbiont_tree, Host_to_Symbiont_df,
                             n_permutations = 100,
                             statistic = NULL, stat_threshold = 0.75,
                             pvalue = NULL, p_threshold = 0.01,
                             observed = NULL, seed = 8675309, verbose = TRUE,
                             ...) {
  # NULL means "the Hommola statistic", resolved against what the result
  # actually carries, since the uncollapsed columns are opt-in.
  .hc <- .default_hommola_cols(observed)
  if (is.null(statistic)) statistic <- unname(.hc[["stat"]])
  if (is.null(pvalue)) pvalue <- unname(.hc[["pvalue"]])
  dots <- list(...)
  run_scan <- function(ht) {
    args <- utils::modifyList(list(), dots)
    args$Host_tree <- ht
    args$Symbiont_tree <- Symbiont_tree
    args$Host_to_Symbiont_df <- Host_to_Symbiont_df
    args$verbose <- FALSE
    args$Save_fp <- tempfile(fileext = ".tsv")
    args$continue <- FALSE
    do.call(codiv, args)
  }

  count_significant <- function(res) {
    if (!all(c(statistic, pvalue) %in% colnames(res))) {
      stop("Result is missing the '", statistic, "' or '", pvalue,
           "' column; set `statistic`/`pvalue` to match the method(s) used.")
    }
    sum(res[[statistic]] > stat_threshold & res[[pvalue]] < p_threshold,
        na.rm = TRUE)
  }

  # observed scan
  if (is.null(observed)) {
    if (verbose) message("Running the observed scan ...")
    observed <- run_scan(Host_tree)
  }
  observed_count <- count_significant(observed)

  # preserve the caller's global RNG
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    saved_rng <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", saved_rng, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }
  if (!is.na(seed)) set.seed(seed)

  # null: permute host tip labels, rerun the whole scan, count significant
  # clades. Each inner scan draws its own progress bar (the intra-scan bar);
  # after each permutation completes we print a counter with elapsed time and an
  # estimated time remaining, so progress is clear over a long run.
  null_counts <- integer(n_permutations)
  start_time <- Sys.time()
  for (i in seq_len(n_permutations)) {
    perm_host <- Host_tree
    perm_host$tip.label <- sample(Host_tree$tip.label)
    null_counts[i] <- count_significant(run_scan(perm_host))
    if (verbose) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      eta <- elapsed / i * (n_permutations - i)
      message(sprintf("Permutation %d/%d complete | elapsed %s | ETA ~%s",
                      i, n_permutations, .fmt_duration(elapsed),
                      .fmt_duration(eta)))
    }
  }

  empirical_fdr <- if (observed_count > 0) {
    min(mean(null_counts) / observed_count, 1)
  } else {
    NA_real_
  }

  list(
    observed = observed_count,
    null_counts = null_counts,
    global_pvalue = (sum(null_counts >= observed_count) + 1) / (n_permutations + 1),
    empirical_fdr = empirical_fdr,
    thresholds = list(statistic = statistic, stat_threshold = stat_threshold,
                      pvalue = pvalue, p_threshold = p_threshold)
  )
}
