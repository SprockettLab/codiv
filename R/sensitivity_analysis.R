#' Parameter Sensitivity Analysis
#'
#' Reruns [codiv()] across a grid of `span_fraction` and `min_symbiont_tips`
#' values and summarizes how the number of scanned and significant nodes changes,
#' so you can see whether findings are stable or threshold-dependent.
#'
#' @param Host_tree A phylo object; binary phylogenetic tree of hosts
#' @param Symbiont_tree A phylo object; binary phylogenetic tree of symbionts
#' @param Host_to_Symbiont_df A data frame with columns "Host" and "Symbiont"
#' @param span_fraction_range Numeric vector of `span_fraction` values to test;
#'   default c(0.05, 0.1, 0.2, 0.3)
#' @param min_symbiont_tips_range Integer vector of `min_symbiont_tips` values;
#'   default c(5, 7, 10, 15)
#' @param permutations Integer; number of permutations, applied to every run
#' @param min_hosts Integer; minimum hosts, applied to every run (>= 3)
#' @param methods Character vector of methods to use
#' @param alpha Significance threshold for counting significant nodes; default 0.05
#' @param cores Number of cores passed to [codiv()]
#' @param seed Random seed for reproducibility
#' @param verbose If TRUE (default), report progress across the grid
#'
#' @return A list with elements:
#'   - `sensitivity_data`: one row per parameter combination, with the number of
#'     nodes scanned and, per method, the number of significant nodes and the
#'     mean of the primary statistic
#'   - `results_by_params`: named list of the full [codiv()] result for each
#'     parameter combination
#'   - `heatmap_data`: long-format data (one row per combination and method) with
#'     `n_sig`, `n_nodes`, and `prop_sig`, ready for plotting
#'
#' @details
#' Runtime is roughly `baseline_time x length(span_fraction_range) x
#' length(min_symbiont_tips_range)`. Tree-distance metrics are turned off for
#' speed, since sensitivity is judged from node counts and significance.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
#' sens <- sensitivity_analysis(sim$host_tree, sim$symbiont_tree, sim$links,
#'                              span_fraction_range = c(0.3, 0.6),
#'                              min_symbiont_tips_range = c(7, 10),
#'                              permutations = 99, methods = "hommola",
#'                              verbose = FALSE)
#' sens$sensitivity_data
#' }
#'
sensitivity_analysis <- function(Host_tree, Symbiont_tree, Host_to_Symbiont_df,
                                 span_fraction_range = c(0.05, 0.1, 0.2, 0.3),
                                 min_symbiont_tips_range = c(5, 7, 10, 15),
                                 permutations = 99, min_hosts = 3,
                                 methods = c("hommola", "paco", "parafit"),
                                 alpha = 0.05,
                                 cores = 1, seed = 8675309, verbose = TRUE) {

  if (!is.numeric(span_fraction_range) || length(span_fraction_range) < 1 ||
      any(span_fraction_range <= 0 | span_fraction_range > 1)) {
    stop("`span_fraction_range` must be numbers in (0, 1]")
  }
  if (!is.numeric(min_symbiont_tips_range) ||
      length(min_symbiont_tips_range) < 1 ||
      any(min_symbiont_tips_range < 1 |
          min_symbiont_tips_range != round(min_symbiont_tips_range))) {
    stop("`min_symbiont_tips_range` must be positive integers")
  }

  # every combination of the two ranges
  grid <- expand.grid(span_fraction = span_fraction_range,
                      min_symbiont_tips = min_symbiont_tips_range,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  n_combos <- nrow(grid)

  results_by_params <- vector("list", n_combos)
  summary_rows <- vector("list", n_combos)
  heatmap_rows <- vector("list", n_combos)

  for (k in seq_len(n_combos)) {
    sf <- grid$span_fraction[k]
    mt <- grid$min_symbiont_tips[k]
    if (verbose) {
      message("[", k, "/", n_combos, "] span_fraction = ", sf,
              ", min_symbiont_tips = ", mt)
    }

    # each combination runs independently with its own checkpoint file
    res <- codiv(Host_tree, Symbiont_tree, Host_to_Symbiont_df,
                 min_hosts = min_hosts, min_symbiont_tips = mt,
                 span_fraction = sf, permutations = permutations,
                 seed = seed, methods = methods, subtree_features = FALSE,
                 cores = cores, verbose = FALSE,
                 Save_fp = tempfile(fileext = ".tsv"), continue = FALSE)

    n_nodes <- nrow(res)
    results_by_params[[k]] <- res

    # per-method significance and mean statistic
    per_method <- lapply(methods, function(m) {
      sig <- .codiv_significant(res, m, alpha)
      mc <- .codiv_method_cols(m, colnames(res))
      mean_stat <- if (n_nodes > 0 && "stat" %in% names(mc)) {
        mean(res[[mc[["stat"]]]], na.rm = TRUE)
      } else {
        NA_real_
      }
      n_sig <- if (is.null(sig)) NA_integer_ else sig$n_sig
      list(method = m, n_sig = n_sig, mean_stat = mean_stat)
    })

    # wide summary row
    wide <- data.frame(span_fraction = sf, min_symbiont_tips = mt,
                       n_nodes = n_nodes, stringsAsFactors = FALSE)
    for (pm in per_method) {
      wide[[paste0(pm$method, "_n_sig")]] <- pm$n_sig
      wide[[paste0(pm$method, "_mean_stat")]] <- pm$mean_stat
    }
    summary_rows[[k]] <- wide

    # long heatmap rows
    heatmap_rows[[k]] <- do.call(rbind, lapply(per_method, function(pm) {
      data.frame(span_fraction = sf, min_symbiont_tips = mt,
                 method = pm$method, n_sig = pm$n_sig, n_nodes = n_nodes,
                 prop_sig = if (n_nodes > 0 && !is.na(pm$n_sig)) {
                   pm$n_sig / n_nodes
                 } else {
                   NA_real_
                 },
                 stringsAsFactors = FALSE)
    }))
  }

  names(results_by_params) <- paste0("span", grid$span_fraction,
                                     "_tips", grid$min_symbiont_tips)

  list(
    sensitivity_data = do.call(rbind, summary_rows),
    results_by_params = results_by_params,
    heatmap_data = do.call(rbind, heatmap_rows)
  )
}
