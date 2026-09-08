#' Leave-One-Host-Out Analysis
#'
#' Reruns [codiv()] with each host removed in turn and compares each run to the
#' full analysis, showing how much each host contributes to the
#' codiversification signal.
#'
#' @param Host_tree A phylo object; binary phylogenetic tree of hosts
#' @param Symbiont_tree A phylo object; binary phylogenetic tree of symbionts
#' @param Host_to_Symbiont_df A data frame with columns "Host" and "Symbiont"
#'   linking symbionts to hosts
#' @param hosts_to_test Character vector of host names to test. If NULL (default),
#'   tests all hosts. Use a subset to reduce computation time.
#' @param alpha Significance threshold for counting significant nodes; default 0.05
#' @param verbose If TRUE (default), report progress across hosts
#' @param ... Additional arguments passed to [codiv()] (e.g. min_hosts,
#'   min_symbiont_tips, span_fraction, permutations, methods, cores, seed)
#'
#' @return A list with elements:
#'   - `full_results`: the [codiv()] result on the complete data
#'   - `loo_results`: named list of the [codiv()] result with each host removed
#'   - `host_importance`: one row per tested host (ordered by `impact_score`),
#'     with the number of nodes scanned, per-method significant-node counts and
#'     their change from the full analysis (`<method>_delta_n_sig`), and an
#'     `impact_score` (mean absolute change in significant-node count across
#'     methods when the host is removed)
#'
#' @details
#' Each run uses its own isolated checkpoint file so removing a host genuinely
#' re-computes the affected nodes. Tree-distance metrics are disabled for speed.
#' Runtime is roughly `baseline_time x length(hosts_to_test)`, so consider
#' `hosts_to_test`, fewer permutations, or a single method for a first pass, and
#' `cores` to parallelize each run.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sim <- simulate_codiv_data(n_hosts = 8, n_clades = 2, seed = 1)
#' loo <- loo_host_analysis(sim$host_tree, sim$symbiont_tree, sim$links,
#'                          permutations = 99, methods = "hommola",
#'                          verbose = FALSE)
#' loo$host_importance
#' }
#'
loo_host_analysis <- function(Host_tree, Symbiont_tree, Host_to_Symbiont_df,
                              hosts_to_test = NULL, alpha = 0.05,
                              verbose = TRUE, ...) {
  Host_to_Symbiont_df$Host <- as.character(Host_to_Symbiont_df$Host)
  Host_to_Symbiont_df$Symbiont <- as.character(Host_to_Symbiont_df$Symbiont)

  all_hosts <- unique(Host_to_Symbiont_df$Host)
  if (is.null(hosts_to_test)) {
    hosts_to_test <- all_hosts
  }
  missing_hosts <- setdiff(hosts_to_test, all_hosts)
  if (length(missing_hosts) > 0) {
    stop("hosts_to_test not present in Host_to_Symbiont_df: ",
         paste(missing_hosts, collapse = ", "))
  }

  dots <- list(...)
  # each run gets an isolated checkpoint and does not resume, so removing a host
  # genuinely recomputes the affected nodes rather than reusing full-data values
  run_codiv <- function(ht, st, df) {
    args <- utils::modifyList(list(subtree_features = FALSE), dots)
    args$Host_tree <- ht
    args$Symbiont_tree <- st
    args$Host_to_Symbiont_df <- df
    args$verbose <- FALSE
    args$Save_fp <- tempfile(fileext = ".tsv")
    args$continue <- FALSE
    do.call(codiv, args)
  }

  if (verbose) message("Running full analysis ...")
  full_results <- run_codiv(Host_tree, Symbiont_tree, Host_to_Symbiont_df)
  methods <- attr(full_results, "codiv_params")$methods
  base_stats <- .loo_method_stats(full_results, methods, alpha)

  loo_results <- vector("list", length(hosts_to_test))
  names(loo_results) <- hosts_to_test
  rows <- vector("list", length(hosts_to_test))

  for (i in seq_along(hosts_to_test)) {
    host_i <- hosts_to_test[i]
    if (verbose) {
      message("[", i, "/", length(hosts_to_test), "] removing host ", host_i)
    }
    df_i <- Host_to_Symbiont_df[Host_to_Symbiont_df$Host != host_i, ]
    tree_i <- ape::drop.tip(Host_tree, host_i)
    res_i <- run_codiv(tree_i, Symbiont_tree, df_i)
    loo_results[[host_i]] <- res_i

    stats_i <- .loo_method_stats(res_i, methods, alpha)
    row <- data.frame(host = host_i, n_nodes = nrow(res_i),
                      stringsAsFactors = FALSE)
    for (m in methods) {
      row[[paste0(m, "_n_sig")]] <- stats_i[[m]]$n_sig
      row[[paste0(m, "_delta_n_sig")]] <- stats_i[[m]]$n_sig -
        base_stats[[m]]$n_sig
      row[[paste0(m, "_mean_stat")]] <- stats_i[[m]]$mean_stat
    }
    deltas <- vapply(methods,
                     function(m) abs(stats_i[[m]]$n_sig - base_stats[[m]]$n_sig),
                     numeric(1))
    row$impact_score <- mean(deltas, na.rm = TRUE)
    rows[[i]] <- row
  }

  host_importance <- do.call(rbind, rows)
  host_importance <- host_importance[order(-host_importance$impact_score), ]
  rownames(host_importance) <- NULL

  list(full_results = full_results,
       loo_results = loo_results,
       host_importance = host_importance)
}

# Per-method significant-node count and mean primary statistic for one result.
# inputs: res (codiv object); methods (character); alpha
# output: named list (by method) of list(n_sig, mean_stat)
.loo_method_stats <- function(res, methods, alpha) {
  stats::setNames(lapply(methods, function(m) {
    sig <- .codiv_significant(res, m, alpha)
    mc <- .codiv_method_cols(m, colnames(res))
    mean_stat <- if (nrow(res) > 0 && "stat" %in% names(mc)) {
      mean(res[[mc[["stat"]]]], na.rm = TRUE)
    } else {
      NA_real_
    }
    list(n_sig = if (is.null(sig)) NA_integer_ else sig$n_sig,
         mean_stat = mean_stat)
  }), methods)
}
