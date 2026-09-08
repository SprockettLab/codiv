#' Summarize which hosts have co-diversifying symbionts
#'
#' Reports, for each host, how many co-diversifying clades it participates in.
#' A node counts as co-diversifying when its statistic exceeds `stat_threshold`
#' and its p-value is below `p_threshold` (the defaults, `Hommola_r > 0.75` and
#' `p < 0.01`, follow Sanders et al. 2023 and Sprockett et al. 2025). The hosts
#' in each node are read from the per-node host tree stored by [codiv()].
#'
#' @param codiv_results A `codiv` result (or data frame) from [codiv()]
#' @param statistic Name of the test-statistic column. NULL (default) uses
#'   `Hommola_r`, the collapsed Hommola correlation that [codiv()] always
#'   writes when the Hommola method runs
#' @param stat_threshold A node is co-diversifying when its statistic exceeds
#'   this value; default 0.75
#' @param pvalue Name of the p-value column. NULL (default) resolves
#'   alongside `statistic`
#' @param p_threshold A node is co-diversifying when its p-value is below this
#'   value; default 0.01
#'
#' @return A data frame with one row per host (ordered by
#'   `n_codiversifying_nodes`, descending):
#'   - `Host`: host tip label
#'   - `n_codiversifying_nodes`: co-diversifying nodes that include this host
#'   - `n_scanned_nodes`: all scanned nodes that include this host
#'   - `prop_codiversifying`: `n_codiversifying_nodes / n_scanned_nodes`
#'   - `max_statistic`: the largest statistic value among this host's
#'     co-diversifying nodes (NA if none)
#'
#' @details
#' Scanned nodes are nested, so a host appears in several overlapping nodes; the
#' counts are a descriptive ranking of which hosts carry the most co-diversifying
#' symbiont lineages, not independent tallies.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sim <- simulate_codiv_data(n_hosts = 10, clades = list(
#'   list(codiversifying = TRUE,  congruence = 1, n_per_host = 2),
#'   list(codiversifying = FALSE, n_per_host = 2)), seed = 3)
#' results <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
#'                  methods = "hommola", permutations = 99, verbose = FALSE)
#' host_codiv_summary(results)
#' }
#'
host_codiv_summary <- function(codiv_results, statistic = NULL,
                               stat_threshold = 0.75, pvalue = NULL,
                               p_threshold = 0.01) {
  # NULL means "the Hommola statistic", resolved against what the result
  # actually carries, since the uncollapsed columns are opt-in.
  .hc <- .default_hommola_cols(codiv_results)
  if (is.null(statistic)) statistic <- unname(.hc[["stat"]])
  if (is.null(pvalue)) pvalue <- unname(.hc[["pvalue"]])
  df <- as.data.frame(codiv_results)
  if (!all(c(statistic, pvalue, "Host_Tree") %in% colnames(df))) {
    stop("Result is missing the '", statistic, "', '", pvalue,
         "', or 'Host_Tree' column; set `statistic`/`pvalue` to match the ",
         "method(s) used.")
  }

  # hosts present in each scanned node (from the stored per-node host tree)
  hosts_per_node <- lapply(df$Host_Tree, function(nw) {
    if (is.na(nw) || !nzchar(nw)) return(character(0))
    tryCatch(ape::read.tree(text = nw)$tip.label, error = function(e) character(0))
  })

  is_codiv <- df[[statistic]] > stat_threshold & df[[pvalue]] < p_threshold
  is_codiv[is.na(is_codiv)] <- FALSE

  all_hosts <- sort(unique(unlist(hosts_per_node)))
  if (length(all_hosts) == 0) {
    return(data.frame(Host = character(0), n_codiversifying_nodes = integer(0),
                      n_scanned_nodes = integer(0),
                      prop_codiversifying = numeric(0),
                      max_statistic = numeric(0)))
  }

  in_node <- vapply(all_hosts, function(h)
    vapply(hosts_per_node, function(x) h %in% x, logical(1)),
    logical(length(hosts_per_node)))
  # in_node is nodes x hosts
  if (is.null(dim(in_node))) in_node <- matrix(in_node, nrow = nrow(df))

  out <- data.frame(
    Host = all_hosts,
    n_codiversifying_nodes = colSums(in_node & is_codiv),
    n_scanned_nodes = colSums(in_node),
    stringsAsFactors = FALSE)
  out$prop_codiversifying <- out$n_codiversifying_nodes / out$n_scanned_nodes
  out$max_statistic <- vapply(all_hosts, function(h) {
    vals <- df[[statistic]][in_node[, h] & is_codiv]
    if (length(vals) == 0) NA_real_ else max(vals, na.rm = TRUE)
  }, numeric(1))

  out <- out[order(-out$n_codiversifying_nodes, -out$prop_codiversifying), ]
  rownames(out) <- NULL
  out
}
