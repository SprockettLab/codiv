# S3 methods for the `codiv` object returned by codiv().
# The object is a data frame subclass, so all data-frame operations still work;
# these methods provide readable summaries in place of dumping every column.

# The Hommola statistic column a result should be read from by default.
# The uncollapsed test is opt-in (see codiv()'s uncollapsed_hommola), so
# Uncollapsed_Hommola_r is frequently absent while Hommola_r is written whenever
# the method ran. Functions that default to "the Hommola statistic" resolve it
# here rather than hard-coding a name that may not exist.
# inputs: df (a codiv result or data frame)
# output: named character vector with stat and pvalue column names
.default_hommola_cols <- function(df) {
  c(stat = "Hommola_r", pvalue = "Hommola_pvalue")
}

# Which methods a result carries, inferred from its columns.
# Splitting on the first underscore does not work: the collapsed columns are
# named Hommola_r and so on, so that yields "Collapsed" for every
# method. Match on the method token wherever it appears instead.
# inputs: cols (character vector of column names)
# output: character vector of method names, in canonical order
.codiv_methods_from_cols <- function(cols) {
  known <- c(hommola = "Hommola", paco = "PACo",
             parafit = "ParaFit", topology = "Topology")
  names(known)[vapply(known, function(tok) any(grepl(tok, cols, fixed = TRUE)),
                      logical(1))]
}

# Primary statistic and p-value column names for a method, keeping only those
# actually present in the results.
# inputs: method (one of "hommola", "paco", "parafit", "topology"); cols
# output: named character vector with any of stat/pvalue that are present
.codiv_method_cols <- function(method, cols) {
  # The uncollapsed Hommola test is opt-in, so Uncollapsed_Hommola_r is often absent while
  # Hommola_r is always written when the method ran. Fall back to the
  # collapsed pair rather than reporting the method as missing.
  spec <- switch(
    method,
    hommola = c(stat = "Hommola_r", pvalue = "Hommola_pvalue"),
    paco = c(stat = "PACo_ss", pvalue = "PACo_pvalue"),
    parafit = c(stat = "ParaFitGlobal",
                pvalue = "ParaFit_pvalue"),
    topology = c(stat = "Topology_congruence", pvalue = "Topology_pvalue"),
    character(0)
  )
  spec[spec %in% cols]
}

# Newick and TreeDist columns are for downstream use, not for reading; hide them
# in the compact preview.
.codiv_hidden_cols <- function() {
  c("Symbiont_Tree", "Host_Tree",
    "TreeDistance", "SharedPhylogeneticInfo",
    "DifferentPhylogeneticInfo", "NyeSimilarity",
    "JaccardRobinsonFoulds", "MatchingSplitDistance",
    "MatchingSplitInfoDistance", "MutualClusteringInfo")
}

# Count nodes with a per-node permutation p-value below `alpha`. This is a
# descriptive count, not a multiple-testing-corrected result (see
# codiv_null_scans for a scan-level false-discovery assessment).
# inputs: df; method; alpha
# output: list(n_sig, total, col) or NULL when no p-value column is present
.codiv_significant <- function(df, method, alpha) {
  mc <- .codiv_method_cols(method, colnames(df))
  col <- mc[["pvalue"]]
  if (is.null(col) || is.na(col)) return(NULL)
  vals <- df[[col]]
  list(n_sig = sum(vals < alpha, na.rm = TRUE),
       total = sum(!is.na(vals)),
       col = col)
}

#' Print a codiv result
#'
#' Compact overview of a codiversification scan: the run settings, how many
#' nodes are significant per method, and a short preview of the strongest nodes
#' (Newick and tree-distance columns are hidden).
#'
#' @param x A `codiv` object from [codiv()].
#' @param alpha Significance threshold for the summary counts; default 0.05.
#' @param n Number of preview rows to show; default 6.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' \donttest{
#' sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
#' res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
#'              methods = "hommola", permutations = 99)
#' res                 # calls print.codiv
#' print(res, n = 10)  # show more preview rows
#' }
#'
#' @exportS3Method print codiv
print.codiv <- function(x, alpha = 0.05, n = 6, ...) {
  params <- attr(x, "codiv_params")
  methods <- if (!is.null(params)) params$methods else
    .codiv_methods_from_cols(colnames(x))

  cat("<codiv> codiversification scan\n")
  cat("  Nodes scanned:", nrow(x), "\n")
  if (!is.null(params)) {
    cat("  Methods:", paste(methods, collapse = ", "), "\n")
    cat("  Permutations:", params$permutations, "\n")
  }

  # per-node significance counts (raw permutation p-value; not corrected across
  # nodes -- see codiv_null_scans for a scan-level FDR)
  sig_bits <- character(0)
  for (m in methods) {
    s <- .codiv_significant(x, m, alpha)
    if (!is.null(s)) {
      sig_bits <- c(sig_bits, sprintf("%s %d/%d", m, s$n_sig, s$total))
    }
  }
  if (length(sig_bits) > 0) {
    cat("  Nodes with p < ", alpha, ": ", paste(sig_bits, collapse = ", "),
        "\n", sep = "")
  }

  if (nrow(x) == 0) {
    cat("\n(no nodes met the filtering thresholds)\n")
    return(invisible(x))
  }

  # compact preview: id/size columns plus each method's statistic and p-value
  base_cols <- intersect(c("Node_ID", "N_Symbionts", "N_Hosts"), colnames(x))
  stat_cols <- unlist(lapply(methods, function(m) {
    mc <- .codiv_method_cols(m, colnames(x))
    c(mc[["stat"]], mc[["pvalue"]])
  }), use.names = FALSE)
  preview_cols <- c(base_cols, stat_cols)
  preview_cols <- preview_cols[!preview_cols %in% .codiv_hidden_cols()]

  # order by the first available primary statistic (descending |value|)
  ord <- .codiv_preview_order(x, methods)
  preview <- as.data.frame(x)[ord, preview_cols, drop = FALSE]
  cat("\nTop nodes:\n")
  print(utils::head(preview, n), row.names = FALSE)
  if (nrow(x) > n) cat("... ", nrow(x) - n, " more nodes\n", sep = "")

  invisible(x)
}

# Row order for the preview: descending absolute value of the first method's
# primary statistic; falls back to original order.
.codiv_preview_order <- function(x, methods) {
  for (m in methods) {
    mc <- .codiv_method_cols(m, colnames(x))
    if ("stat" %in% names(mc)) {
      return(order(abs(x[[mc[["stat"]]]]), decreasing = TRUE))
    }
  }
  seq_len(nrow(x))
}

#' Summarize a codiv result
#'
#' Per-method summary of a codiversification scan: the range of the primary
#' statistic and the number of nodes with a per-node permutation p-value below
#' `alpha`. This p-value count is descriptive and not corrected across nodes;
#' use [codiv_null_scans()] for a scan-level false-discovery assessment.
#'
#' @param object A `codiv` object from [codiv()].
#' @param alpha Significance threshold; default 0.05.
#' @param ... Ignored.
#'
#' @return A data frame with one row per method (invisibly returned by print),
#'   containing the statistic range and the per-node significant count.
#'
#' @examples
#' \donttest{
#' sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
#' res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
#'              methods = c("hommola", "paco"), permutations = 99)
#' summary(res)
#' summary(res, alpha = 0.01)   # stricter significance threshold
#' }
#'
#' @exportS3Method summary codiv
summary.codiv <- function(object, alpha = 0.05, ...) {
  params <- attr(object, "codiv_params")
  methods <- if (!is.null(params)) params$methods else
    .codiv_methods_from_cols(colnames(object))

  # build one summary row from a statistic/p-value column pair
  make_row <- function(label, stat_col, pvalue_col) {
    has_stat <- !is.null(stat_col) && stat_col %in% colnames(object)
    stat <- if (has_stat) object[[stat_col]] else NA_real_
    n_raw <- if (!is.null(pvalue_col) && pvalue_col %in% colnames(object))
      sum(object[[pvalue_col]] < alpha, na.rm = TRUE) else NA_integer_
    data.frame(
      method = label,
      statistic = if (has_stat) stat_col else NA_character_,
      stat_min = suppressWarnings(min(stat, na.rm = TRUE)),
      stat_median = stats::median(stat, na.rm = TRUE),
      stat_max = suppressWarnings(max(stat, na.rm = TRUE)),
      n_sig_p = n_raw,
      n_nodes = nrow(object),
      stringsAsFactors = FALSE
    )
  }

  rows <- list()
  for (m in methods) {
    mc <- .codiv_method_cols(m, colnames(object))
    rows[[length(rows) + 1L]] <- make_row(
      m,
      if ("stat" %in% names(mc)) mc[["stat"]] else NA_character_,
      if ("pvalue" %in% names(mc)) mc[["pvalue"]] else NA_character_)
    # Hommola is run on both the full and collapsed subtree; the primary row
    # reports the full statistic, so surface the collapsed result as its own row
    # Hommola_r is the collapsed test and is the primary row. When the
    # opt-in uncollapsed test was also run, surface it as its own row so the
    # two are never confused for each other.
    if (m == "hommola" &&
        all(c("Uncollapsed_Hommola_r", "Uncollapsed_Hommola_pvalue") %in%
            colnames(object))) {
      rows[[length(rows) + 1L]] <- make_row(
        "hommola_uncollapsed", "Uncollapsed_Hommola_r",
        "Uncollapsed_Hommola_pvalue")
    }
  }
  out <- do.call(rbind, rows)
  cat("codiv summary:", nrow(object), "nodes,",
      "per-node p <", alpha, "\n\n")
  print(out, row.names = FALSE)
  invisible(out)
}

#' Coerce a codiv result to a plain data frame
#'
#' @param x A `codiv` object.
#' @param ... Ignored.
#'
#' @return A plain data frame with the `codiv` class and parameter attribute
#'   removed.
#'
#' @examples
#' \donttest{
#' sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
#' res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
#'              methods = "hommola", permutations = 99)
#' df <- as.data.frame(res)   # plain data frame for filtering / export
#' head(df)
#' }
#'
#' @exportS3Method as.data.frame codiv
as.data.frame.codiv <- function(x, ...) {
  attr(x, "codiv_params") <- NULL
  class(x) <- "data.frame"
  x
}
