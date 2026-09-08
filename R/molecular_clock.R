#' Molecular-clock corroboration of co-diversification
#'
#' For the most strongly co-diversifying clades, regresses symbiont divergence
#' against host clade age. A positive relationship corroborates that symbiont
#' clades diversified contemporaneously with their hosts, and the slope is an
#' estimate of the symbiont molecular rate (substitutions per unit host time).
#' This follows the approach of Sanders et al. (2023) and Sprockett et al.
#' (2025).
#'
#' @param codiv_results A `codiv` result (or data frame) from [codiv()]
#' @param host_time_tree A time-calibrated host tree (phylo, ultrametric) whose
#'   branch lengths are in time units (e.g. millions of years)
#' @param statistic Name of the test-statistic column. NULL (default) uses
#'   `Hommola_r`, the collapsed Hommola correlation that [codiv()] always
#'   writes when the Hommola method runs
#' @param stat_threshold Only clades with a statistic above this value are used;
#'   default 0.95 (strong co-diversification)
#' @param pvalue Name of the p-value column. NULL (default) resolves
#'   alongside `statistic`
#' @param p_threshold Only clades with a p-value below this value are used;
#'   default 0.01
#' @param symbiont_divergence Optional override for the symbiont-divergence
#'   summary. Either a named numeric vector keyed by `Node_ID`, or a function
#'   taking a symbiont subtree (phylo) and returning a single number. If NULL
#'   (default), the clade's crown depth is used (mean root-to-tip distance in
#'   the symbiont tree's branch-length units)
#'
#' @return A list with:
#'   - `model`: the fitted `lm` of symbiont divergence on host age
#'   - `slope`, `intercept`, `r_squared`, `p_value`: regression summaries
#'   - `slope_ci`: 95% confidence interval for the slope
#'   - `per_clade`: data frame with one row per clade used (`Node_ID`,
#'     `n_hosts`, `host_age`, `symbiont_divergence`)
#'
#' @details
#' Host clade age is the crown age of the clade's hosts (the age of their most
#' recent common ancestor in `host_time_tree`). The default symbiont divergence
#' is a simplification: it is the clade crown depth in substitution units, not a
#' dedicated absolute date. For rigorous absolute dating supply your own
#' estimates via `symbiont_divergence`. Nested clades are non-independent; the
#' regression treats each passing node as a point, so consider restricting to
#' non-nested clades for formal inference.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sim <- simulate_codiv_data(n_hosts = 14, clades = list(
#'   list(codiversifying = TRUE, congruence = 1, n_per_host = 3)), seed = 1)
#' res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, span_fraction = 0.9,
#'              methods = "hommola", permutations = 199, verbose = FALSE)
#' # sim$host_tree is ultrametric, so it stands in for a time-calibrated tree
#' mc <- molecular_clock(res, sim$host_tree, stat_threshold = 0.75)
#' mc$slope
#' mc$per_clade
#' }
#'
molecular_clock <- function(codiv_results, host_time_tree,
                            statistic = NULL, stat_threshold = 0.95,
                            pvalue = NULL, p_threshold = 0.01,
                            symbiont_divergence = NULL) {
  df <- as.data.frame(codiv_results)
  # NULL means "the Hommola statistic", resolved against what the result
  # actually carries, since the uncollapsed columns are opt-in.
  .hc <- .default_hommola_cols(df)
  if (is.null(statistic)) statistic <- unname(.hc[["stat"]])
  if (is.null(pvalue)) pvalue <- unname(.hc[["pvalue"]])
  if (!all(c(statistic, pvalue, "Host_Tree", "Symbiont_Tree", "Node_ID")
           %in% colnames(df))) {
    stop("Result is missing a required column (", statistic, ", ", pvalue,
         ", Host_Tree, Symbiont_Tree, or Node_ID).")
  }
  if (!inherits(host_time_tree, "phylo")) {
    stop("`host_time_tree` must be a phylo object")
  }
  if (!ape::is.ultrametric(host_time_tree, tol = 1e-4)) {
    warning("`host_time_tree` is not ultrametric; host ages assume a ",
            "time-calibrated (ultrametric) tree.")
  }

  # host node ages (time before present) for the calibrated tree
  ages <- ape::branching.times(host_time_tree)

  # default symbiont divergence = clade crown depth (mean root-to-tip distance)
  crown_depth <- function(subtree) {
    n_tip <- length(subtree$tip.label)
    if (is.null(subtree$edge.length)) return(NA_real_)
    mean(ape::node.depth.edgelength(subtree)[seq_len(n_tip)])
  }

  keep <- df[[statistic]] > stat_threshold & df[[pvalue]] < p_threshold
  keep[is.na(keep)] <- FALSE
  rows <- which(keep)
  if (length(rows) < 3) {
    stop("Fewer than 3 clades pass the thresholds; cannot fit a regression. ",
         "Loosen `stat_threshold`/`p_threshold`.")
  }

  per <- lapply(rows, function(i) {
    hosts <- tryCatch(ape::read.tree(text = df$Host_Tree[i])$tip.label,
                      error = function(e) character(0))
    hosts <- intersect(hosts, host_time_tree$tip.label)
    if (length(hosts) < 2) return(NULL)
    mrca <- ape::getMRCA(host_time_tree, hosts)
    host_age <- unname(ages[as.character(mrca)])

    sub <- tryCatch(ape::read.tree(text = df$Symbiont_Tree[i]),
                    error = function(e) NULL)
    if (is.null(sub)) return(NULL)
    div <- if (is.null(symbiont_divergence)) {
      crown_depth(sub)
    } else if (is.function(symbiont_divergence)) {
      symbiont_divergence(sub)
    } else {
      unname(symbiont_divergence[as.character(df$Node_ID[i])])
    }
    data.frame(Node_ID = df$Node_ID[i], n_hosts = length(hosts),
               host_age = host_age, symbiont_divergence = div,
               stringsAsFactors = FALSE)
  })
  per_clade <- do.call(rbind, per)
  per_clade <- per_clade[stats::complete.cases(per_clade), ]
  if (nrow(per_clade) < 3) {
    stop("Fewer than 3 clades had a resolvable host age and symbiont ",
         "divergence; cannot fit a regression.")
  }

  fit <- stats::lm(symbiont_divergence ~ host_age, data = per_clade)
  s <- summary(fit)
  list(
    model = fit,
    slope = unname(stats::coef(fit)[2]),
    intercept = unname(stats::coef(fit)[1]),
    r_squared = s$r.squared,
    p_value = stats::pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                        lower.tail = FALSE),
    slope_ci = stats::confint(fit)[2, ],
    per_clade = per_clade
  )
}
