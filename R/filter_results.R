# Identify which calls are nested inside which others.
# Clades of a single rooted tree form a laminar family: any two are either
# disjoint or one contains the other. So when clade B is smaller than clade A,
# a single shared tip is enough to prove containment, and no full set
# comparison is needed.
# inputs: tip_sets (list of character vectors, one per call)
# output: integer vector giving, for each call, the index of its immediate
#   parent call, or NA when nothing contains it
.nesting_parents <- function(tip_sets) {
  n <- length(tip_sets)
  if (n == 0) return(integer(0))
  sizes <- lengths(tip_sets)
  first_tip <- vapply(tip_sets, `[`, character(1), 1)

  # hashed membership per call, so the containment test is a single lookup
  members <- lapply(tip_sets, function(x) {
    e <- new.env(hash = TRUE, size = max(2L, length(x) * 2L))
    for (k in x) assign(k, TRUE, envir = e)
    e
  })

  # visiting candidates smallest first means the first strict container found
  # is the immediate parent
  parent <- rep(NA_integer_, n)
  ord <- order(sizes)
  for (j in seq_len(n)) {
    for (i in ord) {
      if (sizes[i] <= sizes[j]) next
      if (exists(first_tip[j], envir = members[[i]], inherits = FALSE)) {
        parent[j] <- i
        break
      }
    }
  }
  parent
}

# Walk each call up to the outermost call that contains it, which labels the
# nested family a call belongs to.
# inputs: parent (integer vector from .nesting_parents())
# output: integer vector of family roots, one per call
.nesting_root <- function(parent) {
  vapply(seq_along(parent), function(i) {
    cur <- i
    while (!is.na(parent[cur])) cur <- parent[cur]
    cur
  }, integer(1))
}

# Read a sort order out of the filtering expressions.
# A comparison says which direction is "better": `p < 0.01` wants the smallest
# p first, `r > 0.5` wants the largest r first. Only simple
# `column op constant` comparisons carry that meaning, so anything else is
# skipped.
# inputs: conditions (list of unevaluated expressions), nms (column names)
# output: list with `columns` (character) and `decreasing` (logical), in the
#   order the conditions were given, with repeated columns kept once
.sort_keys_from_conditions <- function(conditions, nms) {
  ops_up <- c(">", ">=")     # larger is better
  ops_dn <- c("<", "<=")     # smaller is better
  cols <- character(0)
  dec <- logical(0)
  for (cond in conditions) {
    if (!is.call(cond) || length(cond) != 3L) next
    op <- as.character(cond[[1]])
    if (!op %in% c(ops_up, ops_dn)) next
    lhs <- cond[[2]]
    rhs <- cond[[3]]
    if (is.name(lhs) && as.character(lhs) %in% nms && is.numeric(tryCatch(eval(rhs), error = function(e) NULL))) {
      col <- as.character(lhs)
      decreasing <- op %in% ops_up
    } else if (is.name(rhs) && as.character(rhs) %in% nms && is.numeric(tryCatch(eval(lhs), error = function(e) NULL))) {
      # the comparison is written backwards, so the direction flips
      col <- as.character(rhs)
      decreasing <- op %in% ops_dn
    } else next
    if (col %in% cols) next
    cols <- c(cols, col)
    dec <- c(dec, decreasing)
  }
  list(columns = cols, decreasing = dec)
}

#' Filter a codiversification scan, optionally collapsing nested calls
#'
#' Subsets the output of [codiv()] on any combination of its columns, on which
#' hosts a clade contains, and on nesting. Scanned nodes are nested by
#' construction, so a single codiversifying lineage is usually reported many
#' times over: once for the whole clade and once for most of its subclades.
#' Counting those rows directly overstates how many distinct lineages were
#' found. The `nested` argument reduces each nested series to one call.
#'
#' @param codiv_results A `codiv` result (or data frame) from [codiv()]
#' @param ... Filtering expressions evaluated inside `codiv_results`, in the
#'   style of [subset()], for example `Hommola_pvalue < 0.05` and
#'   `Hommola_r > 0.7`. Multiple expressions are combined with AND.
#'   Rows where an expression is NA are dropped, matching [subset()].
#' @param hosts Optional character vector of host names to require. [codiv()]
#'   writes a `<host>_PRESENT` column for every name passed to its
#'   `focus_hosts` argument, and those columns are what this reads, so a host
#'   can only be filtered on if it was tracked during the scan.
#' @param hosts_match Whether a clade must contain all of `hosts` ("all", the
#'   default) or at least one of them ("any")
#' @param nested How to treat nested calls, applied after every other filter:
#'   * `"keep_all"` (default) leaves them alone, so the result still counts one
#'     lineage many times
#'   * `"outermost"` keeps the largest call in each nested series
#'   * `"innermost"` keeps the calls that contain no other surviving call.
#'     A call can hold two disjoint subclades, so this can return more rows
#'     than `"outermost"`, not fewer
#'   * `"best"` keeps whichever call in each series has the highest
#'     `statistic`. The strongest signal often sits inside a series rather
#'     than at its outermost call, so this frequently differs from
#'     `"outermost"`
#' @param statistic Column used to choose a winner when `nested = "best"`;
#'   default "Hommola_r"
#' @param sort_by_filters If TRUE, order the surviving rows by the columns used
#'   in `...`, in the order those conditions were written, taking the direction
#'   from each comparison: `Hommola_pvalue < 0.01` sorts that column ascending,
#'   `Hommola_r > 0.5` sorts it descending. Only simple `column op constant`
#'   comparisons give a direction, so compound conditions and `==` are ignored
#'   for sorting. Default FALSE, which leaves the scan's own order alone.
#'
#' @return A `codiv` object holding the surviving rows, carrying the same
#'   `codiv_params` attribute as the input plus a `codiv_filter` attribute
#'   recording what was applied.
#'
#' @details
#' Filtering always runs before de-nesting. The other order would pick a
#' representative for a series and then possibly discard it, silently losing
#' the whole series even when other members passed the filter.
#'
#' De-nesting reads the per-node symbiont trees stored by [codiv()], so
#' `Symbiont_Tree` must still be present.
#'
#' @seealso [codiv()], [host_codiv_summary()]
#' @export
#'
#' @examples
#' \donttest{
#' sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
#' res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
#'              methods = "hommola", permutations = 99, verbose = FALSE)
#'
#' # every significant call, including nested duplicates
#' filter_results(res, Hommola_pvalue < 0.05)
#'
#' # one row per distinct lineage
#' filter_results(res, Hommola_pvalue < 0.05, nested = "outermost")
#' }
#'
filter_results <- function(codiv_results, ..., hosts = NULL,
                           hosts_match = c("all", "any"),
                           nested = c("keep_all", "outermost", "innermost", "best"),
                           statistic = "Hommola_r",
                           sort_by_filters = FALSE) {
  hosts_match <- match.arg(hosts_match)
  nested <- match.arg(nested)

  df <- as.data.frame(codiv_results)
  params <- attr(codiv_results, "codiv_params")
  keep <- rep(TRUE, nrow(df))

  # column expressions, evaluated the way subset() does: NA drops the row
  conditions <- as.list(substitute(list(...)))[-1L]
  for (cond in conditions) {
    val <- eval(cond, df, parent.frame())
    if (!is.logical(val)) {
      stop("Filtering expression `", deparse(cond),
           "` must give a logical vector, not ", class(val)[1])
    }
    keep <- keep & !is.na(val) & val
  }

  # host presence, read off the <host>_PRESENT columns codiv() writes
  if (!is.null(hosts)) {
    wanted <- paste0(hosts, "_PRESENT")
    missing_cols <- wanted[!wanted %in% colnames(df)]
    if (length(missing_cols) > 0) {
      available <- sub("_PRESENT$", "", grep("_PRESENT$", colnames(df),
                                             value = TRUE))
      stop("No presence column for host(s): ",
           paste(sub("_PRESENT$", "", missing_cols), collapse = ", "),
           ". Hosts must be passed to codiv(focus_hosts = ) during the scan. ",
           if (length(available) > 0) {
             paste0("Tracked hosts are: ", paste(available, collapse = ", "))
           } else "No hosts were tracked in this scan.")
    }
    present_mat <- as.matrix(df[, wanted, drop = FALSE]) == TRUE
    host_ok <- if (hosts_match == "all") {
      rowSums(present_mat, na.rm = TRUE) == length(wanted)
    } else {
      rowSums(present_mat, na.rm = TRUE) > 0
    }
    keep <- keep & host_ok
  }

  df <- df[keep, , drop = FALSE]

  # de-nesting runs last, over the calls that survived every other filter
  if (nested != "keep_all" && nrow(df) > 1) {
    if (!"Symbiont_Tree" %in% colnames(df)) {
      stop("De-nesting needs the 'Symbiont_Tree' column, which is not present")
    }
    tip_sets <- lapply(df$Symbiont_Tree,
                       function(x) ape::read.tree(text = x)$tip.label)
    parent <- .nesting_parents(tip_sets)

    selected <- switch(
      nested,
      outermost = is.na(parent),
      innermost = !(seq_len(nrow(df)) %in% parent[!is.na(parent)]),
      best = {
        if (!statistic %in% colnames(df)) {
          stop("`statistic` column '", statistic, "' is not in the results")
        }
        root <- .nesting_root(parent)
        stat <- df[[statistic]]
        winner <- vapply(split(seq_len(nrow(df)), root), function(ix) {
          s <- stat[ix]
          # a family whose statistic is entirely NA falls back to its outermost
          if (all(is.na(s))) ix[which.max(lengths(tip_sets[ix]))]
          else ix[which.max(replace(s, is.na(s), -Inf))]
        }, integer(1))
        seq_len(nrow(df)) %in% winner
      }
    )
    df <- df[selected, , drop = FALSE]
  }

  # sorting runs last, so it orders the rows that actually survived
  sort_keys <- list(columns = character(0), decreasing = logical(0))
  if (sort_by_filters && nrow(df) > 1) {
    sort_keys <- .sort_keys_from_conditions(conditions, colnames(df))
    if (length(sort_keys$columns) == 0) {
      warning("`sort_by_filters` is TRUE but none of the conditions is a ",
              "simple `column op constant` comparison; order left unchanged")
    } else {
      df <- df[do.call(order, c(unname(as.list(df[, sort_keys$columns,
                                                  drop = FALSE])),
                                list(decreasing = sort_keys$decreasing,
                                     method = "radix"))), , drop = FALSE]
    }
  }

  attr(df, "codiv_params") <- params
  attr(df, "codiv_filter") <- list(
    conditions = vapply(conditions, function(x) paste(deparse(x), collapse = " "),
                      character(1)),
    hosts = hosts, hosts_match = hosts_match,
    nested = nested, statistic = if (nested == "best") statistic else NA,
    sorted_by = if (length(sort_keys$columns) > 0) {
      paste0(sort_keys$columns, ifelse(sort_keys$decreasing, " (desc)", " (asc)"))
    } else NA
  )
  class(df) <- unique(c("codiv", class(df)))
  df
}
