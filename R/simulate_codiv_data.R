# --- internal simulator helpers -------------------------------------------

# a one-tip phylo, for grafting (pure ape; avoids a phytools dependency)
.sim_single_tip <- function(lab, el) {
  structure(list(edge = matrix(c(2L, 1L), 1, 2), tip.label = lab,
                 edge.length = el, Nnode = 1L), class = "phylo")
}

# graft a new tip as sister to `sister`
.sim_graft_sister <- function(tree, sister, newlab, frac = 0.3) {
  tip <- which(tree$tip.label == sister)
  el <- tree$edge.length[which(tree$edge[, 2] == tip)]
  if (length(el) == 0 || is.na(el)) el <- 1
  ape::bind.tree(tree, .sim_single_tip(newlab, el * frac),
                 where = tip, position = el * frac)
}

# apply `n` random subtree-prune-regraft moves to scramble topology
.sim_spr <- function(tree, n) {
  for (i in seq_len(n)) {
    if (length(tree$tip.label) < 4) break
    tip <- sample(tree$tip.label, 1)
    pruned <- ape::drop.tip(tree, tip)
    e <- sample(seq_len(nrow(pruned$edge)), 1)
    tree <- ape::bind.tree(
      pruned, .sim_single_tip(tip, pruned$edge.length[e] / 2),
      where = pruned$edge[e, 2], position = pruned$edge.length[e] / 2)
  }
  tree
}

# build one clade's symbiont subtree on a fixed host set
# output: list(tree, df, truth)
.sim_build_clade <- function(host_tree, codiversifying, congruence, n_per_host,
                             bl_noise, prefix) {
  n_hosts <- length(host_tree$tip.label)
  if (codiversifying) {
    backbone <- host_tree                       # topology tracks the host tree
    host_of <- host_tree$tip.label              # tip i -> host i
    n_moves <- round((1 - congruence) * n_hosts)
    if (n_moves > 0) {
      backbone$tip.label <- paste0(prefix, "_", seq_len(n_hosts))
      relabel <- stats::setNames(host_of, backbone$tip.label)
      backbone <- .sim_spr(backbone, n_moves)   # scramble to degrade congruence
      host_of <- unname(relabel[backbone$tip.label])
    }
  } else {
    backbone <- ape::rcoal(n_hosts)             # independent topology
    host_of <- sample(host_tree$tip.label)      # arbitrary host assignment
  }
  backbone$tip.label <- paste0(prefix, "_h", seq_len(n_hosts))
  backbone$edge.length <- backbone$edge.length *
    stats::runif(length(backbone$edge.length), 1 - bl_noise, 1 + bl_noise)

  sym <- backbone
  df <- data.frame(Host = host_of, Symbiont = sym$tip.label,
                   stringsAsFactors = FALSE)
  if (n_per_host > 1) {
    base <- sym$tip.label
    for (k in seq_len(n_per_host - 1)) {
      for (i in seq_along(base)) {
        newlab <- paste0(base[i], "_", k)
        sym <- .sim_graft_sister(sym, base[i], newlab)
        df <- rbind(df, data.frame(Host = host_of[i], Symbiont = newlab,
                                   stringsAsFactors = FALSE))
      }
    }
  }
  truth <- data.frame(Symbiont = sym$tip.label, clade = prefix,
                      codiversifying = codiversifying, congruence = congruence,
                      n_per_host = n_per_host, stringsAsFactors = FALSE)
  list(tree = sym, df = df, truth = truth)
}

# combine clade subtrees into one symbiont tree under a common root
.sim_combine <- function(clades, stem = 1) {
  strs <- vapply(clades, function(cl) {
    paste0(sub(";\\s*$", "", ape::write.tree(cl$tree)), ":", stem)
  }, character(1))
  big <- ape::read.tree(text = paste0("(", paste(strs, collapse = ","), ");"))
  big <- ape::multi2di(big)
  big$edge.length[is.na(big$edge.length) | big$edge.length == 0] <- stem
  big
}

#' Simulate host and symbiont trees with known co-diversification
#'
#' Builds a host tree and a symbiont tree assembled from several clades, each
#' with known ground truth, all mapping to the same set of hosts. Use it for a
#' fixed sanity-check or to generate many random datasets to study how tree
#' features affect the scan.
#'
#' Each clade is one monophyletic block of the symbiont tree, controlled by:
#'   - `codiversifying`: whether its topology tracks the host tree
#'   - `congruence`: for co-diversifying clades, how closely it tracks the host
#'     (1 = identical branching order; lower = more scrambled)
#'   - `n_per_host`: symbionts sampled per host (1 = no within-host
#'     pseudoreplication; > 1 introduces it)
#'
#' @param n_hosts Number of host tips (ignored if `host_tree` is supplied)
#' @param clades Optional list of clade specifications, each a list with
#'   `codiversifying` (logical) and optionally `congruence` and `n_per_host`.
#'   If NULL, `n_clades` random clades are generated.
#' @param n_clades Number of random clades to generate when `clades` is NULL
#' @param prop_codiversifying Probability a random clade is co-diversifying
#' @param congruence_range Range to draw random co-diversifying `congruence` from
#' @param n_per_host_range Integer range to draw random `n_per_host` from
#' @param bl_noise Multiplicative branch-length jitter (0 = none)
#' @param host_tree Optional host tree (phylo); if NULL a random coalescent tree
#'   of `n_hosts` tips is used
#' @param seed Optional random seed
#'
#' @return A list with:
#'   - `host_tree`: the host tree (phylo)
#'   - `symbiont_tree`: the combined symbiont tree (phylo)
#'   - `links`: data frame of Host-Symbiont associations
#'   - `truth`: per-symbiont ground truth (clade, codiversifying, congruence,
#'     n_per_host)
#'
#' @export
#'
#' @examples
#' # fixed four-clade sanity check
#' sim <- simulate_codiv_data(
#'   n_hosts = 10,
#'   clades = list(
#'     list(codiversifying = TRUE,  congruence = 1.0, n_per_host = 1),
#'     list(codiversifying = TRUE,  congruence = 0.6, n_per_host = 4),
#'     list(codiversifying = FALSE, n_per_host = 1),
#'     list(codiversifying = FALSE, n_per_host = 4)
#'   ),
#'   seed = 1
#' )
#' sim$host_tree
#' head(sim$links)
#'
#' # many random datasets (vary seed) to study tree-feature effects
#' rand <- simulate_codiv_data(n_hosts = 15, n_clades = 6, seed = 42)
#'
simulate_codiv_data <- function(n_hosts = 12, clades = NULL, n_clades = NULL,
                                prop_codiversifying = 0.5,
                                congruence_range = c(0.5, 1),
                                n_per_host_range = c(1, 4), bl_noise = 0.1,
                                host_tree = NULL, seed = NULL) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)),
              add = TRUE)
    }
    set.seed(seed)
  }

  if (is.null(host_tree)) {
    host_tree <- ape::rcoal(n_hosts)
  }
  host_tree$tip.label <- paste0("H", seq_along(host_tree$tip.label))
  n_hosts <- length(host_tree$tip.label)

  # build the clade specifications
  if (is.null(clades)) {
    if (is.null(n_clades)) {
      stop("Supply either `clades` or `n_clades`.")
    }
    clades <- lapply(seq_len(n_clades), function(i) {
      list(codiversifying = stats::runif(1) < prop_codiversifying,
           congruence = stats::runif(1, congruence_range[1],
                                     congruence_range[2]),
           n_per_host = sample(seq(n_per_host_range[1],
                                   n_per_host_range[2]), 1))
    })
  }

  # build each clade and combine
  built <- lapply(seq_along(clades), function(i) {
    spec <- clades[[i]]
    .sim_build_clade(
      host_tree,
      codiversifying = isTRUE(spec$codiversifying),
      congruence = if (is.null(spec$congruence)) 1 else spec$congruence,
      n_per_host = if (is.null(spec$n_per_host)) 1L else spec$n_per_host,
      bl_noise = bl_noise,
      prefix = paste0("c", i))
  })

  list(
    host_tree = host_tree,
    symbiont_tree = .sim_combine(built),
    links = do.call(rbind, lapply(built, function(b) b$df)),
    truth = do.call(rbind, lapply(built, function(b) b$truth))
  )
}
