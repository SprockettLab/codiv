
#' Convert host-symbiont associations to binary presence-absence matrix
#'
#' Transforms a long-format data frame of host-symbiont associations into
#' a wide-format binary matrix suitable for co-phylogenetic analysis methods
#' like PACo and ParaFit.
#'
#' @param Host_to_Symbiont_df A data frame with columns "Host" and "Symbiont"
#'   where each row represents one association. Symbionts must be unique
#'   (one symbiont per host per row).
#'
#' @return A binary matrix with hosts as rows and symbionts as columns, filled
#'   with 0s and 1s where 1 indicates the symbiont was isolated from that host.
#'
#' @export
#'
#' @examples
#' hs_df <- data.frame(
#'   Host = c("H1", "H1", "H2", "H2", "H3"),
#'   Symbiont = c("S1", "S2", "S3", "S4", "S5")
#' )
#' mat <- host_symbiont_links(hs_df)
#'
host_symbiont_links <- function(Host_to_Symbiont_df) {
  if (sum(colnames(Host_to_Symbiont_df) %in% c("Host", "Symbiont")) != 2) {
    stop("`Host_to_Symbiont_df` must include both 'Host' and 'Symbiont' columns")
  }

  if (anyDuplicated(Host_to_Symbiont_df$Symbiont)) {
    warning("One or more 'Symbiont' values are duplicated in Host_to_Symbiont_df; ",
            "each symbiont should map to a single host.")
  }

  hosts <- unique(Host_to_Symbiont_df$Host)
  symbionts <- unique(Host_to_Symbiont_df$Symbiont)

  # fill a hosts-by-symbionts matrix in one vectorized indexing step
  host_symb_mat <- matrix(0L, nrow = length(hosts), ncol = length(symbionts),
                          dimnames = list(hosts, symbionts))
  row_idx <- match(Host_to_Symbiont_df$Host, hosts)
  col_idx <- match(Host_to_Symbiont_df$Symbiont, symbionts)
  host_symb_mat[cbind(row_idx, col_idx)] <- 1L

  host_symb_mat
}

