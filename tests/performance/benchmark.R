#!/usr/bin/env Rscript
#
# Performance benchmark for codiv() function
# Tests on trees of increasing size to identify scaling characteristics
#
# Usage: Rscript benchmark.R
#

library(codiv)
library(ape)
library(tibble)

# Create benchmark scenarios
benchmark_scenarios <- list(
  list(name = "tiny", n_hosts = 3, n_symbionts = 20),
  list(name = "small", n_hosts = 5, n_symbionts = 50),
  list(name = "medium", n_hosts = 10, n_symbionts = 200),
  list(name = "large", n_hosts = 20, n_symbionts = 500)
)

results <- list()

for (scenario in benchmark_scenarios) {
  cat("\n", "="*60, "\n")
  cat("Scenario:", scenario$name, "\n")
  cat("Hosts:", scenario$n_hosts, "| Symbionts:", scenario$n_symbionts, "\n")
  cat("="*60, "\n")

  # Generate test data
  h_tree <- ape::rtree(scenario$n_hosts)
  s_tree <- ape::rtree(scenario$n_symbionts)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, length.out = scenario$n_symbionts),
    Symbiont = s_tree$tip.label
  )

  # Benchmark with different parameter sets
  benchmarks <- list(
    list(
      name = "Quick scan (permutations=9, span=0.2, min_tips=5)",
      permutations = 9,
      span_fraction = 0.2,
      min_symbiont_tips = 5,
      methods = "hommola"
    ),
    list(
      name = "Standard scan (permutations=99, span=0.1, min_tips=7)",
      permutations = 99,
      span_fraction = 0.1,
      min_symbiont_tips = 7,
      methods = c("hommola")
    )
  )

  for (bench in benchmarks) {
    cat("\nBenchmark:", bench$name, "\n")

    temp_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(temp_file))

    # Time the codiv() call
    start_time <- Sys.time()
    result <- codiv(
      Host_tree = h_tree,
      Symbiont_tree = s_tree,
      Host_to_Symbiont_df = hs_df,
      min_hosts = 3,
      min_symbiont_tips = bench$min_symbiont_tips,
      span_fraction = bench$span_fraction,
      permutations = bench$permutations,
      seed = 8675309,
      verbose = FALSE,
      Save_fp = temp_file,
      methods = bench$methods,
      subtree_features = FALSE
    )
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "secs")

    n_nodes <- nrow(result)
    time_per_node <- as.numeric(elapsed) / max(n_nodes, 1)

    cat("  Nodes scanned:", n_nodes, "\n")
    cat("  Total time:", format(elapsed), "\n")
    cat("  Time per node:", round(time_per_node, 3), "seconds\n")

    results[[paste(scenario$name, bench$name, sep = " | ")]] <- list(
      scenario = scenario$name,
      benchmark = bench$name,
      n_hosts = scenario$n_hosts,
      n_symbionts = scenario$n_symbionts,
      n_nodes_scanned = n_nodes,
      total_time_sec = as.numeric(elapsed),
      time_per_node_sec = time_per_node
    )
  }
}

# Create summary
cat("\n\n", "="*80, "\n")
cat("PERFORMANCE SUMMARY\n")
cat("="*80, "\n\n")

summary_df <- as.data.frame(do.call(rbind, results))
summary_df <- tibble::as_tibble(summary_df)

print(summary_df)

# Save summary to file
write.csv(summary_df, "benchmark_results.csv", row.names = FALSE)
cat("\n✓ Results saved to benchmark_results.csv\n")
