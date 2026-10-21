#Group members:
#Ruiqi Wang s2930743
#Zirui Guo s2948295
#Tong Shen s2939003
#
#Contributions:
#Tong Shen (33.3%): Implemented Tasks 1-3.
#Ruiqi Wang (33.3%): Implemented Task 4-5.
#Zirui Guo (33.3%): Implemented Task 6-7.
#All members contributed to the GitHub repository setup, code review, and debugging.
#
# Sequence duels and the perils of independence
#
# Heads is 1 and tails is 0. The chance of heads is p.
# Each toss is independent of the other tosses.
# Two players watch the same sequence and choose different patterns.
# Their patterns have the same length. The first completed pattern wins.
#
# This program finds patterns, simulates duels, and calculates waiting times.
# It compares all patterns of length three and finds best replies and cycles.
# It also makes two plots and repeats the experiments for biased coins.
#
# Only base R is needed. All functions are defined before the experiments.
# The plots are saved in sequence_duels.pdf in the current working directory.

# In several runuings, total running time is around 7s.

check_pattern <- function(a) {
  # Purpose: check that a pattern is a non-empty binary vector.
  # Inputs: a, the pattern to check.
  # Output: invisible NULL, or an error message.
  # Approach: check the type, length, missing values, and allowed values.
  
  if (!(is.numeric(a) || is.logical(a)) || !is.null(dim(a)) ||
      length(a) == 0 || anyNA(a) || any(!(a %in% c(0, 1)))) {
    stop("A pattern must be a non-empty vector containing only 0 and 1.")
  }
  invisible(NULL)
}


firsthit <- function(x, a) {
  # Purpose: find the first complete match of a in x.
  # Inputs: x, a binary vector which may be empty; a, a non-empty pattern.
  # Output: the end position of the first match, or NA_integer_.
  # Approach: move a window through x and stop at its first full match.
  
  check_pattern(a)
  if (!(is.numeric(x) || is.logical(x)) || !is.null(dim(x)) ||
      anyNA(x) || any(!(x %in% c(0, 1)))) {
    stop("x must be a vector containing only 0 and 1.")
  }
  k <- length(a)
  
  # A short sequence cannot contain a full pattern.
  if (length(x) < k) return(NA_integer_)
  
  # last is the end of the window, not its start.
  for (last in seq.int(k, length(x))) {
    if (all(x[(last - k + 1):last] == a)) return(last)
  }
  NA_integer_
}


check_simulation <- function(p, nreps) {
  # Purpose: check the probability and the number of experiments.
  # Inputs: p, the heads probability; nreps, the number of experiments.
  # Output: invisible NULL, or an error message.
  # Approach: allow a finite probability and a positive integer count.
  
  if (!is.numeric(p) || length(p) != 1 || !is.finite(p) ||
      p <= 0 || p >= 1) {
    stop("p must be a single probability strictly between 0 and 1.")
  }
  if (!is.numeric(nreps) || length(nreps) != 1 || !is.finite(nreps) ||
      nreps < 1 || nreps != floor(nreps)) {
    stop("nreps must be a positive integer.")
  }
  invisible(NULL)
}


simulate_hits <- function(a, b = NULL, p = 0.5, nreps = 1) {
  # Purpose: simulate pattern matches for many independent experiments.
  # Inputs: pattern a, optional pattern b, probability p, and nreps.
  # Output: a matrix with columns winner and tosses.
  # Approach: generate blocks of tosses and continue unfinished rows.
  #
  # Each row is one experiment. Both patterns use the same row.
  # Winner is 1 for a and 0 for b.
  # With no b, the function waits only for a and winner is always 1.
  
  check_pattern(a)
  check_simulation(p, nreps)
  if (!is.null(b)) {
    check_pattern(b)
    if (length(a) != length(b) || all(a == b)) {
      stop("Duel patterns must be distinct and have the same length.")
    }
  }
  
  k <- length(a)
  block_size <- 32L
  result <- matrix(NA_integer_, nrow = nreps, ncol = 2,
                   dimnames = list(NULL, c("winner", "tosses")))
  
  # Keep the original row numbers when finished rows are removed.
  active <- seq_len(nreps)
  
  # The first history is padding, not real tosses.
  # When k is 1, this matrix has zero columns.
  history <- matrix(0L, nrow = nreps, ncol = k - 1)
  elapsed <- 0L
  
  # There is no fixed total length.
  # Each experiment continues until a pattern really appears.
  while (length(active) > 0) {
    tosses <- matrix(rbinom(length(active) * block_size, size = 1, prob = p),
                     nrow = length(active), ncol = block_size)
    
    # A pattern may start in an old block and finish in the new block.
    # Keeping k-1 previous tosses allows these matches.
    buffer <- cbind(history, tosses)
    match_a <- matrix(TRUE, nrow = length(active), ncol = block_size)
    match_b <- matrix(!is.null(b), nrow = length(active), ncol = block_size)
    
    # Each column describes a possible pattern ending in the new block.
    # All positions must match before a full window is accepted.
    for (j in seq_len(k)) {
      window_toss <- buffer[, j:(j + block_size - 1), drop = FALSE]
      match_a <- match_a & window_toss == a[j]
      if (!is.null(b)) match_b <- match_b & window_toss == b[j]
    }
    matches <- match_a | match_b
    
    # Do not count a window that still uses the initial padding.
    # This also works when the pattern is longer than one block.
    incomplete <- min(block_size, k - 1 - elapsed)
    if (incomplete > 0) matches[, seq_len(incomplete)] <- FALSE
    
    # TRUE is larger than FALSE, so its first maximum is the first hit.
    # A row with no TRUE still has a maximum. Check the selected value.
    first <- max.col(matches, ties.method = "first")
    finished <- matches[cbind(seq_len(length(active)), first)]
    hit <- which(finished)
    if (length(hit) > 0) {
      result[active[hit], "winner"] <-
        as.integer(match_a[cbind(hit, first[hit])])
      result[active[hit], "tosses"] <- elapsed + first[hit]
    }
    
    # Keep the last k-1 values only for unfinished experiments.
    # drop = FALSE keeps a matrix when only one row remains.
    # seq_len(0) gives zero history columns when k is 1.
    history <- buffer[!finished, block_size + seq_len(k - 1), drop = FALSE]
    active <- active[!finished]
    elapsed <- elapsed + block_size
  }
  result
}


duel <- function(a, b, p = 0.5) {
  # Purpose: run one duel between two patterns.
  # Inputs: different equal-length patterns a and b, and probability p.
  # Output: a named vector containing winner and tosses.
  # Approach: take one row from the shared simulation function.
  
  simulate_hits(a, b, p, nreps = 1)[1, ]
}


pwin <- function(a, b, p = 0.5, nreps = 10000) {
  # Purpose: estimate the chance that a wins against b.
  # Inputs: different equal-length patterns a and b, p, and nreps.
  # Output: the estimated probability that a wins.
  # Approach: average the winner values from the shared simulator.
  
  results <- simulate_hits(a, b, p, nreps)
  mean(results[, "winner"])
}


ewait <- function(a, p = 0.5) {
  # Purpose: calculate the exact mean waiting time for one pattern.
  # Inputs: pattern a and heads probability p.
  # Output: the expected number of tosses until a first appears.
  # Approach: compare each prefix with the suffix of the same length.
  # For every match, add the inverse of the prefix probability.
  # This function uses no random numbers.
  
  check_pattern(a)
  check_simulation(p, 1)
  k <- length(a)
  expected <- 0
  
  # An overlap allows the end of one copy to help start another copy.
  for (j in seq_len(k)) {
    prefix <- a[seq_len(j)]
    suffix <- a[(k - j + 1):k]
    if (all(prefix == suffix)) {
      heads <- sum(prefix)
      probability <- p^heads * (1 - p)^(j - heads)
      expected <- expected + 1 / probability
    }
  }
  
  # The full pattern always matches itself, so j = k is included.
  expected
}


all_patterns <- function(k = 3) {
  # Purpose: make all binary patterns of length k.
  # Inputs: k, a positive integer.
  # Output: a matrix with one pattern per row and H/T row names.
  # Approach: read binary digits from the integers 0 to 2^k - 1.
  #
  # For k = 3, the row order is:
  # TTT, TTH, THT, THH, HTT, HTH, HHT, HHH.
  
  if (!is.numeric(k) || length(k) != 1 || !is.finite(k) ||
      k < 1 || k != floor(k)) {
    stop("k must be a positive integer.")
  }
  patterns <- matrix(0L, nrow = 2^k, ncol = k)
  ids <- 0:(2^k - 1)
  labels <- rep("", nrow(patterns))
  
  # Read the largest binary place first to keep the required order.
  for (j in seq_len(k)) {
    patterns[, j] <- (ids %/% 2^(k - j)) %% 2
    labels <- paste0(labels, ifelse(patterns[, j] == 1, "H", "T"))
  }
  rownames(patterns) <- labels
  patterns
}


pmat <- function(k = 3, p = 0.5, nreps = 10000) {
  # Purpose: estimate the win probability for every pair of patterns.
  # Inputs: pattern length k, heads probability p, and nreps per pair.
  # Output: a win matrix P with H/T row names and column names.
  # Approach: simulate each pair once and use the opposite result.
  #
  # P[i, j] is the chance that row pattern i beats column pattern j.
  # A pattern does not play itself, so the diagonal stays NA.
  
  check_simulation(p, nreps)
  patterns <- all_patterns(k)
  n <- nrow(patterns)
  P <- matrix(NA_real_, nrow = n, ncol = n,
              dimnames = list(rownames(patterns), rownames(patterns)))
  
  # Different equal-length patterns cannot finish on the same toss.
  # Their two win probabilities therefore add to one.
  for (i in seq_len(n - 1)) {
    for (j in seq.int(i + 1, n)) {
      P[i, j] <- pwin(patterns[i, ], patterns[j, ], p, nreps)
      P[j, i] <- 1 - P[i, j]
    }
  }
  P
}


best_replies <- function(P) {
  # Purpose: find a best reply to each opponent.
  # Inputs: a win probability matrix P.
  # Output: a table of opponents, replies, and estimated probabilities.
  # Approach: find the largest probability in each opponent's column.
  #
  # which.max ignores the diagonal NA.
  # If estimates tie, the first row with that value is selected.
  
  reply <- integer(ncol(P))
  for (j in seq_len(ncol(P))) reply[j] <- which.max(P[, j])
  data.frame(
    opponent = colnames(P),
    best_reply = rownames(P)[reply],
    win_probability = round(P[cbind(reply, seq_len(ncol(P)))], 4),
    row.names = NULL
  )
}


balanced_pairs <- function(P, tolerance = 0.02) {
  # Purpose: find pairs with estimated chances close to one half.
  # Inputs: win matrix P and the allowed distance from 0.5.
  # Output: a table containing each selected pair once.
  # Approach: check only the upper triangle of P.
  #
  # This rule finds approximate balance in the simulation.
  # It does not prove that the exact probability is 0.5.
  
  pairs <- which(upper.tri(P) & abs(P - 0.5) <= tolerance, arr.ind = TRUE)
  data.frame(
    pattern_a = rownames(P)[pairs[, 1]],
    pattern_b = colnames(P)[pairs[, 2]],
    probability_a_wins = round(P[pairs], 4),
    row.names = NULL
  )
}


find_cycle <- function(P, margin = 0.05) {
  # Purpose: find a shortest cycle of clear wins.
  # Inputs: win matrix P and a margin above 0.5.
  # Output: node numbers for a closed cycle, or integer(0).
  # Approach: use breadth-first search from every start node.
  #
  # An arrow goes from a winner to the opponent it beats.
  # The default rule needs an estimated probability above 0.55.
  # The cycle length is found by the search, not chosen in advance.
  
  n <- nrow(P)
  shortest <- integer(0)
  for (start in seq_len(n)) {
    queue <- start
    
    # An NA parent means that the node has not been visited.
    # The start is already visited and must not enter the queue again.
    parent <- rep(NA_integer_, n)
    parent[start] <- 0L
    position <- 1L
    
    while (position <= length(queue)) {
      current <- queue[position]
      position <- position + 1L
      neighbours <- which(P[current, ] > 0.5 + margin)
      
      # The queue visits shorter paths first.
      # An edge back to start closes a shortest cycle for this start.
      if (current != start && start %in% neighbours) {
        path <- current
        while (path[1] != start) path <- c(parent[path[1]], path)
        path <- c(path, start)
        if (length(shortest) == 0 || length(path) < length(shortest)) {
          shortest <- path
        }
        break
      }
      
      # Each node gets a parent only once.
      # Internal nodes in the final path therefore do not repeat.
      for (next_node in neighbours) {
        if (is.na(parent[next_node])) {
          parent[next_node] <- current
          queue <- c(queue, next_node)
        }
      }
    }
  }
  shortest
}


show_cycle <- function(P, margin = 0.05) {
  # Purpose: print a shortest cycle and its edge probabilities.
  # Inputs: win matrix P and a margin above 0.5.
  # Output: printed results and invisible NULL.
  # Approach: call find_cycle and read the selected edges from P.
  
  cycle <- find_cycle(P, margin)
  if (length(cycle) == 0) {
    cat("No cycle detected with every win probability above", 0.5 + margin, "\n")
    return(invisible(NULL))
  }
  cat("Shortest detected cycle:",
      paste(rownames(P)[cycle], collapse = " -> "), "\n")
  
  # The last node repeats the first node to close the path.
  from <- head(cycle, -1)
  to <- tail(cycle, -1)
  print(data.frame(
    winner = rownames(P)[from], opponent = colnames(P)[to],
    win_probability = round(P[cbind(from, to)], 4), row.names = NULL
  ))
  invisible(NULL)
}


same_wait_example <- function(P, waits) {
  # Purpose: find equal waiting times with different duel chances.
  # Inputs: win matrix P and exact waiting times in matching row order.
  # Output: one comparison row, or NULL if no equal waits exist.
  # Approach: find equal waits and select the largest win gap.
  #
  # A small tolerance allows for numerical rounding.
  
  difference <- abs(outer(waits, waits, "-"))
  pairs <- which(upper.tri(P) & difference < 1e-10, arr.ind = TRUE)
  if (nrow(pairs) == 0) return(NULL)
  pair <- pairs[which.max(abs(P[pairs] - 0.5)), ]
  
  # Put the estimated winner first in the table.
  if (P[pair[1], pair[2]] < 0.5) pair <- rev(pair)
  data.frame(
    winner = rownames(P)[pair[1]], opponent = colnames(P)[pair[2]],
    expected_wait_both = waits[pair[1]],
    win_probability = round(P[pair[1], pair[2]], 4), row.names = NULL
  )
}


coin_summary <- function(P, p, reference) {
  # Purpose: compare one coin with the fair-coin results.
  # Inputs: win matrix P, probability p, and a reference win matrix.
  # Output: one row with reply changes, win range, and cycle length.
  # Approach: use the best reply and cycle functions on the results.
  
  replies <- best_replies(P)
  reference_replies <- best_replies(reference)
  cycle <- find_cycle(P)
  
  # Use unrounded probabilities to calculate the range.
  best_probabilities <- numeric(ncol(P))
  for (j in seq_len(ncol(P))) {
    best_probabilities[j] <- max(P[, j], na.rm = TRUE)
  }
  
  # NA means no cycle was found at the chosen cutoff.
  cycle_length <- if (length(cycle) == 0) NA_integer_ else length(cycle) - 1L
  
  data.frame(
    p = p,
    changed_replies = sum(replies$best_reply != reference_replies$best_reply),
    lowest_best_win = round(min(best_probabilities), 4),
    highest_best_win = round(max(best_probabilities), 4),
    shortest_cycle_edges = cycle_length,
    row.names = NULL
  )
}


# All experiments start here.
# Set the seed once. Later experiments use new random numbers.
run_start <- proc.time()[["elapsed"]]
set.seed(27)

# These examples can be checked by hand.
# The last sequence has two matches, but only the first end is returned.
examples <- data.frame(
  stream = c("HTTHTHH", "TTTT", "HT", "HTHTH"),
  pattern = rep("HTH", 4), expected = c(6, NA, NA, 3),
  observed = c(firsthit(c(1, 0, 0, 1, 0, 1, 1), c(1, 0, 1)),
               firsthit(c(0, 0, 0, 0), c(1, 0, 1)),
               firsthit(c(1, 0), c(1, 0, 1)),
               firsthit(c(1, 0, 1, 0, 1), c(1, 0, 1)))
)
cat("\nFirst occurrence examples\n")
print(examples)

# Both patterns use the same tosses in each experiment.
cat("\nOne shared‑stream duel: HTH versus THH\n")
print(duel(c(1, 0, 1), c(0, 1, 1)))
cat("Estimated P(HTH beats THH):", pwin(c(1, 0, 1), c(0, 1, 1)), "\n")

# Compare exact waiting times with simulated waiting times.
# The last example uses a biased coin.
wait_patterns <- list(c(1, 1), c(0, 0), c(1, 0), c(0, 1),
                      c(1, 0, 1), c(1, 0, 1))
wait_p <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.7)[c(1, 2, 3, 4, 5, 7)]
wait_comparison <- data.frame(
  pattern = c("HH", "TT", "HT", "TH", "HTH", "HTH"),
  p = wait_p, exact = NA_real_, simulated_mean = NA_real_, mc_se = NA_real_
)
for (i in seq_along(wait_patterns)) {
  waiting <- simulate_hits(wait_patterns[[i]], p = wait_p[i],
                           nreps = 20000)[, "tosses"]
  wait_comparison$exact[i] <- ewait(wait_patterns[[i]], wait_p[i])
  wait_comparison$simulated_mean[i] <- mean(waiting)
  
  # This is the Monte Carlo standard error of the estimated mean.
  # It is different from the standard deviation of one waiting time.
  wait_comparison$mc_se[i] <- sd(waiting) / sqrt(length(waiting))
}
cat("\nExact isolated waits and Monte Carlo estimates (20000 streams each)\n")
print(wait_comparison, digits = 4)

# Simulate each unordered pair once for the fair coin.
patterns <- all_patterns(3)
P_fair <- pmat(k = 3, p = 0.5, nreps = 10000)

cat("\n5(a) Best replies: fair coin, 10000 duels per pair\n")
print(best_replies(P_fair))

# The selected pairs have estimates between 0.48 and 0.52.
# Another simulation may move a pair across this boundary.
cat("\n5(b) Approximately balanced pairs (0.48 ~ 0.52 win prob)\n")
print(balanced_pairs(P_fair))
cat("\n5(b) Shortest non‑transitive cycle:\n")
show_cycle(P_fair)

# Waiting time describes each pattern alone.
# Mean win probability describes a uniformly chosen other opponent.
waits <- numeric(nrow(patterns))
for (i in seq_len(nrow(patterns))) waits[i] <- ewait(patterns[i, ])
ranking <- data.frame(
  pattern = rownames(patterns), expected_wait = waits,
  mean_win_probability = rowMeans(P_fair, na.rm = TRUE), row.names = NULL
)
cat("\n5(c) Isolated waiting‑time order versus mean head‑to‑head win probability\n")
print(ranking[order(ranking$expected_wait), ], digits = 4, row.names = FALSE)
cat("\n5(c) Example: same expected wait but unequal duel performance\n")
print(same_wait_example(P_fair, waits))

# An empty search sequence is allowed.
cat("\nEmpty sequence example: firsthit(integer(0), c(1, 0, 1)) = ")
print(firsthit(integer(0), c(1, 0, 1)))

# Add position names to the pattern matrix used in the experiments.
colnames(patterns) <- paste0("position_", seq_len(ncol(patterns)))

cat("\nInterpretation\n")
cat("P[i, j] estimates the chance that row i beats column j.\n")
cat("Mean win probability uses a uniformly chosen opponent from the other patterns.\n")
cat("It does not give a ranking that holds against every opponent.\n")

# Independent tosses do not mean independent pattern waiting times.
# Both waiting times use the same sequence.
# A pattern's ending may also help another pattern start.
# These links affect winning chances, even when mean waits are equal.
cat("Independent tosses do not make the two pattern waiting times independent.\n")
cat("Shared tosses and prefix-suffix links affect the duel result.\n")

# A clear win needs an estimate above 0.55.
# Near-balanced pairs give no edge under this rule.
# Missing edges can block short cycles.
# Near ties alone do not prove that every short cycle is absent.
cat("Cycle cutoff: every arrow needs an estimated win probability above 0.55.\n")
cat("Near-balanced pairs give no edge at this cutoff and can block short cycles.\n")
cat("The shortest cycle is found in the estimated graph; it is not a mathematical proof.\n")

# Opposite probabilities add to one.
# Therefore, two opposite edges cannot both have probability above 0.55.
cat("Two-node cycles are impossible because opposite probabilities add to one.\n")


# Use the same pattern pair in both plots.
plot_a <- c(0, 1, 1)
plot_b <- c(1, 1, 1)
plot_nreps <- 20000L

# tosses measures the time until either pattern wins.
# It is not the waiting time for THH alone.
plot_results <- simulate_hits(
  plot_a, plot_b, p = 0.5, nreps = plot_nreps
)
duel_tosses <- plot_results[, "tosses"]

# Keep values with zero counts so the time axis stays complete.
toss_values <- seq_len(max(duel_tosses))
frequencies <- tabulate(
  duel_tosses, nbins = max(duel_tosses)
) / plot_nreps
median_tosses <- median(duel_tosses)

# Run 50 separate estimates without resetting the seed.
estimate_runs <- 50L
estimate_nreps <- 2000L
win_estimates <- numeric(estimate_runs)
for (i in seq_len(estimate_runs)) {
  win_estimates[i] <- pwin(
    plot_a, plot_b, p = 0.5, nreps = estimate_nreps
  )
}

# HHH wins only if the first three tosses are all H.
# After any T, three H complete THH before they complete HHH.
# The exact fair-coin chance that THH wins is therefore 1 - 0.5^3.
exact_reference <- 1 - 0.5^3


# Leave clear space above the data for the median label.
frequency_top <- max(frequencies) * 1.25
plot(
  toss_values, frequencies,
  type = "h",
  lwd = 2,
  col = "steelblue",
  ylim = c(0, frequency_top),
  xlab = "Tosses until either pattern appears",
  ylab = "Observed relative frequency",
  main = "Duel duration: THH versus HHH",
  sub = paste("Fair coin;", plot_nreps, "independent duels")
)

# Small points help show the lower frequencies.
lines(
  toss_values, frequencies,
  type = "p",
  pch = 16,
  cex = 0.55,
  col = "steelblue"
)

# Mark the median with a vertical line.
lines(
  c(median_tosses, median_tosses),
  c(0, frequency_top * 0.90),
  col = "firebrick",
  lty = 2,
  lwd = 2
)
text(
  median_tosses,
  frequency_top * 0.95,
  labels = paste("Median =", median_tosses),
  pos = 4,
  cex = 0.9,
  col = "firebrick"
)

# Leave space above the estimates for the exact reference label.
estimate_range <- range(c(win_estimates, exact_reference))
estimate_limits <- c(
  estimate_range[1] - 0.005,
  estimate_range[2] + 0.015
)

plot(
  seq_len(estimate_runs), win_estimates,
  type = "p",
  pch = 16,
  cex = 0.8,
  col = "steelblue",
  ylim = estimate_limits,
  xlab = "Independent simulation run",
  ylab = "Estimated probability that THH wins",
  main = "Win probability estimates: THH versus HHH",
  sub = paste(
    estimate_runs, "runs;",
    estimate_nreps, "duels per run; fair coin"
  )
)
lines(
  seq_len(estimate_runs), win_estimates,
  col = "steelblue",
  lwd = 1
)
lines(
  c(1, estimate_runs),
  c(exact_reference, exact_reference),
  col = "firebrick",
  lty = 2,
  lwd = 2
)
text(
  1,
  estimate_limits[2] - 0.002,
  labels = paste("Exact reference =", exact_reference),
  adj = c(0, 1),
  col = "firebrick",
  cex = 0.9
)


cat("\nPlot experiment summary\n")
print(data.frame(
  comparison = "THH versus HHH",
  duration_duels = plot_nreps,
  median_tosses = median_tosses,
  estimate_runs = estimate_runs,
  duels_per_run = estimate_nreps,
  mean_estimated_win = round(mean(win_estimates), 4),
  sd_of_estimates = round(sd(win_estimates), 4),
  exact_win_probability = exact_reference,
  row.names = NULL
))
cat("Exact reference: HHH wins only when the first three tosses are HHH.\n")
cat("Two plot pages saved to sequence_duels.pdf.\n")


# Repeat the pair comparisons for two biased coins.
# Keep the same sample size and cycle cutoff.
P_low <- pmat(k = 3, p = 0.25, nreps = 10000)
P_high <- pmat(k = 3, p = 0.7, nreps = 10000)

biased_matrices <- list(P_low, P_high)
biased_probabilities <- c(0.25, 0.7)

for (i in seq_along(biased_matrices)) {
  cat(
    "\nBest replies: p =", biased_probabilities[i],
    "; 10000 duels per pair\n"
  )
  print(best_replies(biased_matrices[[i]]))
  
  # Reuse the same cycle search for each coin.
  cat("Cycle search with estimated win probabilities above 0.55\n")
  show_cycle(biased_matrices[[i]])
}

# Compare actual results after all three matrices have been simulated.
# The table shows changed replies, best win chances, and cycle lengths.
# No result is assumed before running the experiments.
bias_summary <- rbind(
  coin_summary(P_fair, 0.5, P_fair),
  coin_summary(P_low, 0.25, P_fair),
  coin_summary(P_high, 0.7, P_fair)
)

cat("\nComparison of fair and biased coins\n")
print(bias_summary, row.names = FALSE)

# Describe reply changes using the values from this run.
# Close estimated probabilities may select different replies in another run.
for (i in 2:nrow(bias_summary)) {
  cat(
    "At p =", bias_summary$p[i], ",",
    bias_summary$changed_replies[i],
    "of 8 estimated best replies differ from the fair-coin run.\n"
  )
}

# Compare cycle lengths only when the searches found cycles.
# A shorter cycle here means a shorter cycle in the estimated graph.
fair_cycle_length <- bias_summary$shortest_cycle_edges[1]

for (i in 2:nrow(bias_summary)) {
  biased_cycle_length <- bias_summary$shortest_cycle_edges[i]
  
  if (is.na(biased_cycle_length)) {
    cat(
      "At p =", bias_summary$p[i],
      ", no cycle was found above the 0.55 cutoff.\n"
    )
  } else if (is.na(fair_cycle_length)) {
    cat(
      "At p =", bias_summary$p[i],
      ", the shortest detected cycle has",
      biased_cycle_length,
      "edges; no fair-coin cycle was found.\n"
    )
  } else {
    cycle_change <- if (biased_cycle_length < fair_cycle_length) {
      "shorter than"
    } else if (biased_cycle_length > fair_cycle_length) {
      "longer than"
    } else {
      "the same length as"
    }
    
    cat(
      "At p =", bias_summary$p[i],
      ", the shortest detected cycle has",
      biased_cycle_length,
      "edges and is", cycle_change,
      "the fair-coin result.\n"
    )
  }
}

# Coin bias changes the probabilities of pattern prefixes.
# The tables show how reply choices and win chances change in this run.
# These simulated comparisons do not prove exact best replies or cycles.
cat("Reply and cycle comparisons use simulated probabilities and the same cutoff.\n")

# Include simulation, printing, and plotting in the total running time.
run_elapsed <- proc.time()[["elapsed"]] - run_start
cat(sprintf("\nTotal running time: %.2f seconds\n", run_elapsed))
