# Group contributions: replace the placeholders with actual work and shares.
# Ruiqi Wang: [contribution and approximate percentage].
# Zirui Guo: [contribution and approximate percentage].
# Tong Shen: [contribution and approximate percentage].
# Repository: https://github.com/chantalshen1213-dev/statistical_programming_group27
# AI-generated reference draft for the group's study and revision.
#
# Sequence duels and the perils of independence
# Two players watch the SAME stream of independent coin tosses. Each chooses a
# different binary pattern of equal length. The first completed pattern wins.
# Heads is 1, tails is 0, and p is the probability of heads on each toss.
#
# The functions below locate patterns, simulate shared-stream duels and compute
# exact isolated waiting times from self-overlaps. The experiments then compare
# all length-three patterns, find best replies and non-transitive cycles, plot
# simulation results, and repeat the analysis for biased coins. The comparison
# asks whether isolated waiting times can explain head-to-head performance.
#
# Only base R is used. All function definitions precede the experiments.
# Running the file creates a two-page PDF in the current working directory.


# Step 1
check_pattern <- function(a) {
  # Purpose: reject patterns that cannot represent a non-empty coin sequence.
  # Input: a, a candidate pattern vector.
  # Output: invisible NULL, or an informative error for an invalid pattern.
  # Approach: check the vector structure and membership in the set {0, 1}.
  if (!(is.numeric(a) || is.logical(a)) || !is.null(dim(a)) ||
      length(a) == 0 || anyNA(a) || any(!(a %in% c(0, 1)))) {
    stop("A pattern must be a non-empty vector containing only 0 and 1.")
  }
  invisible(NULL)
}

firsthit <- function(x, a) {
  # Purpose: locate the first completed occurrence of a in an existing stream.
  # Inputs: x, a binary toss vector; a, a non-empty binary pattern.
  # Output: the ending position of the first match, or NA if no match exists.
  # Approach: inspect consecutive length-k windows in chronological order.
  check_pattern(a)
  if (!(is.numeric(x) || is.logical(x)) || !is.null(dim(x)) ||
      anyNA(x) || any(!(x %in% c(0, 1)))) {
    stop("x must be a vector containing only 0 and 1.")
  }
  k <- length(a)
  if (length(x) < k) return(NA_integer_)

  for (last in seq.int(k, length(x))) {
    if (all(x[(last - k + 1):last] == a)) return(last)
  }
  NA_integer_
}


# Step 2
check_simulation <- function(p, nreps) {
  # Purpose: enforce the coin model and a valid number of simulated streams.
  # Inputs: p, a heads probability; nreps, a replication count.
  # Output: invisible NULL, or an informative error.
  # Approach: require 0 < p < 1 and a finite positive integer count.
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
  # Purpose: provide the shared simulation engine for duels and waiting times.
  # Inputs: patterns a and optionally b, heads probability p, and nreps streams.
  # Output: a matrix with winner (1 for a, 0 for b) and tosses (completion time).
  # Approach: generate independent rows in blocks, compare all windows in each
  # block, and extend only unfinished rows. With b=NULL, wait for a alone.
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
  active <- seq_len(nreps)
  history <- matrix(0L, nrow = nreps, ncol = k - 1)
  elapsed <- 0L

  while (length(active) > 0) {
    tosses <- matrix(rbinom(length(active) * block_size, size = 1, prob = p),
                     nrow = length(active), ncol = block_size)

    # Retaining k-1 previous tosses preserves matches crossing a block boundary.
    # Within a row, BOTH patterns are compared with this same buffer.
    buffer <- cbind(history, tosses)
    match_a <- matrix(TRUE, nrow = length(active), ncol = block_size)
    match_b <- matrix(!is.null(b), nrow = length(active), ncol = block_size)
    for (j in seq_len(k)) {
      window_toss <- buffer[, j:(j + block_size - 1), drop = FALSE]
      match_a <- match_a & window_toss == a[j]
      if (!is.null(b)) match_b <- match_b & window_toss == b[j]
    }
    matches <- match_a | match_b

    # Initial history contains padding, so a window cannot count until k real
    # tosses exist. This also handles patterns longer than a single block.
    incomplete <- min(block_size, k - 1 - elapsed)
    if (incomplete > 0) matches[, seq_len(incomplete)] <- FALSE

    # In a logical row the first TRUE is its first maximum. An all-FALSE row
    # also has a maximum, so check the selected entry before declaring a hit.
    first <- max.col(matches, ties.method = "first")
    finished <- matches[cbind(seq_len(length(active)), first)]
    hit <- which(finished)
    if (length(hit) > 0) {
      result[active[hit], "winner"] <-
        as.integer(match_a[cbind(hit, first[hit])])
      result[active[hit], "tosses"] <- elapsed + first[hit]
    }

    history <- buffer[!finished, block_size + seq_len(k - 1), drop = FALSE]
    active <- active[!finished]
    elapsed <- elapsed + block_size
  }
  result
}

duel <- function(a, b, p = 0.5) {
  # Purpose: simulate one duel until either pattern first appears.
  # Inputs: distinct equal-length binary patterns a and b; heads probability p.
  # Output: named vector c(winner, tosses), with winner=1 for a and 0 for b.
  # Approach: use one row of the same block simulator used by pwin.
  simulate_hits(a, b, p, nreps = 1)[1, ]
}


# Step 3
pwin <- function(a, b, p = 0.5, nreps = 10000) {
# Purpose: estimate the probability that a beats b on a shared toss stream.
# Inputs: distinct equal-length patterns a and b, p, and a replication count.
# Output: a single Monte Carlo win-probability estimate.
# Approach: average a's win indicators across independent simulated duels.
  results <- simulate_hits(a, b, p, nreps)
  mean(results[, "winner"])
}


# check_pattern()
# check_simulation()
# firsthit()
# simulate_hits()
# duel()
# pwin()
#----------------------------------------------------------------------
# Task 4 ewait: exact expected waiting time (self‑overlap formula)
#----------------------------------------------------------------------
ewait <- function(a, p = 0.5) {
  # Purpose: calculate the exact expected wait for one pattern in isolation.
  # Inputs: binary pattern a and heads probability p.
  # Output: its expected completion time, measured from the start of a stream.
  # Approach: for each matching prefix and suffix of length j, add the inverse
  # prefix probability. Full‑length overlap always contributes to the sum.
  check_pattern(a)
  check_simulation(p, 1)
  k <- length(a)
  expected <- 0
  for (j in seq_len(k)) {
    prefix <- a[seq_len(j)]
    suffix <- a[(k - j + 1):k]
    if (all(prefix == suffix)) {
      heads <- sum(prefix)
      probability <- p^heads * (1 - p)^(j - heads)
      expected <- expected + 1 / probability
    }
  }
  expected
}

#----------------------------------------------------------------------
# Helper: generate all k‑length binary patterns, fixed documented ordering
# Output matrix: each row = pattern vector 0/1; rownames H/T labels
# Order(k=3): TTT,TTH,THT,THH,HTT,HTH,HHT,HHH
#----------------------------------------------------------------------
all_patterns <- function(k = 3) {
  if (!is.numeric(k) || length(k) != 1 || !is.finite(k) ||
      k < 1 || k != floor(k)) {
    stop("k must be a positive integer.")
  }
  patterns <- matrix(0L, nrow = 2^k, ncol = k)
  ids <- 0:(2^k - 1)
  labels <- rep("", nrow(patterns))
  for (j in seq_len(k)) {
    patterns[, j] <- (ids %/% 2^(k - j)) %% 2
    labels <- paste0(labels, ifelse(patterns[, j] == 1, "H", "T"))
  }
  rownames(patterns) <- labels
  patterns
}

#----------------------------------------------------------------------
# Task5 pmat: full 2^k × 2^k duel win‑probability matrix
# P[i,j] = probability pattern i beats pattern j; diagonal NA
# Exploit P[j,i] = 1‑P[i,j] to halve simulation workload
#----------------------------------------------------------------------
pmat <- function(k = 3, p = 0.5, nreps = 10000) {
  check_simulation(p, nreps)
  patterns <- all_patterns(k)
  n <- nrow(patterns)
  P <- matrix(NA_real_, nrow = n, ncol = n,
              dimnames = list(rownames(patterns), rownames(patterns)))
  for (i in seq_len(n - 1)) {
    for (j in seq.int(i + 1, n)) {
      P[i, j] <- pwin(patterns[i, ], patterns[j, ], p, nreps)
      P[j, i] <- 1 - P[i, j]
    }
  }
  P
}

#----------------------------------------------------------------------
# 5(a) helper: get best‑reply for each opponent pattern
#----------------------------------------------------------------------
best_replies <- function(P) {
  reply <- integer(ncol(P))
  for (j in seq_len(ncol(P))) reply[j] <- which.max(P[, j])
  data.frame(
    opponent = colnames(P),
    best_reply = rownames(P)[reply],
    win_probability = round(P[cbind(reply, seq_len(ncol(P)))], 4),
    row.names = NULL
  )
}

#----------------------------------------------------------------------
# 5(b) helper: find approximately evenly‑matched pairs
#----------------------------------------------------------------------
balanced_pairs <- function(P, tolerance = 0.02) {
  pairs <- which(upper.tri(P) & abs(P - 0.5) <= tolerance, arr.ind = TRUE)
  data.frame(
    pattern_a = rownames(P)[pairs[, 1]],
    pattern_b = colnames(P)[pairs[, 2]],
    probability_a_wins = round(P[pairs], 4),
    row.names = NULL
  )
}

#----------------------------------------------------------------------
# 5(b) helper: BFS search shortest non‑transitive win cycle
#----------------------------------------------------------------------
find_cycle <- function(P, margin = 0.05) {
  n <- nrow(P)
  shortest <- integer(0)
  for (start in seq_len(n)) {
    queue <- start
    parent <- rep(NA_integer_, n)
    parent[start] <- 0L
    position <- 1L
    while (position <= length(queue)) {
      current <- queue[position]
      position <- position + 1L
      neighbours <- which(P[current, ] > 0.5 + margin)
      if (current != start && start %in% neighbours) {
        path <- current
        while (path[1] != start) path <- c(parent[path[1]], path)
        path <- c(path, start)
        if (length(shortest) == 0 || length(path) < length(shortest)) {
          shortest <- path
        }
        break
      }
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
  cycle <- find_cycle(P, margin)
  if (length(cycle) == 0) {
    cat("No cycle detected with every win probability above", 0.5 + margin, "\n")
    return(invisible(NULL))
  }
  cat("Shortest detected cycle:",
      paste(rownames(P)[cycle], collapse = " -> "), "\n")
  from <- head(cycle, -1)
  to <- tail(cycle, -1)
  print(data.frame(
    winner = rownames(P)[from], opponent = colnames(P)[to],
    win_probability = round(P[cbind(from, to)], 4), row.names = NULL
  ))
  invisible(NULL)
}

#----------------------------------------------------------------------
# 5(c) helper: find example pair with identical ewait but asymmetric duel
#----------------------------------------------------------------------
same_wait_example <- function(P, waits) {
  difference <- abs(outer(waits, waits, "-"))
  pairs <- which(upper.tri(P) & difference < 1e-10, arr.ind = TRUE)
  if (nrow(pairs) == 0) return(NULL)
  pair <- pairs[which.max(abs(P[pairs] - 0.5)), ]
  if (P[pair[1], pair[2]] < 0.5) pair <- rev(pair)
  data.frame(
    winner = rownames(P)[pair[1]], opponent = colnames(P)[pair[2]],
    expected_wait_both = waits[pair[1]],
    win_probability = round(P[pair[1], pair[2]], 4), row.names = NULL
  )
}

#----------------------------------------------------------------------
# Experiment execution: Task4 validation, Task5 a/b/c, Task6 plot PDF, Task7 biased‑coin analysis
#----------------------------------------------------------------------
run_start <- proc.time()[["elapsed"]]
set.seed(27)

# Task1 manual test example
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

# Task2, Task3 quick demo
cat("\nOne shared‑stream duel: HTH versus THH\n")
print(duel(c(1, 0, 1), c(0, 1, 1)))
cat("Estimated P(HTH beats THH):", pwin(c(1, 0, 1), c(0, 1, 1)), "\n")

# Task4 validation: exact ewait vs Monte‑Carlo simulation
wait_patterns <- list(c(1, 1), c(0, 0), c(1, 0), c(0, 1),
                      c(1, 0, 1), c(1, 0, 1))
wait_p <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.7)
wait_comparison <- data.frame(
  pattern = c("HH", "TT", "HT", "TH", "HTH", "HTH"),
  p = wait_p, exact = NA_real_, simulated_mean = NA_real_, mc_se = NA_real_
)
for (i in seq_along(wait_patterns)) {
  waiting <- simulate_hits(wait_patterns[[i]], p = wait_p[i],
                           nreps = 20000)[, "tosses"]
  wait_comparison$exact[i] <- ewait(wait_patterns[[i]], wait_p[i])
  wait_comparison$simulated_mean[i] <- mean(waiting)
  wait_comparison$mc_se[i] <- sd(waiting) / sqrt(length(waiting))
}
cat("\nExact isolated waits and Monte Carlo estimates (20000 streams each)\n")
print(wait_comparison, digits = 4)

# Task5 k=3 fair coin analysis
patterns <- all_patterns(3)
P_fair <- pmat(k = 3, p = 0.5, nreps = 10000)

cat("\n5(a) Best replies: fair coin, 10000 duels per pair\n")
print(best_replies(P_fair))

cat("\n5(b) Approximately balanced pairs (0.48 ~ 0.52 win prob)\n")
print(balanced_pairs(P_fair))
cat("\n5(b) Shortest non‑transitive cycle:\n")
show_cycle(P_fair)

#5(c) compare isolated waiting‑time ranking vs duel head‑to‑head performance
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

