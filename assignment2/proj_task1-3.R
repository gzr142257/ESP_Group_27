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
