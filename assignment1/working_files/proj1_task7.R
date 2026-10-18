# Task 7. Source this file after the matrix M has been made.

next_note <- function(key, M, w = rep(1, ncol(M) - 1)) {
  if (!is.matrix(M) || ncol(M) < 2 || nrow(M) < 1) {
    stop("M must be a non-empty sequence matrix")
  }

  mlag <- ncol(M) - 1
  if (!is.numeric(w) || length(w) != mlag ||
      any(!is.finite(w)) || any(w < 0)) {
    stop("w must have one non-negative number for each lag")
  }
  if (anyNA(key)) stop("key cannot contain NA")

  key <- tail(key, mlag)
  next_tokens <- integer(0)
  chance <- numeric(0)
  matching_rows <- seq_len(nrow(M))

  # Start with the last token. Add one earlier token each time.
  for (lag in seq_len(length(key))) {
    column <- mlag - lag + 1
    key_position <- length(key) - lag + 1
    matching_rows <- matching_rows[
      M[matching_rows, column] == key[key_position]
    ]

    if (length(matching_rows) == 0) break

    if (w[column] > 0) {
      next_tokens <- c(next_tokens, M[matching_rows, mlag + 1])
      chance <- c(chance, rep(w[column] / length(matching_rows),
                              length(matching_rows)))
    }
  }

  # If no weighted match was found, use all next tokens in M.
  if (length(next_tokens) == 0) {
    next_tokens <- M[, mlag + 1]
    chance <- NULL
  }

  if (length(next_tokens) == 1) return(next_tokens[1])
  sample(next_tokens, size = 1, prob = chance)
}
