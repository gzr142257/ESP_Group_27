#Group members:
#Ruiqi Wang s2930743
#Zirui Guo s2948295
#Tong Shen s2939003

#Contributions:
#Ruiqi Wang (33.3%): Implemented Tasks 3-6, including data import, preprocessing, note-to-token mapping, and the sequence matrix.
#Zirui Guo (33.3%): Implemented Task 7 (next_note()) and contributed validation and documentation for the preprocessing and matrix code.
#Tong Shen (33.3%): Implemented Task 8 (generate_music()), compared the three weight settings, and integrated and checked the full workflow.
#All members contributed to the GitHub repository setup, dataset download, code review, and debugging.

#Task 1
#Create a GitHub repository for the project and clone it

#Task 2
#Download palestrina.txt from Learn and place it in the

#Task 3
#Read the music data into R

# setwd("path/to/your/local/repository")

music_raw <- readLines("palestrina.txt")

#Inspect the first 10-20 lines of the raw data
head(music_raw, 20)

# Task 4
# Preprocess the raw music data

preprocess <- function(music) {
  #Remove the first three lines containing information
  music <- music[-(1:3)]

  #Remove metadata/comment lines beginning with "#".
  music <- music[!grepl("^#", music)]

  #Identify blank lines separating different musical pieces.
  #trimws() removes spaces before checking whether a line is blank.
  blank <- trimws(music) == ""

  #Give each piece a group number.
  piece_id <- cumsum(blank)

  #Keep only non-blank lines.
  music_nonblank <- music[!blank]
  piece_id <- piece_id[!blank]

  #Split the music into separate pieces.
  pieces <- split(music_nonblank, piece_id)

  #Process each piece separately.
  clean_pieces <- lapply(pieces, function(piece) {

    #Remove leading/trailing spaces from each line.
    piece <- trimws(piece)

    #Remove any accidentally empty lines.
    piece <- piece[nzchar(piece)]

    #Add "|" to the end of every measure except the final measure of the piece.
    if (length(piece) > 1) {
      piece[1:(length(piece) - 1)] <-
        paste0(piece[1:(length(piece) - 1)], " |")
    }

    # Add "||" to the final measure of the piece.
    piece[length(piece)] <-
      paste0(piece[length(piece)], " ||")

    #Combine all measures into one character string.
    piece <- paste(piece, collapse = " ")

    #Split the string into individual tokens.
    strsplit(piece, "\\s+")[[1]]
  })

  #Combine all pieces into one long vector.
  music_clean <- unlist(clean_pieces, use.names = FALSE)

  return(music_clean)
}

# Run the preprocessing function
music_clean <- preprocess(music_raw)


# Check the final ten elements.
tail(music_clean, n = 10)

# Task 5
# Create the unique note vector and token vector
notes <- unique(music_clean)

# Convert every element of music_clean into the position of that note in the notes vector.

tokens <- match(music_clean, notes)
# Check that tokens and music_clean have the same length.
length(music_clean)
length(tokens)

# The following should return TRUE.
length(music_clean) == length(tokens)

#Task 6
#Construct the sequence matrix M
#make_matrix creates a matrix containing all consecutive sequences of mlag + 1 tokens.

#Inputs:
#tokens - integer vector containing the token IDs
#end    - token ID corresponding to "||"
#mlag   - maximum lag/order of the model

#Output:
#M - matrix with (n - mlag) rows and (mlag + 1)
#columns

make_matrix <- function(tokens, end, mlag = 4) {

  #Number of tokens in the music data.
  n <- length(tokens)

  #There must be more tokens than the requested lag.
  if (n <= mlag) {
    stop("tokens must contain more elements than mlag")
  }

  #Create an empty matrix with the required dimensions.
  M <- matrix(
    NA_integer_,
    nrow = n - mlag,
    ncol = mlag + 1
  )

  # Fill each column with a shifted version of tokens.
  # Column 1:
  # tokens[1] to tokens[n-mlag]
  # Column 2:
  # tokens[2] to tokens[n-mlag+1]
  # Final column:
  # tokens[mlag+1] to tokens[n]
  for (j in seq_len(mlag + 1)) {
    M[, j] <- tokens[
      j:(n - mlag + j - 1)
    ]
  }
  # Remove rows where "||" occurs within the first
  # mlag columns.
  keep <- rowSums(
    M[, 1:mlag, drop = FALSE] == end
  ) == 0
  M <- M[keep, , drop = FALSE]
  return(M)
}

# Create the matrix M required for the model
# Find the token ID corresponding to "||".
end <- which(notes == "||")

# Construct the sequence matrix using a maximum lag of 4.
M <- make_matrix(
  tokens = tokens,
  end = end,
  mlag = 4
)

# Basic checks for Task 6
# Display the dimensions of M.
dim(M)
# Display the first few rows.
head(M)
# Display the first few rows converted back to note names.
head(notes[M])

#test
#test task3
#length(music_raw)
#head(music_raw, 10)
#tail(music_raw, 10)
#music_clean <- preprocess(music_raw)

#test task4
#length(music_clean)
#head(music_clean, 20)
#tail(music_clean, 10)

#test task5
#notes <- unique(music_clean)
#tokens <- match(music_clean, notes)

#length(music_clean)
#length(tokens)
#length(music_clean) == length(tokens)

#head(notes, 20)
#head(tokens, 20)
#which(notes == "||")
#length(music_clean) == length(tokens)

#test task6
#end <- which(notes == "||")

#M <- make_matrix(
#  tokens = tokens,
#  end = end,
#  mlag = 4
#)

#dim(M)
#head(M)
#sum(M[, 1:4] == end)
#head(notes[M], 5)

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

  # Keep matching rows when one earlier token is added.
  for (lag in seq_len(length(key))) {
    mc <- mlag - lag + 1
    current_key <- tail(key, lag)
    ii <- colSums(!(t(M[matching_rows, mc:mlag, drop = FALSE]) == current_key))
    matching_rows <- matching_rows[which(ii == 0)]

    if (length(matching_rows) == 0) break

    if (w[mc] > 0) {
      next_tokens <- c(next_tokens, M[matching_rows, mlag + 1])
      chance <- c(chance, rep(w[mc] / length(matching_rows),
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

#Task 8: Generate music until one of the stopping conditions is reached.

generate_music <- function(start, M, w = rep(1, ncol(M) - 1), notes) {

  #Find the token IDs for the piece-ending symbol and measure bar.
  end_token <- which(notes == "||")
  bar_token <- which(notes == "|")

  #Initialise the generated sequence with the starting token.
  music_tokens <- start

  #Continue while the sequence is under 200 tokens,
  #has not ended with ||, and has fewer than 8 measure bars.
  while (
    length(music_tokens) < 200 &&
    tail(music_tokens, 1) != end_token &&
    sum(music_tokens == bar_token) < 8
  ) {

    #Ask next_note() to predict one token from the sequence so far.
    new_token <- next_note(
      key = music_tokens,
      M = M,
      w = w
    )

    #Append the predicted token to the end of the sequence.
    music_tokens <- c(music_tokens, new_token)
  }

  #Convert token IDs back to their score symbols.
  music_notes <- notes[music_tokens]

  #Join the symbols into one space-separated string.
  music_string <- paste(music_notes, collapse = " ")

  #Return the generated music string.
  return(music_string)
}


set.seed(1)

#w1 combines context lengths 4, 3, 2, and 1.
sim_music <- generate_music(
  start = 20,
  M = M,
  notes = notes,
  w = c(1, 1, 1, 1)
)

#w2 uses only the longest context: the last four tokens.
sim_2 <- generate_music(
  start = 20,
  M = M,
  notes = notes,
  w = c(1, 0, 0, 0)
)

#w3 uses only the shortest context: the last one token.
sim_3 <- generate_music(
  start = 20,
  M = M,
  notes = notes,
  w = c(0, 0, 0, 1)
)

#The model predicts one token at a time and does not enforce the same total
#note duration in every measure. Therefore, measure durations may vary.

#With w1 = c(1, 1, 1, 1), all four history lengths are used.
#With seed 1, it produced 51 symbols and several runs of nearby notes.
#The eight bars had different numbers of notes.

#With w2 = c(1, 0, 0, 0), only the last four tokens are used.
#With seed 1, it produced 27 symbols. Four bars in a row were "Rb |".
#When a four-token match is missing, the function uses the fallback sample.

#With w3 = c(0, 0, 0, 1), only the last token is used.
#With seed 1, it produced 26 symbols, with short bars and some large
#jumps between notes. Earlier notes do not affect its next prediction.

#source("play_music.r")
#play(str_to_waveform(sim_music, bpm = 200))