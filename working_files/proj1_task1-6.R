#Practical 1: Music Simulator
#Group members:
#1.YOUR NAME - YOUR UNIVERSITY USERNAME
#2.YOUR NAME - YOUR UNIVERSITY USERNAME
#3.YOUR NAME - YOUR UNIVERSITY USERNAME
#Contributions:

#Task 1
#Create a GitHub repository for the project and clone it

#Task 2
#Download palestrina.txt from Learn and place it in the

#Task 3
#Read the music data into R

#setwd("D:/R-4.4.1/statistical_programming_group27")

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
length(music_raw)
head(music_raw, 10)
tail(music_raw, 10)
music_clean <- preprocess(music_raw)
#test task4
length(music_clean)
head(music_clean, 20)
tail(music_clean, 10)
#tset task5
notes <- unique(music_clean)
tokens <- match(music_clean, notes)

length(music_clean)
length(tokens)
length(music_clean) == length(tokens)

head(notes, 20)
head(tokens, 20)
which(notes == "||")
length(music_clean) == length(tokens)
#test task6
end <- which(notes == "||")

M <- make_matrix(
  tokens = tokens,
  end = end,
  mlag = 4
)

dim(M)
head(M)
sum(M[, 1:4] == end)
head(notes[M], 5)
