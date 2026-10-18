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

source("play_music.r")
waveform <- str_to_waveform(sim_music, bpm = 200)
play(waveform)
