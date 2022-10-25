# fibWord.R
#
# Construct a Fibonacci Word of a given length
#
# cf. https://en.wikipedia.org/wiki/Fibonacci_word

fibWord <- function(N = 8*13) {

  s <- integer(N+2)
  i <- 1
  iLast <- i+1

  while (iLast < N+1) {
    if (s[i] == 1) { s[iLast]             <- 0;   iLast <- iLast+1
    } else         { s[iLast:(iLast + 1)] <- 1:0; iLast <- iLast+2
    }
    i <- i+1
  }
  return(paste(s[1:N], collapse = ""))
}

if (FALSE) {

  fibWord(1)
  fibWord(3)
  fibWord(5)
  fibWord(8)
  fibWord(13)
  fibWord(21)
  fibWord(1000)

  s <- fibWord()
  iStart <- ((0:7)*13)+1
  iStop  <- iStart+12
  cat(sprintf("\n%s", substring(s, iStart, iStop)))

}

# [END]
