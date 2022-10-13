# shannonEntropy.R
#
# Demo code for CSB195-2022
#
# 2022-10-12
# Boris Steipe (boris.steipe@utoronto.ca)
#
# Version:  1.0
#
# ==============================================================================

# Load the standard genetic code from the Biostrings package.
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

stdGC <- Biostrings::GENETIC_CODE

# define a function to compute the Shannon entropy of a Genetic Code
HGC <- function(GC) {
  # GC   char   a genetic code with 64 character elements
  if (length(GC)         != 64) {stop("Input does not have 64 codons.")}
  if (length(unique(GC)) != 21) {stop("Input does not have 21 amino acids.")}

  # convert to PMF  (Probability Mas Function). Sum is 1.0.
  V <- table(GC) / 64
  # compute entropy base 2
  H <- -sum(V * (log(V) / log(2)))
  return(H)
}

# Analyse
HGC(stdGC) # 4.218139 bits

# maximum entropyy code
( maxEcode <- c(rep(unique(stdGC), 3), "A") )
HGC(maxEcode) # 4.389098

# minimum entropyy code
( minEcode <- c(unique(stdGC), rep("A", 43)) )
HGC(minEcode) # 2.246641

# We see that the natural code has close to maximal entropy! Only 4%
# less then the maximum. Remarkable.

# How does the standard code compare to random codes? Random codes for this
# purpose are constructed by appending 43 random amino acids to the "alphabet"
# of 21 amino acids. That way we gurantee that all amino acids are present,
# and the code is 64 characters long. For calculating entropy, order does not
# matter, it makes no difference that the first 21 characters are always the
# same.
#
N <- 10000             # run N trials
rndC <- numeric(N)     # prepare vector to store N results
alf <- unique(stdGC)   # extract the alphabet from the stdCode
for (i in 1:N) {
  x <- c(alf, sample(alf, 43, replace = TRUE)) # concatenate 21 fixed and
                                               # 43 random characters
  rndC[i] <- HGC(x)    # store the entropy
}

hist(rndC,
     breaks = 100,
     xlim = c(4.0, 4.4),
     main = "Shannon entropy of random genetic codes",
     xlab = "H (bits)")
abline(v = HGC(stdGC), col = "#BB0000")       # standard code
abline(v = c(HGC(minEcode),
             HGC(maxEcode)), col = "#00DD00") # max- and min-entropy codes

# Here we see that the standard code is virtually indistinguishable from random.


# [END]
