# tocID <- "integrateGeneticCodeOptimality.R"
#
# Purpose:  Computational Biology Foundations:
#              R code demonstrating the integration of observation
#              and simulation in assessing the optimality of the genetic
#              code.
#
# Version:  1.0
#
# Date:     2022-10
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           0.1    First version during lecture.
#           1.0    Finalized and posted.
#
#
#
# TODO:
#    ...
# ==============================================================================
#


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                                  Line
#TOC> ----------------------------------------------------------------------
#TOC>   1        Exploratory thoughts                                     42
#TOC>   2        First pass: lay out the workflow                         77
#TOC>   3        Second pass: functions and key datastructures            89
#TOC>   3.1        str2vec()                                             114
#TOC>   3.2        uCode                                                 129
#TOC>   3.3        randCode()                                            162
#TOC>   3.4        mapAdjacent()                                         211
#TOC>   3.5        A feature table                                       255
#TOC>   3.6        A similarity function                                 308
#TOC>   4        Third pass: implement the workflow                      333
#TOC>   4.1        Compare universal and random codes                    379
#TOC>
#TOC> ==========================================================================


# =    1  Exploratory thoughts  ================================================

# Is the genetic code the best it can be?
#
# What are the fixed properties of the code that we want to consider?
#   - length of a codon: 3
#       (nb. Does this need to be fixed. Cf. Huffman codes
#            ... but we need to remember that the code has _evolved_ .)
#
#   - number of amino acids: (20 + 1)
#


# Do elements that are adjacent in the code table have similar properties?
# Do elements that are not adjacent in the code table have dissimilar properties?
#
#
# Write  a function that takes two amino acids as an input and returns
# a value for their similarity.
#   Similarity:
#     - build a vector for each amino acid and take the dot product.
#     - distance between amino acids in feature space
#        - appropriately scale the feature values
#
# What makes a code _good_?
#
# The "best" code is "more good" than alternatives.
#
# Create random codes and compare to observed code.
#
#  A random mapping must:
#  - contain each amino acid at least once
#  - have the same level of redundancy
#

# =    2  First pass: lay out the workflow  ====================================
# Pseudocode:
#   Represent the universal genetic code (UGC) as a datastructure
#   for N trials
#     Construct a random genetic code (RGC)
#     for all adjacent positions in the code
#       - sum similarity values for the two amino acids that are encoded
#   Tabulate distribution of sums for all simulated RGC
#   Report where the UGC fall on this distribution
#


# =    3  Second pass: functions and key datastructures  =======================
#

# Represent the universal genetic code (UGC) as a data frame:
#     one column for amino acids (one letter code), one column each for
#     first, second, and third nucleobase
#   - copy the genetic code from a repository (e.g. NCBI)
#   - take e string of nucleotides or amino acids and strsplit() them
#     into a vector of one letter elements
#   - needs a function for that ... call it: str2vec()
#   - add each vector into the dataframe
# for N trials
#   Construct a random genetic code (RGC)
#     - needs a function ... call it randCode()
#   for all adjacent positions in the code
#     - needs a function to enumerate all adjacent positions, and a smart
#       datastructure that will make it easy to retrieve amino acid pairs
#       call it mapAdjacent()
#     - sum similarity values for the two amino acids that are encoded
#       - needs a function to calculate similarity values,
#         call it calcSimilarity()
# Tabulate distribution of summed similarity values
# Compare to universal genetic code
#

# ==   3.1  str2vec()  =========================================================

str2vec <- function(s) {
  # split a string into a vector of characters
  # s:     a string
  # value: a vector of characters
  return(strsplit(s, "")[[1]])
}

# verify:
str2vec("undr")           # [1] "u" "n" "d" "r"
str2vec("")               # character(0)
str2vec(c("tweedledum",
          "tweedledee"))  # [1] "t" "w" "e" "e" "d" "l" "e" "d" "u" "m"

# ==   3.2  uCode  =============================================================

# Define the Universal Genetic Code
#
# Genetic code: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
# AAs    = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
# Starts = ---M------**--*----M---------------M----------------------------
# Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
# Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
# Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

AAs   <- str2vec("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
Base1 <- str2vec("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG")
Base2 <- str2vec("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG")
Base3 <- str2vec("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG")

uCode <- data.frame(AA = AAs,
                    B1 = Base1,
                    B2 = Base2,
                    B3 = Base3)
uCode$codon <- paste0(uCode$B1, uCode$B2, uCode$B3) # use paste0() to combine
                                                    # the three nucleotides
                                                    # into a single codon string

rownames(uCode) <- uCode$codon   # this gives us a convenient way to retrieve
                                 # specific rows from the code

# Verify:
head(uCode)
uCode["ATG", ]  # Methionine (M)


# ==   3.3  randCode()  ========================================================

randCode <- function(inCode, doSwap = TRUE, doPermute = TRUE) {
  # inCode:    data frame with one column called "AA"
  # doSwap:    if TRUE (default) randomly swap characters. Note: swapping keeps
  #            the order, but changes what each character is - eg. ABBCCC might
  #            be converted into CAABBB
  # doPermute: if TRUE (default) permute the AA column. Note: permuting changes
  #            the order, but keeps each character's identity - eg. CAABBB might
  #            be converted to BBACBA
  # value:     a code data frame in which the "AA" column has been modified as
  #            directed.
  #

  # randomly swap characters
  if (doSwap) {
    alphabet <- unique(inCode$AA)
    alphrand <- sample(alphabet)
    # collapse these to strings and use chartr() to swap characters:
    inCode$AA <- chartr(paste(alphabet, collapse = ""),
                        paste(alphrand, collapse = ""),
                        inCode$AA)
  }

  # randomly shuffle translated codons
  if (doPermute) {
    inCode$AA <- sample(inCode$AA)
  }


  return(inCode)
}

# verify
randCode(data.frame(AA = c("A", "B", "B", "C", "C", "C")),
         doSwap = FALSE,
         doPermute = FALSE) # should return the same as the input

randCode(data.frame(AA = c("A", "B", "B", "C", "C", "C")),
         doPermute = FALSE) # swap only: keeps the pattern x y y z z z

randCode(data.frame(AA = c("A", "B", "B", "C", "C", "C")),
         doSwap = FALSE)    # permute only: positions shuffled but we still have
                            # one A, two B and three C
                            #
randCode(data.frame(AA = c("A", "B", "B", "C", "C", "C")))
                            # completely randomized positions and frequencies


# ==   3.4  mapAdjacent()  =====================================================

mapAdjacent <- function(inCode) {
  # inCode:    data frame with one column named "codon" and three columns
  #            named B1, B2, and B3 which contain the first, second and third
  #            codon respectively.
  # value:     a data frame with columns c1 and c2 in which each row contains
  #            two codons that are adjacent in the sense that one can be
  #            generated from the other with a single change (point mutation).
  adj <- data.frame(character(),
                    character())

  # to test all elements, think of a pair as a cell in a large square matrix
  # with codons as columns and codons as rows. But we don't need to look at
  # the diagonal (those are identical), and we only need to consider the
  # triangle above the diagonal, since the triangle below the diagonal
  # is symmetric to it.)
  for (i in 1:(nrow(inCode) - 1)) { # for all rows
    for(j in (i+1):nrow(inCode)) {
       # print(sprintf("[%2d, %2d]", i, j)) # verify the correct range of the
                                            # loop variables i, and j.

       x <- inCode[i, c("B1", "B2", "B3")]  # extract two vectors with the three
       y <- inCode[j, c("B1", "B2", "B3")]  # nucleotides of rows i and j

       if ( sum(x != y) == 1) { # if x and y have exactly one difference...
         # append the codon from row i column j to the output
         adj <- rbind(adj, c(inCode$codon[i], inCode$codon[j]))
       }
    }
  }

  colnames(adj) <- c("C1", "C2") # define the column names

  return(adj)

}

# Verify:
myAdj <- mapAdjacent(uCode)
nrow(myAdj)  # 288. Why is this correct? What is the expected number?
myAdj        # inspect ...


# ==   3.5  A feature table  ===================================================

# fetch three feature vectors from seqinr::aaindex
library(seqinr)
data(aaindex)

myF    <- data.frame(f1 = aaindex[[  2]]$I) # hydrophobicity
myF$f2 <-                 aaindex[[112]]$I # volume
myF$f3 <-                 aaindex[[  8]]$I  # flexibility

# Verify that these are largely independent: calculate
# coefficients of correlation ...
cor(myF$f1, myF$f2)
cor(myF$f1, myF$f3)
cor(myF$f2, myF$f3)

# To avoid problems with different absolute values of our features, we use
# scale() which centers the columns (each column's mean becomes 0) and scales
# the valuse so that each vector's standard deviation is one. This type of
# "normalization" is called "standardization". Scale can work on all columns
# of a matrix-like object, like a data frame.


myF <- scale(myF)


# Verify:
colMeans(f2)

# Change the rownames of myF to single-letter code
rownames(myF)                              # what are the rownames now?
seqinr::a(rownames(myF))                   # what do we want them to be ?
rownames(myF) <- seqinr::a(rownames(myF))  # reassign them
rownames(myF)                              # Verify

# This is useful, because we can extract the values of a data frame through
# their respective row- and column-names, for example
myF["C", ]  # ... to retrieve a three-element feature vector from myF


# Finally: what should the similarity of an amino acid and the stop codon be?
# We need this because "*" appears in the table, and stop codons are adjacent
# to other codons.
# We can define this in many ways. For now, I just set the value to
# zero in each column. This at least does not change the mean.
myF <- rbind(myF, c(0.0, 0.0, 0.0))
rownames(myF)[nrow(myF)] <- "*"      # ... the last element of rownames()

# Verify...
tail(myF)



# ==   3.6  A similarity function  =============================================

sim <- function(a, b, F = myF) {
  # a, b:  amino acids in one-letter code
  # F:     an amino acid feature table (default: myF)
  # value: the Euclidian distance between the two vectors F[a] and F[b]
  #
  # Note: since the basic mathematical operations in R are all "vectorized"
  #       we can use this function for any number of columns in our input
  #       table of features, i.e. if we were to add a fourth, fifth etc.
  #       column

  d <- sqrt(sum( (F[a, ] - F[b, ])^2 ))  # distance: square root of the sum of
                                         # squared differences
  return(d)
}

# Verify:
myF["K", ]  #  0.1846380   0.8142622   0.4724765
myF["E", ]  # -0.63866590 -0.02386329  0.84611870

sqrt((0.185 + 0.639)^2 + (0.814 + 0.024)^2 + (0.472 - 0.846)^2) # 1.233327
sim("K", "E") # 1.232839


# =    4  Third pass: implement the workflow  ==================================
#

# prepare datasets
# we have prepared above: uCode    - genetic code
#                         myAdj    - all adjacent codon pairs
#                         myF      - a table of amino acid features


# for N trials
# N <- 10     # Always try with a small number first and get a sense of how long
# N <- 100    # it will take. Immediate? Time for a coffee? A beer? Overnight?
# N <- 1000
N <- 100000

sumRndSim <- numeric(N)  # prepare a vector to collect similarity sums from
                         # random codes

startTime <- Sys.time()
for (i in 1:N) {
   # show a progress bar ... this may take a while
   pBar(i, N)

   #  Construct a random code
   rndCode <- randCode(uCode)
   sims <- numeric(nrow(myAdj))    # prepare  a vector to collect similarities

   for (j in 1:nrow(myAdj)) {      # for all adjacent positions ...
     codon1 <- myAdj$C1[j]         # first codon
     codon2 <- myAdj$C2[j]         # second codon
     aa1 <- rndCode[codon1, "AA"]  # translate to first amino acid
     aa2 <- rndCode[codon2, "AA"]  # translate to second amino acid
     sims[j] <- sim(aa1, aa2)      # compute and store the similarity
   }
   sumRndSim[i] <- sum(sims)       # store the sum for this random code, and
                                   # repeat ...
}
endTime <- Sys.time()

endTime - startTime                # 25 minutes for me with 100,000 trials

saveRDS(sumRndSim, "sumRndSim.Rds") # ... back it up, just in case

summary(sumRndSim)


# ==   4.1  Compare universal and random codes  ================================

# compute a similarity score for the natural code
sims <- numeric(nrow(myAdj))    # prepare  a vector to collect similarities

for (j in 1:nrow(myAdj)) {      # for all adjacent positions ...
  codon1 <- myAdj$C1[j]         # first codon
  codon2 <- myAdj$C2[j]         # second codon
  aa1 <- uCode[codon1, "AA"]    # translate to first amino acid
  aa2 <- uCode[codon2, "AA"]    # translate to second amino acid
  sims[j] <- sim(aa1, aa2)      # compute and store the similarity
}

uCodeSim <- sum(sims)

# plot this value on the histogram of simulated values
hist(sumRndSim,
     xlim = c(400, 800),
     breaks = 50,
     col = "#b3e0dc",
     main = "Similarity of adjacent codons in random codes",
     xlab = "Sum of adjacent pairwise similarities")

# add the natural code's value with a red line.
abline(v = uCodeSim, col = "#AA0000")


# Hm.

# How many standard deviations is this away from the mean?
mySigma <- sd(sumRndSim)

(mean(sumRndSim) - uCodeSim) / mySigma

# What is the density of a normal distribution with these parameters at the
# value of uCodeSim ?

dnorm(uCodeSim, mean = mean(sumRndSim), sd = mySigma)
# What does this mean?


# Just some thoughts:
#   - Could we use this code to _evolve_ a better genetic code?
#   - Can we distinguish between an "exploring" (i.e.: optimize change) and a
#     "defensive" (i.e. minimize disruption) code?
#   - Are we making unjustified assumptions? For example: we are assuming that
#     all adjacent codons can be reached with the same probability...
#   - More?



# [END]
