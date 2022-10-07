# tocID <- "benchmarkAndRcppDemos.R"
#
# Tutorial and Demos for code benchmarking and RCPP speedup
#
# Version:  2.0
#
# Versions: 1.0  2018 Biochemistry bioinformatics graduate course
#           2.0  Extensively commented and package updates for CSB195
#
# 2017-03  - 2022-10
#
# Boris Steipe
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                                Line
#TOC> --------------------------------------------------------------------
#TOC>   1        Preparations                                           67
#TOC>   1.1        Install packages from CRAN and Bioconductor          69
#TOC>   1.2        Update existing packages                            162
#TOC>   1.3        Check that you have all files                       173
#TOC>   2        The simplest Rcpp program                             194
#TOC>   3        Scenario: Sequence to Codons                          205
#TOC>   4        Benchmarking                                          378
#TOC>   4.1        Three ways to benchmark                             390
#TOC>   4.2        Comparing our approaches                            424
#TOC>   5        Three RCPP versions                                   448
#TOC>   5.1        Compile C++ code and bind to an R function          451
#TOC>   5.2        Final comparison ...                                483
#TOC>
#TOC> ==========================================================================


# Speeding up R-code may be valuable, in particular if you are running lengthy
# simulations. Of course, you always have to balance the added development time
#  against the reduced execution time. Thus speedup generally is important only
#  for functions or tools that you will use again and again. Also, these days
#  we would almost never write an entire program in a compiled language
#  like C / C++ / Rust - but we would typically handle
#   - data preprocessing,
#   - bringing the program configuration parameters into a sane state,
#   - statistical analysis of results,
#   - plotting figures,
#   - and generating reports
#  with scripted languages like R or Python. These aspects
#  change frequently, and prototyping code in R and Python is generally faster
#  and easier, AND the large number of packages that are available allows for
#  efficient code reuse and profiting from state-of-the-art implementations
#  that are contributed by a large expert community in practically any
#  application domain. That said, if the "heavy-lifting" parts of your code
#  is a significant burden on your resources, you may be able to achieve
#  speed-ups of some two orders of magnitude if you run these in a compiled
#  language.
#
#  However: you need to SHOW that that's where the bottleneck actually is,
#  not just assume it, and therefore we need to talk about benchmarking before
#  we talk about mixed language code.
#
#  The scenario we use is quite straightforward: write a function that takes
#  a string as input and return a vector of substrings of length 3. Or, in
#  "biospeak": split a nucleotide sequence into its codons.
#


# =    1  Preparations  ========================================================

# ==   1.1  Install packages from CRAN and Bioconductor  =======================


# ====
# RCPP - run C++ code from within R
# ====

if (! requireNamespace("Rcpp", quietly=TRUE)) {
  install.packages("Rcpp")
}
# Package information:
#  library(help = Rcpp)       # basic information
#  browseVignettes("Rcpp")    # available vignettes
#  data(package = "Rcpp")     # available datasets

vignette("Rcpp-introduction") # Basic information on the package


# ==============
# microbenchmark - fine-grained measurement of code execution time
# ==============

if (! requireNamespace("microbenchmark", quietly=TRUE)) {
  install.packages("microbenchmark")
}

library(microbenchmark)
# Package information:
#  library(help = microbenchmark)       # basic information
#  browseVignettes("microbenchmark")    # available vignettes
#  data(package = "microbenchmark")     # available datasets


# =======
# stringi - fast string-processing functions
# =======

if (! requireNamespace("stringi", quietly=TRUE)) {
  install.packages("stringi")
}
#  library(help = stringi)       # basic information


# =======
# stringr - versatile string-processing functions built on top of stringi::
# =======

if (! requireNamespace("stringr", quietly=TRUE)) {
  install.packages("stringr")
}
#  library(help = stringr)       # basic information


# ===
# ore - alternatives to R's inbuilt regular expression libraries
# ===

if (! requireNamespace("ore", quietly=TRUE)) {
  install.packages("ore")
}
#  library(help = ore)       # basic information


# ==========
# Biostrings - from the bioinformatics community's BioConductor world
# ==========

# First: Bioconductor packages have their own package system, not CRAN. We need
# to install that first ...
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Biostrings deals with "biological" string classes ...
if (! requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
# Examine the package information:
# library(help = Biostrings)       # basic information
# browseVignettes("Biostrings")    # available vignettes
# data(package = "Biostrings")     # available datasets


# =======
# ggplot2 - a popular alternative to R's plot functions
# =======

if (! requireNamespace("ggplot2", quietly=TRUE)) {
  install.packages("ggplot2")
}
#  library(help = ggplot2)       # basic information


# ==   1.2  Update existing packages  ==========================================
#

# From time to time, update all packages. But remember, alwas answer "no" when
# asked if you want to install from source, except if you have a pressing reason
# to do that and know what you are doing  :-)
#

utils::update.packages(ask = FALSE)


# ==   1.3  Check that you have all files  =====================================

# The following files are assumed to exist in the current working directory:
#  - hello.cpp
#  - codonSplit1.cpp
#  - codonSplit2.cpp
#  - codonSplit3.cpp
# You may need to copy them there ( or temporarily change your working directory
# with  setwd("<new.working.directory>")  )

# Let's check:
if (! file.exists("hello.cpp") ||
    ! file.exists("codonSplit1.cpp") ||
    ! file.exists("codonSplit2.cpp") ||
    ! file.exists("codonSplit3.cpp")) {
  message("PANIC: missing one or more files.\n:-(\n")
} else {
  message("All good.\n:-)\n")
}


# =    2  The simplest Rcpp program  ===========================================


Rcpp::sourceCpp("hello.cpp")  # It takes a few seconds to compile the code and
                              # register the function ...
hello()                       # Once the function exists, you can call it like
                              # any other R function.

hello  # execute this to see what it looks like on the inside


# =    3  Scenario: Sequence to Codons  ========================================


# Assume we want to split a nucleotide sequence string into codons.
# This could be an ORF like ...
x1 <- "ATGCGGCATTGCAGCTGA"
x2 <- "ATGCGGCATTGCAGCTGAT"  # Note: one truncated codon at the end

# ...and it should be converted into a vector like
# c("ATG", "CGG", "CAT", "TGC", "AGC", "TGA") ... respectively
# c("ATG", "CGG", "CAT", "TGC", "AGC", "TGA", "T")


# Let's try some base-R approaches and define them as functions...
#

# ======================
# strsplit() and paste()  - after subsetting with a logical vector
# ======================
splitCodons1 <- function(x){
    codons <- character()
    y <- unlist(strsplit(x, ""))
    codons <- paste0(y[c(TRUE, FALSE, FALSE)],  # subset each nucleotide
                     y[c(FALSE, TRUE, FALSE)],  # position ...
                     y[c(FALSE, FALSE, TRUE)])  # and paste it together (1)
    return(codons)
}
splitCodons1(x1)
splitCodons1(x2)   # oops ... bitten by "vector recycling"!

# demo of how paste0 works here in principle:
paste0(c("A", "D", "H"),
       c("B", "E", "I"),
       c("C", "F", "J"))


# =======================
# strsplit() and paste0()  - after subsetting with index vectors
# =======================
splitCodons2 <- function(x){  # split and paste with seq of indices
    codons <- character()
    y <- unlist(strsplit(x, ""))
    a <- y[seq(1, length(y), by = 3)]
    b <- y[seq(2, length(y), by = 3)]
    c <- y[seq(3, length(y), by = 3)]
    codons <- paste0(a, b, c)
    return(codons)
}
splitCodons2(x1)
splitCodons2(x2)   # oops ... again: recycling appends a fantasy codon!



# ======
# gsub()  - pattern matching and substitution
# ======
splitCodons3 <- function(x) {
    codons <- character()
    y <- unlist(strsplit(x, ""))
    temp <- gsub(pattern = "([ATCG]{3})",
                 replacement = "\\1,",
                 x = x,
                 fixed = FALSE)
    codons <- unlist(strsplit(temp, split = ","))
    return(codons)
}
splitCodons3(x1)
splitCodons3(x2)   # This one gets it right.


# More ...
#
# ... what other options do we have, and how do they compare in time and memory
# use? In principle we can use methods that search for patterns, and methods
# that extract substrings.


# =======================================================
# use a "positive lookbehind regular expression" directly
# =======================================================
# Solution 1: a regular expression in base R:
strsplit(x1, "(?<=...)", perl=TRUE)[[1]]
strsplit(x2, "(?<=...)", perl=TRUE)[[1]]


# =====================================
# use stringi::stri_extract_all_regex()
# =====================================
# Solution 2: stringi offers a function to extract all substrings that match
# a pattern. The regular expression pattern "..." matches with any three
# consecutive characters:
stringi::stri_extract_all_regex(x1, '...')[[1]]
stringi::stri_extract_all_regex(x2, '...')[[1]]  # Drops the last nucleotide

# Note that we are not getting ALL substrings of length three, but ALL
# CONSECUTIVE triplets, i.e. they do not overlap.



# =====================
# use ore::ore.search()
# =====================
# Solution 2.1:
# The Oniguruma regular expression library allows to precompile regular
# expressions, and we can expect that this speeds up repeated queries
NNN <- ore::ore("...") # regular expression: three dots means "match exactly
                       #                     three characters, whatever they
                       #                     are."
ore::matches(ore::ore.search(NNN, x1, all = TRUE))
ore::matches(ore::ore.search(NNN, x2, all = TRUE))  # Drops the last nucleotide



# ======================================================================
# use substring() or stringr::str_sub() with a vector of start positions
# ======================================================================
# Solution 3:
# We could precalculate a vector of indices and then use the
# vectorized substring() or stringr::str_sub() to extract substrings
# from the indices
pos <- seq(1, nchar(x1), by = 3)
substring(x1, first=pos, last=pos+2)
stringr::str_sub(x1, start=pos, end=pos+2)

pos <- seq(1, nchar(x2), by = 3)
substring(x2, first=pos, last=pos+2)
stringr::str_sub(x2, start=pos, end=pos+2)
# both implement the requirement correctly - but require an additional vector
# that is as long as the sequence, in memory.


# using substring() ...
splitCodons4 <- function(x) {
    pos <- seq(1, nchar(x), by = 3)
    return(substring(x, first=pos, last=pos+2))
}
splitCodons4(x1)
splitCodons4(x2)

# using stringr::str_sub()
str_splitCodons5 <- function(x) {
  pos <- seq(1, nchar(x), by = 3)
  return(stringr::str_sub(x, start=pos, end=pos+2))
}
str_splitCodons5(x1)
str_splitCodons5(x2)

# ======================================================================
# use strsplit() and paste0() with a vector of start positions
# ======================================================================
# Solution 4:
# We could strsplit() into single characters and then use our pos
# vector to reassemble the characters:
pos <- seq(1, nchar(x1), by = 3)
y <- strsplit(x1, "")[[1]]
paste0(y[pos], y[pos+1], y[pos+2])

pos <- seq(1, nchar(x2), by = 3)
y <- strsplit(x2, "")[[1]]
paste0(y[pos], y[pos+1], y[pos+2]) # Two NA in the last codon - not wrong.


# ========================
# use Biostrings::codons()
# ========================
# Solution 5:
# The Biostrings package has a codons() function, but it needs to work
# on a DNAString object
as.character(Biostrings::codons(Biostrings::DNAString(x1)))
as.character(Biostrings::codons(Biostrings::DNAString(x2))) # ... with warning



# =    4  Benchmarking  ========================================================


# How do these solutions compare? How fast are they?
#

# Make a random sequence with 111 codons
cDNA <- paste(sample(c("A", "T", "C", "G"), 999, replace = TRUE), collapse="")
#pos <- seq(1, nchar(cDNA), by = 3)
NNN <- ore::ore("...")  # the ore pattern


# ==   4.1  Three ways to benchmark  ===========================================

# ================
# using Sys.time()
# ================
tStart <- Sys.time()
splitCodons1(cDNA)
tEnd <- Sys.time()
tEnd - tStart
# ... this depends on how fast you hit the keyboard ...


# ===================
# using system.time()
# ===================
#  - this measures the time spent in a function call, or executing a
#  block of code. (Nb. such a block is delimited with curly braces.)
system.time(splitCodons1(cDNA)) # too fast. Do it 10,000 times
system.time({ for (i in 1:10000) { x <- splitCodons1(cDNA) } })

# ====================
# using microbenchmark
# ====================
microbenchmark(y <- splitCodons1(cDNA))
# Note: microbenchmark() by default repeats the function call neval -times and
# reports statistics how long it took

# You can assign the results ... these are the recorded function call durations
# and then plot the distribution
bm <- microbenchmark(y <- splitCodons1(cDNA), splitCodons2(cDNA), times = 1000)

ggplot2::autoplot(bm)   # (This is a so-called "violin plot" - like a smoothed
                        # histogram that is reflected on the x-axis)

# ==   4.2  Comparing our approaches  ==========================================

microbenchmark(y <- splitCodons2(cDNA))
microbenchmark(y <- splitCodons3(cDNA))

microbenchmark(y <- strsplit(cDNA, "(?<=...)", perl=TRUE)[[1]])

microbenchmark(y <- stringi::stri_extract_all_regex(cDNA, '...')[[1]])
microbenchmark(y <- ore::matches(ore::ore.search(NNN, cDNA, all = TRUE)))

microbenchmark(y <- splitCodons4(cDNA))
microbenchmark(y <- str_splitCodons4(cDNA))

microbenchmark(y <- as.character(Biostrings::codons(Biostrings::DNAString(cDNA))))
# Note: this takes milliseconds ! Not microseconds like the others.


# In my experience, stringi::stri_extract_all_regex() is consistently the
# fastest,  it is closely followed by the str_splitCondons4() version using
# stringr::str_sub(), which requires building an additional vector of positions.

# How much speedup do we get with Rcpp?


# =    5  Three RCPP versions  =================================================


# ==   5.1  Compile C++ code and bind to an R function  ========================

file.show("codonSplit1.cpp")
Rcpp::sourceCpp("codonSplit1.cpp")  # patience ... compiling takes time
cpp_codonSplit1(x1)

file.show("codonSplit2.cpp")
Rcpp::sourceCpp("codonSplit2.cpp")
cpp_codonSplit2(x1)

file.show("codonSplit3.cpp")
Rcpp::sourceCpp("codonSplit3.cpp")
cpp_codonSplit3(x1)


# How do they compare
microbenchmark(y <- cpp_codonSplit1(cDNA))  # Slower than some base-R versions
microbenchmark(y <- cpp_codonSplit2(cDNA))  # about 15 times faster !!!
microbenchmark(y <- cpp_codonSplit3(cDNA))  # even faster

# Summary: our fastest Rcpp function is faster than any other by a factor of 3 -
# but we have to be careful how we write the C++ code. It's not a panacea - a
# naive implementation is quite slow, slower than some base R versions. Rcpp
# allows you to code closer "to the metal", but many of the std library
# functions trade convenience (and safety) for speed. Careful considertaion
# of memory use will do a lot. Whatever you do, you must profile carefully.

# Overall, we seem to have achieved a speedup of 3 over anything we can do in R,
# and a 250-fold speedup over Biostrings. Can we be even faster? Probably yes,
# working with pointers and strcpy() should be even faster.


# ==   5.2  Final comparison ...  ==============================================

bm <- microbenchmark(
     as.character(Biostrings::codons(Biostrings::DNAString(cDNA))),
     cpp_codonSplit1(cDNA),
     splitCodons3(cDNA),
     strsplit(cDNA, "(?<=...)", perl=TRUE)[[1]],
     splitCodons2(cDNA),
     ore::matches(ore::ore.search(NNN, cDNA, all = TRUE)),
     splitCodons1(cDNA),
     str_splitCodons4(cDNA),
     stringi::stri_extract_all_regex(cDNA, '...')[[1]],
     splitCodons4(cDNA),
     cpp_codonSplit2(cDNA),
     cpp_codonSplit3(cDNA),
     times = 1000)              # takes about 15 seconds

     # plot all results into a new, detached plot window
     dev.new()
     ggplot2::autoplot(bm)

     # close the new window
     dev.off()



# [END]
