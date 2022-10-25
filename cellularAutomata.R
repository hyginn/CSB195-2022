# tocID <- "cellularAutomata.R"
#
# Demo code for CSB195-2022
#
# 2022-10-12
# Boris Steipe (boris.steipe@utoronto.ca)
#
# Version:  1.1
#
# Versions:
#           1.1   Updates after demo in course
#           1.0   First course version
#
#  To Do:
#
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                 Line
#TOC> -------------------------------------
#TOC>   1        INTRODUCTION            35
#TOC>   2        FUNCTIONS               45
#TOC>   2.1        showRule()            55
#TOC>   2.2        applyRule()           93
#TOC>   2.3        CA()                 115
#TOC>   3        FIRST STEPS            215
#TOC>   4        EXPLORATIONS           252
#TOC> 
#TOC> ==========================================================================


# =    1  INTRODUCTION  ========================================================
#


# The required functions are defined below, you are welcome to look at them in
# detail, but there is a lot of bookkeeping going on which may obscure that they
# are really very simple. Just be sure to understand what they do in principle.
# You can always use the help-function to recall the meaning of the parameters
# and the function details.

# =    2  FUNCTIONS  ===========================================================
#

# These function are sourced from outside this script because they are needed in
# other scripts later on.

source("imPlot.R")  # plot the contents of a matrix
source("fibWord.R") # compute a Fibonacci word as initialization vector


# ==   2.1  showRule()  ========================================================

showRule <- function(R) {
  #' @title showRule()
  #' @description Print out the rule-set for Wolfram code R
  #' @param R  int.  The Wolfram code for a 1D cellular automaton rule-set;
  #'                 an integer in [0, 255].
  #'
  #' @return   int.  The rule set as a named integer vector. The rule set is
  #'                 printed as a side-effect.
  #'
  #' @details `showRule()` prints out the construction of a cellular automaton
  #'          rule-set computed from its Wolfram code. This code is derived from
  #'          the state-configurations of a cell's neighborhood, ordered in
  #'          descending order - in our 1D case: [111][110][101]...[000].
  #'          Each of the eight configurations can be seen as a digit in
  #'          an eight-digit binary number, and thus each of the 255 numbers
  #'          specifies a unique rule set.
  #'
  #' @examples   showRule(18)
  #'             myR <- showRule(110)
  #'             myR["010"]
  if (R < 0 || R > 255) {
    stop("R is not in [0,255].")
  }

  v <- as.integer(intToBits(R)[8:1] == "01")  # convert integer to binary array
  names(v) <- c("111", "110", "101", "100", "011", "010", "001", "000")

  cat(sprintf("\n  Rule: %d\n  %s\n   %s  ",
              R,
              paste(names(v), collapse = " "),
              paste0(v, collapse = "   ")))

  return(invisible(v))
}


# ==   2.2  applyRule()  =======================================================

applyRule <- function(R, s) {
  #' @title applyRule()
  #' @description The new state of a cellular automaton is computed by
  #' applying the Wolfram code rule-set R to the neighbourhood state s.
  #'
  #' @param  R   int. The Wolfram code for a 1D cellular automaton
  #'                      rule set; an integer in [0, 255].
  #' @param  s   int. A binary vector of length three.
  #' @return      0 or 1
  #' @examples
  #'   applyRule(110,  c(1, 0, 1))
  #'   applyRule( 30, c(1, 0, 0))

  v <- as.integer(intToBits(R)[8:1] == "01")  # convert integer to binary array
  idx <-  8 - (s[1]*4) - (s[2]*2) - (s[3]*1)  # compute index of v from s

  return(v[idx])
}


# ==   2.3  CA()  ==============================================================

CA <- function(R, nx = 100, ny = 100, vInit, wrap = TRUE, N = NULL) {
  #' @title CA()
  #' @description Evolve a cellular automaton with Wolfram code R in a nx * ny
  #' matrix.
  #'
  #' @param  R   int. The Wolfram code for a 1D cellular automaton
  #'                      rule set; an integer in [0, 255].
  #' @param  nx  the width of the space in which the automaton evolves
  #' @param  ny  the number of steps to compute
  #' @param  vInit   the initialization vector. See details.
  #' @param  wrap    if TRUE, wrap the vertical edges of the matrix (periodic
  #' boundary).
  #' @param  N    If this is not NULL, the function plots every updated row
  #' until the matrix is filled, then scrolls the matrix for a total of N steps
  #' (or until the process is aborted.).
  #' @return      the evolved matrix as 0s and 1s.
  #' @details If vInit is missing, a single 1 is placed at the centre of the
  #' first row. If vInit is in [0, 1] the vector is filled with random 1s of
  #' p == vInit. If vInit is > 1, round(vInit) 1s are equally spaced into
  #' vInit. If vInit is the keyword "fib", the first row is initialized with
  #' a Fibonacci word. If length(vInit) > 1 the vector is recycled into the
  #' length of vInit. In particular, this allows to pass in a custom
  #' initialization vector of length nx.
  #' @examples
  #'

  if (missing(vInit)) { # place a single 1 in the middle
    vInit <- integer(nx)
    vInit[round((nx) / 2)] <- 1
  } else if (length(vInit) > 1){  # recycle this vector
    vInit <- rep(vInit, ceiling(nx/length(vInit)))[1:nx]
  } else if (is.numeric(vInit) && vInit >= 0 && vInit <= 1) { # probability
    vInit <- sample(0:1, nx, replace= TRUE, prob = c(1-vInit, vInit))
  } else if (is.numeric(vInit) && vInit > 1) { # place 1s
    idx <- round(seq(1, nx, length.out = round(vInit)))
    vInit <- integer(nx)
    vInit[idx] <- 1
  } else if (vInit == "fib") { # Fibonacci word
    vInit <- as.integer(unlist(strsplit(fibWord(nx), "")))
  } else {
    stop("PANIC: I don't know what to do with this vInit.")
  }

  if (is.null(N)) {  # No plotting of progress
    N <- 0
  } else {
    Nmax <- N
    cat("\n")
  }

  stopifnot(nx == length(vInit))

  m <- matrix(integer(nx * ny), ncol = nx)
  m[1, ] <- vInit

  if (wrap) {
    m <- cbind(m[ ,nx], m, m[ ,1]) # wrap first and last column
    nx <- nx + 2                  # increase nx accordingly
  }


  for (i in 1:(ny-1)) {
    for (j in 2:(nx-1)) {
      m[(i+1), j] <- applyRule(R, m[i, (j-1):(j+1)])
    }
    if (wrap) {
      m[i+1, 1] <- m[i+1, j]     # wrap last cell to beginning
      m[i+1, nx] <- m[i+1, 2]    # wrap first cell to end
    }
    if (N > 0) {
      imPlot(m)
      dev.flush()
      N <- N-1
      cat(sprintf("\rStep %d", Nmax - N))
    }
  }

  for (k in seq_len(N)) { # scroll if N is > 0
    m[1:(ny-1), ] <- m[2:ny, ] # shift the matrix one row up
    for (j in 2:(nx-1)) {
      m[(i+1), j] <- applyRule(R, m[i, (j-1):(j+1)])
    }
    m[i+1, 1] <- m[i+1, nx-1]  # wrap last cell to beginning
    m[i+1, nx] <- m[i+1, 2]    # wrap first cell to end
    imPlot(m)
    dev.flush()
    cat(sprintf("\rStep %d", Nmax - (N - k)))

  }

  if (wrap) {
    m <- m[ , -c(1, nx)] # remove wrapped columns
  }

  return(m)
}


# =    3  FIRST STEPS  =========================================================
#

# Let's evolve a few cellular automata.

nx <- 11  # width
ny <- 9   # height
m <- matrix(integer(nx * ny), nrow = ny)  # integer() fills the matrix with 0s

m[1,6] <- 1  # Initialize the first row

mTest <- m  # save this so we can eaily recreate it.

image(m)                                  # R's basic plotting function image()

plot(c(1, nx), c(1,ny), type = "n")       # R's plotting function rasterImage()
rasterImage(m, 1, 1, nx, ny, interpolate = FALSE)

imPlot(m, drawGrid = TRUE)                # Our own imPlot()

# Now look at the cells row by row and evolve a CA. The readlines() function
# pauses every row and waits for your input. Just hit <return>  in the
# console.

R <- 213         # Wofram Code 18
showRule(R)
for (i in 1:(ny-1)) {
  for (j in 2:(nx-1)) {
    m[(i+1), j] <- applyRule(R, m[i, (j-1):(j+1)])
  }
  imPlot(m, drawGrid = TRUE)
  readline("Hit <Return> for next row >")
}

m <- mStart # Reset m and try a different rule - perhaps 213 ...


# =    4  EXPLORATIONS  ========================================================
#

# == First exploration.
# Set i to a particular value. Select the entire if()
# statement. When you hit <cmd>-<return>, i is incremented by one, and
# the CA is run. You can step through the different rules one-by-one simply by
# hitting <cmd>-<return> again and again.
i <- 0
if ((i <- i+1)) {
  m <- CA(i, nx = 150, ny=300)
  imPlot(m, drawGrid = FALSE)
  showRule(i)
}



# == Second exploration.
# Look at the effect of varying initializations.

R <- 18
showRule(R)

# A single starting cell
m <- CA(R, nx = 200, ny=400); imPlot(m, drawGrid = FALSE)

# Three starting cells
m <- CA(R, nx = 200, ny=400, vInit = 3); imPlot(m, drawGrid = FALSE)

# Seven starting cells
m <- CA(R, nx = 200, ny=400, vInit = 7); imPlot(m, drawGrid = FALSE)

# random 1s with p == 0.33
m <- CA(R, nx = 200, ny=400, vInit = 0.33); imPlot(m, drawGrid = FALSE)

# random 1s with p == 0.67
m <- CA(R, nx = 200, ny=400, vInit = 0.67); imPlot(m, drawGrid = FALSE)

# repeating a vector
v <- c(rep(0, 13), rep(1, 13), rep(0, 21), rep(1, 21))
m <- CA(R, nx = 200, ny=400, vInit = v); imPlot(m, drawGrid = FALSE)

# Fibonacci word
m <- CA(R, nx = 200, ny=400, vInit = "fib"); imPlot(m, drawGrid = FALSE)



# == Third exploration.
# Look at every CA.

# This is a loop that plots A LOT. The default plot window doesn't like to be
# updated so frequently. It "buffers" output untiul the process ends. I am sure
# this has a use-case - and perhaps it makes the display more responsive in some
# cases, but for what we are doing here, it breaks the process. That's bad.

# So we'll open a new graphics window. This is reasonable responsive on the Mac,
# but I don't know how well it works on Windows. Let me know.
dev.new()

myCols <- sample(hcl.colors(256, "Temps"))  # random colors

# Also, there is a problem with keeping a little pause here. We would like to be
# able to abort the loop by pressing the stop-sign. However if we use
# Sys.sleep() while running a CPU-intensive process, our interrupts get lost.
# Here is a safe alternative. (It's not really important though ... just ignore
# if you find this confusing..):

safeSleep <- function(tSleep) {
  # alternative to Sys.sleep().
  # cf. https://stackoverflow.com/questions/1174799/how-to-make-execution-pause-sleep-wait-for-x-seconds-in-r
  then <-Sys.time()
  while((as.numeric(Sys.time()) - as.numeric(then)) < tSleep){} #dummy while loop
}



# I'll run this here with "fib" initialization but you can choose any of the
# other modes ...

# ... and just so that we'll not be looking at the same rules again and again,
# lets count down from 255.

for (i in 255:0) {
  m <- CA(i, nx=300, ny=600, vInit = sample(0:1, 300, replace = TRUE))
  showRule(i)
  cat("\n--------------------------------------\n")
  imPlot(m, colMap = c("0" = "#f3f3f3", "1"=myCols[i+1]))
  safeSleep(0.5)  # short pause
}

# == Fourth exploration:

# Scrolling the CA. We may be interested in letingt this run for a long time, to
# see whether the CA ultimately settles into a repeating state. I added a
# scrolling mechanism to the CA function ...

dev.off() # remove the alternative device ...
dev.new() # and make a fresh one.

# Try this ...

m <- CA(18, nx=200, ny=400, vInit =7, N = 4000)

# Remember: you can hit the red "Stop Sign" to stop the scrolling process.

# So. Have fun exploring. And note down if you find anything remarkable, or if
# you think this code doesn't support something you would like to try.



# [END]
