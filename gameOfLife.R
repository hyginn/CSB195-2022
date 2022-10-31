# tocID <- "gameOfLife.R"
#
# Demo code for CSB195-2022
#
# 2022-10-23
# Boris Steipe (boris.steipe@utoronto.ca)
#
# Version:  1.0
#
# ==============================================================================
#

# == 1. Introduction

# == 1. Structures and Functions

source("imPlot.R")
source("seedGOL.R")


# == wrap() ====================================================================

wrap <- function(w) {
# wrap a 2D matrix of size x*y into a torus (x-2) * (y-2)
#
  w[1, ]       <- w[nrow(w)-1, ]    # wrap last actual row to row 1
  w[nrow(w), ] <- w[2,         ]    # wrap row 2 to last row
  w[ ,1]       <- w[ ,ncol(w)-1]    # wrap last actual column to column 1
  w[ ,ncol(w)] <- w[ ,2        ]    # wrap column 2 to last column

  return(w)
}


# == runGOL() ==================================================================
runGOL <- function(w0, rL, rD, nTick = 100) {

  if (missing(rL)) {
    # Neighbours:    0  1  2  3  4  5  6  7  8
    rL <-          c(0, 0, 1, 1, 0, 0, 0, 0, 0)
  }
  if (missing(rD)) {
    # Neighbours:    0  1  2  3  4  5  6  7  8
    rD <-          c(0, 0, 0, 1, 0, 0, 0, 0, 0)
  }

  w0 <- wrap(w0)
  w1 <- w0

  for (tick in seq_len(nTick)) {    # for each tick
    pBar(tick, nTick)               # progress bar

    for (i in 2:(ncol(w0)-1)) {     # for each actual row
      for (j in 2:(nrow(w0)-1)) {   # for each actual column

        # sum of neighbours in 3 x 3 neighbourhood - self
        nNeigh <- sum(w0[(i-1):(i+1), (j-1):(j+1)]) - w0[i, j]

        # Rules
        if (w0[i, j] == 1) { w1[i, j] <- rL[nNeigh + 1]
        } else             { w1[i, j] <- rD[nNeigh + 1] }
 + 1
      } # end column
    } # end row
    w1 <- wrap(w1)

    imPlot(w1[2:(nrow(w1)-1), 2:(ncol(w1)-1)], drawGrid = TRUE)
    w0 <- w1

    safeSleep(0.1)
    dev.flush()

  } # end tick

  return(invisible(w1))
}




if (FALSE) {

  # initialize
  # myWorld <- matrix(integer(X * Y), ncol = X)
  # myWorld[2,4] <- 1    # edge dot
  # myWorld[4,3:5] <- 1  # oscillator
  # myWorld[5:7,7]<-1;myWorld[6,5]<-1;myWorld[7,6]<-1   # glider

  # myWorld <- matrix(integer(X * Y), ncol = X)
  # myWorld <- placeCell("GosperGun.cell", myWorld, 4, 20)

  # myWorld <- matrix(integer(X * Y), ncol = X)
  # myWorld <- placeLif("maxFill.lif", myWorld, 30, 30)

  # myWorld <- matrix(integer(X * Y), ncol = X)
  # myWorld <- placeCell("pentomino.cell", myWorld, 50, 50)


  # Initialize dimensions and probabilities
  X <- 100
  Y <- 100
  p <- 0.5
  dev.new()

  myWorld <- sample(0:1, X * Y, prob=c(1-p, p), replace = TRUE)
  dim(myWorld) <- c(Y,X)


  imPlot(myWorld, drawGrid = TRUE)

  newWorld <- runGOL(myWorld, nTick = 1000)

}

