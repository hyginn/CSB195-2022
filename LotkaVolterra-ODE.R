# tocID <- "LotkaVolterra-ODE.R"
#
# Purpose:  Demonstrate the Lotka-Volterra equation as an ODE.
#
# Version:  2.2
# Date:     2022-10-18
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           2.2  Major refactoring to properly abstract out the equations
#                  from the code that just runs the trajectory
#           2.1  Move plotting code to functions
#           2.0  Add manual iteration code
#           1.0  First demo in tutorial
#
# ToDo:
#
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                     Line
#TOC> ---------------------------------------------------------
#TOC>   1        PACKAGES                                    52
#TOC>   2        FUNCTIONS AND GLOBALS                       72
#TOC>   3        deSolve::ode()                             216
#TOC>   4        EULER'S METHOD                             241
#TOC>   5        ADDING NOISE                               355
#TOC>   5.1        What does added noise look like?         362
#TOC>   5.2        Modifying our parameters                 373
#TOC> 
#TOC> ==========================================================================


#   This code has its origins in a blog entry in "Mind of a Markov Chain" (2010)
#   that is no longer public. It was preserved in the R-Bloggers aggregation.
#   (Shout out to R-Bloggers: https://www.r-bloggers.com/). The link is
#   https://www.r-bloggers.com/2010/03/lotka-volterra-model%C2%A0%C2%A0intro/
#   The author of the original post is "apeescape". Apparently that post
#   (uncredited) inspired Green and Shou (2014) to a very much more detailed
#   exposition in a Methods in Molecular Biology chapter: "Modeling community
#   population dynamics with the open-source language R."
#   https://pubmed.ncbi.nlm.nih.gov/24838889/  (pdf available.)

# ==============================================================================



# =    1  PACKAGES  ============================================================
# Check that required packages have been installed. Install if needed.

if (! requireNamespace("deSolve", quietly=TRUE)) {
  install.packages("deSolve")
}
# Package information:
#  library(help = deSolve)       # basic information
#  browseVignettes("deSolve")    # available vignettes
#  data(package = "deSolve")     # available datasets

if (! requireNamespace("plotrix", quietly=TRUE)) {
  install.packages("plotrix")
}
# Package information:
#  library(help = plotrix)       # basic information
#  browseVignettes("plotrix")    # available vignettes
#  data(package = "plotrix")     # available datasets


# =    2  FUNCTIONS AND GLOBALS  ===============================================


LVmod <- function(t, state, param) {
  #' LVmod - Lotka-Volterra Model
  #'
  #' @param t      num  A vector of time points
  #' @param state  num  A named vector with starting values for nPrey and nPred.
  #' @param param  num  A named vector containing the parameters of the model:
  #'                    preyGrowthRate, predationRate, predatorGrowthRate, and
  #'                    predatorLossRate.
  #'
  #' @return num        A list of population changes "dPrey" and "dPred".
  #'
  #' @details  This function is meant to be used by the deSolve::ode() ODE
  #'           solver to compute prey and predator populations over time.
  #'

  with(as.list(c(state, param)), {
    dPrey = (preyGrowthRate * nPrey) - (predationRate * nPrey * nPred)
    dPred = (predatorGrowthRate * nPrey * nPred) - (predatorLossRate * nPred)
    return(list(c(dPrey, dPred)))
  })
}

PREYCOL <- "#008080"  # Color for prey in plots: "teal"
PREDCOL <- "#800000"  # Color for predators in plot: "maroon"


timePlot <- function(dat,
                     ix = 1,
                     iy = c(2, 3),
                     col = c(PREYCOL, PREDCOL),
                     main = "") {
  #' timePlot - Time-series Plot
  #'
  #' @param  dat  num  a numeric matrix that contains x-axis values (time) and
  #'           `       population values.
  #' @param  ix   int  column index of time values. Default: 1
  #' @param  iy   int  column indices of population values. Default: c(2, 3)
  #' @param  col  chr  color values. Default c(PREYCOL, PREDCOL).
  #' @param  main chr  plot title. Default "".
  #'
  #' @return           NULL. The function plots a time series.
  #'
  #' @details          This is a standard time-series plot but written so that
  #'                   it can plot a variable number of columns. The column
  #'                   indices for the values to be plotted along the y-axis
  #'                   are passed in the variable iy. Axis labels are taken
  #'                   from the colnames of the matrix that contains the values
  #'                   to be plotted. First, an empty frame is plotted with
  #'                   enough space to accommodate the values, then the
  #'                   lines() function is used to plot the value of each
  #'                   column in turn.

  # empty plot frame
  plot(dat[, ix], dat[, iy[1]],
       type = "n",
       ylim = c(0, max(dat[, iy] * 1.3)),
       xlab = colnames(dat)[ix],
       ylab = paste(colnames(dat)[iy], collapse = ", "),
       main = main)

  # values
  for (i in 1:length(iy)) {
    lines(dat[, ix], dat[, iy[i]],
          type = "l",
          col = col[i])
  }

  # add a legend
  legend("topright",
         legend = colnames(dat)[iy],
         lty = rep(1, length(iy)),
         col = col,
         bty = "n")

  return(invisible(NULL))
}


phasePlot <- function(dat,
                      ix = 2,
                      iy = 3,
                      col = c(PREYCOL, PREDCOL),
                      main = "") {
  #' phasePlot - Phase-space Plot
  #'
  #' @param  dat  num  a numeric matrix that contains population values
  #'                   in columns ix and iy.
  #' @param  ix   int  column index of first population. Default: 2
  #' @param  iy   int  column index of second population. Default: 3
  #' @param  col  chr  color values. Default c(PREYCOL, PREDCOL).
  #' @param  main chr  plot title. Default "".
  #'
  #' @return           NULL. The function plots a phase-space trajectory.
  #'
  #' @details          This is a standard phase-space plot which plots two
  #'                   populations from a Lotka-Volterra plot or similar against
  #'                   each other. However it uses a colour gradient to plot the
  #'                   trajectory. The colour gradient is
  #'                   determined along the log-ratio of the two populations.
  #'                   First, an empty frame is plotted, then the
  #'                   plotrix::color.scale.line() function plots the
  #'                   trajectory.


  # compute a color vector in 50 steps

  # 50 colors from col[1] to col[2]
  makeMyPal <- colorRampPalette(c(col[1], col[2]))
  myPal <- makeMyPal(50)

  # utility function ...
  normalize <- function(v) {
    i <- min(v, na.rm = TRUE)
    a <- max(v, na.rm = TRUE)
    v <- (v - i) / (a - i)
    return(v)
  }

  # Fetch each row's color based on the log ratio of dat[ , ix] and dat[ , iy]
  x <- dat[ , ix] / dat[ , iy]      # ratio
  x <- log(x)                       # log ratio
  x <- round(normalize(x) * 49) + 1 # scale from 1 to 50
  col <- myPal[x]                   # assign one of the 50 colour values

  # empty plot frame
  plot(dat[, ix], dat[, iy],
       type = "n",
       xlab = colnames(dat)[ix],
       ylab = colnames(dat)[iy],
       main = main)

  # plot the phase-space trajectory
  plotrix::color.scale.lines(dat[ , ix],
                             dat[ , iy],
                             col = col,
                             lwd = 3)

  return(invisible(NULL))
}


# =    3  deSolve::ode()  ======================================================

# Define the parameters
Par <- c(preyGrowthRate = 1,
         predationRate = 0.2,
         predatorGrowthRate = 0.08,
         predatorLossRate = 0.5)

# Define the initial state populations
State <- c(nPrey = 5, nPred = 2)

# Initialize the time-axis
Time <- seq(0, 100, by = 0.01)

# Run the ODE solver
out <- as.data.frame(deSolve::ode(func  = LVmod,
                                  y     = State,
                                  parms = Par,
                                  times = Time))

timePlot(out, main = "Lotka-Volterra Model (ODE)")

phasePlot(out, main = "Phase-Space of Lotka-Volterra Model (ODE)")


# =    4  EULER'S METHOD  ======================================================

# See also Green and Shou (2014)

# Let us rewrite our code to perform exactly the same computation of the coupled
# differential equations of the Lotka-Volterra Model above (up to numeric
# error), but step-by-step, from a starting condition. This is a method
# published by Leonhard Euler in 1768, thus it is known as "Euler's Method". The
# ODE solver we used above has a number of modern methods available, which
# generally result in a much smaller error. But calculating the process
# "manually" as we do below, makes it easier to explore different functions for
# predator and prey population interactions.

# To be more flexible in our explorations, we tease apart the computational
# steps more cleanly. For example, our parameters concern values that we
# need for running the function over time, parameters that concern only
# the behaviour of prey, and parameters that concern only the predators. We
# can abstract this out into three different sets of parameters.

myPar <- list(nPreyStart = 5,
              nPredStart = 2,
              Time = 100,
              dt = 0.01)

myParPrey <- list(preyGrowthRate = 1,
                  predationRate = 0.2)

myParPred <- list(predatorGrowthRate = 0.08,
                  predatorLossRate = 0.5)

# Next we define a function that does the actual calculation. It takes
# as parameters the values we defined above, and two functions! A function
# in R is an object as any other, it can be stored in lists, passed as a
# parameter, assigned to names etc. etc. Abstracting the core equations
# out of our computation makes it easy to change the functions later
# and experiment with different ways to compute interactions. We just call our
# computation with different functions and we don't need to copy the rest
# of the code again and again. This ensures that our computations remain
# consistent!

LVEuler <- function(par = myPar,
                    parPrey = myParPrey,
                    parPred = myParPred,
                    fPrey,
                    fPred) {

  N <- par$Time / par$dt                         # Number of steps to compute
  outM <-  matrix(numeric(N * 3), ncol = 3)      # Set up the output matrix
  colnames(outM) <- c("time", "nPrey", "nPred")

  nPrey <- par$nPreyStart                        # Initialize the starting
  nPred <- par$nPredStart                        # values

  outM[1, "time"]  <- 0.0                        # Store the starting values
  outM[1, "nPrey"] <- nPrey
  outM[1, "nPred"] <- nPred

  # Iterate ...
  for (i in 2:N) {                               # for all time points

    # compute the population change
    dPrey = fPrey(nSelf = nPrey, nOther = nPred, par = parPrey) * dt
    dPred = fPred(nSelf = nPred, nOther = nPrey, par = parPred) * dt

    # compute the new populations
    newPrey <- nPrey + dPrey
    newPred <- nPred + dPred

    # store the results
    outM[i, ] <- c((i-1) * dt, newPrey, newPred)

    # prepare for the next iteration
    nPrey <- newPrey
    nPred <- newPred
  }

  # Done: return the result
  return(outM)
}

# The equations that compute the population change
# ... for Prey
fLVPrey <- function(nSelf, nOther, par) {
  dSelf <- (par$preyGrowthRate * nSelf) -
           (par$predationRate  * nSelf * nOther)
  return(dSelf)
}

# ... for Predators
fLVPred <- function(nSelf, nOther, par) {
  dSelf <- - (par$predatorLossRate   * nSelf) +
             (par$predatorGrowthRate * nSelf * nOther)
  return(dSelf)
}

# With this "refactoring" pf our code, we clearly distinguish those parts of the
# computation that implement the actual model, and those parts that we only need
# to run the model along a trajectory. This is an example of an important
# objective in software engineering: "separation of concerns".

# == run this ...
outE <- LVEuler(fPrey = fLVPrey,  # Note: no parentheses! We are passing the
                fPred = fLVPred)  # function itself, not the result of
                                  # running the function. And: all other
                                  # parameters are passed via the default
                                  # values

# plot the results ...
timePlot(outE, main = "Lotka-Volterra Model (Euler)")

phasePlot(outE, main = "Phase-Space of Lotka-Volterra Model (Euler)")



# =    5  ADDING NOISE  ========================================================
# Now, let's reap the benefits of the increased abstraction of our code. Let's
# add a bit of noise to the predation rate. There are many ways to define
# "noise". Let's use a Gaussian model that increases or decreases the value of
# the predation rate by +- 1/3 of its value.


# ==   5.1  What does added noise look like?  ==================================
myParPrey$predationRate   # 0.2

x <- rnorm(10000,          # What does this distribution look like?
           mean = myParPrey$predationRate,
           sd = myParPrey$predationRate / 3)
head(x)
hist(x, breaks = 50)
abline(v = myParPrey$predationRate, col = "#AA0000", lwd = 3)


# ==   5.2  Modifying our parameters  ==========================================

# Here is the modified function:
fLVPreyStoch <- function(nSelf, nOther, par) {
  R <- par$predationRate
  R <- rnorm(1, mean = R, sd = R / 3)
  dSelf <- (par$preyGrowthRate * nSelf) - (R * nSelf * nOther)
  return(dSelf)
}

# plot
outS <- LVEuler(fPrey = fLVPreyStoch, fPred = fLVPred)
timePlot(outS, main = "Lotka-Volterra Model (Euler)")
phasePlot(outS, main = "Stochastic MOdel I")

# Add more stochasticity: 33% of noise on the predator loss rate
fLVPredStoch <- function(nSelf, nOther, par) {
  R <- par$predatorLossRate
  R <- rnorm(1, mean = R, sd = R / 3)
  dSelf <- - ( R  * nSelf) + (par$predatorGrowthRate * nSelf * nOther)
  return(dSelf)
}

outS <- LVEuler(fPrey = fLVPreyStoch, fPred = fLVPredStoch)
timePlot(outS, main = "Lotka-Volterra Model (Euler)")
phasePlot(outS, main = "Stochastic Model II")


# 50% of noise on all rates
fLVPreyStoch <- function(nSelf, nOther, par) {
  R1 <- par$predationRate
  R1 <- rnorm(1, mean = R1, sd = R1 / 2)
  R2 <- par$preyGrowthRate
  R2 <- rnorm(1, mean = R2, sd = R2 / 2)
  dSelf <- (R2 * nSelf) - (R2 * nSelf * nOther)
  return(dSelf)
}

fLVPredStoch <- function(nSelf, nOther, par) {
  R1 <- par$predatorLossRate
  R1 <- rnorm(1, mean = R1, sd = R1 / 2)
  R2 <- par$predatorGrowthRate
  R2 <- rnorm(1, mean = R2, sd = R2 / 2)
  dSelf <- - ( R1  * nSelf) + (R2 * nSelf * nOther)
  return(dSelf)
}

outS <- LVEuler(fPrey = fLVPreyStoch, fPred = fLVPredStoch)
timePlot(outS, main = "Lotka-Volterra Model (Euler)")
phasePlot(outS, main = "Stochastic Model II")


# Now we can start getting more creative with our model...



# ====  TESTS  =================================================================

if (FALSE) {





}

# [END]
