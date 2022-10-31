# tocID <- "LotkaVolterra-Model.R"
#
# Purpose:  Explorations of the Lotka-Volterra equation as coupled ODEs.
#
# Version:  3.0
# Date:     2022-10-21
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           3.0  3D relations. Add selection mechanism to plot functions.
#                  Demonstrate strange attractor.
#           2.2  Major refactoring to properly abstract out the equations
#                  from the code that just runs the trajectory
#           2.1  Move plotting code to functions
#           2.0  Add manual iteration code
#           1.0  First demo in tutorial
#
# ToDo:
#   Change variable plot density to second(!) derivative
#   Time-series analysis?
#
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                            Line
#TOC> ----------------------------------------------------------------
#TOC>   1        INTRODUCTION                                       63
#TOC>   2        PACKAGES                                          103
#TOC>   3        FUNCTIONS AND GLOBALS                             137
#TOC>   3.1        LVmod()                                         149
#TOC>   3.2        timePlot()                                      173
#TOC>   3.3        selP()                                          238
#TOC>   3.4        phasePlot()                                     270
#TOC>   4        Lotka-Volterra model via deSolve::ode()           354
#TOC>   5        Lotka-Volterra model via Euler's Method           381
#TOC>   6        ADDING NOISE                                      499
#TOC>   6.1        What does added noise look like?                506
#TOC>   6.2        Modifying our parameters                        517
#TOC>   7        Perturbations via deSolve::ode                    584
#TOC>   8        Going 3D: Chaos                                   626
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
#   https://pubmed.ncbi.nlm.nih.gov/24838889/  (pdf available.). The
#   three-species model with a strange attractor was published by Michael Gilpin
#   (1979) Spiral Chaos in a Predator Prey Model. The American Naturalist.
#   113(2): 306-308. https://www.jstor.org/stable/pdf/2460209.pdf

# ==============================================================================

# =    1  INTRODUCTION  ========================================================

# Combining change and resistance to change are hallmarks of life, as it evolves
# at the boundary between static order and chaotic dissipation. Models of change
# are usually expressed as differential equations: mathematical expressions that
# calculate the change of a state over time. Such "state" can be a feature like
# length or weight, the size of a population, or the mean of a distance -
# depending on the question at hand. Such equations are typically non-linear and
# incorporate geometric growth terms (i.e. something increases by a percentage
# of the current value at each time step, rather than by a constant amount).
# Whenever two or more different entities interact, and their respective values
# appear in the same equation, we speak of a "coupled" differential equation,
# and even though these are simple systems, they are generall not solvable by
# analytic methods but need numeric methods (i.e. we can't solve them by
# applying a formula but need to compute change over time), and they can display
# intriguingly complex behaviour.

# One particular model of two coupled non-linear differential equations are the
# Lotka-Volterra equations, named after Alexander Lotka and Vito Volterra who
# independently discovered their use to describe prdator-prey interactions (cf.
# https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations). The basic
# equations are:
#    dx/dt = ax - bxy
#    dy/dt = gxy - hy
# with x and y as the prey and predator populations and a,b,g,h positive,
# real parameters.

# In this exploration:
#  - we solve the equations and analyze the population numbers over time;
#  - we introduce time-series plots and phase space plots;
#  - we perform the calculation in a stepwise approximation (Euler's method)
#      and explore the effect of random perturbations;
#  - we show that a systematic perturbation puts the system into an altered,
#      but still stable configuration.
#  - we show how a simple three-species model evolves in a "strange" attractor
#      basin and has chaotic behaviour;
#  - and we plot the trajectory in 3D.



# =    2  PACKAGES  ============================================================
# Check that required packages have been installed. Install if needed. Get
# information about the package contents


# A package to compute ODEs ...
if (! requireNamespace("deSolve", quietly=TRUE)) {
  install.packages("deSolve")
}
# Package information:
#  library(help   = deSolve)     # basic information
#  browseVignettes("deSolve")    # available vignettes
#  data(package  = "deSolve")    # available datasets


# A package to support plotting. We need this later to plot lines with a colour
# gradient which we can't do in base R.
if (! requireNamespace("plotrix", quietly=TRUE)) {
  install.packages("plotrix")
}
# Package information:
#  library(help   = plotrix)     # basic information
#  data(package  = "plotrix")    # available datasets


# A package to support plotting 3D plots ...
if (! requireNamespace("plot3D", quietly = TRUE)) {
  install.packages("plot3D")
}
#  library(help   = plot3D)     # basic information
#  data(package  = "plot3D")    # available datasets



# =    3  FUNCTIONS AND GLOBALS  ===============================================


# Color definitions for a consistent appearance of plots
PREYCOL <-  "#04ccc2"  # Color for prey in plots: "oceanic"
PREY2COL <- "#427deb"  # A second prey species: "seljuk blue"
PREDCOL <-  "#c41b65"  # Color for predators in plot: "bright rose"

# check:
# barplot(c(1, 1, 1), axes = FALSE, col = c(PREYCOL, PREY2COL, PREDCOL))


# ==   3.1  LVmod()  ===========================================================
LVmod <- function(t, state, param) {
  #' LVmod - Lotka-Volterra Model for ODE solver
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

# A function to plot time-series
# ==   3.2  timePlot()  ========================================================
timePlot <- function(dat,
                     ix = 1,
                     iy = c(2, 3),
                     col = c(PREYCOL, PREDCOL),
                     main = "",
                     pMax = min(2000, nrow(dat))) {
  #' timePlot - Time-series Plot
  #'
  #' @param  dat  num  a numeric matrix that contains x-axis values (time) and
  #'           `       population values.
  #' @param  ix   int  column index of time values. Default: 1
  #' @param  iy   int  column indices of population values. Default: c(2, 3)
  #' @param  col  chr  color values. Default c(PREYCOL, PREDCOL).
  #' @param  main chr  plot title. Default "".
  #' @param  pMax num  The maximum number of points to plot. Some time-series
  #'                   can be quite long and the plots can take a long time
  #'                   to plot and update. But we generally do not need to
  #'                   plot at the full resolution of the time-series. However
  #'                   just plotting every n-th point or so creates plotting
  #'                   artefacts. Instead, we space points by distance along the
  #'                   plotted line. Default is the min() of 2000 points or
  #'                   the number of rows in the dataset.
  #'
  #' @return           NULL. The function plots a time series.
  #'
  #' @details          This is a standard time-series plot but written so that
  #'                   it can plot a variable number of columns, and use only
  #'                   a subset of equally spaced points. The column
  #'                   indices for the values to be plotted along the y-axis
  #'                   are passed in the variable iy. Axis labels are taken
  #'                   from the colnames of the matrix that contains the values
  #'                   to be plotted. First, an empty frame is plotted with
  #'                   enough space to accommodate the values, then the
  #'                   lines() function is used to plot the value of each
  #'                   column in turn. Finally, a legend is added.

  # empty plot frame
  plot(dat[, ix], dat[, iy[1]],
       type = "n",
       ylim = c(0, max(dat[, iy] * 1.3)),
       xlab = colnames(dat)[ix],
       ylab = paste(colnames(dat)[iy], collapse = ", "),
       main = main)

  # values
  for (i in 1:length(iy)) {
    sel <- selP(dat[,c(ix, iy[i])], pMax)
    lines(dat[sel, ix], dat[sel, iy[i]],
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


# A function to subset large datasets for plotting
# ==   3.3  selP()  ============================================================

selP <- function(x, N) {
  #' selP - select points for a line-plot
  #'
  #' @param x  num
  #' @param N  Number of points to be selected.
  #'
  #' @details  A vector of indices is produced that selects from a vector or
  #'           matrix, according to the Euclidian distance between two
  #'           consecutive elements.
  #'
  #' @return  a vector of N row-indices for subsetting x according to the
  #'          distance between consecutive rows. The first and last points
  #'          are always included.

  if (is.null(nrow(x))) { x <- as.matrix(x, byrow = TRUE) }
  l <- nrow(x)                                    # number of rows
  d <- sqrt(rowSums((x[2:l, ] - x[1:(l-1), ])^2)) # distances
  d <- c(0, d)                                    # add back first element

  # compute the cumulative sums, cut them into (N-1) intervals, get the
  # row indices and add the last element too.
  sel <- cumsum(d)
  sel <- cut(sel, breaks = (N-1), include.lowest=TRUE, labels = FALSE)
  sel <- c(which(! duplicated(sel)), l)

  return(sel)
}


# A function to plot trajectories in phase-space
# ==   3.4  phasePlot()  =======================================================

phasePlot <- function(dat,
                      ix = 2,
                      iy = 3,
                      colMode = "ratio",
                      col = c(PREYCOL, PREDCOL),
                      main = "",
                      pMax = min(2000, nrow(dat))) {
  #' phasePlot - Phase-space Plot
  #'
  #' @param  dat  num  a numeric matrix that contains population values
  #'                   in columns ix and iy.
  #' @param  ix   int  column index of first population. Default: 2
  #' @param  iy   int  column index of second population. Default: 3
  #' @param  colmode chr  "ratio": color by ratio of species1/species2
  #'                      (default); "length": color by length of the
  #'                      trajectory.
  #' @param  col  chr  color values. Default c(PREYCOL, PREDCOL).
  #' @param  main chr  plot title. Default "".
  #' @param  pMax num  The maximum number of points to plot. Some trajectories
  #'                   can be quite long and the plots can take a long time
  #'                   to plot and update. But we generally do not need to
  #'                   plot at the full resolution of the trajectory. However
  #'                   just plotting every n-th point or so creates plotting
  #'                   artefacts. Instead, we space points by distance along the
  #'                   plotted line. Default is the min() of 2000 points or
  #'                   the number of rows in the dataset.
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

  # subset the plotMatrix
  sel <-selP(dat[ , c(ix, iy)], pMax)
  x <- dat[sel, ix]
  y <- dat[sel, iy]

  # construct a color vector for the trajectory, either by
  # log-ratio of x and y, or by length of the trajectory
  if (colMode == "ratio") {

    # 50 colors from col[1] to col[2]
    makeMyPal <- colorRampPalette(c(col[1], col[2]))
    myPal <- makeMyPal(50)


    # Fetch each row's color based on the log ratio of x and y
    lr <- log(x/y)                             # log ratio
    lr <- (lr - min(lr)) / (max(lr) - min(lr)) # normalize
    lr <- round(lr * 49) + 1                   # scale from 1 to 50
    col <- myPal[lr]                           # assign one of the 50 colours


  } else if (colMode == "length") {

    col <- hcl.colors(length(x), "Spectral")    # compute palette directly

  } else {
    stop(sprintf("Unknown or missing colMode \"%s\".", colMode))
  }


  # empty plot frame
  plot(x, y,
       type = "n",
       xlab = colnames(dat)[ix],
       ylab = colnames(dat)[iy],
       main = main)

  # plot the phase-space trajectory
  plotrix::color.scale.lines(x, y, col = col, lwd = 3)

  return(invisible(NULL))
}


# =    4  Lotka-Volterra model via deSolve::ode()  =============================

# Compute and plot a Lotka Volterra model with an ODE solver

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


# =    5  Lotka-Volterra model via Euler's Method  =============================

# See also Green and Shou (2014).

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
    dPrey = fPrey(nSelf = nPrey, nOther = nPred, par = parPrey) * par$dt
    dPred = fPred(nSelf = nPred, nOther = nPrey, par = parPred) * par$dt

    # compute the new populations
    newPrey <- nPrey + dPrey
    newPred <- nPred + dPred

    # store the results
    outM[i, ] <- c((i-1) * par$dt, newPrey, newPred)

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

# With this "refactoring" of our code, we clearly distinguish those parts of the
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

phasePlot(outE,
          main = "Phase-Space of Lotka-Volterra Model (Euler)",
          colMode = "length")



# =    6  ADDING NOISE  ========================================================
# Now, let's reap the benefits of the increased abstraction of our code. Let's
# add a bit of noise to the predation rate. There are many ways to define
# "noise". Let's use a Gaussian model that increases or decreases the value of
# the predation rate by +- 1/3 of its value.


# ==   6.1  What does added noise look like?  ==================================
myParPrey$predationRate   # 0.2

x <- rnorm(10000,          # What does this distribution look like?
           mean = myParPrey$predationRate,
           sd = myParPrey$predationRate / 3)
head(x)
hist(x, breaks = 50)
abline(v = myParPrey$predationRate, col = "#AA0000", lwd = 3)


# ==   6.2  Modifying our parameters  ==========================================

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

# ... or, coloring the trajectory by time
phasePlot(outS, main = "Stochastic Model II", colMode = "length")

# Let's see what happens if we make this trajectory longer

myPar$Time <- 500  # was 100
outS <- LVEuler(fPrey = fLVPreyStoch, fPred = fLVPredStoch)
phasePlot(outS, main = "Stochastic Model II", colMode = "length")

# The trajectory becomes noisy, but on average, the behaviour is no different from the model without perturbations....

myPar$Time <- 500  # was 100
out <- LVEuler(fPrey = fLVPrey, fPred = fLVPred)
phasePlot(out, main = "Unperturbed Model", colMode = "length")


# =    7  Perturbations via deSolve::ode  ======================================


# If we apply a perturbation over a longer window of time however, we see a
# different behaviour.

LVperturb <- function(t, state, param) {

  with(as.list(c(state, param)), {
#    doPerturbation <- runif(1) < 0.005
    dPrey = (preyGrowthRate * nPrey) - (predationRate * nPrey * nPred)
    dPred = (predatorGrowthRate * nPrey * nPred) - (predatorLossRate * nPred)
    if (t > 29.5 && t < 33) {  # for a given time-window
      dPrey <- dPrey * 2       # multiply prey-chenge by two
      dPred =  dPred / 2       # divide predator-change by 2
    }
    return(list(c(dPrey, dPred)))
  })
}

Par <- c(preyGrowthRate = 1,
         predationRate = 0.2,
         predatorGrowthRate = 0.08,
         predatorLossRate = 0.5)
State <- c(nPrey = 5, nPred = 2)
Time <- seq(0, 100, by = 0.01)

# Run the ODE solver
out <- as.data.frame(deSolve::ode(func  = LVperturb,
                                  y     = State,
                                  parms = Par,
                                  times = Time))

timePlot(out, iy = c(2,3), main = "Lotka-Volterra Model (ODE)")

phasePlot(out, main = "Perturbation trajectory")

phasePlot(out, main = "Perturbation trajectory", colMode = "length")


# The perturbation shifts the system out of its equilibrium, but it quickly
# swings back to a new, stable(!) regime.

# =    8  Going 3D: Chaos  =====================================================

# Let's add another species into our model.

# First we abstract out the basic equation. All our equations have an intrinsic
# term (growth rate or loss rate), and an interaction term (predation rate, and
# predatorGrowth rate). But they are fundamentally the same type of equation:
#
#    dx/dt = a*nX - b*nX*nY.
#
# We had switched the order of the two components of the equation for
# predators - but that is just the same as keeping the order and changing
# the sign on the parameters. We can abstract out the equation and just
# change our parameters accordingly.
#

# == fDE3(): difference equation for a three-species model.
fDE3 <- function(p, n1, n2, n3) {
  dx = (p[1]*n1) - (p[2]*n1*n1) - (p[3]*n1*n2) - (p[4]*n1*n3)
  return(dx)
}

# == LV3D() - a three-species LV model that computes differences per time-step
# with fDE3()
LV3D <- function(t, state, param) {

  with(as.list(c(state, param)), {
    dA = fDE3(param[1, ], nA, nB, nC)
    dB = fDE3(param[2, ], nB, nA, nC)
    dC = fDE3(param[3, ], nC, nA, nB)
    return(list(c(dA, dB, dC)))
  })
}

# The initial parameters for starting conditions: species A and B are prey
# species, species C is a predator.
State <- c(nA = 5, nB = 5, nC = 2)
Time <- seq(0, 100, by = 0.01)
col3 <- c(PREYCOL, PREY2COL, PREDCOL)   # colours

# The individual growth rates, and the pairwise interaction terms are as
# follows. The values were taken from Gilpin (1979).

Par <- matrix(numeric(12), nrow = 3)
Par[1,1] <-  1.0      # r1      Prey 1
Par[1,2] <-  0.001    # nA n1  n1n1 nAnA: a11
Par[1,3] <-  0.001    # nB n2  n1n2 nAnB: a12
Par[1,4] <-  0.01     # nC n3  n1n3 nAnC: a13

Par[2,1] <-  1.0      # r2      Prey 2
Par[2,2] <-  0.001    # nB n1  n1n1 nBnB: a22
Par[2,3] <-  0.0015   # nA n2  n1n2 nBnA: a21
Par[2,4] <-  0.001    # nC n3  n1n3 nBnC: a23

Par[3,1] <- -1.0      # r3      Predator
Par[3,2] <-  0.0      # nC n1  n1n1 nCnC: a33
Par[3,3] <- -0.005    # nA n2  n1n2 nCnA: a31
Par[3,4] <- -0.0005   # nB n3  n1n3 nCnB: a32

# Running the  model works just like for the 2D model ...


out <- as.data.frame(deSolve::ode(func=LV3D, y=State, parms=Par, times=Time))
timePlot(out, iy = c(2,3,4), main = "LV3D", col = col3 )

# Interesting ... Let us increase the simulation time

Time <- seq(0, 500, by = 0.01)
out <- as.data.frame(deSolve::ode(func=LV3D, y=State, parms=Par, times=Time))
timePlot(out, iy = c(2,3,4), main = "LV3D", col = col3 )

# What does the pase-space look like ? There are three 2D trajectories:

# Prey A (column 2) vs. the Predator (column 4)
phasePlot(out, ix=2, iy=4, main = "LV 3D", colMode = "length")

# Prey B (column 3) vs. the Predator (column 4)
phasePlot(out, ix=3, iy=4, main = "LV 3D", colMode = "length")

# Prey A (column 2) vs. Prey B (column 3)
phasePlot(out, ix=2, iy=3, main = "LV 3D", colMode = "length")

# These are 2D projections of a 3D dataset. To get a better sense what the data
# actually looks like, we should plot this dataset in 3D:

# The plotting functions we wrote above, have code that samples the number of
# points along the trajectory. For the lines3D function of the plot3D:: package,
# we should prepare the dataset and subset it. R has no problems to handle
# datasets with millions of rows - but one would not want to plot them directly.

N <- 3000                          # limit to N points along the trajectory

sel <- selP(out[ , c(2, 3, 4)], N) # subset according to the distance between
                                   # points on the trajectory

plot3D::lines3D(out[sel,2], out[sel,3], out[sel,4],
                colkey = FALSE,
                phi=20, theta=100)

# Note that the color gradient that lines3D applies colors along the Z axis
# value, not along the time axis. Too bad. Let's increase the time even more,
# and we can also reduce the resolution somewhat.

# was: Time <- seq(0, 500, by = 0.01)

Time <- seq(0, 5000, by = 0.05)      # 100,000 values

out <- as.data.frame(deSolve::ode(func=LV3D, y=State, parms=Par, times=Time))
timePlot(out, iy = c(2,3,4), main = "LV3D", col = col3, pMax = 4000)

# There is a certain periodicity apparent, but it is not exact.

N <- 6000
sel <- selP(out[ , c(2, 3, 4)], N)

plot3D::lines3D(out[sel,2], out[sel,3], out[sel,4],
                colkey = FALSE,
                phi=20, theta=100)


# We can plot this as a 3D Stereoplot

dev.new()
opar <- par(mfcol = c(1,2),      # two images side-by-side
            mar   = c(1,1,1,1))  # reduce plot-margin widths

plot3D::lines3D(out[sel,2], out[sel,3], out[sel,4], # left-hand plot
                colkey = FALSE,
                phi=-3, theta=98)

plot3D::lines3D(out[sel,2], out[sel,3], out[sel,4], # right-hand plot
                colkey = FALSE,
                phi=-3, theta=102)   # rotate CCW by 5 degrees

par(opar)  #reset the plot parameters

# That's it for now.


# ====  TESTS  =================================================================

if (FALSE) {





}

# [END]
