# LotkaVolterra-ODE.R
#
# Purpose:  Demonstrate the Lotka-Volterra equation as an ODE.
#
# Version:  2.0
# Date:     2022-10-12
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           2.0  Add manual iteration code
#           1.0  First demo in tutorial
#
# ToDo:
#
# Notes:

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


# ====  PACKAGES  ==============================================================
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


# ====  FUNCTIONS  =============================================================


LVmod <- function(t, state, param) {
  # Purpose:
	#     Describe ...
	# Parameters:
	#     a: ...
	#     b: ...
	# Value:
	#     result: ...
	#
	#     preyGrowthRate
	#     predationRate
	#     predatorGrowthRate
	#     predatorLossRate
	#

  with(as.list(c(state, param)), {
    dPrey = (preyGrowthRate * nPrey) - (predationRate * nPrey * nPred)
    dPred = (predatorGrowthRate * nPrey * nPred) - (predatorLossRate * nPred)
    return(list(c(dPrey, dPred)))
  })
}



# ====  PROCESS  ===============================================================

Par <- c(preyGrowthRate = 1,
         predationRate = 0.2,
         predatorGrowthRate = 0.08,
         predatorLossRate = 0.5)

State <- c(nPrey = 5, nPred = 2)

Time <- seq(0, 100, by = 0.01)

out <- as.data.frame(deSolve::ode(func  = LVmod,
                                  y     = State,
                                  parms = Par,
                                  times = Time))
preyCol <- "#008080"
predCol <- "#d91c8a"

plot(out[, "time"], out[, "nPrey"],
     type = "l",
     col = preyCol,
     ylim = c(0, max(out[, "nPrey"] * 1.3)),
     xlab = "time", ylab = "population",
     main = "Lotka-Volterra Model")
lines(out[, "time"], out[, "nPred"],
      type = "l",
      col = predCol)

legend("topright", c("Prey", "Predator"),
       lty = c(1, 1), col = c(preyCol, predCol))

plot(out[ , "nPrey"], out[ , "nPred"],
     type="l",
     main = "Phase-space of Lotka Volterra Model",
     xlab = "Prey population",
     ylab = "Predator population")



# color vector: color by predator prey ratio in 50 steps

# 50 colors from preyCol to predCol
makeMyPal <- colorRampPalette(c(predCol, "#778fa8", preyCol))
myPal <- makeMyPal(50)

# utility function ...
normalize <- function(v) {
  i <- min(v, na.rm = TRUE)
  a <- max(v, na.rm = TRUE)
  v <- (v - i) / (a - i)
  return(v)
}

# fetch each row's color based on the prey / predator ratio
x <- round(normalize(out[, "nPrey"] / out[, "nPred"]) * 49) + 1
ppCol <- myPal[x]


# plot the ratio
plot(out[ , "time"], x,
     type="n",
     main = "Prey / Predator ratio",
     xlab = "t",
     ylab = "nPrey / nPred")
plotrix::color.scale.lines(out[ , "time"],
                           x,
                           col = ppCol,
                           lwd = 4)

# plot the phase-space
plot(out[ , "nPrey"], out[ , "nPred"],
     type="n",
     main = "Phase-space of Lotka Volterra Model",
     xlab = "Prey population",
     ylab = "Predator population")

plotrix::color.scale.lines(out[ , "nPrey"],
                           out[ , "nPred"],
                           col = ppCol,
                           lwd = 4)



# ====  Manual computation  =====================================================
# Cf. Green and Shou (2014)

LVmanual <- function(nPrey = 5,
                     nPred = 2,
                     preyGrowthRate = 1,
                     predationRate = 0.2,
                     predatorGrowthRate = 0.08,
                     predatorLossRate = 0.5,
                     Time = 100,
                     dt = 0.01) {
  N <- Time/dt
  outM <-  matrix(numeric(N * 3), ncol = 3)
  colnames(outM) <- c("time", "nPrey", "nPred")

  outM[1, ] <- c(0.0, nPrey, nPred)

  for (i in 2:N) {
    dPrey = ((preyGrowthRate * nPrey) -
             (predationRate * nPrey * nPred)) * dt
    dPred = ((predatorGrowthRate * nPrey * nPred) -
             (predatorLossRate * nPred)) * dt

    newPrey <- nPrey + dPrey
    newPred <- nPred + dPred

    outM[i, ] <- c((i-1) * dt, newPrey, newPred)

    nPrey <- newPrey
    nPred <- newPred
  }

  return(outM)
}

outM <- LVmanual()

plot(outM[, "time"], outM[, "nPrey"],
     type = "l",
     col = preyCol,
     ylim = c(0, max(outM[, "nPrey"] * 1.3)),
     xlab = "time", ylab = "population",
     main = "Lotka-Volterra Model (manual)")
lines(outM[, "time"], outM[, "nPred"],
      type = "l",
      col = predCol)

legend("topright", c("Prey", "Predator"),
       lty = c(1, 1), col = c(preyCol, predCol))

plot(outM[ , "nPrey"], outM[ , "nPred"],
     type="l",
     main = "Phase-space of Lotka Volterra Model (manual)",
     xlab = "Prey population",
     ylab = "Predator population")


# ====  Stochastic model  ======================================================






# ====  TESTS  =================================================================

if (FALSE) {





}

# [END]
