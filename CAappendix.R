# tocID <- "CAappendix.R"
#
#
#
# Demo code
#
# 2022 -10  -  2022 - 11
#
# Boris  Steipe (boris.steipe@utoronto.ca)
#
#
# Version:  1.0
#
# Versions:
#           1.0 First version
#
# ==============================================================================
#
#

graphics.off()    # close all open graphics devices
dev.new()         # re-open the standard device
dev.new()         # open a second, floating graphics window


# ==  Recap: Rule 18 in "Standard mode

# Standard mode, Rule 18, all defaults
m <- CA(18); imPlot(m)

# == Larger ...

X <- 250
Y <- 404
m <- CA(18, nx = X, ny = Y); imPlot(m)


# == Two versions of a probabilistic CA

# pI - "Input" ... misreading an input state ...
m <- CA(18, nx=X, ny=Y,   ruleMode="pI",   pR = 0.0002);  imPlot(m)

# pO - "Output" ... miswriting an output state ...
m <- CA(18, nx=X, ny=Y,   ruleMode="pO",   pR = 0.0002);  imPlot(m)




# == A time-reversible CA!

# Rule 37R
mI <- matrix(integer(2*X), nrow = 2)    # two empty rows
mI[2, sample(1:X, 17)] <- 1             # seventeen random 1s
m <- CA(37, nx=X, ny=Y, ruleMode="tRev", vInit=mI);   imPlot(m)

# Rule 22R (Try this a few times...)
mI <- matrix(integer(2*X), nrow = 2)   # two empty rows
mI[2, sample(1:X, 21)] <- 1             # twenty one random 1s
m <- CA(22, nx=X, ny=Y, ruleMode="tRev", vInit=mI);   imPlot(m)

# Evolutions from the same, sparse staring state, for comparison
mI <- matrix(integer(2*X), nrow = 2)   # two empty rows
mI[2, sample(1:X, 11)] <- 1             # eleven random 1s

# Try a few dozen ... select, run, run again ...
m <- CA((i <- sample(1:255, 1)), nx=X, ny=Y, ruleMode="tRev", vInit=mI)
imPlot(m, main = sprintf("Rule%dR", i))



# === More time-reversibility: Evolving 122R ==============================
#

# A simple pattern evolving through Rule 122R

X <- 80
Y <- 500
mI <- matrix(integer(2*X), nrow = 2)   # two empty rows

mI[1,] <- c(rep(0,31),        rep(1,18),        rep(0,31))
mI[2,] <- c(rep(0,30), 0, 0,  rep(1,16), 0, 0,  rep(0,30))

m <- CA(122, nx=X, ny=Y, ruleMode="tRev", vInit=mI, N = 500)
# Not much going on ... you can stop this. Now we will add just
# TWO additional pixels to the second row... things get chaotic

#                      +------  HERE -------+
#                      |                    |
#                      v                    v
mI[2,] <- c(rep(0,30), 1, 0,  rep(1,16), 0, 1,  rep(0,30))

m <- CA(122, nx=X, ny=Y, ruleMode="tRev", vInit=mI, N = 1200)
# Let this run to the end! Otherwise the matrix we need below won't exist.
# Beautiful, disordered pattern.

# Now we swap the last two rows around, and use them as the starting rows.
# Again: 1,200 steps ... where does this take us?
mI[1, ] <- m[nrow(m),  ]
mI[2, ] <- m[nrow(m)-1,]
m <- CA(122,
        nx=ncol(m), ny=nrow(m),
        ruleMode="tRev",
        vInit=mI,
        N=1200)



# [END]
