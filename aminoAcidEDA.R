# tocID <- "aminoAcidEDA.R"
#
# Purpose:  Computational Biology Foundations:
#              R code to demonstrate exploratory data analysis (EDA)
#              of amino acid features: PCA and t-SNE
#
# Version:  1.1
#
# Date:     2022-10
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    Prepared for tutorial
#           1.1    Refactored bestCor() - finding good correlations
#
#
#
# TODO:
#    ...
# ==============================================================================
#


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                             Line
#TOC> -----------------------------------------------------------------
#TOC>   1        Dimensionality reduction                            45
#TOC>   2        PCA introduction                                    56
#TOC>   2.1        minimal PCA example                               87
#TOC>   2.2        Calculate a PCA of x1 and y2                     125
#TOC>   2.3        prcomp() and princomp()                          156
#TOC>   2.4        Scaling                                          173
#TOC>   3        EDA with PCA                                       190
#TOC>   3.1        Load the aaindex data                            192
#TOC>   3.2        Run the analysis                                 216
#TOC>   3.3        Interpreting the principal components            246
#TOC>   4        t-Stochastic neighbour embedding  (tsne)           314
#TOC>   4.1        tsne of the aaindex                              336
#TOC>   4.2        Final thoughts                                   379
#TOC> 
#TOC> ==========================================================================


# =    1  Dimensionality reduction  ============================================
#

# The seqinr dataset aaindex is a list that contains 544 vectors of amino acid
# properties. But many of them are redundant and highly correlated with each
# other. How many dimensions does a non-redundant dataset of amino acid
# properties have, and what are they?

# A common approach to such dimensionality reduction is "principal component
# analysis", or PCA.

# =    2  PCA introduction  ====================================================


# The goal of Principal Component Analysis (PCA) is to transform a number of
# possibly correlated variables into a smaller number of uncorrelated variables
# called principal components.

# The (possibly) smaller number of variables can be used for data reduction and
# visualization.

# Principal component analysis (PCA) converts a set of observations of possibly
# correlated variables into a set of values of uncorrelated variables called
# principal components. The first principal component is the projection of the
# data into a single dimension that has as high a variance as possible (that is,
# accounts for as much of the variability in the data as possible); each
# succeeding component in turn has the highest variance possible under the
# constraint that it be orthogonal to (uncorrelated with) the preceding
# components. Therefore the PCs provide a view on the structure of the data that
# best explains its variance. This is especially useful for EDA of
# high-dimensional data that can't be intuitively visualized. Given a set of
# points in Euclidean space, the first principal component corresponds to a line
# that passes through the multidimensional mean and minimizes the sum of squares
# of the distances of the points from the line.

# The second principal component is calculated in the same way, after all
# correlation with the first principal component has been subtracted out from
# the points.

# Let's illustrate this with a simple 2D example:


# ==   2.1  minimal PCA example  ===============================================

# 500 normally distributed samples each: uncorrelated data
set.seed(112358)
x1 <- rnorm(500,0,1)
y1 <- rnorm(500,0,1)
plot(x1, y1)
cor(x1, y1)


# generate y2 correlated with (i.e. "dependent on") x1
y2 <- 2*x1 + y1
mean(y2)
y2 <- y2-mean(y2)
mean(y2)
sd(y2)
y2 <- y2 / sd(y2)
sd(y2)

# Now y2 _looks_ like a random distribution ...
hist(y2, breaks = 10)

# .. but its values are actually highly correlated with x1:
plot(x1, y2)
cor(x1, y2)

# Compare:
oPar <- par(mfrow = c(2,2)) # set new and save old graphics state

# four plots ...
hist(x1)
hist(y2)
plot(x1, y1)
plot(x1, y2)

par(oPar) # restore graphics state parameters


# ==   2.2  Calculate a PCA of x1 and y2  ======================================

?prcomp

pcaSample <- prcomp(cbind(x1,y2))

# here are the information items from the returned list of results
pcaSample
str(pcaSample)
pcaSample$sdev
pcaSample$rotation
summary(pcaSample)
head(pcaSample$x)
plot(pcaSample$x, xlim=c(-5,5), ylim=c(-5,5))

# After PCA, most of thew variance is now contained in PC1 - and we could decide
# that PC1 is sufficient to describe the relationships in our data. Thus, we
# have reduced a two-dimensional dataset to a one-dimensional dataset.

# Compare the histograms before and after the rotation:
oPar <- par(mfrow = c(2,2))
hist(x1, xlim=c(-4,4), ylim=c(0,150), main="", breaks=12)
hist(y2, xlim=c(-4,4), ylim=c(0,150), main="", breaks=12)
hist(pcaSample$x[,1], xlim=c(-4,4), ylim=c(0,150),
     main="", col=rgb(0.86,0,0,0.5), , breaks=12)
hist(pcaSample$x[,2], xlim=c(-4,4), ylim=c(0,150),
     main="", col=rgb(0.31, 0.5, 0.74, 0.5), , breaks=12)
par(oPar) # restore graphics state parameters



# ==   2.3  prcomp() and princomp()  ===========================================

# R has two different functions for PCA: prcomp() and princomp(). They use
# different mathematical approaches but the results are virtually identical.
# prcomp() is numerically more stable. However, they also use different names
# for the elements of their result lists.

#  prcomp()  princomp()
#
#  center    center	  The vector that was subtracted to center the data
#  sdev      sdev     Standard deviations for each dimension of the rotated data
#  rotation  loadings The actual principal components
#  x         scores   The rotated data, i.e. after projection along each PC

# e.g. use data$x for the rotated results of a prcomp() call, but use
# data$scores if the result came from princomp()

# ==   2.4  Scaling  ===========================================================
#
# PCA is sensitive to the scaling of the variables.

# If we have just two variables and they have the same sample variance and are
# positively correlated, then the PCA will entail a rotation by 45° and the
# "loadings" for the two variables with respect to the principal component will
# be equal. But if we multiply all values of the first variable by 100, then the
# principal component will be almost the same as that variable, with a small
# contribution from the other variable, whereas the second component will be
# almost aligned with the second original variable. This means that whenever the
# different variables have different units (like temperature and mass), PCA is a
# somewhat arbitrary method of analysis. (Different results would be obtained if
# one used Fahrenheit rather than Celsius for example.) One way to address this
# is to scale variables to have unit variance.


# =    3  EDA with PCA  ========================================================

# ==   3.1  Load the aaindex data  =============================================

library(seqinr)
data(aaindex)

aaFeatures <- data.frame(aaindex[[1]]$I)

for (i in 2:length(aaindex)) {
  aaFeatures <- cbind(aaFeatures, aaindex[[i]]$I)
}

colnames(aaFeatures) <- 1:ncol(aaFeatures)
aaFeatures <- scale(aaFeatures)
colMeans(aaFeatures)  # Note: there are NA values! PCA will not work
                      # if NA values are present.


aaindex[[524]]$D
aaFeatures[ , 524]

sel <- colMeans(aaFeatures)
sel <- ifelse(is.na(sel), FALSE, TRUE)


# ==   3.2  Run the analysis  ==================================================

aaDat  <- t(aaFeatures[, sel])   # Note: prcomp() expects categories in columns
                                 # these are our "things", and values in rows,
                                 # these are our things' propertiers of many
                                 # dimensions

pcaAA <- prcomp(aaDat, scale. = TRUE)

plot(pcaAA)
summary(pcaAA)
str(pcaAA)

# Plot projections along the components into a scatterplot.
# Axes for points are scaled as values, for vectors as variance
# Default for biplot() is the first and second component.

biplot(pcaAA)


# Examine the actual principal components in a parallel-coordinates
# plot.
N <- 4
matplot(pcaAA$rotation[, 1:N],
        type="b", lwd=1,
        xlab = "amino acid", ylab="PCs")
# We could (on a rainy day) print the amino acid codes onto the x-axis.



# ==   3.3  Interpreting the principal components  =============================

# ... what actual tables do they correlate with?
myCorPC <- matrix(numeric(ncol(aaFeatures) * 4), ncol = 4)
PC1 <- pcaAA$rotation[ , 1]
PC2 <- pcaAA$rotation[ , 2]
PC3 <- pcaAA$rotation[ , 3]
PC4 <- pcaAA$rotation[ , 4]

for (i in 1:ncol(aaFeatures)) {
  myCorPC[i, 1] <- cor(aaFeatures[ , i], PC1)
  myCorPC[i, 2] <- cor(aaFeatures[ , i], PC2)
  myCorPC[i, 3] <- cor(aaFeatures[ , i], PC3)
  myCorPC[i, 4] <- cor(aaFeatures[ , i], PC4)
}
summary(myCorPC)

hist(myCorPC[ , 1])

# A function to find the best correlation:

bestCor <- function(V, dat = aaFeatures) {
  # Identify that column in dataframe dat that has the highest correlation
  # or anticorrelation with vector V. Print some information about it.
  # Value: The index of the column (invisibly)
  #
  # ToDo: identify and return n-best correlations
  #
  cors <- numeric(ncol(dat))
  for (i in 1:ncol(dat)) {
    cors[i] <- cor(V, dat[ , i]) # correlation between input vector and each
  }                              # column of the feature set in turn

  myBest <- max(abs(cors), na.rm = TRUE)
  sel <- abs(cors) == myBest   # Is TRUE for the highest absolute value in the
                               # vector of correlation coefficients.
  idx <- which(sel)[1]         # which() finds the index of the TRUE value(s).
                               # Pick only the first in case there are ties.
  cat(sprintf("\n Highest correlation for index %d (%f): %s ",
              idx,
              cors[idx],
              aaindex[[idx]]$D))
  return(invisible(idx))      # return the index (invisibly).
}

bestCor(pcaAA$rotation[ , 1])  # This is PC1 ...
bestCor(pcaAA$rotation[ , 2])  # PC2
bestCor(pcaAA$rotation[ , 3])  # etc.
bestCor(pcaAA$rotation[ , 4])


plot(PC1, PC2)  # uninformative

# ... we need a better way to plot this ... we need to identify the amino acids
plot(PC1, PC2, type = n)
text(PC1, PC2, labels = names(PC1))

# ... and if the one letter codes do not tell us a lot, we need to use
# other features so the plot can tell a story. Study the function code
# of plotAA() in the .util.R script.
#
plotAA(PC1, PC2)
plotAA(PC2, PC3)
plotAA(PC3, PC4)

# What have we learned?


# =    4  t-Stochastic neighbour embedding  (tsne)  ============================

# t-Stochastic Neighbour Embedding is a powerful dimension re-
# duction algorithm developed in the lab of Geoff Hinton at UofT.
#
# see: https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding
#
# Let's try the t-SNE algorithm to explore amino acid similarity

if (! requireNamespace("tsne", quietly=TRUE)) {
  install.packages("tsne")
}
library(tsne)

# Package information:
#  library(help =   tsne)     # basic information
#  browseVignettes("tsne")    # available vignettes
#  data(package =  "tsne")    # available datasets

?tsne


# ==   4.1  tsne of the aaindex  ===============================================


# get a vector to mask all columns in aaFeatures with NA values
sel <- ifelse(is.na(colMeans(aaFeatures)), FALSE, TRUE)

tsnePlot <- function(x) {
  # tsne() plots a 2D matrix of numbers - x, but does not remember rownames
  rownames(x) <- rownames(aaFeatures)
  plotAA(x[ , 1],x[ , 2])
}

# Make the run reproducible
set.seed(4184542)   # Iteration #700 error is: 0.409304304796215

# run tsne
tsneAA <- tsne(aaFeatures[ , sel],
               epoch_callback = tsnePlot,
               epoch = 50,
               perplexity = 160,
               max_iter = 1000)

# Here too we might ask: Can these reduced dimensions be interpreted? tsne()
# returns the "collapsed" dimensions and we have assigned them the tsneAA,
# which is just a matrix with 20 rows and 2 columns. So we can ask again: what
# was it that tsne() found?


bestCor(tsneAA[ , 1])
bestCor(tsneAA[ , 2])

# Are these interpretable?
# Or are the poor correlations actually spurious? E.g.

x1 <- tsneAA[ , 1]
names(x1) <- names(aaindex[[326]]$I)
plotAA(x1, aaindex[[326]]$I)

x2 <- tsneAA[ , 2]
names(x2) <- names(aaindex[[75]]$I)
plotAA(x2, aaindex[[75]]$I)


# ==   4.2  Final thoughts  ====================================================

# Are we actually improving our quantitative measure of similarity when we add
# lots and lots of data, without regard for what the data means? Can the
# algorithm figure things out on its own?
#
# How would we assess that?
#
# The two distributions appear quite random, there are no
# clusters of amino acids that are mutually similar. Everything seems to be
# pretty evenly spaced. Does this mean that biology has found amino acids that
# are quite non-redundant in their properties?



# [END]
