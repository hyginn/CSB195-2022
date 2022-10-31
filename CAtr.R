# tr cellular automata


applyRRule <- function(R, sNow, sPast) {
  #' @title applyRule()
  #' @description Time reversible CA rule computed from Wolfram code
  #'
  #' @param  R   int. The Wolfram code for a 1D cellular automaton
  #'                      rule set; an integer in [0, 255].
  #' @param  sNow   int The current state of the cell and its immediate
  #' neighbourhood. A binary vector of length three.
  #' @param  sNow   int The past state of the cell. Either 0 or 1.
  #' @return      0 or 1
  #' @examples
  #'   applyRRule(110,  c(1, 0, 1), 1)
  #'   applyRule( 30, c(1, 0, 0), 0)

  v <- as.integer(intToBits(R)[8:1] == "01")           # output integers of R
  idx <-  8 - (sNow[1]*4) - (sNow[2]*2) - (sNow[3]*1)  # index in v from sNow
  S <- v[idx]
  S <- xor(S, sPast)

  return(S)
}




CAtr <- function(R, nx = 100, ny = 100,
               m0, m1,
               wrap = TRUE, N = NULL) {
  #' @title CAtr()
  #' @description Evolve a cellular automaton with a time-reversible Wolfram
  #' code R in a nx * ny matrix.
  #'
  #' @param  R   int. The Wolfram code for a 1D cellular automaton rule set; an
  #'   integer in [0, 255].
  #' @param  nx  the width of the space in which the automaton evolves
  #' @param  ny  the number of steps to compute
  #' @param  vInit   the initialization vector. See details.
  #' @param  seed    Initializaton of the RNG. Default NULL.
  #' @param  wrap    if TRUE, wrap the vertical edges of the matrix (periodic
  #'   boundary).
  #' @param  N    If this is not NULL, the function plots every updated row
  #'   until the matrix is filled, then scrolls the matrix for a total of N
  #'   steps (or until the process is aborted.).
  #' @return      the evolved matrix as 0s and 1s.
  #' @details If vInit is missing or NULL, a single 1 is placed at the centre of
  #'   the first row. If vInit is in [0, 1[ the vector is filled with random 1s
  #'   of p == vInit. If vInit is >= 1, round(vInit) 1s are equally spaced into
  #'   vInit. If vInit is the keyword "fib", the first row is initialized with a
  #'   Fibonacci word. If length(vInit) > 1 the vector is recycled into the
  #'   length of vInit. In particular, this allows to pass in a custom
  #'   initialization vector of length nx - this vector must be a numeric vector
  #'   of 0s and 1s. If vInit is a character-string of 0s and 1s, it is
  #'   converted into such an initialization vector.
  #' @examples
  #'



  if (is.null(N)) {  # No scrolling
    N <- 0
  } else {
    Nmax <- N
    cat("\n")
  }


  m <- matrix(integer(nx * ny), ncol = nx)
  m[1,] <- m0
  m[2,] <- m1

  if (wrap) {
    m <- cbind(m[ ,nx], m, m[ ,1]) # wrap first and last column
    nx <- nx + 2                  # increase nx accordingly
  }


  for (i in 2:(ny-1)) {
    for (j in 2:(nx-1)) {
      m[(i+1), j] <- applyRRule(R, m[i, (j-1):(j+1)], m[i-1, j])
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
      m[(i+1), j] <- applyRRule(R, m[i, (j-1):(j+1)], m[i-1, j])
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

  return(invisible(m))
}

if (FALSE) {

  # == Samples
  dev.new()

  X <- 250
  Y <- round(X * 1.618)

  m0 <- rep(0, X)
  m1 <- m0; m1[sample(1:X, 7)] <- 1

  m <- CAtr(1, nx=X, ny=Y, m0=m0, m1=m1); imPlot(m)
  m <- CAtr(37, nx=X, ny=Y, m0=m0, m1=m1); imPlot(m)
  m <- CAtr(73, nx=X, ny=Y, m0=m0, m1=m1); imPlot(m)

  m <- CAtr((i<- sample(1:255, 1)), nx=X, ny=Y, m0=m0, m1=m1)
  imPlot(m, main = sprintf("Rule%dR", i))


  # ===Evolving 122R ==============================
  #
  # A simple pattern evolving through Rule 122R

  dev.new()
  m0 <- c(rep(0,31),        rep(1,18),        rep(0,31))
  m1 <- c(rep(0,30), 0, 0,  rep(1,16), 0, 0,  rep(0,30))
  m <- CAtr(122, nx=80, ny=500, m0 = m0, m1 = m1, N=1200)

  # Add a perturbation ... things get chaotic
  #                  +------  HERE -------+
  #                  |                    |
  #                  v                    v
  m1 <- c(rep(0,30), 1, 0,  rep(1,16), 0, 1,  rep(0,30))
  m <- CAtr(122, nx=80, ny=500, m0 = m0, m1 = m1, N=1200)


  # reverse
  m <- CAtr(122, nx=ncol(m), ny=nrow(m),
            m0=m[nrow(m),  ],
            m1=m[nrow(m)-1,],
            N=1200)



}

