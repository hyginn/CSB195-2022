# imPlot.R
#
# Plot a matrix of 0 and 1 as a rasterImage()
#
# Boris Steipe (boris.steipe@utoronto.ca)
# 2022-10-23
# ==============================================================================

imPlot <- function(m, colMap, main = "", drawGrid = FALSE) {
  #' @title imPlot()
  #' @description Plot the contents of a matrix, color coded according to
  #' colMap, such that the matrix cell [1.1] is at the top left of the
  #' image and the cells are drawn with an aspect ratio of 1.
  #'
  #' @param  m   The input matrix.
  #' @param  colMap   A named vector of valid colour specifications. If missing,
  #' a default map is produced.
  #' @param main character. A title for the plot.
  #' @param drawGrid  logical. Whether to draw gridlines on the output.
  #' @return      NULL, invisible.
  #' @details Input matrix m is converted to a character matrix and the values
  #' are replaced by the colors found in colMap. No checks are made that all
  #' values actually appear in colmap
  #' @examples
  #'   # pseudo-order in the non-repeating Fibonacci word of length 8*13
  #'   s <- c("0100101001001010010100100101001001010010100100101001",
  #'          "0100100101001001010010100100101001001010010100100101")
  #'   s <- unlist(strsplit(s, ""))
  #'   s <- matrix(s, nrow=8, byrow=TRUE)
  #'   print(s)
  #'   imPlot(s, colMap = c("0"="#88B8CE", "1"="#1978A5"), drawGrid = TRUE)


  nx <- dim(m)[2]
  ny <- dim(m)[1]

  if (missing(colMap)) {
    val <- as.character(sort(unique(as.vector(m))))
    colMap <- colorRampPalette(c(
      "#f8f8f8",
      "#1978a5",
      NULL))(length(val))
    names(colMap) <- val
  }

  im <- matrix(colMap[as.character(m)], ncol = nx)

  opar <- par("mar" = c(0.5, 0.5, 0.5, 0.5))

  plot(c(0, nx),c(0,ny),           # empty plot to define the frame
       type = "n",
       axes = FALSE,
       xaxs = "i", yaxs = "i",
       xlab = "", ylab = "",
       asp = 1)

  rasterImage(im,
              xleft = 0,
              ybottom = 0,
              xright = nx,
              ytop = ny,
              main = main,
              interpolate = FALSE)

  if (drawGrid) {
    segments(0:nx, rep(0, nx+1),               # draw vertical gridlines
             0:nx, rep(ny, nx+1),
             col = "#ddddff")
    segments(rep(0, ny+1), 0:ny,               # draw horizontal gridlines
             rep(nx, ny+1), 0:ny,
             col = "#ddddff")
  }

  par(opar)                                # reset the graphics state

  return(invisible(NULL))                  # return nothing
}


if (FALSE) {
  # pseudo-order in the non-repeating Fibonacci word of length 8*13
  s <- c("0100101001001010010100100101001001010010100100101001",
         "0100100101001001010010100100101001001010010100100101")
  s <- unlist(strsplit(s, ""))
  s <- matrix(s, nrow=8, byrow=TRUE)
  print(s)
  imPlot(s, colMap = c("0"="#88B8CE", "1"="#1978A5"), drawGrid = TRUE)
}

# [END]
