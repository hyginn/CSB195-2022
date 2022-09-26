# tocID <- "aminoAcidProperties.R"
#
# Purpose:  Computational Biology Foundations:
#              R code demonstrating amino acid properties with ChimeraX
#              visualization
#
# Version:  1.0
#
# Date:     2022-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#
#           1.0    First version
#
#
# PREREQUISITES:
#   Installation of a recent version of ChimeraX
#   Some familiarity with ChimeraX command syntax
#   Learning Units: BIN-SX-Concepts
#                   BIN-SX-Chimera
#
# TODO:
#    ...
#
#
# == DO NOT SIMPLY  source()  THIS FILE! =======================================
#
# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
# going on. That's not how it works ...
#
# ==============================================================================
#
#
#
# =    1  ChimeraX REMOTE SCRIPTING  ===========================================


# One of the cool features of ChimeraX is that it can be driven by Python code,
# both within a running session and through Python scripts. What I find even
# cooler though is that ChimeraX can be driven from any programming language via
# its remote control function that can listen to commands sent from any other
# application. The interface that is used here is the standard REST (method) -
# the GET and POST verbs that ubiquitously underly the communication of clients
# and servers on the Web.

# In order to establish the communication between this script and ChimeraX, all
# we need to do is:
#  - open ChimeraX;
#  - tell it to listen on a specific "port";
#  - send commands to that port via httr::


# ==   1.1  Defining a Port  ===================================================

# The httr:: package needs to be available

if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
# Package information:
#  library(help = httr)       # basic information
#  browseVignettes("httr")    # available vignettes
#  data(package = "httr")     # available datasets

# We need to think of a port. Any available port number between 49152-65535 is
# fine. We'll choose 61803.

CXPORT <- 61803

# Check that our current version of R supports sockets (default since V 3.3)
capabilities("sockets")   # MUST be TRUE. If not, don't continue.


# ==   1.2  Open ChimeraX  =====================================================

#  - Open a fresh, new session of recently updated version of ChimeraX
#  - type:
#
#       remotecontrol rest start port 61803
#
#    ... or whatever the value of CXPORT is.

# Now watch what happens in ChimeraX when you execute the following line:
( x <- httr::GET("http://127.0.0.1:61803/run?command=open+7ACN") )

# (7ACN is a structure of Aconitase - an enzyme of the citric acid cycle. Cf.
#   - https://pdb101.rcsb.org/motm/154
#   - )


# To make this convenient, we need a function -
#   - that accepts ChimeraX commands as a parameter;
#   - that sends the command to Chimera, correctly formatted
#   - that captures output which we may receive.
#
# Here is such a custom function:

CX <- function(cmd, port = CXPORT, quietly = FALSE) {
  # send a command to ChimeraX listening on port CXPORT via its REST
  # interface.
  # Parameters:
  #   cmd      char     a ChimeraX commandline command
  #   port     int      the portnumber on which ChimeraX is listening
  #   quietly  logical  if FALSE, cat() the contents of the response
  #
  # Value:  the reply by ChimeraX, invisibly.

  CXREST <- sprintf("http://127.0.0.1:%s/run?", CXPORT)

  cmd <- gsub("(^\\s+)|(\\s+$)", "", cmd)  # trim whitespace
  # percent encode reserved characters
  cmd <- gsub("%",   "%25", cmd)          #   %
  cmd <- gsub("#",   "%23", cmd)          #   #
  cmd <- gsub("/",   "%2F", cmd)          #   /
  cmd <- gsub(":",   "%3A", cmd)          #   :
  cmd <- gsub("@",   "%40", cmd)          #   @
  cmd <- gsub(",",   "%2C", cmd)          #   ,
  cmd <- gsub("\\*", "%2A", cmd)          #   *
  cmd <- gsub("\\?", "%3F", cmd)          #   ?
  cmd <- gsub("!",   "%21", cmd)          #   !
  cmd <- gsub("=",   "%3D", cmd)          #   =
  cmd <- gsub("\\(", "%28", cmd)          #   (
  cmd <- gsub("\\)", "%29", cmd)          #   )
  cmd <- gsub("\\[", "%5B", cmd)          #   [
  cmd <- gsub("\\]", "%5D", cmd)          #   ]
  cmd <- gsub("&",   "%26", cmd)          #   &
  cmd <- gsub("\\+", "%2B", cmd)          #   +

  cmd <- gsub("\\s+", "+", cmd)            # whitespace to "+"
  cmd <- URLencode(cmd)                    # encode special characters
  cmd <- paste0(CXREST, "command=", cmd, collapse = "")

  r <- httr::GET(cmd)

  if (! r$status_code == 200) {
    stop("ChimeraX returned status code %d", r$status_code)
  }

  if (length(r$content) == 0) {
    reply <- ""
  } else {
    reply <- rawToChar(r$content)
  }

  if (quietly == FALSE) {
    cat(reply)
  }

  return(invisible(reply))

}

# Use the function to set up a visualization

# CX("camera sbs")
CX("camera mono")
CX("lighting soft")
CX("lighting shadows true intensity 0.8")
# define a color gradient to be able to follow the fold of the protein
# cf. https://www.cgl.ucsf.edu/chimerax/docs/user/commands/colornames.html
CX("color sequential #1 & protein target abc palette #B0C4DE:rgb(44,50,56):orchid:maroon")
CX("hide cartoons")
CX("show atoms")
CX("style sphere")
CX("hide ::name='HOH'")

CX("color #1 & protein slate gray")

# Select sidechains of hydrophobic amino acids ...
CX("select :SER, ASN, ASP, GLN, GLU, ARG, LYS, HIS & sidechain")
CX("color sel royal blue")

# ... and hydrophilic amino acids ...
CX("select :VAL, LEU, ILE, MET, PHE, TRP & sidechain")
CX("color sel gold")
CX("select clear")

# set clipping planes to show a slice of the protein
CX("clip near -2.0 far 2.0 position #1")

# move forward and back
CX("move z 0.5")
CX("move z 0.5")
CX("move z 0.5")
CX("move z -0.5")
CX("move z -0.5")
CX("move z -0.5")

CX("clip near -28.0 far -24.0 position #1")
delta <- 0.025
N     <- 2000
for (i in c(1:N, -(1:N))) {
#  print(sprintf("move z %f", sign(i) * delta))
  CX(sprintf("move z %f", sign(i) * delta))
  CX("wait 1")
}

CX("clip near 0.0 far 20.0 position #1")
CX("roll y 0.1")
CX("stop")

# Reset our view
CX("clip off")
CX("hide atoms")
CX("show cartoons")
CX("color sequential #1 & protein target abc palette #B0C4DE:rgb(44,50,56):orchid:maroon")
CX("hide cartoons")

# create a second copy to work with a single residue only
CX("open 7acn")
CX("sequence chain #2/A")

# ... select a residue
CX("select ~sel")
CX("delete sel")
CX("show atoms")
CX("surface")
CX("style stick")

=
CX("remotecontrol rest stop")  # release the socket
# Done.



# [END]
