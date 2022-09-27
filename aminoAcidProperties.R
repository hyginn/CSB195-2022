# tocID <- "aminoAcidProperties.R"
#
# Purpose:  Computational Biology Foundations:
#              R code demonstrating amino acid properties with ChimeraX
#              visualization
#
# Version:  1.1
#
# Date:     2022-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.1    Clean up post lecture.
#           1.0    First version for lecture.
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


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                               Line
#TOC> -------------------------------------------------------------------
#TOC>   1        ChimeraX REMOTE SCRIPTING                             55
#TOC>   1.1        Defining a Port                                     73
#TOC>   1.2        Prepare ChimeraX                                    95
#TOC>   1.3        A function to communicate with ChimeraX            115
#TOC>   2        Visualize the 7ACN fold.                             193
#TOC>   2.1        Composition of the protein's structure             210
#TOC>   2.1.1          Visualize a "slice" through the protein        241
#TOC>   2.1.2          Continuous translation  via a for-loop ...     255
#TOC>   2.2        Study a single amino acid in context:              289
#TOC> 
#TOC> ==========================================================================


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
#  - send commands to that port via httr:: package functions.


# ==   1.1  Defining a Port  ===================================================

# The httr:: package needs to be installed ...

if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
# Package information:
#  library(help = httr)       # basic information
#  browseVignettes("httr")    # available vignettes
#  data(package = "httr")     # available datasets

# We need to define a port - that's something like a mailbox number which
# processes agree on to share information. Any available port number
# between 49152-65535 is fine. We'll choose 61803.

CXPORT <- 61803

# Check that our current version of R supports sockets (default since V 3.3):
capabilities("sockets")   # MUST be TRUE. If not, don't continue.


# ==   1.2  Prepare ChimeraX  ==================================================

#  - Open a fresh, new session of a recently updated version of ChimeraX
#  - type:
#       remotecontrol rest start port 61803
#
#    (... or whatever value you have given CXPORT.)

# Now watch what happens in ChimeraX when you execute the following line:
( x <- httr::GET("http://127.0.0.1:61803/run?command=open+7ACN") )

# (127.0.0.1 is the "localhost" IP address. Only this one is supported for
#  security reasons, but ChimeraX can be opened: cf.
#  https://www.rbvi.ucsf.edu/pipermail/chimerax-users/2017-July/000113.html .)

# (7ACN is a structure of Aconitase - an enzyme of the citric acid cycle. Cf.
#   - https://pdb101.rcsb.org/motm/154
#   - https://pdb101.rcsb.org/motm/89
#   - https://en.wikipedia.org/wiki/Citric_acid_cycle )

# ==   1.3  A function to communicate with ChimeraX  ===========================
#
# To make it convenient to send commands to Chimera, we need a function -
#   - that accepts ChimeraX commands as a parameter;
#   - that sends the command to Chimera, correctly formatted;
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

  # (A) construct the base address, port, and command:
  CXREST <- sprintf("http://127.0.0.1:%s/run?", CXPORT)

  # (B) sanitize the user-enterd variable cmd to be properly encoded in
  # a http message. (No blanks, some specially handled characters, ...)
  # Here we use gsub(), which seeks for a pattern defined in its first
  # argument, substitutes the contents of the second argument, in the
  # string that is identified with the third argument.
  #
  # Patterns are defined as "regular expressions", we'll encounter those later.
  #
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
  cmd <- URLencode(cmd)                    # encode other special characters

  # Combine the base-address and the current contents of cmd ...
  cmd <- paste0(CXREST, "command=", cmd, collapse = "")

  # send the command to ChimeraX, and capture the response ...
  r <- httr::GET(cmd)
  # ... the response is a list-object and we can analyze it.

  if (! r$status_code == 200) {  # If we did NOT receive a "success" status code
    stop("ChimeraX returned status code %d", r$status_code)
  }

  if (length(r$content) == 0) {
    reply <- ""
  } else {
    reply <- rawToChar(r$content)  # reformat
  }

  if (quietly == FALSE) {
    cat(reply)                     # print the response
  }

  return(invisible(reply))         # return the reply  but do not also
                                   # print it.

}

# =    2  Visualize the 7ACN fold.  ============================================
#
# Use the CX() function to set up a visualization of the 7ACN fold.

CX("camera sbs")   # practice your stereo vision
# CX("camera mono")
CX("lighting soft")
CX("lighting shadows true intensity 0.8")

# Define a color gradient to be able to follow the fold of the protein
# cf. https://www.cgl.ucsf.edu/chimerax/docs/user/commands/colornames.html
CX("color sequential #1 & protein target abc palette #B0C4DE:rgb(44,50,56):orchid:maroon")

# Study the fold. Note where it begins, how different elements of "secondary
# structure" (helices, strands) form compact domains, and how the domains
# assemble to a structural whole.

# ==   2.1  Composition of the protein's structure  ============================

# Now change the view to show the actual atoms (as little spheres.)
#
CX("hide cartoons")
CX("show atoms")
CX("style sphere")

# Note the little red dots. These are water molecules (H2O); what you see
# is the position of each oxygen atom. We hide them so we don't clutter
# the screen too much.
CX("hide ::name='HOH'")

# First we color all atoms gray ...
CX("color #1 & protein slate gray")

# Then we select two categories of amino acids.
# Select sidechains of hydrophilic amino acids ...
# (these have atoms in their side-chains that make favourable
# interactions with water)
CX("select :SER, ASN, ASP, GLN, GLU, ARG, LYS, HIS & sidechain")
CX("color sel royal blue")

# ... and hydrophobic amino acids ...
# (these have mostly only carbon atoms in their side chains, which gives
#  them an "oily" nature)
CX("select :VAL, LEU, ILE, MET, PHE, TRP & sidechain")
CX("color sel gold")
CX("select clear")  # clear the selection so it doe not clutter our view


# ===   2.1.1  Visualize a "slice" through the protein   
#
# set clipping planes to show a slice of the protein
CX("clip near -2.0 far 2.0 position #1")

# move the structure away from us and towards us...
CX("move z 1.0")
CX("move z 1.0")
CX("move z 1.0")
CX("move z -1.0")
CX("move z -1.0")
CX("move z -1.0")


# ===   2.1.2  Continuous translation  via a for-loop ...
# ... instead of manually moving along a a view step by step:

CX("clip near -32.0 far -28.0 position #1")

delta <- 0.06     # step-size
N     <- 1000     # number of steps
for (i in c(1:N, -(1:N))) { # Do this N times forward and N times back
#  print(sprintf("move z %f", sign(i) * delta))
  CX(sprintf("move z %f", sign(i) * delta))
  CX("wait 1")
}

# Watch the protein slowly move through the slice-zone. Pay particular
# attention to where the blue and yellow sidechains are located. Inside or
# on the surface of the protein? Clustered, or distributed? Random, or non-
# random?

# (If this takes too long and you want to abort it, just click on the red
# "Stopsign" in the menu bar of the Console pane.)

# As an alternative to a sliding motion, we can rotate the protein slowly:
CX("clip near 0.0 far 20.0 position #1")
CX("roll y 0.1")   # (Documentation for all commands is available via the
                   #  ChimeraX help functions. Easiest access is
                   #  by clicking the hyperlinked command in the log window
                   #  of Chimerax. It is either open, or can be opened via
                   #  the mnu: Tools > Log.
                   #
CX("stop") # stop the rotation

# Task: what have you observed? Why is this so? What does this imply?


# ==   2.2  Study a single amino acid in context:  =============================

# Reset our view ...
CX("clip off")
CX("hide atoms")
CX("show cartoons")
CX("color sequential #1 & protein target abc palette #B0C4DE:rgb(44,50,56):orchid:maroon")

# Let's focus on a single amino acid and examine its context.

# Create a second copy of the protein and delete everything but this single
# amino acid, to work with a single residue only. (Isn't there a way to
# duplicate just one residue or other molecule? IDK.)
CX("open 7acn")

# Open the Sequence viewer window:
CX("sequence chain #2/A")

# To choose an amino acid of a certain type hover over it.
# You will get information like "7acn #1 \A ARG 114 CG" - which is:
#
#   7acn: the PDB ID
#   #1:   model number 1
#   \A:   chain A
#   ARG:  amino acid type (ARG -> arginine)
#   114:  position in sequence  ("sequence number")
#   CG:   atom topology, "Cgamma", third atom of the side-chain

# Let's display the residue with a solvent accessible surface.
# If, say, your amino acid is HIS 298, define this as a named selection:
CX("name myAA #2/A:298")              # define
CX("select #2/A")                     # select the entire chain
CX("select subtract myAA")            # subtract the named selection
CX("delete sel")                      # delete what is selected, myAA remains
CX("view myAA")                       # center the view
CX("show myAA atoms")                 # show atoms, and ...
CX("style myAA ball")                 # style as ball-and-stick
CX("color myAA byelement target ab")  # color by element
CX("surface myAA")                    # compute surface
CX("color myAA #AAAADD target s")     # color wispy pale blue
CX("transparency myAA 80 target s")   # set high transparency

# Now lets calculate a surface of the protein _without_ myAA, to see how the
# amino acid is packed against the rest of the structure.
CX("name myAA")                       # list the selection again to remind us
CX("select #1/A:298")                 # select ...
CX("delete sel")                      # ... and delete it
CX("surface #1")                      # compute a surface
CX("color #1 #FFB7C5 target s")       # color it cherry-blossom pink
CX("transparency #1 70 target s")     # set transparency

# More things we could do is to:
#   - view atoms that are a certain distance _around_ mAA (command "zone");
#   - show hydrogen bonds of myAA (command "hbonds");
#   - label the residue (command "label");
#   - export a high-resolution image for our journal (command "save");
#   - ... and more.


# When we are done ...
CX("remotecontrol rest stop")  # ... release the socket.



# [END]
