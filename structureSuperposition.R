# tocID <- "structureSuperposition.R"
#
# Purpose:  Computational Biology Foundations:
#              R code demonstrating the superposition of two protein
#              structures: 7ACN and
#
# Version:  1.0
#
# Date:     2022-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    First demo for tutorial.
#
#
# PREREQUISITES:
#           Work with aminoAcidProperties.R
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


# ==   1  Prepare ChimeraX  ==================================================

#  - Open a fresh, new session of a recently updated version of ChimeraX
#  - type:
#       remotecontrol rest start port 61803
#
#    (... or whatever value you have given CXPORT.)


# ==   2  Load Sus scrofa Aconitase  ===========================================

CX("open 7ACN")
CX("camera sbs")   # practice your stereo vision
# CX("camera mono")
CX("lighting soft")
CX("lighting shadows true intensity 0.8")
CX("color sequential #1 & protein target abc palette #B0C4DE:rgb(44,50,56):orchid:maroon")

# ==   2  Load Homo sapiens Aconitase  =================================

CX("open 2B3X")
CX('color sequential #2 & protein target abc palette #7FFFD4:#AFEEEE:teal:#4682B4')


# ==   2  Superimpose  =========================================================

CX("select #1:1-533")
CX("delete sel")

CX("select #2:1-635")
CX("delete sel")


CX("matchmaker #!2 to #1")

CX("cofr #2/A/792")










# When we are done ...
CX("remotecontrol rest stop")  # ... release the socket.



# [END]
