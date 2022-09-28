# .util.R
#
# CSB195-2022 Class project: utility scripts
#
# 2022-09
# boris.steipe@utoronto.ca
#
# This file is source()d upon startup by .Rprofile
#
# ==============================================================================
#


CXPORT <- 61803
cat(sprintf("The ChimeraX port (CXPORT) is now defined as %d.\n", CXPORT))

cat("Defining CX() ...")
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

  # (B) sanitize the user-entered variable cmd to be properly encoded in
  # a http message. (No blanks, some specially handled characters, ...)
  # Here we use gsub(), which seeks for a pattern defined in its first
  # argument, substitutes the contents of the second argument, in the
  # string that is identified with the third argument.
  #
  # Patterns are defined as "regular expressions".
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
cat("done.\n")



# [END]
