# tocID <- "exploringFunction.R"
#
# Purpose:  Computational Biology Foundations:
#              R code to explore the GO graph of function annotations, the
#              function annotations in GOa, and the STRING database interaction
#              information. ID mapping is done via HGNC and BioMart data.
#
# Version:  0.3
#
# Date:     2022-11
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           0.3    Include ID mapping concerns, and add STRING data.
#                  Class code 2022-11-28
#           0.2    Add GOa data; define up with a datamodel
#           0.1    In-class development
#
# TODO:
#    ...
# ==============================================================================
#


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                Line
#TOC> ----------------------------------------------------
#TOC>   1        Review: R idioms                       46
#TOC>   2        Download Gene Ontology Data            98
#TOC>   2.1        Annotation Data                     106
#TOC>   2.2        GO Term Data                        346
#TOC>   2.2.1          Download GO term data           362
#TOC>   2.2.2          Read the data                   377
#TOC>   3        Validate it                           811
#TOC>   4        STRING data                           817
#TOC>   5        Next steps                            840
#TOC>
#TOC> ==========================================================================


# Large databases store gene data, function annotations and cross-references.
# Here we explore such data, highlight the fact that care needs to be taken to
# cross-reference identifiers, and prepare a data-schema for GO terms, GO
# annotations, and STRING interaction data.


# =    1  Review: R idioms  ====================================================

# As far as R code is concerned, there are no particular difficulties, but one
# really needs to sweat through careful applications of various subsetting and
# filtering idioms. Remember: to subset is to extract some rows from a data set,
# to filter is to subset according to some criterion applied to the data. Most
# frequently, the following syntactic constructions are applied:

dat <- data.frame(a = letters, A = LETTERS, i = 1:26)

# the dollar operator extracts a data frame column as a vector
dat$A

# square brackets define a subset
dat$a[1:3]

# subsetting can also be done from a vector of logicals, applying
(sel <- dat$i < 8)
dat$A[sel]
dat$A[dat$i == 13]

# The which() function turns a vector of logicals into a vector of indices.
# Sometimes one is more convenient, sometimes the other.
which(sel)

# The %in% operator takes two vectors and returns true for each element in the
# first vector that is present in the second.
myVowels <- c("a", "i", "u", "e", "o")
myVowels %in% dat$a
(sel <- dat$a %in% myVowels)
dat$A[sel] # Note: this is the order in dat$A, not the order of myVowels

# match() is used less frequently, but it works the same way as in. However it
# returns indices.
(sel <- match(myVowels, dat$a))

# table()
set.seed(31416)
myChar <- sample(dat$a, 52, replace = TRUE)
table(myChar)

# Are there any missing?

# unique() extracts any element exactly once
unique(myChar)

# duplicated() returns TRUE for every element that re-occurs in data
(sel <- duplicated(myChar))
myChar[sel]



# =    2  Download Gene Ontology Data  =========================================

# The GO project integrates two types of information: a collection of
# attributes, the GO "terms", that capture some aspect of biological function
# (GO), and the actual genes that have these terms as an annotation (GOa) apply
# to. We need both, terms and annotations to work with GO.


# ==   2.1  Annotation Data  ===================================================

# Let us start with gene annotation information. This is somewhat complicated,
# but essentially not complex. That is: there are many annotations, but they are
# independent and the relationships between them are simple.


# Get Gene annotations
#  http://current.geneontology.org/products/pages/downloads.html
#
#  Download goa_human.gaf
#
#  For gaf format 2.2 see:
#    http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/

# Download the 115 mb large file to a directory on your own computer. For
# myself, I keep such files in a directory called "DATA" just under my home
# directory. If you assign the constant DATADIR to the path of this folder on
# your computer, then you can just work with this code; the function file.path()
# takes care of the operating system-specific details.

DATADIR <- "~/DATA"

myFile <- file.path(DATADIR, "goa_human.gaf")
file.exists(myFile)   # Confirm that this is the right path - must return TRUE.

# It is useful to have a first look at the file in a text-editor. On the Mac,
# TextEdit will do well since it handles very large files quite well. (Although
# RStudio has an excellent text editor, these files are too large. There is 5Mb
# file-size limit for the RStudio IDE.) Once you open the file, you will notice
# that there are data elements that are separated by tabs (not commas), that
# there are comments, and that there is no header (A header in a
# spreadsheet-like datafile contains the column names). Since there is no
# header, we need to define our own column names. We find suggested column names
# in the .gaf format specification (see above), but for the attributes that we
# actually use we define shorter names:


# Read the file
# =============
tmp <- read.delim(myFile,
                  header = FALSE,
                  comment.char = "!",
                  col.names = c("DB",
                                "DB_Object_ID",
                                "sym",               # Gene symbol
                                "Qualifier",
                                "GOid",              # GO ID
                                "DB_Reference",
                                "EC",                # Annotation evidence code
                                "With_or_From",
                                "Aspect",            # GO sub-ontology
                                "name",
                                "DB_Object_Synonym",
                                "DB_Object_Type",
                                "Taxon",             # Must be 9606 (human)
                                "Date",
                                "Assigned_By",
                                "Annotation_Extension",
                                "Gene_Product_Form_ID"
                  ))  # 635,739 rows

head(tmp)
tail(tmp)

# Examine the contents
#
# GOid
# ====
length(unique(tmp$GOid))  # 18883

# The three ontologies
# ====================
unique(tmp$Aspect)          # F: molecular function
# C: cellular component
# P: biological process

# Taxonomy
# ========
#
# The taxID for homo sapiens is 9606.
sum(tmp$Taxon == "taxon:9606")  # 635,172 ???
unique(tmp$Taxon) # What are these ?? Are these human genes?

# Look at three examples:

tmp[tmp$Taxon == "taxon:9606|taxon:197911", ]
# https://www.genenames.org/tools/search/#!/?query=FCN1
# https://www.genenames.org/tools/search/#!/?query=FCN3
# https://www.ncbi.nlm.nih.gov/taxonomy/197911
# ... Probably a Alphainfluenzavirus receptor ... but definitely a human gene.

tmp[tmp$Taxon == "taxon:9606|taxon:5807", ]
# https://www.genenames.org/tools/search/#!/?query=MYD88
# https://www.ncbi.nlm.nih.gov/taxonomy/5807
# Toll-like innate immunity effector defects affect immune response.


# What is the story here? Why can there be two taxon IDs? The .gaf format
# specification clarifies:
#  "For cardinality 1, the ID of the species encoding the gene product.
#   For cardinality 2, to be used only in conjunction with terms that have
#   the biological process term "multi-organism process" or the cellular
#   component term "host cell" as an ancestor. The first taxon ID should be
#   that of the organism encoding the gene or gene product, and the taxon ID
#   after the pipe should be that of the other organism in the interaction."

# Ok - these are still human genes though. We will keep them.

# Evidence codes:
# ==============
# GO is exemplary in that it carefully keeps track of WHY an annotation was
# made, what the evidence for te abnnotation was. This is encoded in the "EC"
# column of the data:
unique(tmp$EC)

# ... but what do these codes mean? According to
#   http://geneontology.org/docs/guide-go-evidence-codes/  :

EC <- data.frame(code = character(), def = character())
EC[ 1,"code"] <- "IEA"; EC[ 1,"def"] <- "Inferred from Electronic Annotation"
EC[ 2,"code"] <- "EXP"; EC[ 2,"def"] <- "Inferred from Experiment"
EC[ 3,"code"] <- "IDA"; EC[ 3,"def"] <- "Inferred from Direct Assay"
EC[ 4,"code"] <- "IPI"; EC[ 4,"def"] <- "Inferred from Physical Interaction"
EC[ 5,"code"] <- "IMP"; EC[ 5,"def"] <- "Inferred from Mutant Phenotype"
EC[ 6,"code"] <- "IGI"; EC[ 6,"def"] <- "Inferred from Genetic Interaction"
EC[ 7,"code"] <- "IEP"; EC[ 7,"def"] <- "Inferred from Expression Pattern"
EC[ 8,"code"] <- "HTP"; EC[ 8,"def"] <- "Inferred from High Throughput Experiment"
EC[ 9,"code"] <- "HDA"; EC[ 9,"def"] <- "Inferred from High Throughput Direct Assay"
EC[10,"code"] <- "HMP"; EC[10,"def"] <- "Inferred from High Throughput Mutant Phenotype"
EC[11,"code"] <- "HGI"; EC[11,"def"] <- "Inferred from High Throughput Genetic Interaction"
EC[12,"code"] <- "HEP"; EC[12,"def"] <- "Inferred from High Throughput Expression Pattern"
EC[13,"code"] <- "IBA"; EC[13,"def"] <- "Inferred from Biological aspect of Ancestor"
EC[14,"code"] <- "IBD"; EC[14,"def"] <- "Inferred from Biological aspect of Descendant"
EC[15,"code"] <- "IKR"; EC[15,"def"] <- "Inferred from Key Residues"
EC[16,"code"] <- "IRD"; EC[16,"def"] <- "Inferred from Rapid Divergence"
EC[17,"code"] <- "ISS"; EC[17,"def"] <- "Inferred from Sequence or structural Similarity"
EC[18,"code"] <- "ISO"; EC[18,"def"] <- "Inferred from Sequence Orthology"
EC[19,"code"] <- "ISA"; EC[19,"def"] <- "Inferred from Sequence Alignment"
EC[20,"code"] <- "ISM"; EC[20,"def"] <- "Inferred from Sequence Model"
EC[21,"code"] <- "IGC"; EC[21,"def"] <- "Inferred from Genomic Context"
EC[22,"code"] <- "RCA"; EC[22,"def"] <- "Inferred from Reviewed Computational Analysis"
EC[23,"code"] <- "TAS"; EC[23,"def"] <- "Traceable Author Statement"
EC[24,"code"] <- "NAS"; EC[24,"def"] <- "Non-traceable Author Statement"
EC[25,"code"] <- "IC";  EC[25,"def"] <- "Inferred by Curator"
EC[26,"code"] <- "ND";  EC[26,"def"] <- "No biological Data available (ND)"
rownames(EC) <- EC$code

# Have all been used in our dataset?
EC$code[! (EC$code %in% unique(tmp$EC))]

# How frequent is each evidence code
oPar = par(mar=c(0.5, 0.5, 0.5, 0.5))
N <- length(unique(tmp$EC))
myCols <- hcl.colors(N, "Viridis")
pie(sort(table(tmp$EC)), col = myCols)
par(oPar)

# or ...
barplot(sort(table(tmp$EC), decreasing = TRUE),
        col = myCols[N:1],
        cex.axis = 0.5,
        cex.names = 0.3)

# Does each gene in our GOa dataset have a symbol?
any(is.na(tmp$sym)) # FALSE - no symbols are missing
any(tmp$sym == "")  # TRUE  - some symbols are "" (i.e. empty)
sum(tmp$sym == "")  # 167   - how many?


# Hm. What are these genes - are they real?
sel <- which(tmp$sym == "")
tmp[sel[1:5], ]
# https://www.uniprot.org/uniprotkb/A0A2R8YE69
# https://www.uniprot.org/uniprotkb/Q9H521
# https://www.uniprot.org/uniprotkb/A0A2R8YFR7
# https://www.uniprot.org/uniprotkb/A8MVJ9
# https://www.uniprot.org/uniprotkb/A0A087WW49

# Looking at these genes, they DO appear to be mostly code for bona fide
# proteins. We will keep them, but we really need a non-empty symbol for each
# gene - since this will be our primary ID for proteins. However, UniProt IDs
# exist also for those genes that do not have a symbol:
any(tmp$DB_Object_ID[sel] == "") # FALSE

# Therefore we can use the UniProt ID as a symbol if HGNC did not assign one.
sel <- which(tmp$sym == "")
tmp$sym[sel] <- tmp$DB_Object_ID[sel]

# Gene names
# ==========
# What about gene names? Does each gene have one?
any(is.na(tmp$name)) # FALSE - no names are missing
any(tmp$name == "")  # FALSE - all genes have a name


# Redundant records
# =================
# Many of these records are redundant since they were compiled from different
# sources. We assign every row a key to make sure all symbol / GO ID
# combinations are unique. We are not selecting by evidence codes, but if we
# would be, then at this step we would need to make sure to keep the most
# "trustworthy" evidence code of all duplicate symbol / GO ID combinations.

tmp$key <- sprintf("%s|%s",
                   tmp$sym,
                   tmp$GOid)
head(tmp$key)

# Finally we remove duplicates, and keep only the four data columns we want to
# continue working with.

GOa <- tmp[! duplicated(tmp$key), c("sym",
                                    "name",
                                    "GOid",
                                    "EC")]  # 294,050 rows

length(unique(GOa$sym))   # 19818 genes
length(unique(GOa$GOid))  # 18883 annotations

# Summary of our cleaned data in the GOa data frame:

# -  Column "sym" contains HGNC gene symbols, or UniProt IDs if no gene symbol is
#    available.
# -  Column "name" contains gene names.
# -  Column "GOid" contains GO IDs for terms that were annotated to a gene.
# -  Column "EC" contains evidence codes for all annotations.

# Note that gene symbol and gene name is redundant information in this table.
# The annotation depends only on the symbol. Thus we split the table apart:

sel <- ! duplicated(GOa$sym)
gene <- data.frame(sym  = GOa$sym[sel],
                   name = GOa$name[sel])
rownames(gene) <- gene$sym

# ... and we remove the redundant colum from GOa:
GOa$name <- NULL


# ==   2.2  GO Term Data  ======================================================

# Next we need to look at the actual GO terms. We might hope that the "generic
# GO slim" can be used as  a compact representation of the GO database:

# According to http://geneontology.org/docs/go-subset-guide/ a "slim"
# is a subset of the Gene Ontology that is useful for a broad overview of
# "the range of functions and processes in a given organisms [...] genome".

# Unfortunately, this is just a haphazard collection of terms, and the existing
# connections do not define a single graph. Moreover, there is no guarantee that
# a list of genes (e.g. human genes) are annotated to exactly these terms and
# not to others. We actually need the entire GO graph if we are to ask questions
# about how different partitions of function can be obtained from the tree.


# ===   2.2.1  Download GO term data

# To work with GO, we need to download the basic version of the ontology.
# According to GO: (http://geneontology.org/docs/download-ontology/) this data
# set is called "go-basic.obo", the graph is "guaranteed to be acyclic and
# annotations can be propagated up the graph", in particular, it "excludes
# relationships that cross the 3 GO hierarchies". Relations are "is a", "part
# of", "regulates", "negatively regulates" and "positively regulates".
#
# Download the 30 mb large file to a directory on your own computer.

myFile <- file.path(DATADIR, "go-basic.obo")
file.exists(myFile)


# ===   2.2.2  Read the data
# First read

txt <- readLines(myFile)

# Examine what we have

head(txt, 50)

# [Term] and [typedef]
# ====================

# It seems that "[Term]" is a useful token to identify the GO term descriptions.
# But are there other tokens? Let's look at all lins that begin with an
# open square bracket, and table them.

table(txt[grep("^\\[", txt)]) # Look carefully at this: we are using grep() to
# find the indices of all lines in txt that start
# with a square bracket. The we use those indices
# to subset txt and return a vector of exactly
# those lines. This vector of text elements is
# passed to table() to count the occurances.

# [Term] [Typedef]
# 47397         5

# The bad news is that we need to identify terms and typedefs separately. The
# good news is that our regular expression identifies both. So to identify
# only "Term"s, we
#  - fetch all indices to "[" characters at the beginning of a line
#  - subset txt from each of the  "[" indices to the next
#  - ignore the term if it contains elements we don't want

idx <- grep("^\\[", txt)         # looks for "[" at the beginning of a line
idx <- c(idx, length(txt) + 1)   # add an index one past the last line

# test
for (i in 15:35) {
  first <- idx[i] + 1       # one line after the index
  last  <- idx[i+1] - 1     # one line before the next index
  cat(txt[first:last], sep = "\n")
  cat("=======================================================\n")
}

# Record types
# ============

# Next: what kind of information do we have in this file in the first place?
# Lines of this data are prefixed with a keyword, followed by a colon.
# To table() the keywords, we can use strsplit() to split on colons, and then
# collect the first elements of each line.

tmp <- character(length(txt))

for (i in 1:length(tmp)) {
  tmp[i] <- strsplit(txt[i], ":")[[1]][1]
}

table(tmp)

# namespace: and id:
# ==================

# This is good - but does every [Term] have _exactly_ one such namespace?
# And does it have exactly one id?

for (i in 1:(length(idx) - 1)) {
  pBar(i, length(idx))

  if (txt[idx[i]] != "[Term]") { next }

  first <- idx[i] + 1       # one line after the index
  last  <- idx[i+1] - 1     # one line before the next index
  s <- txt[first:last]

  if (sum(grepl("^namespace:", s)) != 1 ||
      sum(grepl("^id:", s)) != 1) {
    cat(s, sep="\n")
    break
  }
}

# Now, it would be great if we can simply extract IDs as row-labels for a data
# frame - but those pesky typedef records are in the way. Rather than
# special-casing them, we remove them outright. That's sane advice: rather than
# special-casing a workflow for rare exceptions, sometimes it is better to
# sanitize a dataset so it satisfies your assumptions.

# Sanitize
# ========

# In this case, I happen to know that the [typedef] records are all at the end
# of the file. But they might not be. So here is a more generic process that
# blanks all lines of an excluded subset, and then removes all empty lines. We
# get rid of the records that are not a [Term], and while we are at it, we also
# get rid of the "is_obsolete: true" records.

for (i in 1:(length(idx) - 1)) {
  pBar(i, length(idx))

  first <- idx[i]           # the line of the index
  last  <- idx[i+1] - 1     # one line before the next index

  if (txt[idx[i]] != "[Term]" ||                            # Blank non-[Term]s
      any(grepl("^is_obsolete: true", txt[first:last]))) {  # Blank obsoletes
    txt[first:last] <- ""
  }
}

# We can also blank the header lines.

txt[1:(idx[1] - 1)] <- ""

# Then we remove all blank lines:
txt <- txt[txt != ""]
length(txt)   # 453,306 ... was 543,440
head(txt, 20)

# ... validate ...
x <- grep("is_obsolete", txt)

# ... and recompute our index

idx <- grep("^\\[", txt)         # looks for "[" at the beginning of a line
idx <- c(idx, length(txt) + 1)   # add an index one past the last line

# Validate:
all(txt[idx[1:(length(idx)-1)]] == "[Term]")  # Note: if we include the last
# element, we get a NA result.


# Now we need to revisit the contents of our records - check the keywords:
tmp <- character(length(txt))

for (i in 1:length(tmp)) {
  pBar(i, length(tmp))
  tmp[i] <- strsplit(txt[i], ":")[[1]][1]
}

table(tmp)

# Much cleaner.

# Relationships
# =============

# But what are the relationships we need? "is_a" is one of them,
# but what others exist? And, does every record even have an "is_a"
# relationship?

noIs_a <- numeric()

for (i in 1:(length(idx) - 1)) {
  pBar(i, length(idx))

  first <- idx[i] + 1       # one line after the [Term]
  last  <- idx[i+1] - 1     # one line before the next [Term]
  s <- txt[first:last]

  if (sum(grepl("^is_a:", s)) == 0 ) {
    noIs_a <- c(noIs_a, i)
  }
}

# What are these [Term]s ?
#
# Terms without "is_a" relationships
#
for (i in 1:length(noIs_a)) {

  first <- idx[noIs_a[i]]             # [Term]
  last  <- idx[noIs_a[i] + 1] - 1     # one line before the next [Term]
  cat(txt[first:last], sep = "\n")
  cat("\n\n")
}


# How many [Term]s have more than one "is_a" ?

nIs_a <- numeric()

for (i in 1:(length(idx) - 1)) {
  pBar(i, length(idx))

  first <- idx[i] + 1       # one line after the [Term]
  last  <- idx[i+1] - 1     # one line before the next [Term]

  nIs_a <- c(nIs_a, sum(grepl("^is_a:", txt[first:last])))
}

table(nIs_a)


# What is the contents of the 15,467 "relationship:" records?

sel <- grep("^relationship", txt)
x <- table(txt[sel])
head(x, 20)

x <- character(length(sel))
for (i in 1:length(sel)) {
  x[i] <- strsplit(txt[sel[i]], " ")[[1]][2]
}
myRelations <- names(table(x))


# This defines the constraints on the data structure we need, to work with this
# information.

# Therefore we next build a data model: a schema that supports work with GO
# data.

# The database is a list that contains data frames - the tables that store
# the information about each entity.

GO <- list()


# Table: version
# ==============
#  - Defines the version so that functions that operate on the database can
#    validate their assumptions about the schema.
#  - atribute version: (char) the version number (as a a string)
#  - attribute date: (char) date when this version was defined. YYYY-MM-DD

GO$version <- data.frame(version = "1.0", date = "2022-11-28")


# Table: edgeType
# ===============
#  - stores the types of term-term relations in the model
#  - attribute ID:   (char) a unique mnemonic for the relation
#  - attribute type: (char) short description of the relation
#  detailed semantics: see http://geneontology.org/docs/ontology-relations/

GO$edgeType <- data.frame(ID = character(), type = character())
GO$edgeType[1, "ID"] <- "isa"; GO$edgeType[1, "type"] <- "is_a"
GO$edgeType[2, "ID"] <- "reg"; GO$edgeType[2, "type"] <- "regulates"
GO$edgeType[3, "ID"] <- "neg"; GO$edgeType[3, "type"] <- "negatively_regulates"
GO$edgeType[4, "ID"] <- "pos"; GO$edgeType[4, "type"] <- "positively_regulates"
GO$edgeType[5, "ID"] <- "prt"; GO$edgeType[5, "type"] <- "part_of"
rownames(GO$edgeType) <- GO$edgeType$type


# Table: term
# ===========
#  - stores basic information on each GO term
#  - attribute ID:   (char) the unique GOid assigned by the GO database
#  - attribute name: (char) descriptive name of the term
#  - attribute ns:   (char) "namespace", "aspect", or "sub ontology". One of
#                           {P, C, F} for "biological process", "cellular
#                           component", or "molecular function", respectively.
#  detailed semantics: see http://geneontology.org/docs/ontology-documentation/

idx <- grep("^\\[Term\\]", txt)  # index of all [Term] records
nTerm <- length(idx)             # number of [Term]s
idx <- c(idx, length(txt) + 1)   # add a final index past the last line of txt

GO$term <- data.frame(ID   = character(nTerm),
                      name = character(nTerm),
                      ns =   character(nTerm))

# Table: edge
# ===========
#  - stores the term-term relations as an edge-list of a directed graph
#  - attribute ID:       (int) a unique integer
#  - attribute from:     (char) the GOid of the source node
#  - attribute to:       (char) the GOid of the target node
#  - attribute edgeType: (char) the ID of the relationship type stored in the
#                               edgeType table.

nEdge <- sum(grepl("^is_a:", txt)) + sum(grepl("^relationship:", txt))
GO$edge <- data.frame(ID = integer(nEdge),
                      from = character(nEdge),
                      to = character(nEdge),
                      edgeType = character(nEdge))



# Fill the term and edge tables

# Helper functions
getID <- function(s) {
  s <- gsub("id: ", "", s)
  s <- gsub(" !.*$", "", s)
  return(s)
}

getName <- function(s) {
  s <- gsub("name: ", "", s)
  s <- gsub(" !.*$", "", s)
  return(s)
}

getNS <- function (s) {
  if      (s == "namespace: biological_process") { ns <- "P" }
  else if (s == "namespace: molecular_function") { ns <- "F" }
  else if (s == "namespace: cellular_component") { ns <- "C" }
  else { stop("Namespace not recognized.") }
  return(ns)
}

getIsa <- function(s) {
  s <- gsub("is_a: ", "is_a ", s)
  s <- gsub(" !.*$", "", s)
  return(s)  # bring in form <type> <term_ID>
}

getRel <- function(s) {
  s <- gsub("relationship: ", "", s)
  s <- gsub(" !.*$", "", s)
  return(s)  # bring in form <type> <term_ID>
}

splitRel <- function(s) {
  s <- unlist(strsplit(s, " "))
  s[1] <- GO$edgeType[s[1], "ID"]  # set the first element to
  # the ID of this relation
  return(s)
}

# process all [Term]s
for (i in 1:(length(idx) - 1)) {
  pBar(i, length(idx))

  first <- idx[i] + 1       # one line after the [Term]
  last  <- idx[i+1] - 1     # one line before the next [Term]
  s <- txt[first:last]

  # find the indices of the elements we want
  iTermID <- grep("^id:", s)
  iName   <- grep("^name:", s)
  iNS     <- grep("^namespace:", s)
  iIsa    <- grep("^is_a:", s)
  iRel    <- grep("^relationship:", s)

  # process term information
  GO$term$ID[i]   <- getID(s[iTermID])
  GO$term$name[i] <- getName(s[iName])
  GO$term$ns[i]   <- getNS(s[iNS])

  # process is_a information
  for (id in seq_along(iIsa)) {
    iEdge <- max(GO$edge$ID) + 1
    rel <- splitRel(getIsa(s[iIsa[id]]))
    GO$edge$ID[iEdge]       <- iEdge
    GO$edge$from[iEdge]     <- GO$term$ID[i]
    GO$edge$to[iEdge]       <- rel[2]
    GO$edge$edgeType[iEdge] <- rel[1]
  }

  # process relationship information
  for (id in seq_along(iRel)) {
    iEdge <- max(GO$edge$ID) + 1
    rel <- splitRel(getRel(s[iRel[id]]))
    GO$edge$ID[iEdge]       <- iEdge
    GO$edge$from[iEdge]     <- GO$term$ID[i]
    GO$edge$to[iEdge]       <- rel[2]
    GO$edge$edgeType[iEdge] <- rel[1]
  }
}
# This takes about ten minutes ...
rownames(GO$term) <- GO$term$ID

# Validate ...
# (i): Find the three root [Term]s: a root [Term] is not listed in
#      the GO$edge table:

sel <- which(! GO$term$ID %in% GO$edge$from)
GO$term[sel, ]

# (ii): get the distribution of isa relationships
sel <- GO$edge$edgeType == "isa"
x <- table(GO$edge$from[sel])
hist(x)
x[x == max(x)]
GO$term[names(x)[x == max(x)], ]

# Find this in txt
sel <- grep(paste("id:", names(x)[x == max(x)]), txt)
first <- idx[which(idx > sel)[1] - 1]
last <-  idx[which(idx > sel)[1]] - 1
cat(txt[first:last], sep = "\n")

# (iii): get the distribution of rel relationships
sel <- GO$edge$edgeType != "isa"
x <- table(GO$edge$from[sel])
hist(x)
x[x == max(x)]
n <- names(x[x == max(x)])
GO$term[GO$term$ID %in% n, ]

# Find the third one in txt
sel <- grep(paste("id:", n[3]), txt)
first <- idx[which(idx > sel)[1] - 1]
last <-  idx[which(idx > sel)[1]] - 1
cat(txt[first:last], sep = "\n")

# Add the annotation table to the schema
# Table: a
# ========
#  - Records the annotations of GO terms to genes
#  - attribute sym:  (chr) the gene symbol (resp. the UniProt ID if no
#                          symbol is available)
#  - attribute GOid: (chr) [Term] ID annotated to this gene
#  - attribute EC:   (chr) [Term] ID annotated to this gene

GO$a <- GOa

# Add the annotation evidence codes to the schema
# Table: EC
# ===========
#  - Records the evidence codes for the annotations
#  - attribute code: (chr) code for this evidence
#  - attribute def:  (chr) definition for this evidence

GO$EC <- EC

# Add the gene table to the schema
# Table: gene
# ===========
#  - Records gene information
#  - attribute sym:  (chr) the gene symbol (resp. the UniProt ID if no
#                          symbol is available)
#  - attribute name: (chr) canonical name for this gene symbol

GO$gene <- gene

# Save the result:
# ================

saveRDS(GO, file.path(DATADIR, "GO.rds")) # Compressed: only 3.5 Mb

#
# =    3  Validate it  =========================================================
# Are our assumptions correct
#
# Statistics:
#  - number of leafs

# =    4  STRING data  =========================================================

# https://string-db.org/cgi/download?sessionId=bAJZqRUVUTwN&species_text=Homo+sapiens

myFile <- file.path(DATADIR, "9606.protein.aliases.v11.5.txt")
file.exists(myFile)


tmp <- read.delim(myFile,
                  skip = 1,
                  col.names = c("ENSP", "ID", "source"))  # 4,213,390 rows

sel <- match(GO$gene$sym, tmp$ID)
sum(is.na(sel))  # 1,252 genes could not be mapped





# Statistics:
#  - number of leafs


# =    5  Next steps  ==========================================================

# - Information?
# - Depth?
#








# [END]
