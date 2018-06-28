#!/usr/bin/env Rscript
# Super script to run Alakazam distance to nearest tuning
#
# Author:  Jason Anthony Vander Heiden
# Date:    2017.05.26
#
# Arguments:
#   -d  Change-O formatted TSV (TAB) file.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the input filename.
#   -o  Output directory.
#       Defaults to current directory.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help.

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("methods"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("alakazam"))
suppressPackageStartupMessages(library("shazam"))

# Set defaults
NPROC <- alakazam::getnproc()

# Define commmandline arguments
opt_list <- list(make_option(c("-d", "--db"), dest="DB",
                             help="Change-O formatted TSV (TAB) file."),
                 make_option(c("-n", "--name"), dest="NAME",
                             help=paste("Sample name or run identifier which will be used as the output file prefix.",
                                        "\n\t\tDefaults to a truncated version of the input filename.")),
                 make_option(c("-o", "--outdir"), dest="OUTDIR", default=".",
                             help=paste("Output directory.", "Defaults to the sample name.")),
                 make_option(c("-p", "--nproc"), dest="NPROC", default=NPROC,
                             help=paste("Number of subprocesses for multiprocessing tools.",
                                        "\n\t\tDefaults to the available processing units.")))
# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

# Check input file
if (!("DB" %in% names(opt))) {
    stop("You must provide a Change-O database file with the -d option.")
}

# Check and fill sample name
if (!("NAME" %in% names(opt))) {
    n <- basename(opt$DB)
    opt$NAME <- tools::file_path_sans_ext(basename(opt$DB))
}

# Create output directory
if (!(dir.exists(opt$OUTDIR))) {
    dir.create(opt$OUTDIR)
}

# Load and subset data
db <- as.data.frame(readChangeoDb(opt$DB))
db <- subset(db, PRCONS %in% c("IGHM","IGHG","IGHA"))

# Calculate and plot rank-abundance curve
a <- estimateAbundance(db, group="PRCONS")
ggsave("abundance.pdf", plot(a, silent=T))

# Generate Hill diversity curve
d <- rarefyDiversity(db, group="PRCONS")
ggsave("diversity.pdf", plot(d, silent=T))

# Calculate CDR3 amino acid properties
p <- aminoAcidProperties(db,  
      seq="JUNCTION", nt=T, trim=T, 
      label="CDR3")

# V family usage by isotype and clone
v <- countGenes(db, 
      gene="V_CALL_GENOTYPED", 
      groups="PRCONS", clone="CLONE",  
      mode="family")

# Build tree from a single clone
x <- makeChangeoClone(subset(db, CLONE == 10), text_fields="PRCONS")
g <- buildPhylipLineage(x, dnapars="/usr/local/bin/dnapars")
ggsave("tree.pdf", plot(g, layout=layout_as_tree, vertex.label=V(g)$PRCONS))

# Retrieve the most ancestral sequence
getMRCA(g, root="Germline")

# Calculate distance from germline
getPathLengths(g, root="Germline")

# Calculate subtree properties
summarizeSubtrees(g, fields="PRCONS")

# Tabulate isotype edge relationships
tableEdges(g, "PRCONS", 
  exclude=c("Germline", NA))
