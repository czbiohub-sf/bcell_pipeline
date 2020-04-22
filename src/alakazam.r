## ALL SHAZAM & ALAKAZAM COMMANDS

# Import required packages
install.packages("gdata", repos = 'https://mirror.las.iastate.edu/CRAN/')
install.packages("tidyverse", repos = 'https://mirror.las.iastate.edu/CRAN/')
install.packages("hexbin", repos = 'https://mirror.las.iastate.edu/CRAN/')
install.packages("alakazam", repos = 'https://mirror.las.iastate.edu/CRAN/')
install.packages("shazam", repos = 'https://mirror.las.iastate.edu/CRAN/')
install.packages("scales", repos = 'https://mirror.las.iastate.edu/CRAN/')
install.packages("igraph", repos = 'https://mirror.las.iastate.edu/CRAN/')
#install.packages("grid", repos = 'https://mirror.las.iastate.edu/CRAN/')
install.packages("lattice", repos = 'https://mirror.las.iastate.edu/CRAN/')
install.packages("gridExtra", repos = 'https://mirror.las.iastate.edu/CRAN/')
suppressPackageStartupMessages(library(hexbin))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(dplyr))

#################################################################################################
################################ GENERIC COMMANDS FOR ALL ANALYSES ##############################
#################################################################################################
# layouts
layouthkl <- rbind(c(1,NA,2),
                   c(1,1,3))
layouthkl2 <- rbind (c(1,1,1,2,3))
layouthkl3a <- rbind(c(1,1,1),
                     c(NA,2,3))
layouthkl3 <- rbind(c(1,NA),
                    c(1,2),
                    c(1,3))
layoutgmkl <- rbind(c(1,3),
                    c(2,4))

layout3piesa <- rbind(c(1, 2, 3))
layout3pies <- rbind(c(1),
                     c(2),
                     c(3))
layout3piesc <- rbind(c(1,3),
                      c(2,4))
layout2 <- rbind(c(1),
                 c(2))

# Define commmandline arguments
opt_list <- list(make_option(c("-d", "--db"), dest="DB",
                             help="Change-O formatted TSV (TAB) file."),
                 make_option(c("-c", "--clone"), dest="CLONE", default="1",
                             help="CLONE number to make ChangeoClone"),
                 make_option(c("-o", "--outdir"), dest="OUTDIR", default=".",
                             help=paste("Output directory.", "Defaults to the sample name.")))
# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

# Create output directory
if (!(dir.exists(opt$OUTDIR))) {
  dir.create(opt$OUTDIR)
}

setwd(opt$OUTDIR)

## Sets BX as the input change-o tab file
BX <- read.delim(opt$DB)

### NOW UNIVERSAL COMMANDS - UNTIL LINEAGE SECTION
BX$PRCONS <- as.character(BX$PRCONS)
BX$CREGION <- as.character(BX$CREGION)
BX$SEQUENCE_IMGT <- as.character(BX$SEQUENCE_IMGT)
BX$GERMLINE_IMGT_D_MASK <- as.character(BX$GERMLINE_IMGT_D_MASK)
BX$JUNCTION_LENGTH2 <- round(BX$JUNCTION_LENGTH/3)*3
BX$CDRH3KABAT_LENGTH <- ((BX$JUNCTION_LENGTH2) / 3) - 4
BX$CDRH3KABAT_LENGTH[BX$CDRH3KABAT_LENGTH <= 1] <- 1
BX$CDRL3KABAT_LENGTH <- ((BX$JUNCTION_LENGTH2) / 3) - 2
BX$CDRL3KABAT_LENGTH[BX$CDRL3KABAT_LENGTH <= 1] <- 1
## NOTE FOR SOME DATASETS (FIRST NOTICED IN ONE OF LESLIE'S) MU_FREQ IS OFF - CAN ALSO USE MU_FREQ2
BX$MU_FREQ2 <- 1 - BX$V_IDENTITY
BX$FAMILY <- getFamily(BX$V_CALL, first=TRUE, strip_d=TRUE)
BX$GENE <- getGene(BX$V_CALL, first=TRUE, strip_d=TRUE)
BX$FAMILY <- as.character(BX$FAMILY)
BX$GENE <- as.character(BX$GENE)

## CHANGING TO REFLECT NEW CREGION ANNOTATION AND USING THAT TO DETERMINE PRCONS2
#BX.H <- subset(BX, PRCONS %in% c("IgM", "IgG", "IgA"))
#BX.L <- subset(BX, PRCONS %in% c("Kappa", "Lambda"))
#BX.kappa <- subset(BX, PRCONS %in% c("Kappa"))
#BX.lambda <- subset(BX, PRCONS %in% c("Lambda"))
BX.H <- subset(BX, CREGION %in% c("IgM", "IgG", "IgA"))
BX.L <- subset(BX, CREGION %in% c("Kappa", "Lambda"))
BX.kappa <- subset(BX, CREGION %in% c("Kappa"))
BX.lambda <- subset(BX, CREGION %in% c("Lambda"))

BX.H <- BX.H[ grep("IGKV", BX.H$V_CALL, invert = TRUE) , ]
BX.H <- BX.H[ grep("IGLV", BX.H$V_CALL, invert = TRUE) , ]
BX.L <- BX.L[ grep("IGHV", BX.L$V_CALL, invert = TRUE) , ]
BX.kappa <- BX.kappa[ grep("IGHV", BX.kappa$V_CALL, invert = TRUE) , ]
BX.kappa <- BX.kappa[ grep("IGLV", BX.kappa$V_CALL, invert = TRUE) , ]
BX.lambda <- BX.lambda[ grep("IGHV", BX.lambda$V_CALL, invert = TRUE) , ]
BX.lambda <- BX.lambda[ grep("IGKV", BX.lambda$V_CALL, invert = TRUE) , ]

#BX$PRCONS2 <- factor(BX$PRCONS, levels = c("IgM", "IgG", "IgA", "Kappa", "Lambda"))
#BX.H$PRCONS2 <- factor(BX.H$PRCONS, levels = c("IgM", "IgG", "IgA"))
#BX.L$PRCONS2 <- factor(BX.L$PRCONS, levels = c("Kappa", "Lambda"))
#BX.kappa$PRCONS2 <- factor(BX.kappa$PRCONS, levels = c("Kappa"))
#BX.lambda$PRCONS2 <- factor(BX.lambda$PRCONS, levels = c("Lambda"))
BX$PRCONS2 <- factor(BX$CREGION, levels = c("IgM", "IgG", "IgA", "Kappa", "Lambda"))
BX.H$PRCONS2 <- factor(BX.H$CREGION, levels = c("IgM", "IgG", "IgA"))
BX.L$PRCONS2 <- factor(BX.L$CREGION, levels = c("Kappa", "Lambda"))
BX.kappa$PRCONS2 <- factor(BX.kappa$CREGION, levels = c("Kappa"))
BX.lambda$PRCONS2 <- factor(BX.lambda$CREGION, levels = c("Lambda"))

#################################################################################################
################################### GENE AND GENE FAMILY PLOTS ##################################
#################################################################################################

## gene and gene family usage by clone
BX.H.family <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="family")
BX.H.gene <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="gene")
BX.L.family <- countGenes(BX.L, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="family")
BX.L.gene <- countGenes(BX.L, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="gene")
BX.kappa.family <- countGenes(BX.kappa, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="family")
BX.kappa.gene <- countGenes(BX.kappa, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="gene")
BX.lambda.family <- countGenes(BX.lambda, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="family")
BX.lambda.gene <- countGenes(BX.lambda, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="gene")
BX.lambda.gene$GENE <- gsub("LV", "LV0", BX.lambda.gene$GENE)
BX.lambda.gene$GENE <- gsub("LV010", "LV10", BX.lambda.gene$GENE)
BX.lambda.family$GENE <- gsub("LV", "LV0", BX.lambda.family$GENE)
BX.lambda.family$GENE <- gsub("LV010", "LV10", BX.lambda.family$GENE)
BX.H.gene$GF <- substring(BX.H.gene$GENE, 1,5)
BX.kappa.gene$GF <- substring(BX.kappa.gene$GENE, 1,5)
BX.lambda.gene$GF <- substring(BX.lambda.gene$GENE, 1,6)

## gene and gene family usage by read
BXs.H.family <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), mode="family")
BXs.H.gene <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), mode="gene")
BXs.L.family <- countGenes(BX.L, gene="V_CALL", groups=c("PRCONS2"), mode="family")
BXs.L.gene <- countGenes(BX.L, gene="V_CALL", groups=c("PRCONS2"), mode="gene")
BXs.kappa.family <- countGenes(BX.kappa, gene="V_CALL", groups=c("PRCONS2"), mode="family")
BXs.kappa.gene <- countGenes(BX.kappa, gene="V_CALL", groups=c("PRCONS2"), mode="gene")
BXs.lambda.family <- countGenes(BX.lambda, gene="V_CALL", groups=c("PRCONS2"), mode="family")
BXs.lambda.gene <- countGenes(BX.lambda, gene="V_CALL", groups=c("PRCONS2"), mode="gene")
BXs.lambda.gene$GENE <- gsub("LV", "LV0", BXs.lambda.gene$GENE)
BXs.lambda.gene$GENE <- gsub("LV010", "LV10", BXs.lambda.gene$GENE)
BXs.lambda.family$GENE <- gsub("LV", "LV0", BXs.lambda.family$GENE)
BXs.lambda.family$GENE <- gsub("LV010", "LV10", BXs.lambda.family$GENE)
BXs.H.gene$GF <- substring(BXs.H.gene$GENE, 1,5)
BXs.kappa.gene$GF <- substring(BXs.kappa.gene$GENE, 1,5)
BXs.lambda.gene$GF <- substring(BXs.lambda.gene$GENE, 1,6)


#################################################################################################
# Plot V family clonal usage
gfh <- ggplot(BX.H.family, aes(x=GENE, y=CLONE_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  ggtitle("Gene Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="") + scale_colour_brewer(palette = "Paired", name="") +
  geom_point(aes(sample="BX.H"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2, ncol=1)
#plot(gfh)

gfk <- ggplot(BX.kappa.family, aes(x=GENE, y=CLONE_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="") + scale_colour_brewer(palette = "Paired", name="") +
  geom_point(aes(sample="BX.kappa"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gfk)

gfl <- ggplot(BX.lambda.family, aes(x=GENE, y=CLONE_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="") + scale_colour_brewer(palette = "Paired", name="") +
  geom_point(aes(sample="BX.lambda"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gfl)

# Plot V gene clonal usage
gh <- ggplot(BX.H.gene, aes(x=GENE, y=CLONE_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  ggtitle("Gene Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size = rel(0.6))) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") +
  geom_point(aes(sample="BX.H"), size=3, alpha=1) +
  facet_wrap(~ PRCONS2, ncol=1)
#plot(gh)

gk <- ggplot(BX.kappa.gene, aes(x=GENE, y=CLONE_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size = rel(0.6))) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") +
  geom_point(aes(sample="BX.kappa"), size=3, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gk)

gl <- ggplot(BX.lambda.gene, aes(x=GENE, y=CLONE_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size = rel(0.6))) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") +
  geom_point(aes(sample="BX.lambda"), size=3, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gl)

gfplots1 <- grid.arrange(gfh,gfk,gfl, layout_matrix = layouthkl3)
ggsave("genefamilyusage_byclone.png", gfplots1, width = 16, height = 12, units = "in")
ggsave("genefamilyusage_byclone.pdf", gfplots1, width = 16, height = 12, units = "in")

gplots1 <- grid.arrange(gh,gk,gl, layout_matrix = layouthkl3)
ggsave("geneusage_byclone.png", gplots1, width = 16, height = 12, units = "in")
ggsave("geneusage_byclone.pdf", gplots1, width = 16, height = 12, units = "in")

####### BY READ INSTEAD OF BY CLONE
# Plot V gene, gene family usage by read and PRCONS
gfhs <- ggplot(BXs.H.family, aes(x=GENE, y=SEQ_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  ggtitle("Gene Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="") + scale_colour_brewer(palette = "Paired", name="") +
  geom_point(aes(sample="BX.H"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2, ncol=1)
#plot(gfhs)

gfks <- ggplot(BXs.kappa.family, aes(x=GENE, y=SEQ_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="") + scale_colour_brewer(palette = "Paired", name="") +
  geom_point(aes(sample="BX.kappa"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gfks)

gfls <- ggplot(BXs.lambda.family, aes(x=GENE, y=SEQ_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="") + scale_colour_brewer(palette = "Paired", name="") +
  geom_point(aes(sample="BX.lambda"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gfls)

# Plot V gene by sample and PRCONS
ghs <- ggplot(BXs.H.gene, aes(x=GENE, y=SEQ_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  ggtitle("Gene Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size = rel(0.6))) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") +
  geom_point(aes(sample="BX.H"), size=3, alpha=1) +
  facet_wrap(~ PRCONS2, ncol=1)
#plot(ghs)

gks <- ggplot(BXs.kappa.gene, aes(x=GENE, y=SEQ_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size = rel(0.6))) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") +
  geom_point(aes(sample="BX.kappa"), size=3, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gks)

gls <- ggplot(BXs.lambda.gene, aes(x=GENE, y=SEQ_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size = rel(0.6))) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") +
  geom_point(aes(sample="BX.lambda"), size=3, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gls)

gfsplots1 <- grid.arrange(gfhs,gfks,gfls, layout_matrix = layouthkl3)
ggsave("genefamilyusage_bysample.png", gfsplots1, width = 16, height = 12, units = "in")
ggsave("genefamilyusage_bysample.pdf", gfsplots1, width = 16, height = 12, units = "in")
gsplots1 <- grid.arrange(ghs,gks,gls, layout_matrix = layouthkl3)
ggsave("geneusage_bysample.png", gsplots1, width = 16, height = 12, units = "in")
ggsave("geneusage_bysample.pdf", gsplots1, width = 16, height = 12, units = "in")

#################################################################################################
#################################### MUTATION ANALYSIS ##########################################
#################################################################################################

####################################
###### IF RELOADING DATA - DO SO HERE!
###BX_hobs <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutstats_h.tsv")
###BX_kobs <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutstats_k.tsv")
###BX_lobs <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutstats_l.tsv")
###clonestatsh <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestats_h.tsv")
###clonestatsk <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestats_k.tsv")
###clonestatsl <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestats_l.tsv")
###clonestatshf <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestatsfiltered_h.tsv")
###clonestatskf <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestatsfiltered_k.tsv")
###clonestatslf <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestatsfiltered_l.tsv")
###BX_hobs$PRCONS2 <- factor(BX_hobs$PRCONS, levels = c("IgM", "IgD", "IgG", "IgA"))
###clonestatsh$PRCONS2 <- factor(clonestatsh$PRCONS, levels = c("IgM", "IgD", "IgG", "IgA"))
###clonestatshf$PRCONS2 <- factor(clonestatshf$PRCONS, levels = c("IgM", "IgD", "IgG", "IgA"))
###################################

# Calculate combined R and S mutation frequencies
BX_hobs <- observedMutations(BX.H, sequenceColumn="SEQUENCE_IMGT",
                             germlineColumn="GERMLINE_IMGT_D_MASK",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)

BX_kobs <- observedMutations(BX.kappa, sequenceColumn="SEQUENCE_IMGT",
                             germlineColumn="GERMLINE_IMGT_D_MASK",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)

BX_lobs <- observedMutations(BX.lambda, sequenceColumn="SEQUENCE_IMGT",
                             germlineColumn="GERMLINE_IMGT_D_MASK",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)

## RELOADING THESE TWO LAMBDA FILES TO REMOVE 0 FROM GENE NAMES...
BX.lambda.gene <- countGenes(BX.lambda, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="gene")
BXs.lambda.gene <- countGenes(BX.lambda, gene="V_CALL", groups=c("PRCONS2"), mode="gene")

# merge two data frames adding both gene counts by clone and by read...
BX_hobs <- merge(BX_hobs,BX.H.gene,by=c("PRCONS2","GENE"))
BX_kobs <- merge(BX_kobs,BX.kappa.gene,by=c("PRCONS2","GENE"))
BX_lobs <- merge(BX_lobs,BX.lambda.gene,by=c("PRCONS2","GENE"))

BX_hobs <- merge(BX_hobs,BXs.H.gene,by=c("PRCONS2","GENE"))
BX_kobs <- merge(BX_kobs,BXs.kappa.gene,by=c("PRCONS2","GENE"))
BX_lobs <- merge(BX_lobs,BXs.lambda.gene,by=c("PRCONS2","GENE"))

## renaming gene counts/weights
BX_hobs <- rename(BX_hobs, GENEFREQ_BYCLONE = CLONE_FREQ)
BX_hobs <- rename(BX_hobs, GENECOUNT_BYCLONE = CLONE_COUNT)
BX_hobs <- rename(BX_hobs, GENEFREQ_BYREAD = SEQ_FREQ)
BX_hobs <- rename(BX_hobs, GENECOUNT_BYREAD = SEQ_COUNT)
BX_hobs$GENEFREQWEIGHT_BYCLONE <- rescale(BX_hobs$GENEFREQ_BYCLONE)
BX_hobs$GENEFREQWEIGHT_BYREAD <- rescale(BX_hobs$GENEFREQ_BYREAD)
BX_kobs <- rename(BX_kobs, GENEFREQ_BYCLONE = CLONE_FREQ)
BX_kobs <- rename(BX_kobs, GENECOUNT_BYCLONE = CLONE_COUNT)
BX_kobs <- rename(BX_kobs, GENEFREQ_BYREAD = SEQ_FREQ)
BX_kobs <- rename(BX_kobs, GENECOUNT_BYREAD = SEQ_COUNT)
BX_kobs$GENEFREQWEIGHT_BYCLONE <- rescale(BX_kobs$GENEFREQ_BYCLONE)
BX_kobs$GENEFREQWEIGHT_BYREAD <- rescale(BX_kobs$GENEFREQ_BYREAD)
BX_lobs <- rename(BX_lobs, GENEFREQ_BYCLONE = CLONE_FREQ)
BX_lobs <- rename(BX_lobs, GENECOUNT_BYCLONE = CLONE_COUNT)
BX_lobs <- rename(BX_lobs, GENEFREQ_BYREAD = SEQ_FREQ)
BX_lobs <- rename(BX_lobs, GENECOUNT_BYREAD = SEQ_COUNT)
BX_lobs$GENEFREQWEIGHT_BYCLONE <- rescale(BX_lobs$GENEFREQ_BYCLONE)
BX_lobs$GENEFREQWEIGHT_BYREAD <- rescale(BX_lobs$GENEFREQ_BYREAD)

#BX_hobs$MU_WEIGHT <- rescale(BX_hobs$MU_FREQ)
#BX_kobs$MU_WEIGHT <- rescale(BX_kobs$MU_FREQ)
#BX_lobs$MU_WEIGHT <- rescale(BX_lobs$MU_FREQ)

## only after merging should we change LV names for plotting
BX_lobs$GENE <- gsub("LV", "LV0", BX_lobs$GENE)
BX_lobs$GENE <- gsub("LV010", "LV10", BX_lobs$GENE)
BX_lobs$FAMILY <- gsub("LV", "LV0", BX_lobs$FAMILY)
BX_lobs$FAMILY <- gsub("LV010", "LV10", BX_lobs$FAMILY)

BX.lambda.gene$GENE <- gsub("LV", "LV0", BX.lambda.gene$GENE)
BX.lambda.gene$GENE <- gsub("LV010", "LV10", BX.lambda.gene$GENE)
BX.lambda.gene$GF <- substring(BX.lambda.gene$GENE, 1,6)
BXs.lambda.gene$GENE <- gsub("LV", "LV0", BXs.lambda.gene$GENE)
BXs.lambda.gene$GENE <- gsub("LV010", "LV10", BXs.lambda.gene$GENE)
BXs.lambda.gene$GF <- substring(BXs.lambda.gene$GENE, 1,6)

###################################
### moved all subsequent calculations here - then save .tsv files...
BX_hobs <- BX_hobs %>% add_count(CLONE)
BX_kobs <- BX_kobs %>% add_count(CLONE)
BX_lobs <- BX_lobs %>% add_count(CLONE)

## removing GF.x and GF.y
BX_hobs$GF.x <- NULL
BX_hobs$GF.y <- NULL
BX_kobs$GF.x <- NULL
BX_kobs$GF.y <- NULL

## TRY ADDING AS.CHARACTER HERE TO PRCONS2, GANDA_SUBTYPE (PRCONS AND CREGION ALREADY ARE)
BX_hobs$PRCONS2 <- as.character(BX_hobs$PRCONS2)
BX_hobs$GANDA_SUBTYPE <- as.character(BX_hobs$GANDA_SUBTYPE)
BX_kobs$PRCONS2 <- as.character(BX_kobs$PRCONS2)
BX_lobs$PRCONS2 <- as.character(BX_lobs$PRCONS2)

### creating lists per clone not per read 
## NOTE UPDATING 5/1/19 TO USE MODE FOR PRCONS2 AND GANDA_SUBTYPE
#clonestatsh1 <- BX_hobs %>%
#  group_by(CLONE) %>%
#  summarize_at(c("PRCONS2","GENE","V_CALL","D_CALL","J_CALL","JUNCTION_LENGTH","PRCONS","GANDA_SUBTYPE","JUNCTION_LENGTH2","CDRH3KABAT_LENGTH","FAMILY"), first)
#clonestatsh2 <- BX_hobs %>%
#  group_by(CLONE) %>%
#  summarize_if(is.numeric, mean)
#clonestatsh3 <- BX_hobs %>%
#  group_by(CLONE) %>%
#  summarize_at(c("PRCONS2","GANDA_SUBTYPE"), n_distinct)
#clonestatsh3 <- rename(clonestatsh3, PRCONS2_d = PRCONS2)
#clonestatsh3 <- rename(clonestatsh3, GANDA_SUBTYPE_d = GANDA_SUBTYPE)
#clonestatsh4 <- inner_join(clonestatsh1, clonestatsh2)
#clonestatsh <- inner_join(clonestatsh4, clonestatsh3)

clonestatsh1 <- BX_hobs %>%
  group_by(CLONE) %>%
  summarize_at(c("PRCONS2","GENE","V_CALL","D_CALL","J_CALL","JUNCTION_LENGTH","PRCONS","GANDA_SUBTYPE","JUNCTION_LENGTH2","CDRH3KABAT_LENGTH","FAMILY"), first)
clonestatsh2 <- BX_hobs %>%
  group_by(CLONE) %>%
  summarize_if(is.numeric, mean)
clonestatsh3 <- BX_hobs %>%
  group_by(CLONE) %>%
  summarize_at(c("PRCONS2","GANDA_SUBTYPE"), n_distinct)
#clonestatsh3b <- BX_hobs %>%
#  group_by(CLONE) %>%
#  summarize_at(c("DUPCOUNT","CONSCOUNT"), max)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
clonestatsh5 <- BX_hobs %>%
  group_by(CLONE) %>%
  summarize_at(c("PRCONS2","GANDA_SUBTYPE"), Mode)
clonestatsh1 <- rename(clonestatsh1, PRCONS2f = PRCONS2)
clonestatsh1 <- rename(clonestatsh1, GANDA_SUBTYPEf = GANDA_SUBTYPE)
clonestatsh3 <- rename(clonestatsh3, PRCONS2_d = PRCONS2)
clonestatsh3 <- rename(clonestatsh3, GANDA_SUBTYPE_d = GANDA_SUBTYPE)
clonestatsh4 <- inner_join(clonestatsh1, clonestatsh2)
clonestatsh6 <- inner_join(clonestatsh4, clonestatsh3)
clonestatsh <- inner_join(clonestatsh6, clonestatsh5)

clonestatsk1 <- BX_kobs %>%
  group_by(CLONE) %>%
  summarize_at(c("PRCONS2","GENE","V_CALL","J_CALL","JUNCTION_LENGTH","PRCONS","JUNCTION_LENGTH2","CDRL3KABAT_LENGTH","FAMILY"), first)
clonestatsk2 <- BX_kobs %>%
  group_by(CLONE) %>%
  summarize_if(is.numeric, mean)
clonestatsk <- inner_join(clonestatsk1, clonestatsk2)

clonestatsl1 <- BX_lobs %>%
  group_by(CLONE) %>%
  summarize_at(c("PRCONS2","GENE","V_CALL","J_CALL","JUNCTION_LENGTH","PRCONS","JUNCTION_LENGTH2","CDRL3KABAT_LENGTH","FAMILY"), first)
clonestatsl2 <- BX_lobs %>%
  group_by(CLONE) %>%
  summarize_if(is.numeric, mean)
clonestatsl <- inner_join(clonestatsl1, clonestatsl2)

## filter to remove clones found only once!
clonestatshf <- clonestatsh %>% filter(n > 1)
clonestatskf <- clonestatsk %>% filter(n > 1)
clonestatslf <- clonestatsl %>% filter(n > 1)

###################
## getting some summary stats for each isotype...
hmeans1 <- BX_hobs %>%
  group_by(PRCONS2) %>%
  summarize(mutation_mean_reads = mean(MU_FREQ), n1 = n())
hsds1 <- BX_hobs %>%
  group_by(PRCONS2) %>%
  summarize(mutation_sd_reads = sd(MU_FREQ), n1b = n())
hmeans2 <- clonestatsh %>%
  group_by(PRCONS2) %>%
  summarize(mutation_mean_clones = mean(MU_FREQ), n2 = n())
hsds2 <- clonestatsh %>%
  group_by(PRCONS2) %>%
  summarize(mutation_sd_clones = sd(MU_FREQ), n2b = n())
hmeans3 <- clonestatshf %>%
  group_by(PRCONS2) %>%
  summarize(mutation_mean_filteredclones = mean(MU_FREQ), n3 = n())
hsds3 <- clonestatshf %>%
  group_by(PRCONS2) %>%
  summarize(mutation_sd_filteredclones = sd(MU_FREQ), n3b = n())
#hmeans_and_sds <- bind_cols(hmeans1,hsds1,hmeans2,hsds2,hmeans3,hsds3)
hmeans_and_sds <- hmeans1 %>%
  left_join(hsds1, by='PRCONS2') %>%
  left_join(hmeans2, by='PRCONS2') %>%
  left_join(hsds2, by='PRCONS2') %>%
  left_join(hmeans3, by='PRCONS2') %>%
  left_join(hsds3, by='PRCONS2')

kmeans1 <- BX_kobs %>%
  group_by(PRCONS2) %>%
  summarize(mutation_mean_reads = mean(MU_FREQ), n1 = n())
ksds1 <- BX_kobs %>%
  group_by(PRCONS2) %>%
  summarize(mutation_sd_reads = sd(MU_FREQ), n1b = n())
kmeans2 <- clonestatsk %>%
  group_by(PRCONS2) %>%
  summarize(mutation_mean_clones = mean(MU_FREQ), n2 = n())
ksds2 <- clonestatsk %>%
  group_by(PRCONS2) %>%
  summarize(mutation_sd_clones = sd(MU_FREQ), n2b = n())
kmeans3 <- clonestatskf %>%
  group_by(PRCONS2) %>%
  summarize(mutation_mean_filteredclones = mean(MU_FREQ), n3 = n())
ksds3 <- clonestatskf %>%
  group_by(PRCONS2) %>%
  summarize(mutation_sd_filteredclones = sd(MU_FREQ), n3b = n())
#kmeans_and_sds <- bind_cols(kmeans1,ksds1,kmeans2,ksds2,kmeans3,ksds3)
kmeans_and_sds <- kmeans1 %>%
  left_join(ksds1, by='PRCONS2') %>%
  left_join(kmeans2, by='PRCONS2') %>%
  left_join(ksds2, by='PRCONS2') %>%
  left_join(kmeans3, by='PRCONS2') %>%
  left_join(ksds3, by='PRCONS2')

lmeans1 <- BX_lobs %>%
  group_by(PRCONS2) %>%
  summarize(mutation_mean_reads = mean(MU_FREQ), n1 = n())
lsds1 <- BX_lobs %>%
  group_by(PRCONS2) %>%
  summarize(mutation_sd_reads = sd(MU_FREQ), n1b = n())
lmeans2 <- clonestatsl %>%
  group_by(PRCONS2) %>%
  summarize(mutation_mean_clones = mean(MU_FREQ), n2 = n())
lsds2 <- clonestatsl %>%
  group_by(PRCONS2) %>%
  summarize(mutation_sd_clones = sd(MU_FREQ), n2b = n())
lmeans3 <- clonestatslf %>%
  group_by(PRCONS2) %>%
  summarize(mutation_mean_filteredclones = mean(MU_FREQ), n3 = n())
lsds3 <- clonestatslf %>%
  group_by(PRCONS2) %>%
  summarize(mutation_sd_filteredclones = sd(MU_FREQ), n3b = n())
#lmeans_and_sds <- bind_cols(lmeans1,lsds1,lmeans2,lsds2,lmeans3,lsds3)
lmeans_and_sds <- lmeans1 %>%
  left_join(lsds1, by='PRCONS2') %>%
  left_join(lmeans2, by='PRCONS2') %>%
  left_join(lsds2, by='PRCONS2') %>%
  left_join(lmeans3, by='PRCONS2') %>%
  left_join(lsds3, by='PRCONS2')
allmeans_and_sds <- bind_rows(hmeans_and_sds,kmeans_and_sds,lmeans_and_sds)

### do same for subtype
submeans1 <- BX_hobs %>%
  group_by(GANDA_SUBTYPE) %>%
  summarize(mutation_mean_reads = mean(MU_FREQ), n1 = n())
#  summarize(mutation_mean_reads = mean(MU_FREQ), n1 = n()) %>%
#  replace_na(GANDA_SUBTYPE = "OTHER")
subsds1 <- BX_hobs %>%
  group_by(GANDA_SUBTYPE) %>%
  summarize(mutation_sd_reads = sd(MU_FREQ), n1b = n())
submeans2 <- clonestatsh %>%
  group_by(GANDA_SUBTYPE) %>%
  summarize(mutation_mean_clones = mean(MU_FREQ), n2 = n())
subsds2 <- clonestatsh %>%
  group_by(GANDA_SUBTYPE) %>%
  summarize(mutation_sd_clones = sd(MU_FREQ), n2b = n())
submeans3 <- clonestatshf %>%
  group_by(GANDA_SUBTYPE) %>%
  summarize(mutation_mean_filteredclones = mean(MU_FREQ), n3 = n())
subsds3 <- clonestatshf %>%
  group_by(GANDA_SUBTYPE) %>%
  summarize(mutation_sd_filteredclones = sd(MU_FREQ), n3b = n())
## TRYING REPLACE_NA - DOESN'T SEEM TO WORK - PLAN B TO USE DROP_NA, BUT INSTEAD JUST SUBSET EARLIER...AFTER ADDING AS.CHARACTER REPLACE_NA WORKS!
submeans1 <- replace_na(submeans1, list(GANDA_SUBTYPE = "other"))
submeans2 <- replace_na(submeans2, list(GANDA_SUBTYPE = "other"))
submeans3 <- replace_na(submeans3, list(GANDA_SUBTYPE = "other"))
subsds1 <- replace_na(subsds1, list(GANDA_SUBTYPE = "other"))
subsds2 <- replace_na(subsds2, list(GANDA_SUBTYPE = "other"))
subsds3 <- replace_na(subsds3, list(GANDA_SUBTYPE = "other"))

#submeans_and_sds <- bind_cols(submeans1,subsds1,submeans2,subsds2,submeans3,subsds3)
submeans_and_sds <- submeans1 %>%
  left_join(subsds1, by='GANDA_SUBTYPE') %>%
  left_join(submeans2, by='GANDA_SUBTYPE') %>%
  left_join(subsds2, by='GANDA_SUBTYPE') %>%
  left_join(submeans3, by='GANDA_SUBTYPE') %>%
  left_join(subsds3, by='GANDA_SUBTYPE')
## Oct 1 Alakazam errors occurring around here...first replacing cbind with bind_cols
## after this replacement the bind_cols for lmeans_and_sds works, but the second one just here gives error:
## PROBABABLY BECAUSE OF NA IN SUBTYPES (EVERYTHING NOT IGG OR IGA) - TRY SET_NAMES TO RENAME NA TO OTHER...
##        Error in cbind_all(x) : Argument 5 must be length 7, not 6     Calls: bind_cols -> cbind_all -> .Call

## 9p 10/2 - still getting error with bind_cols below in reflow version: think because of NAs still in some files below...

###################################
### stats for pie charts
lcmeans_and_sds0 <- bind_rows(kmeans_and_sds,lmeans_and_sds)
hcmeans_and_sds <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM", "IgG", "IgA"))


## 10/4 - taking subsets from here now...
mumeans_and_sds <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM"))
mulcmeans_and_sds <- bind_rows(mumeans_and_sds,kmeans_and_sds,lmeans_and_sds)
## NOTE WARNING MESSAGE HERE "IN BIND_ROWS_ BINDING FACTOR AND CHARACTER VECTOR, COERCING INTO CHARACTER VECTOR..."
#mulcmeans_and_sds <- rename(mulcmeans_and_sds, subtypeandn1 = isotypeandn1)
#mulcmeans_and_sds <- rename(mulcmeans_and_sds, subtypeandn2 = isotypeandn2)
mulcmeans_and_sds <- rename(mulcmeans_and_sds, GANDA_SUBTYPE = PRCONS2)
mulcmeans_and_sds5 <- select(mulcmeans_and_sds, GANDA_SUBTYPE,n1,n2,n3,mutation_sd_reads,mutation_mean_reads,mutation_sd_clones,mutation_mean_clones,mutation_sd_filteredclones,mutation_mean_filteredclones)
#mumeans_and_sds <- rename(mumeans_and_sds, GANDA_SUBTYPE = PRCONS2)
#rm(isosubmeans_and_sds)
## ONLY FACTOR AFTER ALL CALCULATIONS ARE DONE BEFORE GRAPHING!!
#allmeans_and_sds$PRCONS2 <- factor(allmeans_and_sds$PRCONS2, levels = c("IgM", "IgG", "IgA", "Kappa", "Lambda"))

### MAKING NEW TABLE COMBINING ISOTYPES AND SUBTYPES - step 1 subset before mutating all
allmeans_and_sds0 <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM", "IgG", "IgA", "Kappa", "Lambda"))
allmeans_and_sds <- allmeans_and_sds0 %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(PRCONS2, n1, col="isotypeandn1", sep=", ", remove = FALSE) %>%
  unite(PRCONS2, n2, col="isotypeandn2", sep=", ", remove = FALSE)
lcmeans_and_sds <- lcmeans_and_sds0 %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(PRCONS2, n1, col="isotypeandn1", sep=", ", remove = FALSE) %>%
  unite(PRCONS2, n2, col="isotypeandn2", sep=", ", remove = FALSE)
#hcmeans_and_sds <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM", "IgG", "IgA"))
hcmeans_and_sds <- hcmeans_and_sds %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(PRCONS2, n1, col="isotypeandn1", sep=", ", remove = FALSE) %>%
  unite(PRCONS2, n2, col="isotypeandn2", sep=", ", remove = FALSE)

#####
subtypemeans_and_sds0 <- subset(submeans_and_sds, GANDA_SUBTYPE %in% c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4"))
subtypemeans_and_sds0$GANDA_SUBTYPE <- gsub("IGHA", "IgA", subtypemeans_and_sds0$GANDA_SUBTYPE)
subtypemeans_and_sds0$GANDA_SUBTYPE <- gsub("IGHG", "IgG", subtypemeans_and_sds0$GANDA_SUBTYPE)

subtypemeans_and_sds <- subtypemeans_and_sds0 %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(GANDA_SUBTYPE, n1, col="subtypeandn1", sep=", ", remove = FALSE) %>%
  unite(GANDA_SUBTYPE, n2, col="subtypeandn2", sep=", ", remove = FALSE)
#rm(isotypeandsubtypemeans_and_sds)
### same color scheme for subtypes - ALSO FORCES ADD OF IGG4
target2 <- tibble(GANDA_SUBTYPE = c("IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
subtypemeans_and_sds1 <- left_join(data.frame(GANDA_SUBTYPE=target2),subtypemeans_and_sds0,by="GANDA_SUBTYPE")
subtypemeans_and_sds1 <- subtypemeans_and_sds1 %>% mutate_at(vars(n1,n1b,n2,n2b,n3,n3b,mutation_sd_reads,mutation_mean_reads,mutation_sd_clones,mutation_mean_clones,mutation_sd_filteredclones,mutation_mean_filteredclones), replace_na, 0)
subtypemeans_and_sds1 <- subtypemeans_and_sds1 %>% replace(is.na(.), "IgG4")
#subtypemeans_and_sds1 <- replace_na(subtypemeans_and_sds1, list(n1 = 0))
#subtypemeans_and_sds1 <- replace_na(subtypemeans_and_sds1, list(n2 = 0))
#subtypemeans_and_sds1$subtypeandn1 <- gsub("NA", "0", subtypemeans_and_sds1$subtypeandn1)
#subtypemeans_and_sds1$subtypeandn2 <- gsub("NA", "0", subtypemeans_and_sds1$subtypeandn2)
subtypemeans_and_sds1 <- subtypemeans_and_sds1 %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(GANDA_SUBTYPE, n1, col="subtypeandn1", sep=", ", remove = FALSE) %>%
  unite(GANDA_SUBTYPE, n2, col="subtypeandn2", sep=", ", remove = FALSE)
#subtypemeans_and_sds1$GANDA_SUBTYPE <- factor(subtypemeans_and_sds1$GANDA_SUBTYPE, levels = c("IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
## 10/4 - taking subsets from above now...

### MAKING NEW TABLE COMBINING ISOTYPES AND SUBTYPES
#mumeans_and_sds <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM"))
### NOTE IF NO IgG4, need another set of subytpe plots - alternately try replacing all NAs with 0s

### NEW IDEA - SELECT ONLY ROWS TO BIND, THEN DO NOT HAVE TO REMOVE SO MANY...
subtypemeans_and_sds5 <- select(subtypemeans_and_sds0, GANDA_SUBTYPE,n1,n2,n3,mutation_sd_reads,mutation_mean_reads,mutation_sd_clones,mutation_mean_clones,mutation_sd_filteredclones,mutation_mean_filteredclones)
isosubmeans_and_sds0 <- bind_rows(subtypemeans_and_sds5,mulcmeans_and_sds5)

#isosubmeans_and_sds0 <- isosubmeans_and_sds0 %>% mutate_at(vars(GANDA_SUBTYPE1,GANDA_SUBTYPE2,GANDA_SUBTYPE3,GANDA_SUBTYPE4,GANDA_SUBTYPE5,PRCONS21,PRCONS22,PRCONS23,PRCONS24,PRCONS25), replace_na, "other")
#isosubmeans_and_sds0 <- isosubmeans_and_sds0 %>% replace(is.na(.), "other")

#isosubmeans_and_sds0$GANDA_SUBTYPE <- factor(isosubmeans_and_sds0$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "Kappa", "Lambda"))
#isosubmeans_and_sds1 <- isosubmeans_and_sds0
### this re-orders the rows...
target <- tibble(GANDA_SUBTYPE = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "Kappa", "Lambda"))
isosubmeans_and_sds <- left_join(data.frame(GANDA_SUBTYPE=target),isosubmeans_and_sds0,by="GANDA_SUBTYPE")
## NOTE WARNING MESSAGE: "Column `GANDA_SUBTYPE` joining character vector and factor, coercing into character vector"
#isosubmeans_and_sds <- replace_na(isosubmeans_and_sds, list(n1 = 0))
#isosubmeans_and_sds <- replace_na(isosubmeans_and_sds, list(n2 = 0))
#isosubmeans_and_sds$subtypeandn1 <- gsub("NA", "0", isosubmeans_and_sds$subtypeandn1)
#isosubmeans_and_sds$subtypeandn2 <- gsub("NA", "0", isosubmeans_and_sds$subtypeandn2)
#isosubmeans_and_sds <- isosubmeans_and_sds %>% mutate_at(vars(n1,n2,n3,mutation_sd_reads,mutation_mean_reads,mutation_sd_clones,mutation_mean_clones,mutation_sd_filteredclones,mutation_mean_filteredclones), replace_na, 0)
isosubmeans_and_sds <- isosubmeans_and_sds %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(GANDA_SUBTYPE, n1, col="subtypeandn1", sep=", ", remove = FALSE) %>%
  unite(GANDA_SUBTYPE, n2, col="subtypeandn2", sep=", ", remove = FALSE)
#isosubmeans_and_sds <- isosubmeans_and_sds %>% mutate_at(vars(GANDA_SUBTYPE1,GANDA_SUBTYPE2,GANDA_SUBTYPE3,GANDA_SUBTYPE4,GANDA_SUBTYPE5,PRCONS21,PRCONS22,PRCONS23,PRCONS24,PRCONS25), replace_na, "IgG4")
isosubmeans_and_sds <- isosubmeans_and_sds %>% replace(is.na(.), "IgG4")

#isosubmeans_and_sds$GANDA_SUBTYPE <- factor(isosubmeans_and_sds$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "Kappa", "Lambda"))
### THIS MIGHT BE WHERE RSCRIPT IN REFLOW CRASHES??!?? - do not need if replace_na is just before it??
## practicing na replaces
#subtypemeans_and_sds1 <- subtypemeans_and_sds1 %>% replace(is.na(.), 0)
#isosubmeans_and_sds1 <- isosubmeans_and_sds1 %>% replace(is.na(.), "other")
#isosubmeans_and_sds0 <- isosubmeans_and_sds0 %>% replace(is.na(.), "other")

## now HC isotypes and subtypes
hcisosubmeans_and_sds0 <- subset(isosubmeans_and_sds0, GANDA_SUBTYPE %in% c("IgA1", "IgA2", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

#hcisosubmeans_and_sds0 <- bind_rows(subtypemeans_and_sds0,mumeans_and_sds)
#hcisosubmeans_and_sds0$GANDA_SUBTYPE <- factor(hcisosubmeans_and_sds0$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
#hcisosubmeans_and_sds1 <- hcisosubmeans_and_sds0
### this re-orders the rows...
htarget <- tibble(GANDA_SUBTYPE = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
hcisosubmeans_and_sds <- left_join(data.frame(GANDA_SUBTYPE=htarget),hcisosubmeans_and_sds0,by="GANDA_SUBTYPE")
hcisosubmeans_and_sds <- hcisosubmeans_and_sds %>% mutate_at(vars(n1,n2,n3,mutation_sd_reads,mutation_mean_reads,mutation_sd_clones,mutation_mean_clones,mutation_sd_filteredclones,mutation_mean_filteredclones), replace_na, 0)

#hcisosubmeans_and_sds <- replace_na(hcisosubmeans_and_sds, list(n1 = 0))
#hcisosubmeans_and_sds <- replace_na(hcisosubmeans_and_sds, list(n2 = 0))
#hcisosubmeans_and_sds$subtypeandn1 <- gsub("NA", "0", hcisosubmeans_and_sds$subtypeandn1)
#hcisosubmeans_and_sds$subtypeandn2 <- gsub("NA", "0", hcisosubmeans_and_sds$subtypeandn2)
hcisosubmeans_and_sds <- hcisosubmeans_and_sds %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(GANDA_SUBTYPE, n1, col="subtypeandn1", sep=", ", remove = FALSE) %>%
  unite(GANDA_SUBTYPE, n2, col="subtypeandn2", sep=", ", remove = FALSE)
#hcisosubmeans_and_sds$GANDA_SUBTYPE <- factor(hcisosubmeans_and_sds$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))

## taking this out in reflow version - just save the two files separately...
#isotypeandsubtypemeans_and_sds <- bind_rows(allmeans_and_sds,subtypemeans_and_sds)

## adding tibbles to reorder rows in allmeans and hcmeans
alltarget <- tibble(PRCONS2 = c("IgM", "IgG", "IgA", "Kappa", "Lambda"))
allmeans_and_sds <- left_join(data.frame(PRCONS2=alltarget),allmeans_and_sds,by="PRCONS2")
hctarget <- tibble(PRCONS2 = c("IgM", "IgG", "IgA"))
hcmeans_and_sds <- left_join(data.frame(PRCONS2=hctarget),hcmeans_and_sds,by="PRCONS2")

## ADDING COLLAPSE INTO CLONES FROM SHAZAM SCRIPTS, THEN SAVING AS TSV FILES
clonesh <- collapseClones(BX_hobs, regionDefinition=IMGT_V, 
                          method="thresholdedFreq", minimumFrequency=0.6,
                          includeAmbiguous=FALSE, breakTiesStochastic=FALSE, 
                          nproc=4)
clonesk <- collapseClones(BX_kobs, regionDefinition=IMGT_V, 
                          method="thresholdedFreq", minimumFrequency=0.6,
                          includeAmbiguous=FALSE, breakTiesStochastic=FALSE, 
                          nproc=4)
clonesl <- collapseClones(BX_lobs, regionDefinition=IMGT_V, 
                          method="thresholdedFreq", minimumFrequency=0.6,
                          includeAmbiguous=FALSE, breakTiesStochastic=FALSE, 
                          nproc=4)
write.table(clonesh, "onereadperclone_h.tsv", sep = "\t", row.names = FALSE)
write.table(clonesk, "onereadperclone_k.tsv", sep = "\t", row.names = FALSE)
write.table(clonesl, "onereadperclone_l.tsv", sep = "\t", row.names = FALSE)

## saving all created files as .tsv files..
write.table(BX_hobs, "mutstats_h.tsv", sep = "\t", row.names = FALSE)
write.table(BX_kobs, "mutstats_k.tsv", sep = "\t", row.names = FALSE)
write.table(BX_lobs, "mutstats_l.tsv", sep = "\t", row.names = FALSE)
write.table(clonestatsh, "mutclonestats_h.tsv", sep = "\t", row.names = FALSE)
write.table(clonestatsk, "mutclonestats_k.tsv", sep = "\t", row.names = FALSE)
write.table(clonestatsl, "mutclonestats_l.tsv", sep = "\t", row.names = FALSE)
write.table(clonestatshf, "mutclonestatsfiltered_h.tsv", sep = "\t", row.names = FALSE)
write.table(clonestatskf, "mutclonestatsfiltered_k.tsv", sep = "\t", row.names = FALSE)
write.table(clonestatslf, "mutclonestatsfiltered_l.tsv", sep = "\t", row.names = FALSE)
#write.table(isotypeandsubtypemeans_and_sds, "summary_mutationstats_all.tsv", sep = "\t", row.names = FALSE)
write.table(allmeans_and_sds, "summary_mutationstats_isotypes.tsv", sep = "\t", row.names = FALSE)
write.table(subtypemeans_and_sds, "summary_mutationstats_subtypes.tsv", sep = "\t", row.names = FALSE)

#clonestatsh$N_WEIGHT <- rescale(clonestatsh$n)
#clonestatsk$N_WEIGHT <- rescale(clonestatsk$n)
#clonestatsl$N_WEIGHT <- rescale(clonestatsl$n)
#clonestatshf$N_WEIGHT <- rescale(clonestatshf$n)
#clonestatskf$N_WEIGHT <- rescale(clonestatskf$n)
#clonestatslf$N_WEIGHT <- rescale(clonestatslf$n)

## NOW FACTOR EVERYTHING FOR GRAPHING
allmeans_and_sds$PRCONS2 <- factor(allmeans_and_sds$PRCONS2, levels = c("IgM", "IgG", "IgA", "Kappa", "Lambda"))
hcmeans_and_sds$PRCONS2 <- factor(hcmeans_and_sds$PRCONS2, levels = c("IgM", "IgG", "IgA"))
#isosubmeans_and_sds0$GANDA_SUBTYPE <- factor(isosubmeans_and_sds0$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "Kappa", "Lambda"))
isosubmeans_and_sds$GANDA_SUBTYPE <- factor(isosubmeans_and_sds$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "Kappa", "Lambda"))
subtypemeans_and_sds1$GANDA_SUBTYPE <- factor(subtypemeans_and_sds1$GANDA_SUBTYPE, levels = c("IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
hcisosubmeans_and_sds0$GANDA_SUBTYPE <- factor(hcisosubmeans_and_sds0$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
hcisosubmeans_and_sds$GANDA_SUBTYPE <- factor(hcisosubmeans_and_sds$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))

BX.H$PRCONS2 <- factor(BX.H$PRCONS2, levels = c("IgM", "IgG", "IgA"))
BX_hobs$PRCONS2 <- factor(BX_hobs$PRCONS2, levels = c("IgM", "IgG", "IgA"))
clonestatsh$PRCONS2 <- factor(clonestatsh$PRCONS2, levels = c("IgM", "IgG", "IgA"))
clonestatshf$PRCONS2 <- factor(clonestatshf$PRCONS2, levels = c("IgM", "IgG", "IgA"))

####################################################################
################## NOW REST OF ALL PLOTS ###########################
####################################################################
### Mutation distribution by Gene plots by read
ghmutv <- ggplot(BX_hobs, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(ghmutv)
ghcdr3v <- ggplot(BX_hobs, aes(x=GENE, y=CDRH3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() + ggtitle("CDR3 Length by Gene") +
  xlab("Gene") + ylab("CDRH3 Length, Kabat (aa)") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(ghcdr3v)

gkmutv <- ggplot(BX_kobs, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gkmutv)
gkcdr3v <- ggplot(BX_kobs, aes(x=GENE, y=CDRL3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() +
  xlab("Gene") + ylab("CDRL3 Length, Kabat (aa)") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gkcdr3v)
glmutv <- ggplot(BX_lobs, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(glmutv)
glcdr3v <- ggplot(BX_lobs, aes(x=GENE, y=CDRL3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() +
  xlab("Gene") + ylab("CDRL3 Length, Kabat (aa)") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(glcdr3v)

mutplots1 <- grid.arrange(ghmutv,gkmutv,glmutv, layout_matrix = layouthkl3)
ggsave("mutation_bygene_byread.png", mutplots1, width = 16, height = 12, units = "in")
ggsave("mutation_bygene_byread.pdf", mutplots1, width = 16, height = 12, units = "in")
cdr3plots1 <- grid.arrange(ghcdr3v,gkcdr3v,glcdr3v, layout_matrix = layouthkl3)
ggsave("CDR3_bygene_byread.png", cdr3plots1, width = 16, height = 12, units = "in")
ggsave("CDR3_bygene_byread.pdf", cdr3plots1, width = 16, height = 12, units = "in")


### MUTATION AND CDR3 PLOTS BY GENE - BUT NOW BY CLONE (AND BY FILTERED CLONE)
ghmutvc <- ggplot(clonestatsh, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gmutv) + scale_y_continuous(labels = scales::percent)
ghcdr3vc <- ggplot(clonestatsh, aes(x=GENE, y=CDRH3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() + ggtitle("CDR3 Length by Gene") +
  xlab("Gene") + ylab("CDRH3 Length, Kabat (aa)") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gcdr3v) + scale_y_continuous(labels = scales::percent)
gkmutvc <- ggplot(clonestatsk, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gkmutv)
gkcdr3vc <- ggplot(clonestatsk, aes(x=GENE, y=CDRL3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() +
  xlab("Gene") + ylab("CDRL3 Length, Kabat (aa)") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gkcdr3v) + scale_y_continuous(labels = scales::percent)
glmutvc <- ggplot(clonestatsl, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(glmutv)
glcdr3vc <- ggplot(clonestatsl, aes(x=GENE, y=CDRL3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() +
  xlab("Gene") + ylab("CDRL3 Length, Kabat (aa)") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(glcdr3v) + scale_y_continuous(labels = scales::percent)

mutplots1c <- grid.arrange(ghmutvc,gkmutvc,glmutvc, layout_matrix = layouthkl3)
ggsave("mutation_bygene_byclone.png", mutplots1c, width = 16, height = 12, units = "in")
ggsave("mutation_bygene_byclone.pdf", mutplots1c, width = 16, height = 12, units = "in")
cdr3plots1c <- grid.arrange(ghcdr3vc,gkcdr3vc,glcdr3vc, layout_matrix = layouthkl3)
ggsave("CDR3_bygene_byclone.png", cdr3plots1c, width = 16, height = 12, units = "in")
ggsave("CDR3_bygene_byclone.pdf", cdr3plots1c, width = 16, height = 12, units = "in")
###

## filtered clone steps removed

### force mutation frequency (y-axis here) for all HC to be 0-35%
ghmutvh35 <- ggplot(BX_hobs, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(ghmutv)
mutplots1h35 <- grid.arrange(ghmutvh35,gkmutv,glmutv, layout_matrix = layouthkl3)
ggsave("mutation_bygene_byreadh35.png", mutplots1h35, width = 16, height = 12, units = "in")
ggsave("mutation_bygene_byreadh35.pdf", mutplots1h35, width = 16, height = 12, units = "in")
ghmutvch35 <- ggplot(clonestatsh, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gmutv) + scale_y_continuous(labels = scales::percent)
mutplots1ch35 <- grid.arrange(ghmutvch35,gkmutvc,glmutvc, layout_matrix = layouthkl3)
ggsave("mutation_bygene_bycloneh35.png", mutplots1ch35, width = 16, height = 12, units = "in")
ggsave("mutation_bygene_bycloneh35.pdf", mutplots1ch35, width = 16, height = 12, units = "in")

## WITH SIMPLE V-GENE MUTATION FREQ
ghmut2vh35 <- ggplot(BX_hobs, aes(x=GENE, y=MU_FREQ2, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
gkmut2v <- ggplot(BX_kobs, aes(x=GENE, y=MU_FREQ2, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
glmut2v <- ggplot(BX_lobs, aes(x=GENE, y=MU_FREQ2, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
mut2plots1h35 <- grid.arrange(ghmut2vh35,gkmut2v,glmut2v, layout_matrix = layouthkl3)
ggsave("mutation2_bygene_byreadh35.png", mut2plots1h35, width = 16, height = 12, units = "in")
ggsave("mutation2_bygene_byreadh35.pdf", mut2plots1h35, width = 16, height = 12, units = "in")

ghmut2vch35 <- ggplot(clonestatsh, aes(x=GENE, y=MU_FREQ2, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
gkmut2vc <- ggplot(clonestatsk, aes(x=GENE, y=MU_FREQ2, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
glmut2vc <- ggplot(clonestatsl, aes(x=GENE, y=MU_FREQ2, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
mut2plots1ch35 <- grid.arrange(ghmut2vch35,gkmut2vc,glmut2vc, layout_matrix = layouthkl3)
ggsave("mutation2_bygene_bycloneh35.png", mut2plots1ch35, width = 16, height = 12, units = "in")
ggsave("mutation2_bygene_bycloneh35.pdf", mut2plots1ch35, width = 16, height = 12, units = "in")


###################################
### histograms of mutation
###################################
## by read
ghmuthistogram <- ggplot(BX_hobs, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol = 1) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(ghmuthistogram)
gkmuthistogram <- ggplot(BX_kobs, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(gkmuthistogram)
glmuthistogram <- ggplot(BX_lobs, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(glmuthistogram)
### free y axis for % and not just counts below
ghmuthistogramfree <- ggplot(BX_hobs, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol = 1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(ghmuthistogramfree)

ghmuthistogramcounts <- ggplot(BX_hobs, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(ghmuthistogramcounts)
gkmuthistogramcounts <- ggplot(BX_kobs, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(gkmuthistogramcounts)
glmuthistogramcounts <- ggplot(BX_lobs, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(glmuthistogramcounts)

mutplotshist <- grid.arrange(ghmuthistogram,gkmuthistogram,glmuthistogram, layout_matrix = layouthkl3)
ggsave("mutation_histogram_hlseparatepanels.png", mutplotshist, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_hlseparatepanels.pdf", mutplotshist, width = 16, height = 12, units = "in")

mutplotshistcounts <- grid.arrange(ghmuthistogramcounts,gkmuthistogramcounts,glmuthistogramcounts, layout_matrix = layouthkl3)
ggsave("mutation_histogram_rawcounts_hlseparatepanels.png", mutplotshistcounts, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_rawcounts_hlseparatepanels.pdf", mutplotshistcounts, width = 16, height = 12, units = "in")

mutplotshistfree <- grid.arrange(ghmuthistogramfree,gkmuthistogram,glmuthistogram, layout_matrix = layouthkl3)
ggsave("mutation_histogram_hlseparatefreepanels.png", mutplotshistfree, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_hlseparatefreepanels.pdf", mutplotshistfree, width = 16, height = 12, units = "in")

### force x axis for all HC to be 0-35%
ghmuthistogramh35 <- ggplot(BX_hobs, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol = 1) + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(ghmuthistogramh35)
ghmuthistogramcountsh35 <- ggplot(BX_hobs, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(ghmuthistogramcountsh35)
mutplotshisth35 <- grid.arrange(ghmuthistogramh35,gkmuthistogram,glmuthistogram, layout_matrix = layouthkl3)
ggsave("mutation_histogram_hlseparatepanelsh35.png", mutplotshisth35, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_hlseparatepanelsh35.pdf", mutplotshisth35, width = 16, height = 12, units = "in")
mutplotshistcountsh35 <- grid.arrange(ghmuthistogramcountsh35,gkmuthistogramcounts,glmuthistogramcounts, layout_matrix = layouthkl3)
ggsave("mutation_histogram_rawcounts_hlseparatepanelsh35.png", mutplotshistcountsh35, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_rawcounts_hlseparatepanelsh35.pdf", mutplotshistcountsh35, width = 16, height = 12, units = "in")

## same mutation plots but by clone
ghmuthistogrambyclone0 <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ghmuthistogrambyclone0)
ghmuthistogrambyclone <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ghmuthistogrambyclone)
gkmuthistogrambyclone <- ggplot(clonestatsk, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(gkmuthistogrambyclone)
glmuthistogrambyclone <- ggplot(clonestatsl, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(glmuthistogrambyclone)
### free y axis for % and not just counts below
ghmuthistogrambyclonefree <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ghmuthistogrambyclonefree)


ghmuthistogramcountsbyclone0 <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(ghmuthistogramcountsbyclone0)
ghmuthistogramcountsbyclone <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(ghmuthistogramcountsbyclone)
gkmuthistogramcountsbyclone <- ggplot(clonestatsk, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(gkmuthistogramcountsbyclone)
glmuthistogramcountsbyclone <- ggplot(clonestatsl, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(glmuthistogramcountsbyclone)

mutplotshistbyclone <- grid.arrange(ghmuthistogrambyclone,gkmuthistogrambyclone,glmuthistogrambyclone, layout_matrix = layouthkl3)
ggsave("mutation_histogram_byclone_hlseparatepanels.png", mutplotshistbyclone, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_hlseparatepanels.pdf", mutplotshistbyclone, width = 16, height = 12, units = "in")
mutplotshistcountsbyclone <- grid.arrange(ghmuthistogramcountsbyclone,gkmuthistogramcountsbyclone,glmuthistogramcountsbyclone, layout_matrix = layouthkl3)
ggsave("mutation_histogram_rawcounts_byclone_hlseparatepanels.png", mutplotshistcountsbyclone, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_rawcounts_byclone_hlseparatepanels.pdf", mutplotshistcountsbyclone, width = 16, height = 12, units = "in")

mutplotshistbyclonefree <- grid.arrange(ghmuthistogrambyclonefree,gkmuthistogrambyclone,glmuthistogrambyclone, layout_matrix = layouthkl3)
ggsave("mutation_histogram_byclone_hlseparatefreepanels.png", mutplotshistbyclonefree, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_hlseparatefreepanels.pdf", mutplotshistbyclonefree, width = 16, height = 12, units = "in")

### force x axis for all HC to be 0-35%
ghmuthistogrambycloneh35 <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1) + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ghmuthistogrambycloneh35)
ghmuthistogramcountsbycloneh35 <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(ghmuthistogramcountsbycloneh35)
mutplotshistbycloneh35 <- grid.arrange(ghmuthistogrambycloneh35,gkmuthistogrambyclone,glmuthistogrambyclone, layout_matrix = layouthkl3)
ggsave("mutation_histogram_byclone_hlseparatepanelsh35.png", mutplotshistbycloneh35, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_hlseparatepanelsh35.pdf", mutplotshistbycloneh35, width = 16, height = 12, units = "in")
mutplotshistcountsbycloneh35 <- grid.arrange(ghmuthistogramcountsbycloneh35,gkmuthistogramcountsbyclone,glmuthistogramcountsbyclone, layout_matrix = layouthkl3)
ggsave("mutation_histogram_rawcounts_byclone_hlseparatepanelsh35.png", mutplotshistcountsbycloneh35, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_rawcounts_byclone_hlseparatepanelsh35.pdf", mutplotshistcountsbycloneh35, width = 16, height = 12, units = "in")

## WITH SIMPLE V-GENE MUTATION FREQ
ghmut2histogramh35 <- ggplot(BX_hobs, aes(x = MU_FREQ2)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol = 1) + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(ghmuthistogramh35)
ghmut2histogramcountsh35 <- ggplot(BX_hobs, aes(x = MU_FREQ2)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(ghmuthistogramcountsh35)
gkmut2histogram <- ggplot(BX_kobs, aes(x = MU_FREQ2)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(gkmuthistogram)
glmut2histogram <- ggplot(BX_lobs, aes(x = MU_FREQ2)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(glmuthistogram)
### free y axis for % and not just counts below
ghmut2histogramfree <- ggplot(BX_hobs, aes(x = MU_FREQ2)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol = 1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(ghmuthistogramfree)
gkmut2histogramcounts <- ggplot(BX_kobs, aes(x = MU_FREQ2)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(gkmuthistogramcounts)
glmut2histogramcounts <- ggplot(BX_lobs, aes(x = MU_FREQ2)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(glmuthistogramcounts)
mut2plotshistfree <- grid.arrange(ghmut2histogramfree,gkmut2histogram,glmut2histogram, layout_matrix = layouthkl3)
ggsave("mutation2_histogram_hlseparatefreepanels.png", mut2plotshistfree, width = 16, height = 12, units = "in")
ggsave("mutation2_histogram_hlseparatefreepanels.pdf", mut2plotshistfree, width = 16, height = 12, units = "in")
mut2plotshisth35 <- grid.arrange(ghmut2histogramh35,gkmut2histogram,glmut2histogram, layout_matrix = layouthkl3)
ggsave("mutation2_histogram_hlseparatepanelsh35.png", mut2plotshisth35, width = 16, height = 12, units = "in")
ggsave("mutation2_histogram_hlseparatepanelsh35.pdf", mut2plotshisth35, width = 16, height = 12, units = "in")
mut2plotshistcountsh35 <- grid.arrange(ghmut2histogramcountsh35,gkmut2histogramcounts,glmut2histogramcounts, layout_matrix = layouthkl3)
ggsave("mutation2_histogram_rawcounts_hlseparatepanelsh35.png", mut2plotshistcountsh35, width = 16, height = 12, units = "in")
ggsave("mutation2_histogram_rawcounts_hlseparatepanelsh35.pdf", mut2plotshistcountsh35, width = 16, height = 12, units = "in")

ghmut2histogrambycloneh35 <- ggplot(clonestatsh, aes(x = MU_FREQ2)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1) + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ghmuthistogrambycloneh35)
gkmut2histogrambyclone <- ggplot(clonestatsk, aes(x = MU_FREQ2)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(gkmuthistogrambyclone)
glmut2histogrambyclone <- ggplot(clonestatsl, aes(x = MU_FREQ2)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(glmuthistogrambyclone)
### free y axis for % and not just counts below
ghmut2histogrambyclonefree <- ggplot(clonestatsh, aes(x = MU_FREQ2)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ghmuthistogrambyclonefree)
ghmut2histogramcountsbycloneh35 <- ggplot(clonestatsh, aes(x = MU_FREQ2)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(ghmuthistogramcountsbycloneh35)
gkmut2histogramcountsbyclone <- ggplot(clonestatsk, aes(x = MU_FREQ2)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(gkmuthistogramcountsbyclone)
glmut2histogramcountsbyclone <- ggplot(clonestatsl, aes(x = MU_FREQ2)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(glmuthistogramcountsbyclone)

mut2plotshistbycloneh35 <- grid.arrange(ghmut2histogrambycloneh35,gkmut2histogrambyclone,glmut2histogrambyclone, layout_matrix = layouthkl3)
ggsave("mutation2_histogram_byclone_hlseparatepanelsh35.png", mut2plotshistbycloneh35, width = 16, height = 12, units = "in")
ggsave("mutation2_histogram_byclone_hlseparatepanelsh35.pdf", mut2plotshistbycloneh35, width = 16, height = 12, units = "in")
mut2plotshistcountsbycloneh35 <- grid.arrange(ghmut2histogramcountsbycloneh35,gkmut2histogramcountsbyclone,glmut2histogramcountsbyclone, layout_matrix = layouthkl3)
ggsave("mutation2_histogram_rawcounts_byclone_hlseparatepanelsh35.png", mut2plotshistcountsbycloneh35, width = 16, height = 12, units = "in")
ggsave("mutation2_histogram_rawcounts_byclone_hlseparatepanelsh35.pdf", mut2plotshistcountsbycloneh35, width = 16, height = 12, units = "in")
mut2plotshistbyclonefree <- grid.arrange(ghmut2histogrambyclonefree,gkmut2histogrambyclone,glmut2histogrambyclone, layout_matrix = layouthkl3)
ggsave("mutation2_histogram_byclone_hlseparatefreepanels.png", mut2plotshistbyclonefree, width = 16, height = 12, units = "in")
ggsave("mutation2_histogram_byclone_hlseparatefreepanels.pdf", mut2plotshistbyclonefree, width = 16, height = 12, units = "in")

###################
## now clone plots but filter so leaving out all single clones!
## removed this also for reflow

################
### RBIND ALL HC AND LC IN SINGLE PANELS - DO NOT LOOK AS GOOD THOUGH
## REMOVING THIS FOR REFLOW VERSION


###################################
##### mutation vs CDR3 hex plots
## by read
gmutandcdr3hexhr <- ggplot(BX_hobs, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexhr)
gmutandcdr3hexkr <- ggplot(BX_kobs, aes(x=CDRL3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads",  breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexkr)
gmutandcdr3hexlr <- ggplot(BX_lobs, aes(x=CDRL3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexlr)
mutvsCDR3r <- grid.arrange(gmutandcdr3hexhr,gmutandcdr3hexkr,gmutandcdr3hexlr, layout_matrix = layouthkl3)
ggsave("mutationvsCDR3_byread_hlseparatepanels.png", mutvsCDR3r, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byread_hlseparatepanels.pdf", mutvsCDR3r, width = 16, height = 12, units = "in")

## by clone
gmutandcdr3hexh <- ggplot(clonestatsh, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexh)
gmutandcdr3hexk <- ggplot(clonestatsk, aes(x=CDRL3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexk)
gmutandcdr3hexl <- ggplot(clonestatsl, aes(x=CDRL3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexl)
mutvsCDR3 <- grid.arrange(gmutandcdr3hexh,gmutandcdr3hexk,gmutandcdr3hexl, layout_matrix = layouthkl3)
ggsave("mutationvsCDR3_byclone_hlseparatepanels.png", mutvsCDR3, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byclone_hlseparatepanels.pdf", mutvsCDR3, width = 16, height = 12, units = "in")

### force mutation frequency (y-axis here) for all HC to be 0-35%
gmutandcdr3hexhrh35 <- ggplot(BX_hobs, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexhr)
mutvsCDR3rh35 <- grid.arrange(gmutandcdr3hexhrh35,gmutandcdr3hexkr,gmutandcdr3hexlr, layout_matrix = layouthkl3)
ggsave("mutationvsCDR3_byread_hlseparatepanelsh35.png", mutvsCDR3rh35, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byread_hlseparatepanelsh35.pdf", mutvsCDR3rh35, width = 16, height = 12, units = "in")
gmutandcdr3hexhh35 <- ggplot(clonestatsh, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexh)
mutvsCDR3h35 <- grid.arrange(gmutandcdr3hexhh35,gmutandcdr3hexk,gmutandcdr3hexl, layout_matrix = layouthkl3)
ggsave("mutationvsCDR3_byclone_hlseparatepanelsh35.png", mutvsCDR3h35, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byclone_hlseparatepanelsh35.pdf", mutvsCDR3h35, width = 16, height = 12, units = "in")


## WITH SIMPLE V-GENE MUTATION FREQ
gmut2andcdr3hexhrh35 <- ggplot(BX_hobs, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ2)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexhr)
gmut2andcdr3hexkr <- ggplot(BX_kobs, aes(x=CDRL3KABAT_LENGTH, y=MU_FREQ2)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads",  breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexkr)
gmut2andcdr3hexlr <- ggplot(BX_lobs, aes(x=CDRL3KABAT_LENGTH, y=MU_FREQ2)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
mut2vsCDR3rh35 <- grid.arrange(gmut2andcdr3hexhrh35,gmut2andcdr3hexkr,gmut2andcdr3hexlr, layout_matrix = layouthkl3)
ggsave("mutation2vsCDR3_byread_hlseparatepanelsh35.png", mut2vsCDR3rh35, width = 16, height = 12, units = "in")
ggsave("mutation2vsCDR3_byread_hlseparatepanelsh35.pdf", mut2vsCDR3rh35, width = 16, height = 12, units = "in")

gmut2andcdr3hexhh35 <- ggplot(clonestatsh, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ2)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexh)
gmut2andcdr3hexk <- ggplot(clonestatsk, aes(x=CDRL3KABAT_LENGTH, y=MU_FREQ2)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexk)
gmut2andcdr3hexl <- ggplot(clonestatsl, aes(x=CDRL3KABAT_LENGTH, y=MU_FREQ2)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
mut2vsCDR3h35 <- grid.arrange(gmut2andcdr3hexhh35,gmut2andcdr3hexk,gmut2andcdr3hexl, layout_matrix = layouthkl3)
ggsave("mutation2vsCDR3_byclone_hlseparatepanelsh35.png", mut2vsCDR3h35, width = 16, height = 12, units = "in")
ggsave("mutation2vsCDR3_byclone_hlseparatepanelsh35.pdf", mut2vsCDR3h35, width = 16, height = 12, units = "in")


  
## test violin plot
## removed for reflow version

## by filtered clone
## removed for reflow version

### SINGLE PANELS & GMKL ONLY
## removed for reflow version

#################################################################################################
################################## HISTOGRAMS OF CLONES #########################################
#################################################################################################
ghnreadshistogrambyclone <- ggplot(clonestatsh, aes(x = n)) + geom_histogram(aes(y=0.1*..density..), binwidth = 0.1) + 
  facet_wrap(~ PRCONS2, ncol=1) + scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("# Reads per Clonal Family") + xlab("# Reads per Clonal Family") + ylab("Proportion of clones")
#plot(ghnreadshistogrambyclone)
### free y axis for % and not just counts below
ghnreadshistogrambyclonefree <- ggplot(clonestatsh, aes(x = n)) + geom_histogram(aes(y=0.1*..density..), binwidth = 0.1) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("# Reads per Clonal Family") + xlab("# Reads per Clonal Family") + ylab("Proportion of clones")
#plot(ghnreadshistogrambyclonefree)
ghnreadshistogramcountsbyclone <- ggplot(clonestatsh, aes(x = n)) + geom_histogram(binwidth = 0.1) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_log10(breaks = c(1, 10, 100, 1000)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("# Reads per Clonal Family, Raw Counts") + xlab("# Reads per Clonal Family") + ylab("Clones")
#plot(ghnreadshistogramcountsbyclone)

gknreadshistogrambyclone <- ggplot(clonestatsk, aes(x = n)) + geom_histogram(aes(y=0.1*..density..), binwidth = 0.1) + 
  facet_wrap(~ PRCONS2) + scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("# Reads per Clonal Family") + ylab("Proportion of clones")
#plot(gknreadshistogrambyclone)
glnreadshistogrambyclone <- ggplot(clonestatsl, aes(x = n)) + geom_histogram(aes(y=0.1*..density..), binwidth = 0.1) + 
  facet_wrap(~ PRCONS2) + scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("# Reads per Clonal Family") + ylab("Proportion of clones")
#plot(glnreadshistogrambyclone)
gknreadshistogramcountsbyclone <- ggplot(clonestatsk, aes(x = n)) + geom_histogram(binwidth = 0.1) + 
  facet_wrap(~ PRCONS2) + scale_x_log10(breaks = c(1, 10, 100, 1000)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("# Reads per Clonal Family") + ylab("Clones")
#plot(gknreadshistogramcountsbyclone)
glnreadshistogramcountsbyclone <- ggplot(clonestatsl, aes(x = n)) + geom_histogram(binwidth = 0.1) + 
  facet_wrap(~ PRCONS2) + scale_x_log10(breaks = c(1, 10, 100, 1000)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("# Reads per Clonal Family") + ylab("Clones")
#plot(glnreadshistogramcountsbyclone)

nreadshistogrambyclone <- grid.arrange(ghnreadshistogrambyclone,gknreadshistogrambyclone,glnreadshistogrambyclone, layout_matrix = layouthkl3)
ggsave("clonalfamily_n_histogram_byclone_hlseparatepanels.png", nreadshistogrambyclone, width = 16, height = 12, units = "in")
ggsave("clonalfamily_n_histogram_byclone_hlseparatepanels.pdf", nreadshistogrambyclone, width = 16, height = 12, units = "in")
nreadshistogrambyclonefree <- grid.arrange(ghnreadshistogrambyclonefree,gknreadshistogrambyclone,glnreadshistogrambyclone, layout_matrix = layouthkl3)
ggsave("clonalfamily_n_histogram_byclone_hlfreeseparatepanels.png", nreadshistogrambyclonefree, width = 16, height = 12, units = "in")
ggsave("clonalfamily_n_histogram_byclone_hlfreeseparatepanels.pdf", nreadshistogrambyclonefree, width = 16, height = 12, units = "in")
nreadshistogramcountsbyclone <- grid.arrange(ghnreadshistogramcountsbyclone,gknreadshistogramcountsbyclone,glnreadshistogramcountsbyclone, layout_matrix = layouthkl3)
ggsave("clonalfamily_n_histogram_byclone_hlseparatepanels_counts.png", nreadshistogramcountsbyclone, width = 16, height = 12, units = "in")
ggsave("clonalfamily_n_histogram_byclone_hlseparatepanels_counts.pdf", nreadshistogramcountsbyclone, width = 16, height = 12, units = "in")

## HC plots combining all isotypes
ghnreadshistogrambycloneallh <- ggplot(clonestatsh, aes(x = n)) + geom_histogram(aes(y=0.1*..density..), binwidth = 0.1) + 
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = percent_format()) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("# Reads per Clonal Family") + xlab("# Reads per Clonal Family") + ylab("Proportion of clones")
#plot(ghnreadshistogrambycloneallh)
ghnreadshistogramcountsbycloneallh <- ggplot(clonestatsh, aes(x = n)) + geom_histogram(binwidth = 0.1) + 
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("# Reads per Clonal Family, Raw Counts") + xlab("# Reads per Clonal Family") + ylab("Clones")
#plot(ghnreadshistogramcountsbycloneallh)
ggsave("clonalfamily_n__allHC_reads.png", ghnreadshistogrambycloneallh, width = 8, height = 6, units = "in")
ggsave("clonalfamily_n__allHC_readcounts.png", ghnreadshistogramcountsbycloneallh, width = 8, height = 6, units = "in")
## EXCLUDING IgM
clonestatsga <- subset(clonestatsh, PRCONS2 != "IgM")
ghnreadshistogrambyclonega <- ggplot(clonestatsga, aes(x = n)) + geom_histogram(aes(y=0.1*..density..), binwidth = 0.1) + 
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = percent_format()) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("# Reads per Clonal Family") + xlab("# Reads per Clonal Family") + ylab("Proportion of clones")
#plot(ghnreadshistogrambyclonega)
ghnreadshistogramcountsbyclonega <- ggplot(clonestatsga, aes(x = n)) + geom_histogram(binwidth = 0.1) + 
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("# Reads per Clonal Family, Raw Counts") + xlab("# Reads per Clonal Family") + ylab("Clones")
#plot(ghnreadshistogramcountsbyclonega)
ggsave("clonalfamily_n_IgGIgA_reads.png", ghnreadshistogrambyclonega, width = 8, height = 6, units = "in")
ggsave("clonalfamily_n_IgGIgA_readcounts.png", ghnreadshistogramcountsbyclonega, width = 8, height = 6, units = "in")

## Hex plots of avg mutation vs. number of reads
gmutandnhexhbyisotypeh <- ggplot(clonestatsh, aes(x=n, y=MU_FREQ)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexhbyisotypeh)
gmutandnhexhbyisotypek <- ggplot(clonestatsk, aes(x=n, y=MU_FREQ)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexhbyisotypek)
gmutandnhexhbyisotypel <- ggplot(clonestatsl, aes(x=n, y=MU_FREQ)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexhbyisotypel)
gmutandnhexhbyisotype <- grid.arrange(gmutandnhexhbyisotypeh,gmutandnhexhbyisotypek,gmutandnhexhbyisotypel, layout_matrix = layouthkl3)
ggsave("clonalfamily_n_vs_mutation_hlseparatepanels.png", gmutandnhexhbyisotype, width = 16, height = 12, units = "in")
ggsave("clonalfamily_n_vs_mutation_hlseparatepanels.pdf", gmutandnhexhbyisotype, width = 16, height = 12, units = "in")

### force mutation frequency (y-axis here) for all HC to be 0-35%
gmutandnhexhbyisotypehh35 <- ggplot(clonestatsh, aes(x=n, y=MU_FREQ)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexhbyisotypeh)
gmutandnhexhbyisotypeh35 <- grid.arrange(gmutandnhexhbyisotypehh35,gmutandnhexhbyisotypek,gmutandnhexhbyisotypel, layout_matrix = layouthkl3)
ggsave("clonalfamily_n_vs_mutation_hlseparatepanelsh35.png", gmutandnhexhbyisotypeh35, width = 16, height = 12, units = "in")
ggsave("clonalfamily_n_vs_mutation_hlseparatepanelsh35.pdf", gmutandnhexhbyisotypeh35, width = 16, height = 12, units = "in")

## HC plots combining all isotypes
gmutandnhexallh <- ggplot(clonestatsh, aes(x=n, y=MU_FREQ)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation, HC Reads") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexallh)
ggsave("clonalfamily_n_vs_mutation_allHC_reads.png", gmutandnhexallh, width = 8, height = 6, units = "in")
## EXCLUDING IgM
gmutandnhexga <- ggplot(clonestatsga, aes(x=n, y=MU_FREQ)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation, IgG & IgA Reads") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexga)
ggsave("clonalfamily_n_vs_mutation_IgGIgA_reads.png", gmutandnhexga, width = 8, height = 6, units = "in")

gmutandnhexallga <- grid.arrange(gmutandnhexallh,gmutandnhexga, layout_matrix = layout2)
ggsave("clonalfamily_n_vs_mutation_2HCformats.png", gmutandnhexallga, width = 16, height = 12, units = "in")
ggsave("clonalfamily_n_vs_mutation_2HCformats.pdf", gmutandnhexallga, width = 16, height = 12, units = "in")


## WITH SIMPLE V-GENE MUTATION FREQ
gmut2andnhexhbyisotypehh35 <- ggplot(clonestatsh, aes(x=n, y=MU_FREQ2)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexhbyisotypeh)
gmut2andnhexhbyisotypek <- ggplot(clonestatsk, aes(x=n, y=MU_FREQ2)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexhbyisotypek)
gmut2andnhexhbyisotypel <- ggplot(clonestatsl, aes(x=n, y=MU_FREQ2)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
gmut2andnhexhbyisotypeh35 <- grid.arrange(gmut2andnhexhbyisotypehh35,gmut2andnhexhbyisotypek,gmut2andnhexhbyisotypel, layout_matrix = layouthkl3)
ggsave("clonalfamily_n_vs_mutation2_hlseparatepanelsh35.png", gmut2andnhexhbyisotypeh35, width = 16, height = 12, units = "in")
ggsave("clonalfamily_n_vs_mutation2_hlseparatepanelsh35.pdf", gmut2andnhexhbyisotypeh35, width = 16, height = 12, units = "in")

gmut2andnhexallh <- ggplot(clonestatsh, aes(x=n, y=MU_FREQ2)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation, HC Reads") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexallh)
ggsave("clonalfamily_n_vs_mutation_allHC_reads.png", gmutandnhexallh, width = 8, height = 6, units = "in")
## EXCLUDING IgM
gmut2andnhexga <- ggplot(clonestatsga, aes(x=n, y=MU_FREQ2)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation, IgG & IgA Reads") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexga)
ggsave("clonalfamily_n_vs_mutation2_IgGIgA_reads.png", gmut2andnhexga, width = 8, height = 6, units = "in")
gmut2andnhexallga <- grid.arrange(gmut2andnhexallh,gmut2andnhexga, layout_matrix = layout2)
ggsave("clonalfamily_n_vs_mutation2_2HCformats.png", gmut2andnhexallga, width = 16, height = 12, units = "in")
ggsave("clonalfamily_n_vs_mutation2_2HCformats.pdf", gmut2andnhexallga, width = 16, height = 12, units = "in")

#################################################################################################
################################## PIE CHARTS OF ISOTYPES #######################################
#################################################################################################
piechartallcounts <- ggplot(allmeans_and_sds, aes(x = "", y = n1, fill = PRCONS2, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("salmon", "purple","blue", "skyblue", "lightgreen"), labels = allmeans_and_sds$isotypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("Reads by Isotype") +
  xlab("") + ylab("Counts") + geom_text(aes(label = scales::percent(pern1)), position = position_stack(vjust = 0.5))
piechartallcounts2 <- piechartallcounts + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartallcounts2)
piechartallclones <- ggplot(allmeans_and_sds, aes(x = "", y = n2, fill = PRCONS2, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("salmon", "purple","blue", "skyblue", "lightgreen"), labels = allmeans_and_sds$isotypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("Clones by Isotype") +
  xlab("") + ylab("Clones") + geom_text(aes(label = scales::percent(pern2)), position = position_stack(vjust = 0.5))
piechartallclones2 <- piechartallclones + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartallclones2)
ggsave("piechart_all_reads.png", piechartallcounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_all_reads.pdf", piechartallcounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_all_clones.png", piechartallclones2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_all_clones.pdf", piechartallclones2, width = 8, height = 6, units = "in", bg = "transparent")

## hc only
piecharthccounts <- ggplot(hcmeans_and_sds, aes(x = "", y = n1, fill = PRCONS2, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("salmon", "purple","blue", "skyblue", "lightgreen"), labels = hcmeans_and_sds$isotypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Reads by Isotype") +
  xlab("") + ylab("Counts") + geom_text(aes(label = scales::percent(pern1)), position = position_stack(vjust = 0.5))
piecharthccounts2 <- piecharthccounts + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piecharthccounts)
piecharthcclones <- ggplot(hcmeans_and_sds, aes(x = "", y = n2, fill = PRCONS2, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("salmon", "purple","blue", "skyblue", "lightgreen"), labels = hcmeans_and_sds$isotypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Clones by Isotype") +
  xlab("") + ylab("Clones") + geom_text(aes(label = scales::percent(pern2)), position = position_stack(vjust = 0.5))
piecharthcclones2 <- piecharthcclones + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piecharthcclones)
ggsave("piechart_hc_reads.png", piecharthccounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hc_reads.pdf", piecharthccounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hc_clones.png", piecharthcclones2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hc_clones.pdf", piecharthcclones2, width = 8, height = 6, units = "in", bg = "transparent")

## lc only
piechartlccounts <- ggplot(lcmeans_and_sds, aes(x = "", y = n1, fill = PRCONS2, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("skyblue", "lightgreen"), labels = lcmeans_and_sds$isotypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("LC Reads by Isotype") +
  xlab("") + ylab("Counts") + geom_text(aes(label = scales::percent(pern1)), position = position_stack(vjust = 0.5))
piechartlccounts2 <- piechartlccounts + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartlccounts)
piechartlcclones <- ggplot(lcmeans_and_sds, aes(x = "", y = n2, fill = PRCONS2, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("skyblue", "lightgreen"), labels = lcmeans_and_sds$isotypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("LC Clones by Isotype") +
  xlab("") + ylab("Clones") + geom_text(aes(label = scales::percent(pern2)), position = position_stack(vjust = 0.5))
piechartlcclones2 <- piechartlcclones + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartlcclones)
ggsave("piechart_lc_reads.png", piechartlccounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_lc_reads.pdf", piechartlccounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_lc_clones.png", piechartlcclones2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_lc_clones.pdf", piechartlcclones2, width = 8, height = 6, units = "in", bg = "transparent")


## SUBTYPES...
piechartsubcounts <- ggplot(subtypemeans_and_sds, aes(x = "", y = n1, fill = GANDA_SUBTYPE, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_brewer(name = "Subtype", palette = "BuPu", labels = subtypemeans_and_sds$subtypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Reads by Subtype") +
  xlab("") + ylab("Counts") + geom_text(aes(label = scales::percent(pern1)), position = position_stack(vjust = 0.5))
piechartsubcounts2 <- piechartsubcounts + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartsubcounts2)
piechartsubclones <- ggplot(subtypemeans_and_sds, aes(x = "", y = n2, fill = GANDA_SUBTYPE, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_brewer(name = "Subtype", palette = "BuPu", labels = subtypemeans_and_sds$subtypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Clones by Subtype") +
  xlab("") + ylab("Clones") + geom_text(aes(label = scales::percent(pern2)), position = position_stack(vjust = 0.5))
piechartsubclones2 <- piechartsubclones + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartsubclones2)
ggsave("piechart_subtype_reads.png", piechartsubcounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype_reads.pdf", piechartsubcounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype_clones.png", piechartsubclones2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype_clones.pdf", piechartsubclones2, width = 8, height = 6, units = "in", bg = "transparent")

# new color scheme
piechartsub1counts <- ggplot(subtypemeans_and_sds1, aes(x = "", y = n1, fill = GANDA_SUBTYPE, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Subtype", values = c("orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3"), labels = subtypemeans_and_sds1$subtypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Reads by Subtype") +
  xlab("") + ylab("Counts") + geom_text(aes(label = scales::percent(pern1)), position = position_stack(vjust = 0.5))
piechartsub1counts2 <- piechartsub1counts + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartsub1counts2)
piechartsub1clones <- ggplot(subtypemeans_and_sds1, aes(x = "", y = n2, fill = GANDA_SUBTYPE, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Subtype", values = c("orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3"), labels = subtypemeans_and_sds1$subtypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Clones by Subtype") +
  xlab("") + ylab("Clones") + geom_text(aes(label = scales::percent(pern2)), position = position_stack(vjust = 0.5))
piechartsub1clones2 <- piechartsub1clones + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartsub1clones2)
ggsave("piechart_subtype1_reads.png", piechartsub1counts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype1_reads.pdf", piechartsub1counts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype1_clones.png", piechartsub1clones2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype1_clones.pdf", piechartsub1clones2, width = 8, height = 6, units = "in", bg = "transparent")


### sets of piecharts in threes...
pies_byread <- grid.arrange(piechartallcounts2,piechartsubcounts2,piecharthccounts2,piechartlccounts2, layout_matrix = layout3piesc)
pies_byclone <- grid.arrange(piechartallclones2,piechartsubclones2,piecharthcclones2,piechartlcclones2, layout_matrix = layout3piesc)

### ISOTYPES + SUBTYPES
piechartisosubcounts <- ggplot(isosubmeans_and_sds, aes(x = "", y = n1, fill = GANDA_SUBTYPE, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype/Subtype", values = c("salmon", "orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3", "skyblue", "lightgreen"), labels = isosubmeans_and_sds$subtypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("Reads by Isotype & Subtype") +
  xlab("") + ylab("Counts") + geom_text(aes(label = scales::percent(pern1)), position = position_stack(vjust = 0.5))
piechartisosubcounts2 <- piechartisosubcounts + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartisosubcounts2)
piechartisosubclones <- ggplot(isosubmeans_and_sds, aes(x = "", y = n2, fill = GANDA_SUBTYPE, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype/Subtype", values = c("salmon", "orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3", "skyblue", "lightgreen"), labels = isosubmeans_and_sds$subtypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("Clones by Isotype & Subtype") +
  xlab("") + ylab("Clones") + geom_text(aes(label = scales::percent(pern2)), position = position_stack(vjust = 0.5))
piechartisosubclones2 <- piechartisosubclones + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartisosubclones2)
ggsave("piechart_isosub_reads.png", piechartisosubcounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_isosub_reads.pdf", piechartisosubcounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_isosub_clones.png", piechartisosubclones2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_isosub_clones.pdf", piechartisosubclones2, width = 8, height = 6, units = "in", bg = "transparent")

## HC ISOTYPES + SUBTYPES
piecharthcisosubcounts <- ggplot(hcisosubmeans_and_sds, aes(x = "", y = n1, fill = GANDA_SUBTYPE, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype/Subtype", values = c("salmon", "orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3"), labels = hcisosubmeans_and_sds$subtypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Reads by Isotype & Subtype") +
  xlab("") + ylab("Counts") + geom_text(aes(label = scales::percent(pern1)), position = position_stack(vjust = 0.5))
piecharthcisosubcounts2 <- piecharthcisosubcounts + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piecharthcisosubcounts2)
piecharthcisosubclones <- ggplot(hcisosubmeans_and_sds, aes(x = "", y = n2, fill = GANDA_SUBTYPE, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype/Subtype", values = c("salmon", "orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3"), labels = hcisosubmeans_and_sds$subtypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Clones by Isotype & Subtype") +
  xlab("") + ylab("Clones") + geom_text(aes(label = scales::percent(pern2)), position = position_stack(vjust = 0.5))
piecharthcisosubclones2 <- piecharthcisosubclones + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piecharthcisosubclones2)
ggsave("piechart_hcisosub_reads.png", piecharthcisosubcounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hcisosub_reads.pdf", piecharthcisosubcounts2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hcisosub_clones.png", piecharthcisosubclones2, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hcisosub_clones.pdf", piecharthcisosubclones2, width = 8, height = 6, units = "in", bg = "transparent")

###############################
### NOTE IF NO IgG4, need another set of subytpe plots - note instead at first calculation I now replace all NAs with 0s


####
## multi-page pdfs...
#pdfset1 <- grid.arrange(gfsplots1, gsplots1, mutplots1, cdr3plots1, mutvsCDR3r, mutplotshist, mutplotshistcounts)
#pdfset2 <- grid.arrange(gfplots1, gplots1, mutplots1c, cdr3plots1c, mutvsCDR3, mutplotshistbyclone, mutplotshistcountsbyclone)
#pdfpies1 <- grid.arrange(piechartallcounts2, piecharthccounts2, piechartlccounts2)
#pdfpies2 <- grid.arrange(piechartallclones2, piecharthcclones2, piechartlcclones2)
#manypdfs1 <- marrangeGrob(pdfset1, nrow=1, ncol=1)
#manypdfs2 <- marrangeGrob(pdfset2, nrow=1, ncol=1)

## old isotype color scheme
#pdfset1 <- c(list(piechartallcounts2, piechartsubcounts2, piecharthccounts2, piechartlccounts2, piechartisosubcounts2, gfsplots1, gsplots1, mutplots1, cdr3plots1, mutvsCDR3r, mutplotshist, mutplotshistcounts))
#manypdfs1 <- marrangeGrob(pdfset1, nrow=1, ncol=1)
#ggsave("results_byreads.pdf", manypdfs1, width = 16, height = 12, units = "in")
#pdfset2 <- c(list(piechartallclones2, piechartsubclones2, piecharthcclones2, piechartlcclones2, piechartisosubclones2, gfplots1, gplots1, mutplots1c, cdr3plots1c, mutvsCDR3, mutplotshistbyclone, mutplotshistcountsbyclone))
#manypdfs2 <- marrangeGrob(pdfset2, nrow=1, ncol=1)
#ggsave("results_byclones.pdf", manypdfs2, width = 16, height = 12, units = "in")

## with new isotype color scheme
pdfset1 <- c(list(piechartallcounts2, piechartsub1counts2, piecharthccounts2, piechartlccounts2, piechartisosubcounts2, piecharthcisosubcounts2, gfsplots1, gsplots1, mutplots1h35, cdr3plots1, mutvsCDR3rh35, mutplotshisth35, mutplotshistcountsh35))
manypdfs1 <- marrangeGrob(pdfset1, nrow=1, ncol=1)
ggsave("results_byreads.pdf", manypdfs1, width = 16, height = 12, units = "in")
pdfset2 <- c(list(piechartallclones2, piechartsub1clones2, piecharthcclones2, piechartlcclones2, piechartisosubclones2, piecharthcisosubclones2, gfplots1, gplots1, mutplots1ch35, cdr3plots1c, mutvsCDR3h35, mutplotshistbycloneh35, mutplotshistcountsbycloneh35, nreadshistogrambyclone, nreadshistogramcountsbyclone, gmutandnhexhbyisotypeh35))
manypdfs2 <- marrangeGrob(pdfset2, nrow=1, ncol=1)
ggsave("results_byclones.pdf", manypdfs2, width = 16, height = 12, units = "in")

###############################
## PIECHARTS WITHOUT LABELS (will be in second set of pdfs)
piechartallcountsblank <- ggplot(allmeans_and_sds, aes(x = "", y = n1, fill = PRCONS2, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("salmon", "purple","blue", "skyblue", "lightgreen"), labels = allmeans_and_sds$isotypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("Reads by Isotype")
piechartallcounts2b <- piechartallcountsblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartallcounts2b)
piechartallclonesblank <- ggplot(allmeans_and_sds, aes(x = "", y = n2, fill = PRCONS2, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("salmon", "purple","blue", "skyblue", "lightgreen"), labels = allmeans_and_sds$isotypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("Clones by Isotype")
piechartallclones2b <- piechartallclonesblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartallclones2b)
ggsave("piechart_all_reads_blank.png", piechartallcounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_all_reads_blank.pdf", piechartallcounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_all_clones_blank.png", piechartallclones2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_all_clones_blank.pdf", piechartallclones2b, width = 8, height = 6, units = "in", bg = "transparent")

## hc only
piecharthccountsblank <- ggplot(hcmeans_and_sds, aes(x = "", y = n1, fill = PRCONS2, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("salmon", "purple","blue", "skyblue", "lightgreen"), labels = hcmeans_and_sds$isotypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Reads by Isotype")
piecharthccounts2b <- piecharthccountsblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piecharthccounts2b)
piecharthcclonesblank <- ggplot(hcmeans_and_sds, aes(x = "", y = n2, fill = PRCONS2, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("salmon", "purple","blue", "skyblue", "lightgreen"), labels = allmeans_and_sds$isotypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Clones by Isotype")
piecharthcclones2b <- piecharthcclonesblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piecharthcclones2b)
ggsave("piechart_hc_reads_blank.png", piecharthccounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hc_reads_blank.pdf", piecharthccounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hc_clones_blank.png", piecharthcclones2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hc_clones_blank.pdf", piecharthcclones2b, width = 8, height = 6, units = "in", bg = "transparent")

## lc only
piechartlccountsblank <- ggplot(lcmeans_and_sds, aes(x = "", y = n1, fill = PRCONS2, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("skyblue", "lightgreen"), labels = lcmeans_and_sds$isotypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("LC Reads by Isotype")
piechartlccounts2b <- piechartlccountsblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartlccounts2b)
piechartlcclonesblank <- ggplot(lcmeans_and_sds, aes(x = "", y = n2, fill = PRCONS2, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype", values = c("skyblue", "lightgreen"), labels = lcmeans_and_sds$isotypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("LC Clones by Isotype")
piechartlcclones2b <- piechartlcclonesblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartlcclones2b)
ggsave("piechart_lc_reads_blank.png", piechartlccounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_lc_reads_blank.pdf", piechartlccounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_lc_clones_blank.png", piechartlcclones2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_lc_clones_blank.pdf", piechartlcclones2b, width = 8, height = 6, units = "in", bg = "transparent")

## SUBTYPES...
piechartsubcountsblank <- ggplot(subtypemeans_and_sds, aes(x = "", y = n1, fill = GANDA_SUBTYPE, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_brewer(name = "Subtype", palette = "BuPu", labels = subtypemeans_and_sds$subtypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Reads by Subtype")
piechartsubcounts2b <- piechartsubcountsblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartsubcounts2b)
piechartsubclonesblank <- ggplot(subtypemeans_and_sds, aes(x = "", y = n2, fill = GANDA_SUBTYPE, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_brewer(name = "Subtype", palette = "BuPu", labels = subtypemeans_and_sds$subtypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Clones by Subtype")
piechartsubclones2b <- piechartsubclonesblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartsubclones2b)
ggsave("piechart_subtype_reads_blank.png", piechartsubcounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype_reads_blank.pdf", piechartsubcounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype_clones_blank.png", piechartsubclones2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype_clones_blank.pdf", piechartsubclones2b, width = 8, height = 6, units = "in", bg = "transparent")

# new color scheme
piechartsub1countsblank <- ggplot(subtypemeans_and_sds1, aes(x = "", y = n1, fill = GANDA_SUBTYPE, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Subtype", values = c("orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3"), labels = subtypemeans_and_sds1$subtypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Reads by Subtype")
piechartsub1counts2b <- piechartsub1countsblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartsub1counts2b)
piechartsub1clonesblank <- ggplot(subtypemeans_and_sds1, aes(x = "", y = n2, fill = GANDA_SUBTYPE, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Subtype", values = c("orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3"), labels = subtypemeans_and_sds1$subtypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Clones by Subtype")
piechartsub1clones2b <- piechartsub1clonesblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartsub1clones2b)
ggsave("piechart_subtype1_reads_blank.png", piechartsub1counts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype1_reads_blank.pdf", piechartsub1counts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype1_clones_blank.png", piechartsub1clones2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_subtype1_clones_blank.pdf", piechartsub1clones2b, width = 8, height = 6, units = "in", bg = "transparent")

### ISOTYPES + SUBTYPES
piechartisosubcountsblank <- ggplot(isosubmeans_and_sds, aes(x = "", y = n1, fill = GANDA_SUBTYPE, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype/Subtype", values = c("salmon", "orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3", "skyblue", "lightgreen"), labels = isosubmeans_and_sds$subtypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("Reads by Isotype & Subtype")
piechartisosubcounts2b <- piechartisosubcountsblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartisosubcounts2b)
piechartisosubclonesblank <- ggplot(isosubmeans_and_sds, aes(x = "", y = n2, fill = GANDA_SUBTYPE, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype/Subtype", values = c("salmon", "orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3", "skyblue", "lightgreen"), labels = isosubmeans_and_sds$subtypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("Clones by Isotype & Subtype")
piechartisosubclones2b <- piechartisosubclonesblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piechartisosubclones2b)
ggsave("piechart_isosub_reads_blank.png", piechartisosubcounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_isosub_reads_blank.pdf", piechartisosubcounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_isosub_clones_blank.png", piechartisosubclones2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_isosub_clones_blank.pdf", piechartisosubclones2b, width = 8, height = 6, units = "in", bg = "transparent")

## HC ISOTYPES + SUBTYPES
piecharthcisosubcountsblank <- ggplot(hcisosubmeans_and_sds, aes(x = "", y = n1, fill = GANDA_SUBTYPE, label = scales::percent(pern1))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype/Subtype", values = c("salmon", "orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3"), labels = hcisosubmeans_and_sds$subtypeandn1) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Reads by Isotype & Subtype")
piecharthcisosubcounts2b <- piecharthcisosubcountsblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piecharthcisosubcounts2b)
piecharthcisosubclonesblank <- ggplot(hcisosubmeans_and_sds, aes(x = "", y = n2, fill = GANDA_SUBTYPE, label = scales::percent(pern2))) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Isotype/Subtype", values = c("salmon", "orchid1", "mediumorchid2", "darkorchid3", "darkorchid4", "blue", "blue3"), labels = hcisosubmeans_and_sds$subtypeandn2) +
  coord_polar("y", start = 0, direction = -1) + theme_void() + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ggtitle("HC Clones by Isotype & Subtype")
piecharthcisosubclones2b <- piecharthcisosubclonesblank + theme(plot.title = element_text(face = "bold", size = rel(1.5)), legend.title = element_text(size = rel(1.5)), legend.text = element_text(lineheight = 0.5, size = rel(1.5)))
#plot(piecharthcisosubclones2b)
ggsave("piechart_hcisosub_reads_blank.png", piecharthcisosubcounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hcisosub_reads_blank.pdf", piecharthcisosubcounts2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hcisosub_clones_blank.png", piecharthcisosubclones2b, width = 8, height = 6, units = "in", bg = "transparent")
ggsave("piechart_hcisosub_clones_blank.pdf", piecharthcisosubclones2b, width = 8, height = 6, units = "in", bg = "transparent")

### PDFS
pdfset1b <- c(list(piechartallcounts2b, piechartsub1counts2b, piecharthccounts2b, piechartlccounts2b, piechartisosubcounts2b, piecharthcisosubcounts2b, gfsplots1, gsplots1, mut2plots1h35, cdr3plots1, mut2vsCDR3rh35, mut2plotshisth35, mut2plotshistcountsh35))
manypdfs1b <- marrangeGrob(pdfset1b, nrow=1, ncol=1)
ggsave("results_byreads_blankpies.pdf", manypdfs1b, width = 16, height = 12, units = "in")
pdfset2b <- c(list(piechartallclones2b, piechartsub1clones2b, piecharthcclones2b, piechartlcclones2b, piechartisosubclones2b, piecharthcisosubclones2b, gfplots1, gplots1, mut2plots1ch35, cdr3plots1c, mut2vsCDR3h35, mut2plotshistbycloneh35, mut2plotshistcountsbycloneh35, nreadshistogrambyclonefree, nreadshistogramcountsbyclone, gmut2andnhexhbyisotypeh35))
manypdfs2b <- marrangeGrob(pdfset2b, nrow=1, ncol=1)
ggsave("results_byclones_blankpies.pdf", manypdfs2b, width = 16, height = 12, units = "in")

#xlab("Median Household Income (thousands)") +
#  labs(title = "Median Household Income by Income Tier Across U.S. Metropolitan Areas",
#       subtitle = "Average median income across 229 metros decreased from $67,863 in 1999 to $62,662 in 2014, representing an 8% loss in \nincome. The lower income class experienced the largest impact with a 11% decrease while the middle and upper class median \nhousehold income decreased by 6% and 8% respectively.",
#       caption = "Source: Pew Research Center analysis of the \n2000 decennial census and 2014 American \nCommunity Survey (IPUMS)")

##########################################################################################################
##########################################################################################################
#############################       LINEAGE RECONSTRUCTION         #######################################
##########################################################################################################
##########################################################################################################
