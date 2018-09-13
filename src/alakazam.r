## ALL SHAZAM & ALAKAZAM COMMANDS

# Import required packages
install.packages("gdata")
install.packages("tidyverse")
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

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
BX$SEQUENCE_IMGT <- as.character(BX$SEQUENCE_IMGT)
BX$GERMLINE_IMGT_D_MASK <- as.character(BX$GERMLINE_IMGT_D_MASK)
BX$JUNCTION_LENGTH2 <- round(BX$JUNCTION_LENGTH/3)*3
BX$CDR3KABAT_LENGTH <- (BX$JUNCTION_LENGTH2) / 3

BX$FAMILY <- getFamily(BX$V_CALL, first=TRUE, strip_d=TRUE)
BX$GENE <- getGene(BX$V_CALL, first=TRUE, strip_d=TRUE)
BX$FAMILY <- as.character(BX$FAMILY)
BX$GENE <- as.character(BX$GENE)

BX.H <- subset(BX, PRCONS %in% c("IgM", "IgG", "IgA"))
BX.L <- subset(BX, PRCONS %in% c("Kappa", "Lambda"))
BX.kappa <- subset(BX, PRCONS %in% c("Kappa"))
BX.lambda <- subset(BX, PRCONS %in% c("Lambda"))

BX.H <- BX.H[ grep("IGKV", BX.H$V_CALL, invert = TRUE) , ]
BX.H <- BX.H[ grep("IGLV", BX.H$V_CALL, invert = TRUE) , ]
BX.L <- BX.L[ grep("IGHV", BX.L$V_CALL, invert = TRUE) , ]
BX.kappa <- BX.kappa[ grep("IGHV", BX.kappa$V_CALL, invert = TRUE) , ]
BX.kappa <- BX.kappa[ grep("IGLV", BX.kappa$V_CALL, invert = TRUE) , ]
BX.lambda <- BX.lambda[ grep("IGHV", BX.lambda$V_CALL, invert = TRUE) , ]
BX.lambda <- BX.lambda[ grep("IGKV", BX.lambda$V_CALL, invert = TRUE) , ]

BX$PRCONS2 <- factor(BX$PRCONS, levels = c("IgM", "IgG", "IgA", "Kappa", "Lambda"))
BX.H$PRCONS2 <- factor(BX.H$PRCONS, levels = c("IgM", "IgG", "IgA"))
BX.L$PRCONS2 <- factor(BX.L$PRCONS, levels = c("Kappa", "Lambda"))
BX.kappa$PRCONS2 <- factor(BX.kappa$PRCONS, levels = c("Kappa"))
BX.lambda$PRCONS2 <- factor(BX.lambda$PRCONS, levels = c("Lambda"))

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

### creating lists per clone not per read
clonestatsh1 <- BX_hobs %>%
  group_by(CLONE) %>%
  summarize_at(c("PRCONS2","GENE","V_CALL","D_CALL","J_CALL","JUNCTION_LENGTH","PRCONS","GANDA_SUBTYPE","JUNCTION_LENGTH2","CDR3KABAT_LENGTH","FAMILY"), first)
clonestatsh2 <- BX_hobs %>%
  group_by(CLONE) %>%
  summarize_if(is.numeric, mean)
clonestatsh3 <- BX_hobs %>%
  group_by(CLONE) %>%
  summarize_at(c("PRCONS2","GANDA_SUBTYPE"), n_distinct)
clonestatsh3 <- rename(clonestatsh3, PRCONS2_d = PRCONS2)
clonestatsh3 <- rename(clonestatsh3, GANDA_SUBTYPE_d = GANDA_SUBTYPE)
clonestatsh4 <- inner_join(clonestatsh1, clonestatsh2)
clonestatsh <- inner_join(clonestatsh4, clonestatsh3)

clonestatsk1 <- BX_kobs %>%
  group_by(CLONE) %>%
  summarize_at(c("PRCONS2","GENE","V_CALL","J_CALL","JUNCTION_LENGTH","PRCONS","JUNCTION_LENGTH2","CDR3KABAT_LENGTH","FAMILY"), first)
clonestatsk2 <- BX_kobs %>%
  group_by(CLONE) %>%
  summarize_if(is.numeric, mean)
clonestatsk <- inner_join(clonestatsk1, clonestatsk2)

clonestatsl1 <- BX_lobs %>%
  group_by(CLONE) %>%
  summarize_at(c("PRCONS2","GENE","V_CALL","J_CALL","JUNCTION_LENGTH","PRCONS","JUNCTION_LENGTH2","CDR3KABAT_LENGTH","FAMILY"), first)
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
hmeans_and_sds <- cbind(hmeans1,hsds1,hmeans2,hsds2,hmeans3,hsds3)

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
kmeans_and_sds <- cbind(kmeans1,ksds1,kmeans2,ksds2,kmeans3,ksds3)

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
lmeans_and_sds <- cbind(lmeans1,lsds1,lmeans2,lsds2,lmeans3,lsds3)
allmeans_and_sds <- rbind.fill(hmeans_and_sds,kmeans_and_sds,lmeans_and_sds)

### do same for subtype
submeans1 <- BX_hobs %>%
  group_by(GANDA_SUBTYPE) %>%
  summarize(mutation_mean_reads = mean(MU_FREQ), n1 = n())
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
submeans_and_sds <- cbind(submeans1,subsds1,submeans2,subsds2,submeans3,subsds3)

###################################
### stats for pie charts
lcmeans_and_sds0 <- rbind.fill(kmeans_and_sds,lmeans_and_sds)
hcmeans_and_sdswithd <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM", "IgD", "IgG", "IgA"))
hcmeans_and_sds <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM", "IgG", "IgA"))


### MAKING NEW TABLE COMBINING ISOTYPES AND SUBTYPES - step 1 subset before mutating all
allmeans_and_sds0 <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM", "IgG", "IgA", "Kappa", "Lambda"))
mumeans_and_sds <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM"))

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

allmeans_and_sdswithd <- allmeans_and_sds %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(PRCONS2, n1, col="isotypeandn1", sep=", ", remove = FALSE) %>%
  unite(PRCONS2, n2, col="isotypeandn2", sep=", ", remove = FALSE)
hcmeans_and_sdswithd <- hcmeans_and_sdswithd %>%
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
isotypeandsubtypemeans_and_sds <- rbind.fill(allmeans_and_sds,subtypemeans_and_sds)

### MAKING NEW TABLE COMBINING ISOTYPES AND SUBTYPES
mumeans_and_sds <- subset(allmeans_and_sds0, PRCONS2 %in% c("IgM"))
mulcmeans_and_sds <- rbind.fill(mumeans_and_sds,lcmeans_and_sds0)
#mulcmeans_and_sds <- rename(mulcmeans_and_sds, subtypeandn1 = isotypeandn1)
#mulcmeans_and_sds <- rename(mulcmeans_and_sds, subtypeandn2 = isotypeandn2)
mulcmeans_and_sds <- rename(mulcmeans_and_sds, GANDA_SUBTYPE = PRCONS2)
mumeans_and_sds <- rename(mumeans_and_sds, GANDA_SUBTYPE = PRCONS2)
#rm(isosubmeans_and_sds)

### NOTE IF NO IgG4, need another set of subytpe plots - alternately try replacing all NAs with 0s
isosubmeans_and_sds0 <- rbind.fill(subtypemeans_and_sds0,mulcmeans_and_sds)
isosubmeans_and_sds0$GANDA_SUBTYPE <- factor(isosubmeans_and_sds0$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "Kappa", "Lambda"))
isosubmeans_and_sds1 <- isosubmeans_and_sds0
### this re-orders the rows...
target <- tibble(GANDA_SUBTYPE = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "Kappa", "Lambda"))
isosubmeans_and_sds <- left_join(data.frame(GANDA_SUBTYPE=target),isosubmeans_and_sds0,by="GANDA_SUBTYPE")
isosubmeans_and_sds <- replace_na(isosubmeans_and_sds, list(n1 = 0))
isosubmeans_and_sds <- replace_na(isosubmeans_and_sds, list(n2 = 0))
isosubmeans_and_sds$subtypeandn1 <- gsub("NA", "0", isosubmeans_and_sds$subtypeandn1)
isosubmeans_and_sds$subtypeandn2 <- gsub("NA", "0", isosubmeans_and_sds$subtypeandn2)
isosubmeans_and_sds <- isosubmeans_and_sds %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(GANDA_SUBTYPE, n1, col="subtypeandn1", sep=", ", remove = FALSE) %>%
  unite(GANDA_SUBTYPE, n2, col="subtypeandn2", sep=", ", remove = FALSE)
isosubmeans_and_sds$GANDA_SUBTYPE <- factor(isosubmeans_and_sds$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "Kappa", "Lambda"))
### same color scheme for subtypes
target2 <- tibble(GANDA_SUBTYPE = c("IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
subtypemeans_and_sds1 <- left_join(data.frame(GANDA_SUBTYPE=target2),subtypemeans_and_sds0,by="GANDA_SUBTYPE")
subtypemeans_and_sds1 <- replace_na(subtypemeans_and_sds1, list(n1 = 0))
subtypemeans_and_sds1 <- replace_na(subtypemeans_and_sds1, list(n2 = 0))
subtypemeans_and_sds1$subtypeandn1 <- gsub("NA", "0", subtypemeans_and_sds1$subtypeandn1)
subtypemeans_and_sds1$subtypeandn2 <- gsub("NA", "0", subtypemeans_and_sds1$subtypeandn2)
subtypemeans_and_sds1 <- subtypemeans_and_sds1 %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(GANDA_SUBTYPE, n1, col="subtypeandn1", sep=", ", remove = FALSE) %>%
  unite(GANDA_SUBTYPE, n2, col="subtypeandn2", sep=", ", remove = FALSE)
subtypemeans_and_sds1$GANDA_SUBTYPE <- factor(subtypemeans_and_sds1$GANDA_SUBTYPE, levels = c("IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))

## now HC isotypes and subtypes
hcisosubmeans_and_sds0 <- rbind.fill(subtypemeans_and_sds0,mumeans_and_sds)
hcisosubmeans_and_sds0$GANDA_SUBTYPE <- factor(hcisosubmeans_and_sds0$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
hcisosubmeans_and_sds1 <- hcisosubmeans_and_sds0
### this re-orders the rows...
htarget <- tibble(GANDA_SUBTYPE = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
hcisosubmeans_and_sds <- left_join(data.frame(GANDA_SUBTYPE=htarget),hcisosubmeans_and_sds0,by="GANDA_SUBTYPE")
hcisosubmeans_and_sds <- replace_na(hcisosubmeans_and_sds, list(n1 = 0))
hcisosubmeans_and_sds <- replace_na(hcisosubmeans_and_sds, list(n2 = 0))
hcisosubmeans_and_sds$subtypeandn1 <- gsub("NA", "0", hcisosubmeans_and_sds$subtypeandn1)
hcisosubmeans_and_sds$subtypeandn2 <- gsub("NA", "0", hcisosubmeans_and_sds$subtypeandn2)
hcisosubmeans_and_sds <- hcisosubmeans_and_sds %>%
  mutate(pern1 = n1/sum(n1)) %>%
  mutate(pern2 = n2/sum(n2)) %>%
  unite(GANDA_SUBTYPE, n1, col="subtypeandn1", sep=", ", remove = FALSE) %>%
  unite(GANDA_SUBTYPE, n2, col="subtypeandn2", sep=", ", remove = FALSE)
hcisosubmeans_and_sds$GANDA_SUBTYPE <- factor(hcisosubmeans_and_sds$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))

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
write.table(isotypeandsubtypemeans_and_sds, "summary_mutationstats_all.tsv", sep = "\t", row.names = FALSE)

#clonestatsh$N_WEIGHT <- rescale(clonestatsh$n)
#clonestatsk$N_WEIGHT <- rescale(clonestatsk$n)
#clonestatsl$N_WEIGHT <- rescale(clonestatsl$n)
#clonestatshf$N_WEIGHT <- rescale(clonestatshf$n)
#clonestatskf$N_WEIGHT <- rescale(clonestatskf$n)
#clonestatslf$N_WEIGHT <- rescale(clonestatslf$n)
#BX.H$PRCONS2 <- factor(BX.H$PRCONS2, levels = c("IgM", "IgG", "IgA"))
#BX_hobs$PRCONS2 <- factor(BX_hobs$PRCONS2, levels = c("IgM", "IgG", "IgA"))
#clonestatsh$PRCONS2 <- factor(clonestatsh$PRCONS2, levels = c("IgM", "IgG", "IgA"))
#clonestatshf$PRCONS2 <- factor(clonestatshf$PRCONS2, levels = c("IgM", "IgG", "IgA"))


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
ghcdr3v <- ggplot(BX_hobs, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
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
gkcdr3v <- ggplot(BX_kobs, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
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
glcdr3v <- ggplot(BX_lobs, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
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
ghcdr3vc <- ggplot(clonestatsh, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
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
gkcdr3vc <- ggplot(clonestatsk, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
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
glcdr3vc <- ggplot(clonestatsl, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
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

ghmutvcf <- ggplot(clonestatshf, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gmutv) + scale_y_continuous(labels = scales::percent)
ghcdr3vcf <- ggplot(clonestatshf, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() + ggtitle("CDR3 Length by gene") +
  xlab("Gene") + ylab("CDRH3 Length, Kabat (aa)") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gcdr3v) + scale_y_continuous(labels = scales::percent)
gkmutvcf <- ggplot(clonestatskf, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gkmutv)
gkcdr3vcf <- ggplot(clonestatskf, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() +
  xlab("Gene") + ylab("CDRL3 Length, Kabat (aa)") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gkcdr3v) + scale_y_continuous(labels = scales::percent)
glmutvcf <- ggplot(clonestatslf, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(glmutv)
glcdr3vcf <- ggplot(clonestatslf, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() +
  xlab("Gene") + ylab("CDRL3 Length, Kabat (aa)") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(glcdr3v) + scale_y_continuous(labels = scales::percent)

mutplots1cf <- grid.arrange(ghmutvcf,gkmutvcf,glmutvcf, layout_matrix = layouthkl3)
ggsave("mutation_bygene_byclone_filtered.png", mutplots1cf, width = 16, height = 12, units = "in")
ggsave("mutation_bygene_byclone_filtered.pdf", mutplots1cf, width = 16, height = 12, units = "in")
cdr3plots1cf <- grid.arrange(ghcdr3vcf,gkcdr3vcf,glcdr3vcf, layout_matrix = layouthkl3)
ggsave("CDR3_bygene_byclone_filtered.png", cdr3plots1cf, width = 16, height = 12, units = "in")
ggsave("CDR3_bygene_byclone_filtered.pdf", cdr3plots1cf, width = 16, height = 12, units = "in")

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

###################
## now clone plots but filter so leaving out all single clones!
ghmuthistogrambyclonef <- ggplot(clonestatshf, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of clones with ≥ 2 Members")
#plot(ghmuthistogrambyclonef)
gkmuthistogrambyclonef <- ggplot(clonestatskf, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Proportion of clones with ≥ 2 Members")
#plot(gkmuthistogrambyclonef)
glmuthistogrambyclonef <- ggplot(clonestatslf, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Proportion of clones with ≥ 2 Members")
#plot(glmuthistogrambyclonef)

ghmuthistogramcountsbyclonef <- ggplot(clonestatshf, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Clones with ≥ 2 Members")
#plot(ghmuthistogramcountsbyclonef)
gkmuthistogramcountsbyclonef <- ggplot(clonestatskf, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Clones with ≥ 2 Members")
#plot(gkmuthistogramcountsbyclonef)
glmuthistogramcountsbyclonef <- ggplot(clonestatslf, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + xlab("% Somatic Hypermutation") + ylab("Clones with ≥ 2 Members")
#plot(glmuthistogramcountsbyclonef)

mutplotshistbyclonef <- grid.arrange(ghmuthistogrambyclonef,gkmuthistogrambyclonef,glmuthistogrambyclonef, layout_matrix = layouthkl3)
ggsave("mutation_histogram_byclone_filtered_hlseparatepanels.png", mutplotshistbyclonef, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_filtered_hlseparatepanels.pdf", mutplotshistbyclonef, width = 16, height = 12, units = "in")
mutplotshistcountsbyclonef <- grid.arrange(ghmuthistogramcountsbyclonef,gkmuthistogramcountsbyclonef,glmuthistogramcountsbyclonef, layout_matrix = layouthkl3)
ggsave("mutation_histogram_rawcounts_byclone_filtered_hlseparatepanels.png", mutplotshistcountsbyclonef, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_rawcounts_byclone_filtered_hlseparatepanels.pdf", mutplotshistcountsbyclonef, width = 16, height = 12, units = "in")

gmutviolinhf <- ggplot(clonestatshf, aes(x=PRCONS2, y=MU_FREQ, stroke = 0.001)) +
  theme_bw() + ggtitle("% Somatic Hypermutation by Clone") +
  xlab("Isotype") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + scale_y_continuous(labels = scales::percent) +
  geom_violin(fill = "gray") + theme(axis.text.x = element_text(angle=45, hjust=1))
#plot(gmutviolinhf)
ggsave("mutation_violinplot_byclone_filtered_h.png", gmutviolinhf, width = 16, height = 8, units = "in")


################
### RBIND ALL HC AND LC IN SINGLE PANELS - DO NOT LOOK AS GOOD THOUGH
BX_allobs <- rbind(BX_hobs, BX_kobs, BX_lobs)
BX_allobs$PRCONS3 <- factor(BX_allobs$PRCONS2, levels = c("IgM", "Kappa", "IgG", "Lambda", "IgA"))
gallmuthistogram <- ggplot(BX_allobs, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS3, ncol=2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of reads")
#plot(gallmuthistogram)
gallmuthistogramcounts <- ggplot(BX_allobs, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS3, ncol=2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(gallmuthistogramcounts)
BX_allobs4 <- BX_allobs[ grep("IgD", BX_allobs$PRCONS2, invert = TRUE) , ]
BX_allobs4 <- BX_allobs[ grep("IgA", BX_allobs$PRCONS2, invert = TRUE) , ]
BX_allobs4$PRCONS4 <- factor(BX_allobs4$PRCONS2, levels = c("IgM", "Kappa", "IgG", "Lambda"))
ggmklmuthistogram <- ggplot(BX_allobs4, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS4, ncol=2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of reads")
#plot(ggmklmuthistogram)
ggmklmuthistogramcounts <- ggplot(BX_allobs4, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS4, ncol=2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(ggmklmuthistogramcounts)
ggsave("mutation_histogram_onepanel.png", gallmuthistogram, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_onepanel.pdf", gallmuthistogram, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_rawcounts_onepanel.png", gallmuthistogramcounts, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_rawcounts_onepanel.pdf", gallmuthistogramcounts, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_gmkl.png", ggmklmuthistogram, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_gmkl.pdf", ggmklmuthistogram, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_rawcounts_gmkl.png", ggmklmuthistogramcounts, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_rawcounts_gmkl.pdf", ggmklmuthistogramcounts, width = 16, height = 12, units = "in")
###
clonestatsall <- rbind.fill(clonestatsh, clonestatsk, clonestatsl)
clonestatsall$PRCONS3 <- factor(clonestatsall$PRCONS2, levels = c("IgM", "Kappa", "IgG", "Lambda", "IgA"))
gallmuthistogrambyclone <- ggplot(clonestatsall, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS3, ncol=2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(gallmuthistogrambyclone)
gallmuthistogramcountsbyclone <- ggplot(clonestatsall, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS3, ncol=2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(gallmuthistogramcountsbyclone)
clonestatsall4 <- clonestatsall[ grep("IgD", clonestatsall$PRCONS2, invert = TRUE) , ]
clonestatsall4 <- clonestatsall[ grep("IgA", clonestatsall$PRCONS2, invert = TRUE) , ]
clonestatsall4$PRCONS4 <- factor(clonestatsall4$PRCONS2, levels = c("IgM", "Kappa", "IgG", "Lambda"))
ggmklmuthistogrambyclone <- ggplot(clonestatsall4, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS4, ncol=2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ggmklmuthistogrambyclone)
ggmklmuthistogramcountsbyclone <- ggplot(clonestatsall4, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS4, ncol=2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(ggmklmuthistogramcountsbyclone)
ggsave("mutation_histogram_byclone_onepanel.png", gallmuthistogrambyclone, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_onepanel.pdf", gallmuthistogrambyclone, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_rawcounts_onepanel.png", gallmuthistogramcountsbyclone, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_rawcounts_onepanel.pdf", gallmuthistogramcountsbyclone, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_gmkl.png", ggmklmuthistogrambyclone, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_gmkl.pdf", ggmklmuthistogrambyclone, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_rawcounts_gmkl.png", ggmklmuthistogramcountsbyclone, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_rawcounts_gmkl.pdf", ggmklmuthistogramcountsbyclone, width = 16, height = 12, units = "in")

## HC plots with only IgM & IgG - THIS MAKES 4 SEPARATE PLOTS EACH WITH OWN RANGE AND LABEL
clonestatsgm <- clonestatsh
clonestatsg <- clonestatsgm[ grep("IgG", clonestatsgm$PRCONS2, invert = FALSE) , ]
clonestatsm <- clonestatsgm[ grep("IgM", clonestatsgm$PRCONS2, invert = FALSE) , ]
ggmuthistogrambyclone <- ggplot(clonestatsg, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ggmuthistogrambyclone)
gmmuthistogrambyclone <- ggplot(clonestatsm, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ggmuthistogrambyclone)

mutplotshistbyclonegmkl <- grid.arrange(gmmuthistogrambyclone,ggmuthistogrambyclone,gkmuthistogrambyclone,glmuthistogrambyclone, layout_matrix = layoutgmkl)
ggsave("mutation_histogram_byclone_gmkl_4panels.png", mutplotshistbyclonegmkl, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_gmkl_4panels.pdf", mutplotshistbyclonegmkl, width = 16, height = 12, units = "in")
#ggmmuthistogrambyclone <- ggplot(clonestatsgm, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
#  facet_wrap(~ PRCONS2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ggmmuthistogrambyclone)

####
clonestatsallf <- rbind.fill(clonestatshf, clonestatskf, clonestatslf)
clonestatsallf$PRCONS3 <- factor(clonestatsallf$PRCONS2, levels = c("IgM", "Kappa", "IgG", "Lambda", "IgA"))
gallmuthistogrambyclonef <- ggplot(clonestatsallf, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS3, ncol=2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of clones with ≥ 2 Members")
#plot(gallfmuthistogrambyclonef)
gallmuthistogramcountsbyclonef <- ggplot(clonestatsallf, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS3, ncol=2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Clones with ≥ 2 Members")
#plot(gallfmuthistogramcountsbyclonef)
clonestatsallf4 <- clonestatsallf[ grep("IgD", clonestatsallf$PRCONS2, invert = TRUE) , ]
clonestatsallf4 <- clonestatsallf[ grep("IgA", clonestatsallf$PRCONS2, invert = TRUE) , ]
clonestatsallf4$PRCONS4 <- factor(clonestatsallf4$PRCONS2, levels = c("IgM", "Kappa", "IgG", "Lambda"))
ggmklmuthistogrambyclonef <- ggplot(clonestatsallf4, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS4, ncol=2) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of clones with ≥ 2 Members")
#plot(ggmklmuthistogrambyclonef)
ggmklmuthistogramcountsbyclonef <- ggplot(clonestatsallf4, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS4, ncol=2) + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Clones with ≥ 2 Members")
#plot(ggmklmuthistogramcountsbyclonef)
ggsave("mutation_histogram_byclone_filtered_onepanel.png", gallmuthistogrambyclonef, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_filtered_onepanel.pdf", gallmuthistogrambyclonef, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_filtered_rawcounts_onepanel.png", gallmuthistogramcountsbyclonef, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_filtered_rawcounts_onepanel.pdf", gallmuthistogramcountsbyclonef, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_filtered_gmkl.png", ggmklmuthistogrambyclonef, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_filtered_gmkl.pdf", ggmklmuthistogrambyclonef, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_filtered_rawcounts_gmkl.png", ggmklmuthistogramcountsbyclonef, width = 16, height = 12, units = "in")
ggsave("mutation_histogram_byclone_filtered_rawcounts_gmkl.pdf", ggmklmuthistogramcountsbyclonef, width = 16, height = 12, units = "in")


###################################
##### mutation vs CDR3 hex plots
## by read
gmutandcdr3hexhr <- ggplot(BX_hobs, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexhr)
gmutandcdr3hexkr <- ggplot(BX_kobs, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads",  breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexkr)
gmutandcdr3hexlr <- ggplot(BX_lobs, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexlr)
mutvsCDR3r <- grid.arrange(gmutandcdr3hexhr,gmutandcdr3hexkr,gmutandcdr3hexlr, layout_matrix = layouthkl3)
ggsave("mutationvsCDR3_byread_hlseparatepanels.png", mutvsCDR3r, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byread_hlseparatepanels.pdf", mutvsCDR3r, width = 16, height = 12, units = "in")

## by clone
gmutandcdr3hexh <- ggplot(clonestatsh, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexh)
gmutandcdr3hexk <- ggplot(clonestatsk, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexk)
gmutandcdr3hexl <- ggplot(clonestatsl, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexl)
mutvsCDR3 <- grid.arrange(gmutandcdr3hexh,gmutandcdr3hexk,gmutandcdr3hexl, layout_matrix = layouthkl3)
ggsave("mutationvsCDR3_byclone_hlseparatepanels.png", mutvsCDR3, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byclone_hlseparatepanels.pdf", mutvsCDR3, width = 16, height = 12, units = "in")

## test violin plot
gmutviolinh <- ggplot(clonestatsh, aes(x=PRCONS2, y=MU_FREQ, stroke = 0.001)) +
  theme_bw() + ggtitle("% Somatic Hypermutation by Clone") +
  xlab("Isotype") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + scale_y_continuous(labels = scales::percent) +
  geom_violin() + theme(axis.text.x = element_text(angle=45, hjust=1))
#plot(gmutviolinh)
ggsave("mutation_violinplot_byclone_h.png", gmutviolinh, width = 16, height = 12, units = "in")

## by filtered clone
gmutandcdr3hexhf <- ggplot(clonestatshf, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexhf)
gmutandcdr3hexkf <- ggplot(clonestatskf, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexkf)
gmutandcdr3hexlf <- ggplot(clonestatslf, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() +
  xlab("CDRL3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexlf)
mutvsCDR3f <- grid.arrange(gmutandcdr3hexhf,gmutandcdr3hexkf,gmutandcdr3hexlf, layout_matrix = layouthkl3)
ggsave("mutationvsCDR3_byclone_filtered_hlseparatepanels.png", mutvsCDR3f, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byclone_filtered_hlseparatepanels.pdf", mutvsCDR3f, width = 16, height = 12, units = "in")

### SINGLE PANELS & GMKL ONLY
clonestatsall <- rbind.fill(clonestatsh, clonestatsk, clonestatsl)
clonestatsall$PRCONS3 <- factor(clonestatsall$PRCONS2, levels = c("IgM", "Kappa", "IgG", "Lambda", "IgA"))
gmutandcdr3hexall <- ggplot(clonestatsall, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDR3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS3, ncol=2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexall)
clonestatsall4 <- clonestatsall[ grep("IgD", clonestatsall$PRCONS2, invert = TRUE) , ]
clonestatsall4 <- clonestatsall[ grep("IgA", clonestatsall$PRCONS2, invert = TRUE) , ]
clonestatsall4$PRCONS4 <- factor(clonestatsall4$PRCONS2, levels = c("IgM", "Kappa", "IgG", "Lambda"))
gmutandcdr3hexgmkl <- ggplot(clonestatsall4, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDR3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS4, ncol=2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexgmkl)
ggsave("mutationvsCDR3_byclone_onepanel.png", gmutandcdr3hexall, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byclone_onepanel.pdf", gmutandcdr3hexall, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byclone_gmkl.png", gmutandcdr3hexgmkl, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byclone_gmkl.pdf", gmutandcdr3hexgmkl, width = 16, height = 12, units = "in")

### 
clonestatsallf <- rbind.fill(clonestatshf, clonestatskf, clonestatslf)
clonestatsallf$PRCONS3 <- factor(clonestatsallf$PRCONS2, levels = c("IgM", "Kappa", "IgG", "Lambda", "IgA"))
gmutandcdr3hexallf <- ggplot(clonestatsallf, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDR3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS3, ncol=2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexallf)
clonestatsallf4 <- clonestatsallf[ grep("IgD", clonestatsallf$PRCONS2, invert = TRUE) , ]
clonestatsallf4 <- clonestatsallf[ grep("IgA", clonestatsallf$PRCONS2, invert = TRUE) , ]
clonestatsallf4$PRCONS4 <- factor(clonestatsallf4$PRCONS2, levels = c("IgM", "Kappa", "IgG", "Lambda"))
gmutandcdr3hexgmkl <- ggplot(clonestatsallf4, aes(x=CDR3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDR3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS4, ncol=2) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones", breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000))
#plot(gmutandcdr3hexgmkl)
ggsave("mutationvsCDR3_byclone_filtered_onepanel.png", gmutandcdr3hexallf, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byclone_filtered_onepanel.pdf", gmutandcdr3hexallf, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byclone_filtered_gmkl.png", gmutandcdr3hexgmkl, width = 16, height = 12, units = "in")
ggsave("mutationvsCDR3_byclone_filtered_gmkl.pdf", gmutandcdr3hexgmkl, width = 16, height = 12, units = "in")

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
pdfset1 <- c(list(piechartallcounts2, piechartsub1counts2, piecharthccounts2, piechartlccounts2, piechartisosubcounts2, piecharthcisosubcounts2, gfsplots1, gsplots1, mutplots1, cdr3plots1, mutvsCDR3r, mutplotshist, mutplotshistcounts))
manypdfs1 <- marrangeGrob(pdfset1, nrow=1, ncol=1)
ggsave("results_byreads.pdf", manypdfs1, width = 16, height = 12, units = "in")
pdfset2 <- c(list(piechartallclones2, piechartsub1clones2, piecharthcclones2, piechartlcclones2, piechartisosubclones2, piecharthcisosubclones2, gfplots1, gplots1, mutplots1c, cdr3plots1c, mutvsCDR3, mutplotshistbyclone, mutplotshistcountsbyclone))
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
pdfset1b <- c(list(piechartallcounts2b, piechartsub1counts2b, piecharthccounts2b, piechartlccounts2b, piechartisosubcounts2b, piecharthcisosubcounts2b, gfsplots1, gsplots1, mutplots1, cdr3plots1, mutvsCDR3r, mutplotshist, mutplotshistcounts))
manypdfs1b <- marrangeGrob(pdfset1b, nrow=1, ncol=1)
ggsave("results_byreads_blankpies.pdf", manypdfs1b, width = 16, height = 12, units = "in")
pdfset2b <- c(list(piechartallclones2b, piechartsub1clones2b, piecharthcclones2b, piechartlcclones2b, piechartisosubclones2b, piecharthcisosubclones2b, gfplots1, gplots1, mutplots1c, cdr3plots1c, mutvsCDR3, mutplotshistbyclone, mutplotshistcountsbyclone))
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


#### FOR PARTICULAR CLONE, USES SAME BX FILE AS INPUT ABOVE
## THE FOLLOWING IS COMMENTED OUT FOR NOW, UNTIL WE KNOW WHAT CLONE NUMBER TO USE

# BX.sub <- subset(BX, CLONE == OPT$CLONE)

# BX.sub$PRCONS <- as.character(BX.sub$PRCONS)
# BX.sub$SEQUENCE_IMGT <- as.character(BX.sub$SEQUENCE_IMGT)
# BX.sub$GERMLINE_IMGT_D_MASK <- as.character(BX.sub$GERMLINE_IMGT_D_MASK)

# BX.sub$JUNCTION_LENGTH2 <- round(BX.sub$JUNCTION_LENGTH/3)*3
# BX.sub$CDR3KABAT_LENGTH <- (BX.sub$JUNCTION_LENGTH2) / 3

# BX.sub.H <- subset(BX.sub, PRCONS %in% c("IgM", "IgD", "IgG", "IgA"))
# BX.sub.L <- subset(BX.sub, PRCONS %in% c("Kappa", "Lambda"))
# BX.sub.kappa <- subset(BX.sub, PRCONS %in% c("Kappa"))
# BX.sub.lambda <- subset(BX.sub, PRCONS %in% c("Lambda"))

# BX.sub.H <- BX.sub.H[ grep("IGKV", BX.sub.H$V_CALL, invert = TRUE) , ]
# BX.sub.H <- BX.sub.H[ grep("IGLV", BX.sub.H$V_CALL, invert = TRUE) , ]
# BX.sub.L <- BX.sub.L[ grep("IGHV", BX.sub.L$V_CALL, invert = TRUE) , ]
# BX.sub.kappa <- BX.sub.kappa[ grep("IGHV", BX.sub.kappa$V_CALL, invert = TRUE) , ]
# BX.sub.lambda <- BX.sub.lambda[ grep("IGHV", BX.sub.lambda$V_CALL, invert = TRUE) , ]



# ## now plotting sizes of clonal members based on barcode count
# ## but using 2 x ln of count, then rounded to nearest five...
# ### if including single sequences

# BX.sub$CONSCOUNT2 <- (2 * log(BX.sub$CONSCOUNT)) + 1
# BX.sub$CONSCOUNT2 <- round(BX.sub$CONSCOUNT2)


# ## ---- eval=TRUE----------------------------------------------------------
# # This example data set does not have ragged ends
# # Preprocess clone without ragged end masking (default)
# clone2 <- makeChangeoClone(BX.sub, max_mask = NULL, pad_end = TRUE, text_fields=c("PRCONS"), 
#                            num_fields="CONSCOUNT2")

# # Show combined annotations
# clone2@data[, c("PRCONS", "CONSCOUNT2")]

# ## ---- eval=FALSE---------------------------------------------------------
# #  # Run PHYLIP and parse output
# #  dnapars_exec <- "~/apps/phylip-3.69/dnapars"
# #  graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE)

# ## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
# # Load data insted of running phylip
# # Clone 3138 is at index 23
# #graph <- ExampleTrees[[23]]
# dnapars_exec <- "~/phylip/exe/dnapars"
# graph <- buildPhylipLineage(clone2, dnapars_exec, rm_temp=TRUE)

# ## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
# # The graph has shared annotations for the clone
# data.frame(CLONE=graph$clone,
#            JUNCTION_LENGTH=graph$junc_len,
#            V_GENE=graph$v_gene,
#            J_GENE=graph$j_gene)

# # The vertices have sequence specific annotations
# data.frame(SEQUENCE_ID=V(graph)$name, 
#            PRCONS=V(graph)$PRCONS,
#            CONSCOUNT2=V(graph)$CONSCOUNT2)

# ## ---- eval=TRUE----------------------------------------------------------
# # Plot graph with defaults
# plot(graph)

# ## ---- eval=TRUE----------------------------------------------------------
# # Modify graph and plot attributes
# V(graph)$color <- "salmon"
# V(graph)$color[V(graph)$name == "Germline"] <- "black"
# V(graph)$color[V(graph)$PRCONS == "IgD"] <- "purple"
# V(graph)$color[V(graph)$PRCONS == "IgD,IgM"] <- "purple"
# V(graph)$color[V(graph)$PRCONS == "IgG"] <- "steelblue"
# V(graph)$color[V(graph)$PRCONS == "IgG1"] <- "cyan"
# V(graph)$color[V(graph)$PRCONS == "IgG1,IgG3"] <- "cyan"
# V(graph)$color[V(graph)$PRCONS == "IgG2"] <- "green"
# V(graph)$color[V(graph)$PRCONS == "IgG2,IgG3"] <- "steelblue"
# V(graph)$color[V(graph)$PRCONS == "IgG3"] <- "steelblue"
# V(graph)$color[V(graph)$PRCONS == "IgG4"] <- "green"
# V(graph)$color[V(graph)$PRCONS == "IgGmil"] <- "darkblue"
# V(graph)$color[V(graph)$PRCONS == "IgM,IgG"] <- "salmon"
# V(graph)$color[V(graph)$PRCONS == "IgM,IgG3"] <- "salmon"
# V(graph)$color[V(graph)$PRCONS == "IgG3,IgM"] <- "salmon"
# V(graph)$color[V(graph)$PRCONS == "IgA"] <- "orange"
# V(graph)$color[V(graph)$PRCONS == "IgA,IgG"] <- "orange"
# V(graph)$color[V(graph)$PRCONS == "IgA,IgM"] <- "orange"
# V(graph)$color[V(graph)$PRCONS == "IgA,IgG3"] <- "orange"
# V(graph)$color[V(graph)$PRCONS == "IgA,IgG4"] <- "orange"
# V(graph)$color[V(graph)$PRCONS == "IgA,IgG1,IgG2,IgG3"] <- "orange"
# #V(graph)$color[V(graph)$PRCONS == "INF17_H"] <- "grey"
# V(graph)$color[V(graph)$PRCONS == "CTB_32A"] <- "grey"
# V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
# V(graph)$color[grepl("INF", V(graph)$name)] <- "grey"
# V(graph)$label <- V(graph)$PRCONS
# E(graph)$label <- ""



# V(graph)$size <- 5
# #V(graph)$size[V(graph)$CONSCOUNT2 == 5] <- 5
# #V(graph)$size[V(graph)$CONSCOUNT2 == 10] <- 10
# #V(graph)$size[V(graph)$CONSCOUNT2 == 15] <- 15
# #V(graph)$size[V(graph)$CONSCOUNT2 == 20] <- 20
# #V(graph)$size[V(graph)$CONSCOUNT2 == 25] <- 25
# #V(graph)$size[V(graph)$CONSCOUNT2 == 30] <- 30
# V(graph)$size[V(graph)$CONSCOUNT2 == 2] <- 6
# V(graph)$size[V(graph)$CONSCOUNT2 == 3] <- 8
# V(graph)$size[V(graph)$CONSCOUNT2 == 4] <- 8
# V(graph)$size[V(graph)$CONSCOUNT2 == 5] <- 8
# V(graph)$size[V(graph)$CONSCOUNT2 == 6] <- 12
# V(graph)$size[V(graph)$CONSCOUNT2 == 7] <- 12
# V(graph)$size[V(graph)$CONSCOUNT2 == 8] <- 12
# V(graph)$size[V(graph)$CONSCOUNT2 == 9] <- 16
# V(graph)$size[V(graph)$CONSCOUNT2 == 10] <- 16
# V(graph)$size[V(graph)$CONSCOUNT2 == 11] <- 16
# V(graph)$size[V(graph)$CONSCOUNT2 == 12] <- 16
# V(graph)$size[V(graph)$CONSCOUNT2 == 13] <- 16
# V(graph)$size[V(graph)$CONSCOUNT2 == 14] <- 16
# V(graph)$size[V(graph)$CONSCOUNT2 == 16] <- 16
# V(graph)$size[V(graph)$CONSCOUNT2 == 17] <- 16
# V(graph)$size[V(graph)$CONSCOUNT2 > 17] <- 20
# #V(graph)$size[V(graph)$CONSCOUNT2 > 40] <- 40

# V(graph)$size[V(graph)$name == "Germline"] <- 4
# V(graph)$size[V(graph)$name == "CTB.32A"] <- 4
# V(graph)$size[grepl("Inferred", V(graph)$name)] <- 2 
# V(graph)$size[grepl("INF", V(graph)$name)] <- 4 

# #V(graph)$size[grepl("TATATGTAGCGTGGCG", V(graph)$name)] <- 2 
# #V(graph)$size[grepl("CCTGGAAGCTTCCCGT", V(graph)$name)] <- 2 


# getPathLengths(graph, root="Germline")
# # Remove large default margins
# par(mar=c(0, 0, 0, 0) + 0.1)
# # Plot graph
# #layout <- layout.reingold.tilford(graph, circular=F)
# plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=1, edge.label.cex=4, vertex.frame.color="grey",
#      vertex.label.color="black", edge.label.color="black")

# V(graph)$label <- ""
# plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=1, edge.label.cex=4, vertex.frame.color="grey",
#      vertex.label.color="white", edge.label.color="white")

# # Add legend
# #legend("topright", c("Germline", "Inferred", "Sample"), 
# #       fill=c("black", "white", "steelblue"), cex=0.75)



# #############################
# ## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
# # Preprocess clones
# clones <- BX %>%
#   group_by(CLONE) %>%
#   do(CHANGEO=makeChangeoClone(., text_fields=c("PRCONS"), 
#                               num_fields="DUPCOUNT"))

# ## ---- eval=FALSE---------------------------------------------------------
# #  # Build lineages
# #  dnapars_exec <- "~/apps/phylip-3.69/dnapars"
# #  graphs <- lapply(clones$CHANGEO, buildPhylipLineage,
# #                   dnapars_exec=dnapars_exec, rm_temp=TRUE)

# ## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
# # Load data insted of running phylip
# graphs <- ExampleTrees

# ## ---- eval=TRUE----------------------------------------------------------
# # Note, clones with only a single sequence will not be processed.
# # A warning will be generated and NULL will be returned by buildPhylipLineage
# # These entries may be removed for clarity
# graphs[sapply(graphs, is.null)] <- NULL

# # The set of tree may then be subset by node count for further 
# # analysis, if desired.
# graphs <- graphs[sapply(graphs, vcount) >= 5]


# #####################################################################
# #####################################################################
# ##################       TOPOLOGY ANALYSIS              #############
# #####################################################################
# #####################################################################

