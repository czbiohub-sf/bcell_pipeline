## ALL SHAZAM & ALAKAZAM COMMANDS

# Import required packages
install.packages("gridExtra")
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(lattice))

#################################################################################################
################################ GENERIC COMMANDS FOR ALL ANALYSES ##############################
#################################################################################################


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

## NOTE - NEED TO SPECIFY INPUT FILE HERE!!! - IN FUTURE HAVE THIS BE AN OUTSIDE CALL??
BX <- read.delim(opt$DB)

#DB4s <- read.delim("/users/eric.waltari/immcantation_pipeline/BX_min5/changeo/d4_1m_clone-pass.tab")
#DB4b <- read.delim("/users/eric.waltari/immcantation_pipeline/BX_min5/changeo/d4_1m_beads_clone-pass.tab")
#DB5b <- read.delim("/users/eric.waltari/immcantation_pipeline/BX_min5/changeo/d5_1m_beads_clone-pass.tab")

## UPDATE DB1
#BX <- DB4s

BX$PRCONS <- as.character(BX$PRCONS)
BX$SEQUENCE_IMGT <- as.character(BX$SEQUENCE_IMGT)
BX$GERMLINE_IMGT_D_MASK <- as.character(BX$GERMLINE_IMGT_D_MASK)
BX$JUNCTION_LENGTH2 <- round(BX$JUNCTION_LENGTH/3)*3
BX$CDR3KABAT_LENGTH <- (BX$JUNCTION_LENGTH2) / 3

BX$FAMILY <- getFamily(BX$V_CALL, first=TRUE, strip_d=TRUE)
BX$GENE <- getGene(BX$V_CALL, first=TRUE, strip_d=TRUE)
BX$FAMILY <- as.character(BX$FAMILY)
BX$GENE <- as.character(BX$GENE)

BX.H <- subset(BX, PRCONS %in% c("IgM", "IgD", "IgG", "IgA"))
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

BX$PRCONS2 <- factor(BX$PRCONS, levels = c("IgM", "IgD", "IgG", "IgA", "Kappa", "Lambda"))
BX.H$PRCONS2 <- factor(BX.H$PRCONS, levels = c("IgM", "IgD", "IgG", "IgA"))
BX.L$PRCONS2 <- factor(BX.L$PRCONS, levels = c("Kappa", "Lambda"))
BX.kappa$PRCONS2 <- factor(BX.kappa$PRCONS, levels = c("Kappa"))
BX.lambda$PRCONS2 <- factor(BX.lambda$PRCONS, levels = c("Lambda"))

#################################################################################################
################################### GENE AND GENE FAMILY PLOTS ##################################
#################################################################################################

#BX$GENE <- getGene(BX$V_CALL, first=TRUE, strip_d=TRUE)
#BX$GENE <- as.character(BX$GENE)


###############
## if using clone-pass
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


#################################################################################################
# Plot V family clonal usage by sample and PRCONS
gfh <- ggplot(BX.H.family, aes(x=GENE, y=CLONE_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  ggtitle("Clonal Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene") + scale_color_discrete(name="Gene") +
  geom_point(aes(sample="BX.H"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gfh)

gfk <- ggplot(BX.kappa.family, aes(x=GENE, y=CLONE_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  ggtitle("Clonal Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene") + scale_color_discrete(name="Gene") +
  geom_point(aes(sample="BX.kappa"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gfk)

gfl <- ggplot(BX.lambda.family, aes(x=GENE, y=CLONE_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  ggtitle("Clonal Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene") + scale_color_discrete(name="Gene") +
  geom_point(aes(sample="BX.lambda"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gfl)

# Plot V gene clonal usage by sample and PRCONS
BX.H.gene$GF <- substring(BX.H.gene$GENE, 1,5)
BX.kappa.gene$GF <- substring(BX.kappa.gene$GENE, 1,5)
BX.lambda.gene$GF <- substring(BX.lambda.gene$GENE, 1,6)

gh <- ggplot(BX.H.gene, aes(x=GENE, y=CLONE_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  ggtitle("Clonal Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_color_discrete(name="Gene Family") +
   geom_point(aes(sample="BX.H"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gh)

gk <- ggplot(BX.kappa.gene, aes(x=GENE, y=CLONE_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  ggtitle("Clonal Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_color_discrete(name="Gene Family") +
  geom_point(aes(sample="BX.kappa"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gk)

gl <- ggplot(BX.lambda.gene, aes(x=GENE, y=CLONE_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  ggtitle("Clonal Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_color_discrete(name="Gene Family") +
  geom_point(aes(sample="BX.lambda"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gl)

### multiple plots - then save as pdf
#ggsave(multiplot(gfh,gfk,gfl, cols=2, layout=matrix(c(1,2,1,3))), file = "test.pdf")

layouthkl <- rbind(c(1,1,2),
             c(1,1,3))
#grid.arrange(gfh,gfk,gfl, layout_matrix = layouthkl)
gfplots1 <- grid.arrange(gfh,gfk,gfl, layout_matrix = layouthkl)
ggsave("genefamily.png", gfplots1, width = 24, height = 8, units = "in")
ggsave("genefamily.pdf", gfplots1, width = 24, height = 8, units = "in")

gplots1 <- grid.arrange(gh,gk,gl, layout_matrix = layouthkl)
ggsave("genes.png", gplots1, width = 24, height = 8, units = "in")
ggsave("genes.pdf", gplots1, width = 24, height = 8, units = "in")

####### BY COPY NUMBER INSTEAD OF CLONAL USAGE
## if not using clone-pass
BX.H.family <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), mode="family")
BX.H.gene <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), mode="gene")
BX.L.family <- countGenes(BX.L, gene="V_CALL", groups=c("PRCONS2"), mode="family")
BX.L.gene <- countGenes(BX.L, gene="V_CALL", groups=c("PRCONS2"), mode="gene")

BX.kappa.family <- countGenes(BX.kappa, gene="V_CALL", groups=c("PRCONS2"), mode="family")
BX.kappa.gene <- countGenes(BX.kappa, gene="V_CALL", groups=c("PRCONS2"), mode="gene")
BX.lambda.family <- countGenes(BX.lambda, gene="V_CALL", groups=c("PRCONS2"), mode="family")
BX.lambda.gene <- countGenes(BX.lambda, gene="V_CALL", groups=c("PRCONS2"), mode="gene")

BX.lambda.gene$GENE <- gsub("LV", "LV0", BX.lambda.gene$GENE)
BX.lambda.gene$GENE <- gsub("LV010", "LV10", BX.lambda.gene$GENE)
BX.lambda.family$GENE <- gsub("LV", "LV0", BX.lambda.family$GENE)
BX.lambda.family$GENE <- gsub("LV010", "LV10", BX.lambda.family$GENE)


# Plot V gene, gene family usage by sample and PRCONS
gfhs <- ggplot(BX.H.family, aes(x=GENE, y=SEQ_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  ggtitle("Copy Number") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene") + scale_color_discrete(name="Gene") +
  geom_point(aes(sample="BX.H"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gfhs)

gfks <- ggplot(BX.kappa.family, aes(x=GENE, y=SEQ_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  ggtitle("Copy Number") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene") + scale_color_discrete(name="Gene") +
  geom_point(aes(sample="BX.kappa"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gfks)

gfls <- ggplot(BX.lambda.family, aes(x=GENE, y=SEQ_FREQ, fill=GENE, color=GENE)) +
  theme_bw() +
  ggtitle("Copy Number") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene Family") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene") + scale_color_discrete(name="Gene") +
  geom_point(aes(sample="BX.lambda"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gfls)

# Plot V gene by sample and PRCONS
BX.H.gene$GF <- substring(BX.H.gene$GENE, 1,5)
BX.kappa.gene$GF <- substring(BX.kappa.gene$GENE, 1,5)
BX.lambda.gene$GF <- substring(BX.lambda.gene$GENE, 1,6)

ghs <- ggplot(BX.H.gene, aes(x=GENE, y=SEQ_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  ggtitle("Copy Number") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_color_discrete(name="Gene Family") +
  geom_point(aes(sample="BX.H"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(ghs)

gks <- ggplot(BX.kappa.gene, aes(x=GENE, y=SEQ_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  ggtitle("Copy Number") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_color_discrete(name="Gene Family") +
  geom_point(aes(sample="BX.kappa"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gks)

gls <- ggplot(BX.lambda.gene, aes(x=GENE, y=SEQ_FREQ, fill=GF, color=GF)) +
  theme_bw() +
  ggtitle("Copy Number") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ylab("Percent of repertoire") +
  xlab("Gene") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name="Gene Family") + scale_color_discrete(name="Gene Family") +
  geom_point(aes(sample="BX.lambda"), size=5, alpha=1) +
  facet_wrap(~ PRCONS2)
#plot(gls)

gfsplots1 <- grid.arrange(gfhs,gfks,gfls, layout_matrix = layouthkl)
ggsave("genefamily_bysample.png", gfsplots1, width = 24, height = 8, units = "in")
ggsave("genefamily_bysample.pdf", gfsplots1, width = 24, height = 8, units = "in")

gsplots1 <- grid.arrange(ghs,gks,gls, layout_matrix = layouthkl)
ggsave("genes_bysample.png", gsplots1, width = 24, height = 8, units = "in")
ggsave("genes_bysample.pdf", gsplots1, width = 24, height = 8, units = "in")

#################################################################################################
#################################### MUTATION ANALYSIS ##########################################
#################################################################################################
## note need to reset all categories first
BX.H <- subset(BX, PRCONS %in% c("IgM", "IgD", "IgG", "IgA"))
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

BX$PRCONS2 <- factor(BX$PRCONS, levels = c("IgM", "IgD", "IgG", "IgA", "Kappa", "Lambda"))
BX.H$PRCONS2 <- factor(BX.H$PRCONS, levels = c("IgM", "IgD", "IgG", "IgA"))
BX.L$PRCONS2 <- factor(BX.L$PRCONS, levels = c("Kappa", "Lambda"))
BX.kappa$PRCONS2 <- factor(BX.kappa$PRCONS, levels = c("Kappa"))
BX.lambda$PRCONS2 <- factor(BX.lambda$PRCONS, levels = c("Lambda"))

## if using clone-pass
BX.H.family <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="family")
BX.H.gene <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="gene")
BX.L.family <- countGenes(BX.L, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="family")
BX.L.gene <- countGenes(BX.L, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="gene")
BX.kappa.family <- countGenes(BX.kappa, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="family")
BX.kappa.gene <- countGenes(BX.kappa, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="gene")
BX.lambda.family <- countGenes(BX.lambda, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="family")
BX.lambda.gene <- countGenes(BX.lambda, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="gene")

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

#BX_hobs <- BX_obs
# merge two data frames
BX_hobs <- merge(BX_hobs,BX.H.gene,by=c("PRCONS2","GENE"))
BX_hobs$MU_WEIGHT <- rescale(BX_hobs$MU_FREQ)
BX_hobs$CLONE_WEIGHT <- rescale(BX_hobs$CLONE_FREQ)

BX_kobs <- merge(BX_kobs,BX.kappa.gene,by=c("PRCONS2","GENE"))
BX_kobs$MU_WEIGHT <- rescale(BX_kobs$MU_FREQ)
BX_kobs$CLONE_WEIGHT <- rescale(BX_kobs$CLONE_FREQ)

BX_lobs <- merge(BX_lobs,BX.lambda.gene,by=c("PRCONS2","GENE"))
BX_lobs$MU_WEIGHT <- rescale(BX_lobs$MU_FREQ)
BX_lobs$CLONE_WEIGHT <- rescale(BX_lobs$CLONE_FREQ)

## only after merging should we change LV names for plotting
BX_lobs$GENE <- gsub("LV", "LV0", BX_lobs$GENE)
BX_lobs$GENE <- gsub("LV010", "LV10", BX_lobs$GENE)
BX.lambda.gene$GENE <- gsub("LV", "LV0", BX.lambda.gene$GENE)
BX.lambda.gene$GENE <- gsub("LV010", "LV10", BX.lambda.gene$GENE)
BX.lambda.family$GENE <- gsub("LV", "LV0", BX.lambda.family$GENE)
BX.lambda.family$GENE <- gsub("LV010", "LV10", BX.lambda.family$GENE)
BX_lobs$FAMILY <- gsub("LV", "LV0", BX_lobs$FAMILY)
BX_lobs$FAMILY <- gsub("LV010", "LV10", BX_lobs$FAMILY)

### plots
ghmutvs <- ggplot(BX_hobs, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=MU_WEIGHT)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Gene") + ylab("Mutation frequency") +
  scale_fill_discrete(name="Family") + scale_color_discrete(name="Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))
#plot(gmutvs)

# ghmutv <- ggplot(BX_obs, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=CLONE_WEIGHT)) +
#   theme_bw() + ggtitle("Total mutations") +
#   xlab("Gene") + ylab("Mutation frequency") +
#   scale_fill_discrete(name="Family") + scale_color_discrete(name="Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
#   geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))
#plot(gmutv) + scale_y_continuous(labels = scales::percent)

## saving as 1800 x 900


# ghcdr3v <- ggplot(BX_obs, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=CLONE_WEIGHT)) +
#   theme_bw() + ggtitle("CDRH3 Length") +
#   xlab("Gene") + ylab("CDRH3 Length, Kabat (aa)") +
#   scale_fill_discrete(name="Family") + scale_color_discrete(name="Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) +
#   geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))
#plot(gcdr3v) + scale_y_continuous(labels = scales::percent)

### LC plots
## saving as 900 x 450

gkmutv <- ggplot(BX_kobs, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=CLONE_WEIGHT)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Gene") + ylab("Mutation frequency") +
  scale_fill_discrete(name="Family") + scale_color_discrete(name="Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))
#plot(gkmutv)

gkcdr3v <- ggplot(BX_kobs, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=CLONE_WEIGHT)) +
  theme_bw() + ggtitle("CDRL3 Length") +
  xlab("Gene") + ylab("CDRL3 Length, Kabat (aa)") +
  scale_fill_discrete(name="Family") + scale_color_discrete(name="Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))
#plot(gkcdr3v) + scale_y_continuous(labels = scales::percent)

glmutv <- ggplot(BX_lobs, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=CLONE_WEIGHT)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Gene") + ylab("Mutation frequency") +
  scale_fill_discrete(name="Family") + scale_color_discrete(name="Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))
#plot(glmutv)

glcdr3v <- ggplot(BX_lobs, aes(x=GENE, y=CDR3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=CLONE_WEIGHT)) +
  theme_bw() + ggtitle("CDRL3 Length") +
  xlab("Gene") + ylab("CDRL3 Length, Kabat (aa)") +
  scale_fill_discrete(name="Family") + scale_color_discrete(name="Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))
#plot(glcdr3v) + scale_y_continuous(labels = scales::percent)

layouthkl <- rbind(c(1,1,2),
                   c(1,1,3))
mutplots1 <- grid.arrange(ghmutvs,gkmutv,glmutv, layout_matrix = layouthkl)
ggsave("mutation_bysample.png", mutplots1, width = 24, height = 8, units = "in")
ggsave("mutation_bysample.pdf", mutplots1, width = 24, height = 8, units = "in")

# cdr3plots1 <- grid.arrange(ghcdr3v,gkcdr3v,glcdr3v, layout_matrix = layouthkl)
# ggsave("CDR3_bysample.png", cdr3plots1, width = 24, height = 8, units = "in")
# ggsave("CDR3_bysample.pdf", cdr3plots1, width = 24, height = 8, units = "in")


##########################################################################################################
##########################################################################################################
#############################       LINEAGE RECONSTRUCTION         #######################################
##########################################################################################################
##########################################################################################################


#### FOR PARTICULAR CLONE
# BX <- read.delim(opt$DB)
# BX <- read.delim("/users/eric.waltari/immcantation_pipeline/BX_min5/changeo/BXjan5_matcheswithss_clone-pass.tab")
#BX2 <- read.delim("/users/eric.waltari/immcantation_pipeline/BX_min5/changeo/BXjan5_matcheswithss_clone-pass2.tab")
#BX2s <- read.delim("/users/eric.waltari/immcantation_pipeline/BX_min5/changeo/BXjan5_matcheswithss_clone-pass2_min2b.tab")
#BX2t <- read.delim("/users/eric.waltari/immcantation_pipeline/BX_min5/changeo/BXjan5_matcheswithss_clone-pass2_min2and1m.tab")
#BX2u <- read.delim("/users/eric.waltari/immcantation_pipeline/BX_min5/changeo/BXjan5_matcheswithss_clone-pass2_min2and1m_b.tab")
#BX2v <- read.delim("/users/eric.waltari/immcantation_pipeline/BX_min5/changeo/BXjan5_matcheswithss_clone-pass2_min2and1m_c.tab")
#BX <- read.delim("/users/eric.waltari/immcantation_pipeline/BX1mil_min2/BXjan2_withss_clone-pass.tab")
#BXt <- read.delim("/users/eric.waltari/immcantation_pipeline/BX_min5/changeo/BX_clone-pass_test.tab")
#BXt.sub <- subset(BXt, CLONE == 15632)

# BX.sub <- subset(BX, CLONE == 27281)
# BX.sub <- subset(BX, CLONE == 114872)
# BX.sub <- subset(BX, CLONE == opt$CLONE)

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

# ## if only >5...
# #BXt.sub$CONSCOUNT2 <- (2 * log(BXt.sub$CONSCOUNT))
# #BXt.sub$CONSCOUNT2 <- round(BXt.sub$CONSCOUNT2/5)*5

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
# # dnapars_exec <- "~/phylip/exe/dnapars"
# dnapars_exec <- "/usr/local/bin/dnapars"
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

# Add legend
#legend("topright", c("Germline", "Inferred", "Sample"), 
#       fill=c("black", "white", "steelblue"), cex=0.75)



#############################
## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
# Preprocess clones
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

