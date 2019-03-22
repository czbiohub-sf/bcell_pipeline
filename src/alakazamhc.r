## ALL SHAZAM & ALAKAZAM COMMANDS

# Import required packages
install.packages("gdata")
install.packages("tidyverse")
install.packages("hexbin")
install.packages("alakazam")
install.packages("shazam")
install.packages("scales")
install.packages("igraph")
#install.packages("grid")
install.packages("lattice")
install.packages("gridExtra")
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

BX$FAMILY <- getFamily(BX$V_CALL, first=TRUE, strip_d=TRUE)
BX$GENE <- getGene(BX$V_CALL, first=TRUE, strip_d=TRUE)
BX$FAMILY <- as.character(BX$FAMILY)
BX$GENE <- as.character(BX$GENE)

## CHANGING TO REFLECT NEW CREGION ANNOTATION AND USING THAT TO DETERMINE PRCONS2
BX.H <- subset(BX, CREGION %in% c("IgM", "IgG", "IgA"))

BX.H <- BX.H[ grep("IGKV", BX.H$V_CALL, invert = TRUE) , ]
BX.H <- BX.H[ grep("IGLV", BX.H$V_CALL, invert = TRUE) , ]

BX$PRCONS2 <- factor(BX$CREGION, levels = c("IgM", "IgG", "IgA", "Kappa", "Lambda"))
BX.H$PRCONS2 <- factor(BX.H$CREGION, levels = c("IgM", "IgG", "IgA"))

#################################################################################################
################################### GENE AND GENE FAMILY PLOTS ##################################
#################################################################################################

## gene and gene family usage by clone
BX.H.family <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="family")
BX.H.gene <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), clone="CLONE", mode="gene")
BX.H.gene$GF <- substring(BX.H.gene$GENE, 1,5)

## gene and gene family usage by read
BXs.H.family <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), mode="family")
BXs.H.gene <- countGenes(BX.H, gene="V_CALL", groups=c("PRCONS2"), mode="gene")
BXs.H.gene$GF <- substring(BXs.H.gene$GENE, 1,5)


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



# merge two data frames adding both gene counts by clone and by read...
BX_hobs <- merge(BX_hobs,BX.H.gene,by=c("PRCONS2","GENE"))

BX_hobs <- merge(BX_hobs,BXs.H.gene,by=c("PRCONS2","GENE"))

## renaming gene counts/weights
BX_hobs <- rename(BX_hobs, GENEFREQ_BYCLONE = CLONE_FREQ)
BX_hobs <- rename(BX_hobs, GENECOUNT_BYCLONE = CLONE_COUNT)
BX_hobs <- rename(BX_hobs, GENEFREQ_BYREAD = SEQ_FREQ)
BX_hobs <- rename(BX_hobs, GENECOUNT_BYREAD = SEQ_COUNT)
BX_hobs$GENEFREQWEIGHT_BYCLONE <- rescale(BX_hobs$GENEFREQ_BYCLONE)
BX_hobs$GENEFREQWEIGHT_BYREAD <- rescale(BX_hobs$GENEFREQ_BYREAD)


###################################
### moved all subsequent calculations here - then save .tsv files...
BX_hobs <- BX_hobs %>% add_count(CLONE)

## removing GF.x and GF.y
BX_hobs$GF.x <- NULL
BX_hobs$GF.y <- NULL

## TRY ADDING AS.CHARACTER HERE TO PRCONS2, GANDA_SUBTYPE (PRCONS AND CREGION ALREADY ARE)
BX_hobs$PRCONS2 <- as.character(BX_hobs$PRCONS2)
BX_hobs$GANDA_SUBTYPE <- as.character(BX_hobs$GANDA_SUBTYPE)

### creating lists per clone not per read
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
clonestatsh3 <- rename(clonestatsh3, PRCONS2_d = PRCONS2)
clonestatsh3 <- rename(clonestatsh3, GANDA_SUBTYPE_d = GANDA_SUBTYPE)
clonestatsh4 <- inner_join(clonestatsh1, clonestatsh2)
clonestatsh <- inner_join(clonestatsh4, clonestatsh3)

## filter to remove clones found only once!
clonestatshf <- clonestatsh %>% filter(n > 1)

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
## for reflow test
hmeans_and_sds <- hmeans1 %>%
  left_join(hsds1, by='PRCONS2') %>%
  left_join(hmeans2, by='PRCONS2') %>%
  left_join(hsds2, by='PRCONS2') %>%
  left_join(hmeans3, by='PRCONS2') %>%
  left_join(hsds3, by='PRCONS2')

allmeans_and_sds <- hmeans_and_sds

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
## for reflow test
submeans_and_sds <- submeans1 %>%
  left_join(subsds1, by='GANDA_SUBTYPE') %>%
  left_join(submeans2, by='GANDA_SUBTYPE') %>%
  left_join(subsds2, by='GANDA_SUBTYPE') %>%
  left_join(submeans3, by='GANDA_SUBTYPE') %>%
  left_join(subsds3, by='GANDA_SUBTYPE')

#submeans_and_sds <- left_join(submeans1,subsds1,submeans2,subsds2,submeans3,subsds3, by = "GANDA_SUBTYPE", suffix = c(".x", ".y"))
#submeans_and_sds <- merge(submeans1,subsds1,submeans2,subsds2,submeans3,subsds3, by = "GANDA_SUBTYPE")


## Oct 1 Alakazam errors occurring around here...first replacing cbind with bind_cols
## after this replacement the bind_cols for lmeans_and_sds works, but the second one just here gives error:
## PROBABABLY BECAUSE OF NA IN SUBTYPES (EVERYTHING NOT IGG OR IGA) - TRY SET_NAMES TO RENAME NA TO OTHER...
##        Error in cbind_all(x) : Argument 5 must be length 7, not 6     Calls: bind_cols -> cbind_all -> .Call

## 9p 10/2 - still getting error with bind_cols below in reflow version: think because of NAs still in some files below...

###################################
### stats for pie charts
hcmeans_and_sds <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM", "IgG", "IgA"))

## 10/4 - taking subsets from here now...
mumeans_and_sds <- subset(allmeans_and_sds, PRCONS2 %in% c("IgM"))
mulcmeans_and_sds <- mumeans_and_sds  
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
isosubmeans_and_sds <- isosubmeans_and_sds %>% mutate_at(vars(n1,n2,n3,mutation_sd_reads,mutation_mean_reads,mutation_sd_clones,mutation_mean_clones,mutation_sd_filteredclones,mutation_mean_filteredclones), replace_na, 0)
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

## adding tibbles to reorder allmeans and hcmeans
hctarget <- tibble(PRCONS2 = c("IgM", "IgG", "IgA"))
hcmeans_and_sds <- left_join(data.frame(PRCONS2=hctarget),hcmeans_and_sds,by="PRCONS2")

## ADDING COLLAPSE INTO CLONES FROM SHAZAM SCRIPTS, THEN SAVING AS TSV FILES
clonesh <- collapseClones(BX_hobs, regionDefinition=IMGT_V, 
                          method="thresholdedFreq", minimumFrequency=0.6,
                          includeAmbiguous=FALSE, breakTiesStochastic=FALSE, 
                          nproc=4)
write.table(clonesh, "onereadperclone_h.tsv", sep = "\t", row.names = FALSE)

## saving all created files as .tsv files..
write.table(BX_hobs, "mutstats_h.tsv", sep = "\t", row.names = FALSE)
write.table(clonestatsh, "mutclonestats_h.tsv", sep = "\t", row.names = FALSE)
write.table(clonestatshf, "mutclonestatsfiltered_h.tsv", sep = "\t", row.names = FALSE)
#write.table(isotypeandsubtypemeans_and_sds, "summary_mutationstats_all.tsv", sep = "\t", row.names = FALSE)

## NOW FACTOR EVERYTHING FOR GRAPHING
hcmeans_and_sds$PRCONS2 <- factor(hcmeans_and_sds$PRCONS2, levels = c("IgM", "IgG", "IgA"))
subtypemeans_and_sds1$GANDA_SUBTYPE <- factor(subtypemeans_and_sds1$GANDA_SUBTYPE, levels = c("IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
hcisosubmeans_and_sds0$GANDA_SUBTYPE <- factor(hcisosubmeans_and_sds0$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
hcisosubmeans_and_sds$GANDA_SUBTYPE <- factor(hcisosubmeans_and_sds$GANDA_SUBTYPE, levels = c("IgM", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))

BX.H$PRCONS2 <- factor(BX.H$PRCONS2, levels = c("IgM", "IgG", "IgA"))
BX_hobs$PRCONS2 <- factor(BX_hobs$PRCONS2, levels = c("IgM", "IgG", "IgA"))
clonestatsh$PRCONS2 <- factor(clonestatsh$PRCONS2, levels = c("IgM", "IgG", "IgA"))
clonestatshf$PRCONS2 <- factor(clonestatshf$PRCONS2, levels = c("IgM", "IgG", "IgA"))


###### testing RELOADING DATA later
#BX_hobs <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutstats_h.tsv")
#BX_kobs <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutstats_k.tsv")
#BX_lobs <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutstats_l.tsv")
#clonestatsh <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestats_h.tsv")
#clonestatsk <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestats_k.tsv")
#clonestatsl <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestats_l.tsv")
###clonestatshf <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestatsfiltered_h.tsv")
###clonestatskf <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestatsfiltered_k.tsv")
###clonestatslf <- read.delim("/users/eric.waltari/changeo/BX_miseq/mutclonestatsfiltered_l.tsv")


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


### force mutation frequency (y-axis here) for all HC to be 0-35%
ghmutvh35 <- ggplot(BX_hobs, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYREAD)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(ghmutv)


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

mutplots1c <- grid.arrange(ghmutvc,gkmutvc,glmutvc, layout_matrix = layouthkl3)
ggsave("mutation_bygene_byclone.png", mutplots1c, width = 16, height = 12, units = "in")
ggsave("mutation_bygene_byclone.pdf", mutplots1c, width = 16, height = 12, units = "in")
cdr3plots1c <- grid.arrange(ghcdr3vc,gkcdr3vc,glcdr3vc, layout_matrix = layouthkl3)
ggsave("CDR3_bygene_byclone.png", cdr3plots1c, width = 16, height = 12, units = "in")
ggsave("CDR3_bygene_byclone.pdf", cdr3plots1c, width = 16, height = 12, units = "in")
###

### force mutation frequency (y-axis here) for all HC to be 0-35%
ghmutvch35 <- ggplot(clonestatsh, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gmutv) + scale_y_continuous(labels = scales::percent)

### by filtered clone
ghmutvcf <- ggplot(clonestatshf, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gmutv) + scale_y_continuous(labels = scales::percent)
ghcdr3vcf <- ggplot(clonestatshf, aes(x=GENE, y=CDRH3KABAT_LENGTH, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() + ggtitle("CDR3 Length by gene") +
  xlab("Gene") + ylab("CDRH3 Length, Kabat (aa)") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) +
  geom_violin(width=1) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gcdr3v) + scale_y_continuous(labels = scales::percent)

#mutplots1cf <- grid.arrange(ghmutvcf,gkmutvcf,glmutvcf, layout_matrix = layouthkl3)
#ggsave("mutation_bygene_byclone_filtered.png", mutplots1cf, width = 16, height = 12, units = "in")
#ggsave("mutation_bygene_byclone_filtered.pdf", mutplots1cf, width = 16, height = 12, units = "in")
#cdr3plots1cf <- grid.arrange(ghcdr3vcf,gkcdr3vcf,glcdr3vcf, layout_matrix = layouthkl3)
#ggsave("CDR3_bygene_byclone_filtered.png", cdr3plots1cf, width = 16, height = 12, units = "in")
#ggsave("CDR3_bygene_byclone_filtered.pdf", cdr3plots1cf, width = 16, height = 12, units = "in")


### force mutation frequency (y-axis here) for all HC to be 0-35%
ghmutvcfh35 <- ggplot(clonestatshf, aes(x=GENE, y=MU_FREQ, fill=FAMILY, color=FAMILY, stroke = 0.001, alpha=GENEFREQ_BYCLONE)) +
  theme_bw() + ggtitle("Mutation distribution by Gene") +
  xlab("Gene") + ylab("% Somatic Hypermutation") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + facet_wrap(~ PRCONS2, ncol=1) + scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_violin(width=1.25) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#plot(gmutv) + scale_y_continuous(labels = scales::percent)

#mutplots1cfh35 <- grid.arrange(ghmutvcfh35,gkmutvcf,glmutvcf, layout_matrix = layouthkl3)
#ggsave("mutation_bygene_byclone_filteredh35.png", mutplots1cfh35, width = 16, height = 12, units = "in")
#ggsave("mutation_bygene_byclone_filteredh35.pdf", mutplots1cfh35, width = 16, height = 12, units = "in")


###################################
### histograms of mutation
###################################
## by read
ghmuthistogram <- ggplot(BX_hobs, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol = 1) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(ghmuthistogram)
### free y axis for % and not just counts below
ghmuthistogramfree <- ggplot(BX_hobs, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol = 1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(ghmuthistogramfree)


ghmuthistogramcounts <- ggplot(BX_hobs, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(ghmuthistogramcounts)


### force x axis for all HC to be 0-35%
ghmuthistogramh35 <- ggplot(BX_hobs, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol = 1) + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Reads")
#plot(ghmuthistogramh35)
ghmuthistogramcountsh35 <- ggplot(BX_hobs, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Reads")
#plot(ghmuthistogramcountsh35)

## same mutation plots but by clone
ghmuthistogrambyclone0 <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ghmuthistogrambyclone0)
ghmuthistogrambyclone <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ghmuthistogrambyclone)
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


### force x axis for all HC to be 0-35%
ghmuthistogrambycloneh35 <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1) + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of Clones")
#plot(ghmuthistogrambycloneh35)
ghmuthistogramcountsbycloneh35 <- ggplot(clonestatsh, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Clones")
#plot(ghmuthistogramcountsbycloneh35)


###################
## now clone plots but filter so leaving out all single clones!
ghmuthistogrambyclonef <- ggplot(clonestatshf, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of clones with ≥ 2 Members")
#plot(ghmuthistogrambyclonef)
### free y axis for % and not just counts below
ghmuthistogrambycloneffree <- ggplot(clonestatshf, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of clones with ≥ 2 Members")
#plot(ghmuthistogrambycloneffree)

ghmuthistogramcountsbyclonef <- ggplot(clonestatshf, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Clones with ≥ 2 Members")
#plot(ghmuthistogramcountsbyclonef)

#mutplotshistbyclonef <- grid.arrange(ghmuthistogrambyclonef,gkmuthistogrambyclonef,glmuthistogrambyclonef, layout_matrix = layouthkl3)
#ggsave("mutation_histogram_byclone_filtered_hlseparatepanels.png", mutplotshistbyclonef, width = 16, height = 12, units = "in")
#ggsave("mutation_histogram_byclone_filtered_hlseparatepanels.pdf", mutplotshistbyclonef, width = 16, height = 12, units = "in")
#mutplotshistcountsbyclonef <- grid.arrange(ghmuthistogramcountsbyclonef,gkmuthistogramcountsbyclonef,glmuthistogramcountsbyclonef, layout_matrix = layouthkl3)
#ggsave("mutation_histogram_rawcounts_byclone_filtered_hlseparatepanels.png", mutplotshistcountsbyclonef, width = 16, height = 12, units = "in")
#ggsave("mutation_histogram_rawcounts_byclone_filtered_hlseparatepanels.pdf", mutplotshistcountsbyclonef, width = 16, height = 12, units = "in")

### force x axis for all HC to be 0-35%
ghmuthistogrambyclonefh35 <- ggplot(clonestatshf, aes(x = MU_FREQ)) + geom_histogram(aes(y=0.01*..density..), binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1) + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + scale_y_continuous(labels = scales::percent) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation") + xlab("% Somatic Hypermutation") + ylab("Proportion of clones with ≥ 2 Members")
#plot(ghmuthistogrambyclonef)
ghmuthistogramcountsbyclonefh35 <- ggplot(clonestatshf, aes(x = MU_FREQ)) + geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ PRCONS2, ncol=1, scales = "free_y") + scale_x_continuous(labels = scales::percent, limits = c(-0.01, .35)) + theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) + ggtitle("Somatic Hypermutation, Raw Counts") + xlab("% Somatic Hypermutation") + ylab("Clones with ≥ 2 Members")
#plot(ghmuthistogramcountsbyclonef)

#mutplotshistbyclonefh35 <- grid.arrange(ghmuthistogrambyclonefh35,gkmuthistogrambyclonef,glmuthistogrambyclonef, layout_matrix = layouthkl3)
#ggsave("mutation_histogram_byclone_filtered_hlseparatepanelsh35.png", mutplotshistbyclonefh35, width = 16, height = 12, units = "in")
#ggsave("mutation_histogram_byclone_filtered_hlseparatepanelsh35.pdf", mutplotshistbyclonefh35, width = 16, height = 12, units = "in")
#mutplotshistcountsbyclonefh35 <- grid.arrange(ghmuthistogramcountsbyclonefh35,gkmuthistogramcountsbyclonef,glmuthistogramcountsbyclonef, layout_matrix = layouthkl3)
#ggsave("mutation_histogram_rawcounts_byclone_filtered_hlseparatepanelsh35.png", mutplotshistcountsbyclonefh35, width = 16, height = 12, units = "in")
#ggsave("mutation_histogram_rawcounts_byclone_filtered_hlseparatepanelsh35.pdf", mutplotshistcountsbyclonefh35, width = 16, height = 12, units = "in")

#mutplotshistbycloneffree <- grid.arrange(ghmuthistogrambycloneffree,gkmuthistogrambyclonef,glmuthistogrambyclonef, layout_matrix = layouthkl3)
#ggsave("mutation_histogram_byclone_filtered_hlseparatefreepanels.png", mutplotshistbycloneffree, width = 16, height = 12, units = "in")
#ggsave("mutation_histogram_byclone_filtered_hlseparatefreepanels.pdf", mutplotshistbycloneffree, width = 16, height = 12, units = "in")

gmutviolinhf <- ggplot(clonestatshf, aes(x=PRCONS2, y=MU_FREQ, stroke = 0.001)) +
  theme_bw() + ggtitle("% Somatic Hypermutation by Clone") +
  xlab("Isotype") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_fill_brewer(palette = "Paired", name="Gene Family") + scale_colour_brewer(palette = "Paired", name="Gene Family") + scale_alpha(guide = "none") + scale_y_continuous(labels = scales::percent) +
  geom_violin(fill = "gray") + theme(axis.text.x = element_text(angle=45, hjust=1))
#plot(gmutviolinhf)
ggsave("mutation_violinplot_byclone_filtered_h.png", gmutviolinhf, width = 16, height = 8, units = "in")


################
### RBIND ALL HC AND LC IN SINGLE PANELS - DO NOT LOOK AS GOOD THOUGH
## removing all for HC
###################################
##### mutation vs CDR3 hex plots
## by read
gmutandcdr3hexhr <- ggplot(BX_hobs, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexhr)
#mutvsCDR3r <- grid.arrange(gmutandcdr3hexhr,gmutandcdr3hexkr,gmutandcdr3hexlr, layout_matrix = layouthkl3)
#ggsave("mutationvsCDR3_byread_hlseparatepanels.png", mutvsCDR3r, width = 16, height = 12, units = "in")
#ggsave("mutationvsCDR3_byread_hlseparatepanels.pdf", mutvsCDR3r, width = 16, height = 12, units = "in")

### force x axis for all HC to be 0-35%
gmutandcdr3hexhrh35 <- ggplot(BX_hobs, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Read") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Reads",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexhr)
#mutvsCDR3rh35 <- grid.arrange(gmutandcdr3hexhrh35,gmutandcdr3hexkr,gmutandcdr3hexlr, layout_matrix = layouthkl3)
#ggsave("mutationvsCDR3_byread_hlseparatepanelsh35.png", mutvsCDR3rh35, width = 16, height = 12, units = "in")
#ggsave("mutationvsCDR3_byread_hlseparatepanelsh35.pdf", mutvsCDR3rh35, width = 16, height = 12, units = "in")


## by clone
gmutandcdr3hexh <- ggplot(clonestatsh, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexh)
#mutvsCDR3 <- grid.arrange(gmutandcdr3hexh,gmutandcdr3hexk,gmutandcdr3hexl, layout_matrix = layouthkl3)
#ggsave("mutationvsCDR3_byclone_hlseparatepanels.png", mutvsCDR3, width = 16, height = 12, units = "in")
#ggsave("mutationvsCDR3_byclone_hlseparatepanels.pdf", mutvsCDR3, width = 16, height = 12, units = "in")

### force x axis for all HC to be 0-35%
gmutandcdr3hexhh35 <- ggplot(clonestatsh, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexh)
#mutvsCDR3h35 <- grid.arrange(gmutandcdr3hexhh35,gmutandcdr3hexk,gmutandcdr3hexl, layout_matrix = layouthkl3)
#ggsave("mutationvsCDR3_byclone_hlseparatepanelsh35.png", mutvsCDR3h35, width = 16, height = 12, units = "in")
#ggsave("mutationvsCDR3_byclone_hlseparatepanelsh35.pdf", mutvsCDR3h35, width = 16, height = 12, units = "in")

## by filtered clone
gmutandcdr3hexhf <- ggplot(clonestatshf, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexhf)
#mutvsCDR3f <- grid.arrange(gmutandcdr3hexhf,gmutandcdr3hexkf,gmutandcdr3hexlf, layout_matrix = layouthkl3)
#ggsave("mutationvsCDR3_byclone_filtered_hlseparatepanels.png", mutvsCDR3f, width = 16, height = 12, units = "in")
#ggsave("mutationvsCDR3_byclone_filtered_hlseparatepanels.pdf", mutvsCDR3f, width = 16, height = 12, units = "in")

### force x axis for all HC to be 0-35%
gmutandcdr3hexhfh35 <- ggplot(clonestatshf, aes(x=CDRH3KABAT_LENGTH, y=MU_FREQ)) +
  theme_bw() + ggtitle("Somatic Hypermutation & CDR3") +
  xlab("CDRH3 Length, Kabat (aa)") + ylab("Average % Somatic Hypermutation per Clone") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandcdr3hexhf)
#mutvsCDR3fh35 <- grid.arrange(gmutandcdr3hexhfh35,gmutandcdr3hexkf,gmutandcdr3hexlf, layout_matrix = layouthkl3)
#ggsave("mutationvsCDR3_byclone_filtered_hlseparatepanelsh35.png", mutvsCDR3fh35, width = 16, height = 12, units = "in")
#ggsave("mutationvsCDR3_byclone_filtered_hlseparatepanelsh35.pdf", mutvsCDR3fh35, width = 16, height = 12, units = "in")


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

#gmutandnhexhbyisotype <- grid.arrange(gmutandnhexhbyisotypeh,gmutandnhexhbyisotypek,gmutandnhexhbyisotypel, layout_matrix = layouthkl3)
#ggsave("clonalfamily_n_vs_mutation_hlseparatepanels.png", gmutandnhexhbyisotype, width = 16, height = 12, units = "in")
#ggsave("clonalfamily_n_vs_mutation_hlseparatepanels.pdf", gmutandnhexhbyisotype, width = 16, height = 12, units = "in")

### force x axis for all HC to be 0-35%
gmutandnhexhbyisotypehh35 <- ggplot(clonestatsh, aes(x=n, y=MU_FREQ)) +
  theme_bw() + ggtitle("# Reads per Clonal Family & Somatic Hypermutation") +
  xlab("# Reads per Clonal Family") + ylab("Average % Somatic Hypermutation") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + scale_y_continuous(labels = scales::percent, limits = c(-0.01, .35)) +
  geom_hex(aes(fill=log10(..count..))) + facet_wrap(~ PRCONS2, ncol=1) + scale_fill_gradient(low = "light blue", high = "magenta", name = "Number of Clones",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))
#plot(gmutandnhexhbyisotypeh)
#gmutandnhexhbyisotypeh35 <- grid.arrange(gmutandnhexhbyisotypehh35,gmutandnhexhbyisotypek,gmutandnhexhbyisotypel, layout_matrix = layouthkl3)
#ggsave("clonalfamily_n_vs_mutation_hlseparatepanelsh35.png", gmutandnhexhbyisotypeh35, width = 16, height = 12, units = "in")
#ggsave("clonalfamily_n_vs_mutation_hlseparatepanelsh35.pdf", gmutandnhexhbyisotypeh35, width = 16, height = 12, units = "in")


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

#gmutandnhexallga <- grid.arrange(gmutandnhexallh,gmutandnhexga, layout_matrix = layout2)
#ggsave("clonalfamily_n_vs_mutation_2HCformats.png", gmutandnhexallga, width = 16, height = 12, units = "in")
#ggsave("clonalfamily_n_vs_mutation_2HCformats.pdf", gmutandnhexallga, width = 16, height = 12, units = "in")

#################################################################################################
################################## PIE CHARTS OF ISOTYPES #######################################
#################################################################################################

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


## SUBTYPES...

### ISOTYPES + SUBTYPES

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


bottom = quote(paste("page", g, "of", npages))
## with new isotype color scheme
pdfset1 <- c(list(piecharthccounts2, piecharthcisosubcounts2, gfhs, ghs, ghmutvh35, ghcdr3v, gmutandcdr3hexhrh35, ghmuthistogramh35, ghmuthistogramcountsh35))
manypdfs1 <- marrangeGrob(pdfset1, nrow=1, ncol=1, bottom = quote(paste("page", g, "of", npages)))
ggsave("results_byreads.pdf", manypdfs1, width = 8, height = 10, units = "in")
pdfset2 <- c(list(piecharthcclones2, piecharthcisosubclones2, gfh, gh, ghmutvch35, ghcdr3vc, gmutandcdr3hexhh35, ghmuthistogrambycloneh35, ghmuthistogramcountsbycloneh35, ghnreadshistogrambyclone, ghnreadshistogramcountsbyclone, gmutandnhexhbyisotypehh35))
manypdfs2 <- marrangeGrob(pdfset2, nrow=1, ncol=1, bottom = quote(paste("page", g, "of", npages)))
ggsave("results_byclones.pdf", manypdfs2, width = 8, height = 10, units = "in")

###############################
## PIECHARTS WITHOUT LABELS (will be in second set of pdfs)

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
pdfset1b <- c(list(piecharthccounts2b, piecharthcisosubcounts2b, gfhs, ghs, ghmutv, ghcdr3v, gmutandcdr3hexhr, ghmuthistogram, ghmuthistogramcounts))
manypdfs1b <- marrangeGrob(pdfset1b, nrow=1, ncol=1, bottom = quote(paste("page", g, "of", npages)))
ggsave("results_byreads_blankpies.pdf", manypdfs1b, width = 8, height = 10, units = "in")
pdfset2b <- c(list(piecharthcclones2b, piecharthcisosubclones2b, gfh, gh, ghmutvc, ghcdr3vc, gmutandcdr3hexh, ghmuthistogrambyclone, ghmuthistogramcountsbyclone, ghnreadshistogrambyclone, ghnreadshistogramcountsbyclone, gmutandnhexallh))
manypdfs2b <- marrangeGrob(pdfset2b, nrow=1, ncol=1, bottom = quote(paste("page", g, "of", npages)))
ggsave("results_byclones_blankpies.pdf", manypdfs2b, width = 8, height = 10, units = "in")
