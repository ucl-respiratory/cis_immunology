################################################################################
# This script creates all figures and supplementary data for this paper.
# As described in README.md, significant data preprocessing is required to run this script.
# It is provided as a reference rather than an easily reproducible pipeline.
################################################################################

version <- "5.0"

# Output directory
figdir <- paste0("~/Dropbox/CIS_Immunology/cis_immunology2/results/v", version, "/figures/")
supdir <- paste0("~/Dropbox/CIS_Immunology/cis_immunology2/results/v", version, "/supdata/")
dir.create(figdir, showWarnings = F, recursive = T)
dir.create(supdir, showWarnings = F, recursive = T)

# Create a text file for other results quoted in the text:
opfile <- paste0(supdir, "results.in.text.txt")
write(c("# Statistical results quoted in the main text:", ""), file = opfile)

# Quick function to add a log to this file:
do.log <- function(x) {
  write(x, file = opfile, append = T)
  write("", file = opfile, append = T)
}
sep <- "###############################################################"

# Load required libraries:
libs <- c("tidyverse", "ggplot2", "ggsignif", "ggrepel", "ggpubr", "cowplot", "magick", "WriteXLS", "foreach", "ChAMP", "pathview", "pheatmap", "grid", "ggplotify", "biomaRt",
          "httr", "jsonlite", "xml2", "sva", "gdata", "BSgenome.Hsapiens.UCSC.hg19", "Homo.sapiens", "limma", "lme4", "lmerTest", "tibble", "dndscv",
          "fgsea", "gage", "emmeans", "RColorBrewer",
          "effsize", "FDRsampsize", "pwr", "WriteXLS", "Rtsne")

for(lib in libs) {
  library(lib, character.only = T)
}

# Export the current versions of required libraries to the results folder
vfile <- paste0(supdir, "package.versions.csv")
version.data <- installed.packages()[libs,]
write.csv(version.data, file = vfile)
file.copy(vfile, '.', overwrite = T)


# Figure width - assume double column = 183mm (see https://www.nature.com/nature/for-authors/formatting-guide)
# -> 7.204" (needed for save_plot base_width method)
fig.width <- 7.204
fig.maxheight <- 9.72

# Format figures in 12 point Times New Roman (suggested by Nature formatting guide)
theme_set(theme_cowplot(font_size=11, font_family = "Times"))

theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank())

# Set default colour palette - for consistency, redefine the ggplot function as follows:
# Use Dark2 from Rcolorbrewer, which is colourblind safe
pal <- 'Set2'
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette=pal, direction = -1) + scale_fill_brewer(palette = pal, direction = -1)
# Define manual colours for progressive, regressive, control in that order (used in a couple of plots)
plotcols <- brewer.pal(3, pal)[c(2,1,3)]

# Red and blu for 'hot' and 'cold' - taken from Set2 on colorbrewer2.org
hotcol <- "#e31a1c"
coldcol <- "#1f78b4"

###############################################################################################
# Load data
###############################################################################################
if(file.exists('data/cached.RData')) {
  message("Loading cached data")
  load('data/cached.RData')
} else {
  message("Refreshing data cache")
  source('load.data.R')
}

message("Loading required functions")
for(f in list.files('utility_functions', pattern=".R$", full.names = T)){
  message(paste("loading", f))
  source(f)
}

# Load some genes used in analysis
message("Loading gene lists")
load('resources/gene.lists.RData')



###############################################################################################
# Immune cell quantification from gene expression data
# Uses danaher method
# Apply this to all GXN datasets
###############################################################################################
message("Running Danaher gene expression deconvolution")
gdata.danaher <- do.danaher(gdata)
# Paired samples for comparison:
gdata.danaher.t <- do.danaher(gdata.pair.t)
gdata.danaher.s <- do.danaher(gdata.pair.s)


###############################################################################################
# Methylation differential analysis
# Done for progressive vs regressive CIS, and for TCGA cancer vs control.
# Results are cached as they are computationally intensive.
# Note that we use all methylation data here rather than a discovery set, hence results differ from those previously published. Conclusions remain very similar.
###############################################################################################
if(file.exists('data/cached.mdiff.RData')) {
  load('data/cached.mdiff.RData')
} else {
  dmps <- champ.DMP(mdata, pheno=mpheno$Sample_Group, compare.group=c('Progressive', 'Regressive'), adjPVal = 1)
  dmps <- dmps$Progressive_to_Regressive
  
  dmrs <- champ.DMR(beta=as.matrix(mdata), pheno=mpheno$Sample_Group, compare.group = c("Progressive", "Regressive"), method="ProbeLasso")
  dmrs <- dmrs$ProbeLassoDMR
  
  # Additionally calculate TCGA DMRs for comparison:
  load('data/mdata.RData')
  tcga.mpheno.tmp <- tcga.mpheno
  tcga.mpheno.tmp$Sample_Group <- make.names(tcga.mpheno.tmp$Sample_Group)
  # Impute to remove NAs
  tcga.mdata.imputed <- champ.impute(beta=as.matrix(tcga.mdata), SampleCutoff = 0.5, ProbeCutoff = 0.5, pd=tcga.mpheno.tmp)
  tcga.mpheno.tmp <- tcga.mdata.imputed$pd
  tcga.mdata.imputed <- tcga.mdata.imputed$beta
  # Strange bug - only use probes in package data(illumina450Gr) for DM (removes very few probes)
  data(illumina450Gr)
  tcga.mdata.imputed <- tcga.mdata.imputed[which(rownames(tcga.mdata.imputed) %in% names(illumina450Gr)),]
  # Find DMPs
  dmps.tcga <- champ.DMP(beta=tcga.mdata.imputed, pheno=tcga.mpheno.tmp$Sample_Group, compare.group = c("TCGA.SqCC", "TCGA.Control"))
  dmps.tcga <- dmps.tcga$TCGA.Control_to_TCGA.SqCC
  # Find DMRs
  dmrs.tcga <- champ.DMR(beta=tcga.mdata.imputed, pheno=tcga.mpheno.tmp$Sample_Group, compare.group=c("TCGA.SqCC", "TCGA.Control"), method="ProbeLasso")
  dmrs.tcga <- dmrs.tcga$ProbeLassoDMR
  
  save(dmps, dmrs, dmps.tcga, dmrs.tcga, file='data/cached.mdiff.RData')
  
  # We need to reload cached data, as loading mdata.RData overwrites mdata and mpheno
  load('data/cached.RData')
}
data('probe.features')

# Add total methylation measures
pheno$total.meth.nosnp <- NA
pheno$median.meth.nosnp <- NA
sel <- which(pheno$SampleID %in% colnames(mdata.nosnp))
pheno$total.meth.nosnp[sel] <- apply(mdata.nosnp[,pheno$SampleID[sel]], 2, function(x) {mean(x, na.rm=T)})
pheno$median.meth.nosnp[sel] <- apply(mdata.nosnp[,pheno$SampleID[sel]], 2, function(x) {median(x, na.rm=T)})


##############################################################################
# Figure 1: IHC and Danaher data
# TODO: adjust for patient ID
##############################################################################
message("Plotting Figure 1")

# Plot IHC data
plotdata <- pheno[which(pheno$IHC_K | pheno$HandE),]
# plotdata <- ihc
plotdata$Outcome <- gsub("ression", ".", plotdata$Outcome)
p.lymphocytes <- compare.fn(lymphocytes_perArea ~ Outcome + (1 | patient), data = plotdata)
p.cd4 <- compare.fn(cd4_perArea ~ Outcome + (1 | patient), data = plotdata)
p.cd8 <- compare.fn(cd8_perArea ~ Outcome + (1 | patient), data = plotdata)
p.foxp3 <- compare.fn(foxp3_perArea ~ Outcome + (1 | patient), data = plotdata)
# Assess _net_ cytotoxic activity (CD8/FOXP3), excluding samples with no FOXP3 cells
p.cd8foxp3 <- compare.fn(cd8_perArea/foxp3_perArea ~ Outcome + (1 | patient), data = plotdata[which(plotdata$foxp3_perArea > 0),])

# Include stromal comparisons
p.lymphocytess <- compare.fn(stroma_lymphocytes_perArea ~ Outcome + (1 | patient), data = plotdata)
p.cd4s <- compare.fn(stroma_cd4_perArea ~ Outcome + (1 | patient), data = plotdata)
p.cd8s <- compare.fn(stroma_cd8_perArea ~ Outcome + (1 | patient), data = plotdata)
p.foxp3s <- compare.fn(stroma_foxp3_perArea ~ Outcome + (1 | patient), data = plotdata)
p.cd8foxp3s <- compare.fn(stroma_cd8_perArea/stroma_foxp3_perArea ~ Outcome + (1 | patient), data = plotdata[which(plotdata$stroma_foxp3_perArea > 0),])

######################################################
# Alternate figure post-review
######################################################

pdata2 <- lapply(c("lymphocytes", "cd4", "cd8", "foxp3"), function(x) {
  ext <- '_perArea'
  y <- plotdata
  y$val <- y[,paste0(x, ext)]
  y$name <- paste0(toupper(x), " CIS")
  y$group <- 'CIS'
  
  y2 <- y
  y2$val <- y2[,paste0("stroma_", x, ext)]
  y2$name <- paste0(toupper(x), " stroma")
  y2$group <- 'Stroma'
  
  z <- rbind(y, y2)
  
  return(z)
})
pdata2 <- do.call('rbind', pdata2)
pdata2$name <- factor(pdata2$name, levels=c("LYMPHOCYTES CIS", "CD4 CIS", "CD8 CIS", "FOXP3 CIS", "LYMPHOCYTES stroma", "CD4 stroma", "CD8 stroma", "FOXP3 stroma"))

m <- max(pdata2$val, na.rm = T)
fig.a <- ggplot(pdata2, aes(x=name, y=val)) + 
  geom_boxplot(aes(fill=Outcome), outlier.size = 0.5) +
  geom_point(position = position_dodge(width=0.75), aes(group=Outcome), size=0.5) +
  ylab('Cell count per unit area') +
  theme(axis.title.x = element_blank()) +
  geom_vline(xintercept = 4.5, color="#CCCCCC") +
  ylim(c(0, m)) +
  annotate("text", x=2.5, y=m, label="Epithelium") +
  annotate("text", x=6.5, y=m, label="Stroma") +
  scale_x_discrete(labels=c("Lymphocytes","CD4", "CD8", "FOXP3", "Lymphocytes","CD4", "CD8", "FOXP3"))

# Add stars iff significant
fig.a <- ann_fun(fig.a, 1, p.lymphocytes$p, m=0.75)
fig.a <- ann_fun(fig.a, 2, p.cd4$p, m=0.75)
fig.a <- ann_fun(fig.a, 3, p.cd8$p, m=0.75)
fig.a <- ann_fun(fig.a, 4, p.cd4$p, m=0.75)
fig.a <- ann_fun(fig.a, 5, p.lymphocytess$p, m=0.75)
fig.a <- ann_fun(fig.a, 6, p.cd4s$p, m=0.75)
fig.a <- ann_fun(fig.a, 7, p.cd8s$p, m=0.75)
fig.a <- ann_fun(fig.a, 8, p.cd4s$p, m=0.75)

fig.b <- ggdraw() + draw_image("resources/fig1_progressive.jpg", scale=0.95) +
  draw_label('Progressive CIS lesion', y=0.13, size = 11)
fig.c <- ggdraw() + draw_image("resources/fig1_regressive.jpg", scale=0.95) +
  draw_label('Regressive CIS lesion', y=0.13, size = 11)
fig.bc <- plot_grid(fig.b, fig.c, labels=c('b','c'))

fig <- plot_grid(fig.a, fig.bc, nrow = 2, labels = c("a", ""))

save_plot(paste0(figdir, "fig1.tiff"), fig, nrow = 2, ncol=1, base_width = fig.width, base_height = 2.9)

# Save all comparisons in a text file
write(paste0("H&E lymphocyte data, n=", length(which(pheno$HandE)), ' (',length(which(pheno$HandE & pheno$Outcome == 'Progression')), 'P / ', length(which(pheno$HandE & pheno$Outcome == 'Regression')), 'R)'), file = opfile, append = T)
write(paste0("IHC data, n=", length(which(pheno$IHC_K)), ' (',length(which(pheno$IHC_K & pheno$Outcome == 'Progression')), 'P / ', length(which(pheno$IHC_K & pheno$Outcome == 'Regression')), 'R)'), file = opfile, append = T)

write(paste0("Lymphocytes P vs R p=", p.lymphocytes$p), file = opfile, append = T)
write(paste0("CD4 P vs R p=", p.cd4$p), file = opfile, append = T)
write(paste0("CD8 P vs R p=", p.cd8$p), file = opfile, append = T)
write(paste0("FOXP3 P vs R p=", p.foxp3$p), file = opfile, append = T)
write(paste0("CD8:FOXP3 P vs R p=", p.cd8foxp3$p), file = opfile, append = T)

# Include stromal comparisons here
write(paste0("Stromal Lymphocytes P vs R p=", p.lymphocytess$p), file = opfile, append = T)
write(paste0("Stromal CD4 P vs R p=", p.cd4s$p), file = opfile, append = T)
write(paste0("Stromal CD8 P vs R p=", p.cd8s$p), file = opfile, append = T)
write(paste0("Stromal FOXP3 P vs R p=", p.foxp3s$p), file = opfile, append = T)
write(paste0("Stromal CD8:FOXP3 P vs R p=", p.cd8foxp3s$p), file = opfile, append = T)

# Lastly include analyses of gradient (CIS concentration - stromal concentration)
p.cd4g <- compare.fn(cd4_perArea - stroma_cd4_perArea ~ Outcome + (1 | patient), data = plotdata)
p.cd8g <- compare.fn(cd8_perArea - stroma_cd8_perArea ~ Outcome + (1 | patient), data = plotdata)
p.foxp3g <- compare.fn(foxp3_perArea - stroma_foxp3_perArea ~ Outcome + (1 | patient), data = plotdata)
write(paste0("CIS - Stromal CD4 P vs R p=", p.cd4s$p), file = opfile, append = T)
write(paste0("CIS - Stromal CD8 P vs R p=", p.cd8s$p), file = opfile, append = T)
write(paste0("CIS - Stromal FOXP3 P vs R p=", p.foxp3s$p), file = opfile, append = T)

# Include legend text
do.log(c(
  sep, "# Figure 1 legend", sep,
  paste0("Combined quantitative immunohistochemistry data of CD4, CD8 and FOXP3 staining (n=",
         length(which(pheno$IHC_K)),"; ",length(which(pheno$IHC_K & pheno$Outcome == 'Progression'))," progressive, ",length(which(pheno$IHC_K & pheno$Outcome == 'Regression')),
         " regressive) with total lymphocyte quantification from H&E images (n=",
         length(which(pheno$HandE)),"; ",length(which(pheno$HandE & pheno$Outcome == 'Progression'))," progressive, ",length(which(pheno$HandE & pheno$Outcome == 'Regression')),
         " regressive) shown. We observe increased lymphocytes (p=",signif(p.lymphocytes$p,2),
         ") and CD8+ cells (p=",signif(p.cd8$p, 2),
         ") per unit area of epithelium within regressive CIS lesions compared to progressive"),
  sep, "# End figure 1 legend", sep, ""
))

####################################################################
# Figure S11 (added here as the same data are used):
# Late Progressors
####################################################################
# Look at the late progressors:
plotdata3 <- plotdata # Data for comparing groups
plotdata3$Outcome[which(plotdata3$exclude.reg == 'TRUE')] <- "LateProg"
p.lymphocytes <- compare.fn(lymphocytes_perArea ~ Outcome + (1 | patient), data = plotdata3)
p.cd4 <- compare.fn(cd4_perArea ~ Outcome + (1 | patient), data = plotdata3)
p.cd8 <- compare.fn(cd8_perArea ~ Outcome + (1 | patient), data = plotdata3)
p.foxp3 <- compare.fn(foxp3_perArea ~ Outcome + (1 | patient), data = plotdata3)

pdata3 <- pdata2 # Expanded data for ggplot plotting
pdata3$Outcome[which(pdata3$exclude.reg == 'TRUE')] <- "LateProg"
fig.lp <- ggplot(pdata3, aes(x=name, y=val)) + 
  geom_boxplot(aes(fill=Outcome), outlier.size = 0.5) +
  geom_point(position = position_dodge(width=0.75), aes(group=Outcome), size=0.5) +
  ylab('Cell count per unit area') +
  theme(axis.title.x = element_blank()) +
  geom_vline(xintercept = 4.5, color="#CCCCCC") +
  ylim(c(0, m)) +
  annotate("text", x=2.5, y=m, label="Epithelium") +
  annotate("text", x=6.5, y=m, label="Stroma") +
  scale_x_discrete(labels=c("Lymphocytes","CD4", "CD8", "FOXP3", "Lymphocytes","CD4", "CD8", "FOXP3"))

fig.lp <- ann_fun(fig.lp, 1, compare.fn(lymphocytes_perArea ~ Outcome + (1 | patient), data=plotdata3)$p, m=0.75)
fig.lp <- ann_fun(fig.lp, 2, compare.fn(cd4_perArea ~ Outcome + (1 | patient), data=plotdata3)$p, m=0.75)
fig.lp <- ann_fun(fig.lp, 3, compare.fn(cd8_perArea ~ Outcome + (1 | patient), data=plotdata3)$p, m=0.75)
fig.lp <- ann_fun(fig.lp, 4, compare.fn(foxp3_perArea ~ Outcome + (1 | patient), data=plotdata3)$p, m=0.75)
fig.lp <- ann_fun(fig.lp, 5, compare.fn(stroma_lymphocytes_perArea ~ Outcome + (1 | patient), data=plotdata3)$p, m=0.75)
fig.lp <- ann_fun(fig.lp, 6, compare.fn(stroma_cd4_perArea ~ Outcome + (1 | patient), data=plotdata3)$p, m=0.75)
fig.lp <- ann_fun(fig.lp, 7, compare.fn(stroma_cd8_perArea ~ Outcome + (1 | patient), data=plotdata3)$p, m=0.75)
fig.lp <- ann_fun(fig.lp, 8, compare.fn(stroma_foxp3_perArea ~ Outcome + (1 | patient), data=plotdata3)$p, m=0.75)

save_plot(paste0(figdir, 'figS11.tiff'), fig.lp, base_width = fig.width)

# Calculate the stats - here we need pair-wise comparisons...
lmm <- lmer(lymphocytes_perArea ~ Outcome + (1 | patient), data=plotdata3, REML = F)
anova(lmm) # 0.17
emmeans(lmm, list(pairwise ~ Outcome), adjust = "tukey") # NS
lmm <- lmer(cd8_perArea ~ Outcome + (1 | patient), data=plotdata3, REML = F)
anova(lmm) # 0.08
emmeans(lmm, list(pairwise ~ Outcome), adjust = "tukey") # NS (P-R 0.11)
lmm <- lmer(cd4_perArea ~ Outcome + (1 | patient), data=plotdata3, REML = F)
anova(lmm) # 0.3
emmeans(lmm, list(pairwise ~ Outcome), adjust = "tukey") # NS
lmm <- lmer(foxp3_perArea ~ Outcome + (1 | patient), data=plotdata3, REML = F)
anova(lmm) # 0.5
emmeans(lmm, list(pairwise ~ Outcome), adjust = "tukey") # NS
# Stroma:
lmm <- lmer(stroma_lymphocytes_perArea ~ Outcome + (1 | patient), data=plotdata3, REML = F)
anova(lmm) # 0.02 **
emmeans(lmm, list(pairwise ~ Outcome), adjust = "tukey") # Lateprog - Reg 0.05
lmm <- lmer(stroma_cd8_perArea ~ Outcome + (1 | patient), data=plotdata3, REML = F)
anova(lmm) # 0.16
emmeans(lmm, list(pairwise ~ Outcome), adjust = "tukey") # NS
lmm <- lmer(stroma_cd4_perArea ~ Outcome + (1 | patient), data=plotdata3, REML = F)
anova(lmm) # 0.4
emmeans(lmm, list(pairwise ~ Outcome), adjust = "tukey") # NS
lmm <- lmer(stroma_foxp3_perArea ~ Outcome + (1 | patient), data=plotdata3, REML = F)
anova(lmm) # 0.8
emmeans(lmm, list(pairwise ~ Outcome), adjust = "tukey") # NS

# Write ANOVA p-values to file:
for(x in c("lymphocytes", "cd4", "cd8", "foxp3")) {
  f <- as.formula(paste0(x, "_perArea ~ Outcome + (1 | patient)"))
  p <- compare.fn(f, data=plotdata3)$p
  write(paste0("ANOVA for P vs R vs late_prog ", x, " (CIS): p=", p), file = opfile, append = T)
  
  f <- as.formula(paste0("stroma_", x, "_perArea ~ Outcome + (1 | patient)"))
  p <- compare.fn(f, data=plotdata3)$p
  write(paste0("ANOVA for P vs R vs late_prog ", x, " (stroma): p=", p), file = opfile, append = T)
}

######################################################
# End
######################################################

##############################################################################
# Figure 2: Clustering analyses
# Cluster based on Danaher data and methylCIBERSORT
##############################################################################
message("Plotting Figure 2")

annot.colors <- list(outcome = c('Progression' = plotcols[1], 'Regression' = plotcols[2], 'Progressive' = plotcols[1], 'Regressive' = plotcols[2], 'LateProg' = plotcols[3]))
# Clustering on affy gxn immune genes:
annot.gxn <- data.frame(outcome = factor(gpheno.v$Outcone, levels=c('Regression', 'Progression')), row.names = gpheno.v$sampleID)

annot.mihc <- data.frame(outcome = factor(gsub('ive', 'ion', pheno.mihc$Outcome)), row.names = pheno.mihc$name)

# Remove CD45 as others are a subset
to.use <- c(1,2,3,4,7,9,10,11,12)
fig.a1 <- pheatmap(mihc.cis.norm2[], annotation_col = annot.mihc,  show_colnames = F, legend = F, annotation_legend = F, cutree_cols = 2, treeheight_row = 0, annotation_colors = annot.colors, fontsize_row = 6)


# Similar using deconvoluted Danaher values:
# Remove columns with empty data.
# Also remove "T-cells", macrophages and Treg as these are not corroborated by IHC (at least one of mIHC and Danaher is significant and results are in the opposite direction)
dan.to.rm <- c(2,11,15,16, 6,13)
fig.a <- pheatmap(t(gdata.danaher.t[,-dan.to.rm]), annotation_col = annot.gxn, scale='row', show_colnames = F, legend = F, annotation_legend = F, cutree_cols = 2, treeheight_row = 0, annotation_colors = annot.colors, fontsize_row = 6)


# Similar using methylCIBERSORT
# Exclude Controls - these won't have TILs as they are brushings
annot.meth <- data.frame(outcome = factor(mpheno$Sample_Group, levels=c('Regressive', 'Progressive')), row.names = mpheno$sampleID)
sel <- which(annot.meth$outcome %in% c("Progressive", "Regressive"))
plotdata <- t(methCS[sel, 1:11])
# Rename some columns
rownames(plotdata)[which(rownames(plotdata) == 'Eos')] <- 'Eosinophils'
rownames(plotdata)[which(rownames(plotdata) == 'CD4_Eff')] <- 'CD4'
rownames(plotdata)[which(rownames(plotdata) == 'Neu')] <- 'Neutrophils'
fig.b <- pheatmap(plotdata, annotation_col = annot.meth, show_colnames = F, legend = F, cutree_cols = 2, annotation_legend = F, treeheight_row = 0, annotation_colors = annot.colors, fontsize_row = 6)

# This shows a cluster with lots of cancer i.e. low immune infiltration

# Quick repeat including "Late Progressors"
annot.meth$outcome <- as.character(annot.meth$outcome)
annot.meth$outcome[which(rownames(annot.meth) %in% pheno$SampleID[which(pheno$exclude.reg == 'TRUE')])] <- 'LateProg'
pheatmap(plotdata, annotation_col = annot.meth, show_colnames = F, legend = F, cutree_cols = 2, annotation_legend = T, treeheight_row = 0, annotation_colors = annot.colors)

########################
# Danaher TIL scores prog vs Reg
# Colour points as 'hot' or 'cold' as per fig A
groups <- cutree(fig.a$tree_col, k = 2)
plotdata <- data.frame(
  til.score = gdata.danaher.t$til.score, 
  Outcome = gpheno.pair$Outcome, 
  patient = factor(gpheno.pair$Patient.Number),
  group = c('hot', 'cold')[groups[gpheno.pair$SampleID]],
  sampleID = gpheno.pair$SampleID
)
plotdata$Outcome <- gsub("ression", ".", plotdata$Outcome)
p.gxn <- compare.fn(til.score ~ Outcome + (1 | patient), data = plotdata)
fig.c <- ggplot(plotdata, aes(x=Outcome, y=til.score)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point(color=ifelse(plotdata$group == 'hot', hotcol, coldcol)) +
  stat_pvalue_manual(p.gxn, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  theme(axis.title.x = element_blank(), legend.direction = 'horizontal', legend.justification = 'center') +
  ylab('Danaher TIL score')
leg.pvr <- get_legend(fig.c)
fig.c <- fig.c + guides(fill=FALSE)

# methylCS scores:
sel <- which(mpheno$Sample_Group %in% c("Progressive", "Regressive"))
methCS.immune.cols <- c(2,3,4,5,6,8,10,11)
groups <- cutree(fig.b$tree_col, k = 2)
plotdata <- data.frame(
  pc.immune = rowSums(methCS[sel,methCS.immune.cols]), 
  Outcome = mpheno$Sample_Group[sel], 
  patient = factor(mpheno$Patient_ID[sel]),
  group = c('cold', 'hot')[groups[mpheno$sampleID[sel]]]
)
plotdata$Outcome <- gsub('ive$', 'ion', plotdata$Outcome)
plotdata$Outcome <- gsub("ression", ".", plotdata$Outcome)
p.methcs <- compare.fn(pc.immune ~ Outcome + (1 | patient), data = plotdata)
fig.d <- ggplot(plotdata, aes(x=Outcome, y=pc.immune)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point(color=ifelse(plotdata$group == 'hot', hotcol, coldcol)) +
  stat_pvalue_manual(p.methcs, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  guides(fill=FALSE) +
  theme(axis.title.x = element_blank()) +
  ylab('MethylCS % immune cells') 

# mIHC scores
# Use all lymphoid cells (T cells + B cells)
groups <- cutree(fig.a1$tree_col, k = 2)
plotdata <- data.frame(
  pc.lym = as.numeric(mihc.cis.norm2["All T cells",]) + as.numeric(mihc.cis.norm2["B cells",]), 
  Outcome = pheno.mihc$Outcome, 
  patient = factor(pheno.mihc$Patient),
  group = c('hot', 'cold')[groups[pheno.mihc$name]],
  sampleID = pheno.mihc$name
)
plotdata$Outcome <- gsub("ressive", ".", plotdata$Outcome)
p.mihc <- compare.fn(pc.lym ~ Outcome + (1 | patient), data = plotdata)
fig.e <- ggplot(plotdata, aes(x=Outcome, y=pc.lym)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point(color=ifelse(plotdata$group == 'hot', hotcol, coldcol)) +
  stat_pvalue_manual(p.mihc, label='p={p}', xmin = "group1", xmax="group2", tip.length = max(plotdata$pc.lym)/15, size=3) +
  guides(fill=FALSE) +
  theme(axis.title.x = element_blank()) +
  ylab('mIHC % lymphoid cells')

fig.main <- plot_grid(fig.e, fig.c, fig.d, as.grob(fig.a1), as.grob(fig.a), as.grob(fig.b), nrow = 2, labels=c('a','b','c','d','e','f'), rel_heights = c(3,5))

# Make a separate legend plot under the others
x <- ggplot(data.frame(x = c(1,2), y = c(1,2), Cluster = c('hot', 'cold')), aes(x=x, y=y, color=Cluster)) + geom_point() +
  scale_color_manual(values = c(coldcol, hotcol)) +
  theme(legend.direction = 'horizontal', legend.justification = 'center')
leg.clust <- get_legend(x)
plot(leg.clust)

fig <- plot_grid(
  fig.main, 
  plot_grid(leg.pvr, leg.clust, nrow = 1), 
  ncol = 1, rel_heights = c(15,1)
)

save_plot(paste0(figdir, "fig2.tiff"), fig, nrow = 2, ncol=2, base_width = fig.width/2)

do.log(c(
  sep, "# Figure 2 legend", sep,
  paste0(
    "Regressive lesions harbored significantly more infiltrating lymphocytes as assessed by multiplex immunohistochemistry (a; p=", p.mihc$p, "). ",
    "This finding was corroborated by molecular data in partially overlapping datasets; regressive lesions had higher gene-expression derived Tumor Infiltrating Lymphocyte (TIL) scores (b; p=", p.gxn$p, ") and ",
    "a higher proportion of immune cells as estimated from methylation data using methylCIBERSORT (c; p=", p.methcs$p, "). ",
    " d) Immune cell quantification from IHC data (n=",dim(pheno.mihc)[1],") shows an ‘immune cold’ cluster (left) in which all lesions progressed to cancer, and an ‘immune hot’ cluster (right) in which the majority regressed. Similar clustering patterns are seen in deconvoluted gene expression data (e; n=",dim(gdata.pair.t)[2],") and on methylation-derived cell subtypes using methylCIBERSORT (n=",length(which(mpheno$Sample_Group %in% c("Progressive", "Regressive"))), ")."
  ),
  sep, "# End figure 2 legend", sep, ""
))


##############################################################################
# Figure 3: Summary of genomic, epigenetic and transcriptomic changes
##############################################################################
message("Plotting Figure 3")

# Consider Sergio's gene lists:
names = list(
  'MHC Class I' = 'MHC I',
  'MHC Class II' = 'MHC II',
  'MHC Non Classical' = 'MHC nc',
  'Immunostimulator' = 'Stim',
  'Immunoinhibitor' = 'Inhib',
  'Antigen Processing' = 'Ag-Proc'
)
gene.groups = lapply(names(names), function(x) {
  gene.lists$somatic.immune.genes$GENES[which(gene.lists$somatic.immune.genes$REGULATION == x)]
})
names(gene.groups) <- as.character(names)
fig <- plot.muts(gene.lists$somatic.immune.genes$GENES, gene.groups = gene.groups)
fig
n = 1.6
save_plot(paste0(figdir, "fig3.tiff"), fig, base_width = fig.width*n*1.45, base_height = fig.width*n)


# Consider also an analysis of Thorson and Wellenstein genes (not included in final paper):
genes <- unique(c(gene.lists$wellenstein.genes, gene.lists$thorsson.drivers, gene.lists$thorsson.pos, gene.lists$thorsson.neg))
fig <- plot.muts(genes)
gene.groups$wellenstein <- gene.lists$wellenstein.genes
gene.groups$thorson.drivers <- gene.lists$thorsson.drivers
gene.groups$thorson.pos <- gene.lists$thorsson.pos
gene.groups$thorson.neg <- gene.lists$thorsson.neg

# Run Prog/Reg comparisons of each group using compareBurden
x <- lapply(1:length(gene.groups), function(i) {
  print(i)
  name = names(gene.groups)[i]
  # gs <- intersect(gene.lists$somatic.immune.genes$GENES, gene.groups[[i]])
  gs <- gene.groups[[i]]
  if(length(gs) == 0) {return(NA)}
  y <- NA
  result = tryCatch({
    y <- compareBurden(gs)
  }, error = function(e) {
    message(paste("FAILED", name))
    message(e)
  })
  y
})
names(x) <- names(gene.groups)
x[['All']] <- compareBurden(gene.lists$somatic.immune.genes$GENES)

# Output this readably:
dfToText <- function(df, n=6) {
  return(as.character(c(
    paste(str_pad(str_trunc(colnames(df),n, ellipsis = ''), n), collapse = '|'),
    apply(df, 1, function(x) {paste(str_pad(str_trunc(x, n, ellipsis = ''), n), collapse = '|')})
  )))
}
write('', file = opfile, append = T)
for(name in names(x)) {
  do.log(c(sep, paste0('*** Comparison of Prog vs Reg for genes classed as ', name), sep))
  comp <- x[[name]]
  if(is.na(comp)) {next}
  if(!is.na(comp$uncorrected)) {
    write(paste0('Direct muts+CNAs PvR comparsion p=', comp$uncorrected[,'Pr(>F)']), file = opfile, append = T)
  }
  if(!is.na(comp$muts.uncor)) {
    write(paste0('Uncorrected muts PvR comparsion p=', comp$muts.uncor['Outcome','Pr(>F)']), file = opfile, append = T)
  }
  if(!is.na(comp$cnas.uncor)) {
    write(paste0('Uncorrected CNAs PvR comparsion p=', comp$cnas.uncor['Outcome','Pr(>F)']), file = opfile, append = T)
  }
  
  if(!is.na(comp$both)) {
    write(paste0('Corrected muts+CNAs PvR comparsion p=', comp$both['Outcome','Pr(>F)']), file = opfile, append = T)
  }
  if(!is.na(comp$muts)) {
    write(paste0('Corrected muts PvR comparsion p=', comp$muts['Outcome','Pr(>F)']), file = opfile, append = T)
  }
  if(!is.na(comp$cnas)) {
    write(paste0('Corrected CNAs PvR comparsion p=', comp$cnas['Outcome','Pr(>F)']), file = opfile, append = T)
  }
  if(!is.na(comp$loh)) {
    write(paste0('Corrected LOH PvR comparsion p=', comp$loh['Outcome','Pr(>F)']), file = opfile, append = T)
  }
  if(!is.na(comp$loss)) {
    write(paste0('Corrected losses PvR comparsion p=', comp$loss['Outcome','Pr(>F)']), file = opfile, append = T)
  }
  if(!is.na(comp$gain)) {
    write(paste0('Corrected gains PvR comparsion p=', comp$gain['Outcome','Pr(>F)']), file = opfile, append = T)
  }
  
  if(!is.na(comp$dnds.global)) {
    write("Global dNdS:", file = opfile, append = T)
    write(dfToText(comp$dnds.global), file = opfile, append = T)
  }
  if(!is.na(comp$dnds)) {
    write("Gene-level dNdS:", file = opfile, append = T)
    write(dfToText(comp$dnds), file = opfile, append = T)
  }
  
  do.log(sep)
}


comp <- x$All
do.log(c(
  sep, "# Figure 3 legend", sep,
  paste0(
    "The mutational status is shown for ", length(gene.lists$somatic.immune.genes$GENES), " genes involved in the immune response, which are expressed by antigen presenting (tumor) cells. ",
    "Genes are categorised as belonging to the Major Histocompatibility Complex (MHC) class I or II; stimulators (Stim) and inhibitors (Inhib) of the immune response, and genes involved in antigen processing (Ag-Proc). ",
    "Mutations and copy number aberrations (CNAs) are shown for each of ", length(which(pheno$Whole.Genome.Sequencing & pheno$Outcome == 'Progression')), " progressive and ", length(which(pheno$Whole.Genome.Sequencing & pheno$Outcome == 'Regression')), " regressive samples. ",
    "Loss of heterozygosity (LOH) events are shown as mutations to avoid confusion with copy number loss, relative to ploidy. ",
    "The GXN PvR column displays the fold-change in expression of each gene between progressive and regressive samples, defined in a partially overlapping set of ", dim(gdata.pair.t)[2], " samples. Significant genes, defined as False Discovery Rate < 0.05, are highlighted in blue. ",
    "The TILcor column displays the Pearson's correlation coefficient between the expression of each gene and the gene-expression based Tumour Infiltrating Lymphocyte (TIL) score, derived by the Danaher method. ",
    "Progressive samples had more mutations (p=",signif(comp$muts.uncorrected[, "Pr(>F)"], 2),
    ") and CNAs (p=",signif(comp$cnas.uncorrected[, "Pr(>F)"], 2),") than regressive in this gene set. ",
    "dN/dS analysis identified ", paste(comp$dnds$gene_name, collapse = ', '), " as showing evidence of selection."
  ),
  sep, "# End figure 3 legend", sep, ""
))



##############################################################################
# Figure 4: Immune escape mechanisms in CIS beyond antigen presentation
##############################################################################
message("Plotting Figure 4")

# Volcano plot
uvv <- limmaCompare(gdata.pair.s, gpheno.pair, fdr_limit = 1)
fig.a <- ggplot(uvv, aes(x=logratio, y=-log(fdr))) +
  geom_point(size = 0.1) +
  ylim(0, -log(0.05)) +
  geom_hline(yintercept=-log(0.05), color='grey') +
  annotate("text", x=1, y=2.8, label="FDR < 0.05") +
  xlab('Log fold change') +
  ylab('-log(FDR)')

# PCA
gdata.pair.combined <- cbind(gdata.pair.t, gdata.pair.s)
p <- data.frame(
  Outcome = rep(gpheno.pair$Outcome, 2),
  Tissue = c(rep('Epithelium', dim(gpheno.pair)[1]), rep('Stroma', dim(gpheno.pair)[1])),
  sampleID = rep(gpheno.pair$SampleID, 2)
)
pca <- prcomp(t(gdata.pair.combined))


plotdata <- p
plotdata$PC1 <- pca$x[,1]
plotdata$PC2 <- pca$x[,2]
fig.b <- ggplot(plotdata, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=Outcome, shape=Tissue)) +
  scale_shape_manual(values = c(4,20))

fig.ab <- plot_grid(fig.a, fig.b, ncol=2, rel_widths = c(1,2), labels = c('a','b'))

# Immunomodulatory signals:
# TNFSF9/CCL27/PDL1
plotdata <- gpheno.pair
plotdata$CD274 <- as.numeric(gdata.pair.t["CD274",])
plotdata$PDCD1 <- as.numeric(gdata.pair.t["PDCD1",])
plotdata$`CD274:PDCD1` <- plotdata$CD274 / plotdata$PDCD1

plotdata$TNFSF9 <- as.numeric(gdata.pair.t["TNFSF9",])
plotdata$TNFRSF9 <- as.numeric(gdata.pair.t["TNFRSF9",])
plotdata$`TNFSF9:TNFRSF9` <- plotdata$TNFSF9 / plotdata$TNFRSF9

plotdata$CCL27 <- as.numeric(gdata.pair.t["CCL27",])
plotdata$CCR10 <- as.numeric(gdata.pair.t["CCR10",])
plotdata$`CCL27:CCR10` <- plotdata$CCL27 / plotdata$CCR10

l <- lapply(c('CD274', 'PDCD1', 'CD274:PDCD1', 'TNFSF9', 'TNFRSF9', 'TNFSF9:TNFRSF9', 'CCL27', 'CCR10', 'CCL27:CCR10'), function(x) {
  df <- data.frame(
    name = x,
    val = plotdata[,x],
    Outcome = plotdata$Outcome,
    patient = plotdata$Patient.Number
  )
  # Include a p-value for each:
  p <- compare.fn(val ~ Outcome + (1 | patient), data = df)$p
  write(paste0("Prog vs Reg p-value for ", x, ": ", p), file = opfile, append = T)
  
  return(df)
})
plotdata2 <- do.call('rbind', l)

plotdata.a <- plotdata2[-grep(':', plotdata2$name, fixed = T),] 
plotdata.a$name <- factor(plotdata.a$name, levels = c("CD274", "PDCD1", "TNFSF9", "TNFRSF9", "CCL27", "CCR10"))
a <- ggplot(plotdata.a, aes(x=name, y=val)) +
  geom_boxplot(aes(fill=Outcome), outlier.size = 0.5) +
  geom_point(position = position_dodge(width=0.75), aes(group=Outcome), size=0.5) +
  guides(fill=F) +
  ylab('Gene Expression') +
  theme(axis.title.x = element_blank())
p.cd274 <- compare.fn(val ~ Outcome + (1 | patient), data = plotdata.a[which(plotdata.a$name == 'CD274'),])$p
a <- ann_fun(a, 1, p.cd274)
p.pdcd1 <- compare.fn(val ~ Outcome + (1 | patient), data = plotdata.a[which(plotdata.a$name == 'PDCD1'),])$p
a <- ann_fun(a, 2, p.pdcd1)
p.tnfsf9 <- compare.fn(val ~ Outcome + (1 | patient), data = plotdata.a[which(plotdata.a$name == 'TNFSF9'),])$p
a <- ann_fun(a, 3, p.tnfsf9)
p.tnfrsf9 <- compare.fn(val ~ Outcome + (1 | patient), data = plotdata.a[which(plotdata.a$name == 'TNFRSF9'),])$p
a <- ann_fun(a, 4, p.tnfrsf9)
p.ccl27 <- compare.fn(val ~ Outcome + (1 | patient), data = plotdata.a[which(plotdata.a$name == 'CCL27'),])$p
a <- ann_fun(a, 5, p.ccl27)
p.ccr10 <- compare.fn(val ~ Outcome + (1 | patient), data = plotdata.a[which(plotdata.a$name == 'CCR10'),])$p
a <- ann_fun(a, 6, p.ccr10)

plotdata.b <- plotdata2[grep(':', plotdata2$name, fixed = T),]
plotdata.b$name <- factor(plotdata.b$name, levels = c("CD274:PDCD1", "TNFSF9:TNFRSF9", "CCL27:CCR10"))
b <- ggplot(plotdata.b, aes(x=name, y=val)) +
  geom_boxplot(aes(fill=Outcome), outlier.size = 0.5) +
  geom_point(position = position_dodge(width=0.75), aes(group=Outcome), size=0.5) +
  guides(fill=F) +
  ylab('Gene Expression Ratio') +
  theme(axis.title.x = element_blank())
p.cd274_pdcd1 <- compare.fn(val ~ Outcome + (1 | patient), data = plotdata.b[which(plotdata.b$name == 'CD274:PDCD1'),])$p
b <- ann_fun(b, 1, p.cd274_pdcd1)
p.tnfsf9_tnfrsf9 <-compare.fn(val ~ Outcome + (1 | patient), data = plotdata.b[which(plotdata.b$name == 'TNFSF9:TNFRSF9'),])$p 
b <- ann_fun(b, 2, p.tnfsf9_tnfrsf9)
p.ccl27_ccr10 <- compare.fn(val ~ Outcome + (1 | patient), data = plotdata.b[which(plotdata.b$name == 'CCL27:CCR10'),])$p
b <- ann_fun(b, 3, p.ccl27_ccr10)

fig.cd <- plot_grid(a,b,rel_widths = c(2,2), labels = c('c','d'))


#### Staining data

# mIHC for PDL1:
# PDL1 from mIHC
plotdata <- pheno.mihc
plotdata$PDL1.pc <- (mihc.cis.norm["CD45-PanCK+Ki67+PDL1+",] + mihc.cis.norm["CD45-PanCK+Ki67-PDL1+",]) / mihc.cis.norm["CD45-PanCK+",]
pdl1.plot <- ggplot(plotdata, aes(x=Outcome, y = 100 * PDL1.pc)) +
  geom_boxplot(aes(fill=Outcome), outlier.size = 0.5) +
  geom_point(position = position_dodge(width=0.75), aes(group=Outcome), size=0.5) +
  guides(fill=F) +
  ylab('% PD-L1 positive cells') +
  theme(axis.title.x = element_blank())

# plot_grid(pdl1.plot, a,b,rel_widths = c(0.8, 2,1.2), labels = c('g','h'), nrow = 1)
load('data/staining_for_plots.RData')
mixdf.tnf <- mixdf.tnf[which(mixdf.tnf$mu1 > 0.15),]
mixdf.ccl27 <- mixdf.ccl27[which(mixdf.ccl27$mu1 > 0.05),]
mixdf.ccr10 <- mixdf.ccr10[which(mixdf.ccr10$mu1 > 0.2),]
# Combine to get a ratio:
samps <- intersect(mixdf.ccl27$sampleID, mixdf.ccr10$sampleID)
samps <- pheno[which(pheno$SampleID %in% samps),]
samps$Patient <- factor(samps$Patient.Number)
samps$ccl27.pospc <- mixdf.ccl27$pc.pos[match(samps$SampleID, mixdf.ccl27$sampleID)]
samps$ccr10.pospc <- mixdf.ccr10$pc.pos[match(samps$SampleID, mixdf.ccr10$sampleID)]
plotdata <- rbind(
  data.frame(
    name = "PD-L1",
    val = (mihc.cis.norm["CD45-PanCK+Ki67+PDL1+",] + mihc.cis.norm["CD45-PanCK+Ki67-PDL1+",]) / mihc.cis.norm["CD45-PanCK+",],
    Outcome = pheno.mihc$Outcome,
    patient = pheno.mihc$Patient
  ),
  data.frame(
    name = 'TNFSF9',
    val = mixdf.tnf$pc.pos,
    Outcome = mixdf.tnf$Outcome,
    patient = mixdf.tnf$Patient
  ),
  data.frame(
    name = 'CCL27',
    val = mixdf.ccl27$pc.pos,
    Outcome = mixdf.ccl27$Outcome,
    patient = mixdf.ccl27$Patient
  ),
  data.frame(
    name = 'CCR10',
    val = mixdf.ccr10$pc.pos,
    Outcome = mixdf.ccr10$Outcome,
    patient = mixdf.ccr10$Patient
  )
)
plotdata$name <- factor(plotdata$name, levels = c('PD-L1', 'TNFSF9', 'CCL27', 'CCR10'))
plots.s <- ggplot(plotdata, aes(x=name, y=val)) +
  geom_boxplot(aes(fill=Outcome), outlier.size = 0.5) +
  geom_point(position = position_dodge(width=0.75), aes(group=Outcome), size=0.5) +
  guides(fill=F) +
  ylab('% positive cells by IHC') +
  theme(axis.title.x = element_blank())
p.pdl1_prot <- compare.fn(val ~ Outcome + (1 | patient), data = plotdata[which(plotdata$name == 'PD-L1'),])$p
plots.s <- ann_fun(plots.s, 1, p.pdl1_prot)
p.tnfsf9_prot <- wilcox.test(val ~ Outcome, data = plotdata[which(plotdata$name == 'TNFSF9'),])$p.value
plots.s <- ann_fun(plots.s, 2, p.tnfsf9_prot)
p.ccl27_prot <- wilcox.test(val ~ Outcome, data = plotdata[which(plotdata$name == 'CCL27'),])$p.value
plots.s <- ann_fun(plots.s, 3, p.ccl27_prot)
p.ccr10_prot <- wilcox.test(val ~ Outcome, data = plotdata[which(plotdata$name == 'CCR10'),])$p.value
plots.s <- ann_fun(plots.s, 4, p.ccr10_prot)
plots.s


plot.cclccrstain <- ggplot(samps, aes(x="CCL27:CCR10", y=ccl27.pospc / ccr10.pospc)) +
  geom_boxplot(aes(fill=Outcome), outlier.size = 0.5) +
  geom_point(position = position_dodge(width=0.75), aes(group=Outcome), size=0.5) +
  guides(fill=F) +
  ylab('Protein expression ratio') +
  theme(axis.title.x = element_blank())
p.ccl27_ccr10_prot <- compare.fn(ccl27.pospc / ccr10.pospc ~ Outcome + (1 | Patient), data = samps)$p
# plot.cclccrstain <- ann_fun(plot.cclccrstain, 1, p.ccl27_ccr10_prot)


# Final plot
# fig <- plot_grid(fig.ab, fig.cd, fig.ef, fig.gh, ncol=1)
fig <- plot_grid(
  fig.ab, 
  fig.cd,
  plot_grid(
    plots.s,
    ggdraw() + draw_image("data/ihc_demo_image.tif") + theme(plot.margin = margin(15,5,15,15)),
    plot.cclccrstain,
    nrow = 1,
    rel_widths = c(2.1, 2.1, 1),
    labels = c("e","f", "g")
  ),
  ncol = 1
)

save_plot(paste0(figdir, "fig4.tiff"), fig, nrow = 4, ncol=4, base_width = fig.width/4, base_height = fig.width/3)


do.log(c(
  sep, "# Figure 4 legend", sep,
  paste0(
    "(a) Volcano plot of gene expression differential analysis of laser-captured stroma comparing progressive (n=",length(which(gpheno.pair$Outcome == 'Progression')),") and regressive (n=",length(which(gpheno.pair$Outcome == 'Regression')),") CIS samples. No genes were significant with FDR < 0.05 following adjustment for multiple testing. (b) Principle Component Analysis plot of the same ",dim(gpheno.pair)[1]," CIS samples, showing laser-captured epithelium and matched stroma. ",
    "(c-f) RNA analysis of immunomodulatory molecules and cytokine:receptor pairs in n=", dim(gdata.pair.t)[2]," CIS samples identified TNFSF9 and CCL27:CCR10 as significantly differentially expressed between progressive and regressive samples (p=", signif(p.tnfsf9, 2), " and p=", signif(p.ccl27_ccr10, 2), " respectively). ",
    "Immunohistochemistry showed that TNFSF9 was similarly differentially expressed at the protein level (p=",signif(p.tnfsf9_prot, 2),"; n=",dim(mixdf.tnf)[1]," with successful staining). ",
    "CCL27 and CCR10 showed a similar trend at the protein level to the RNA level; whilst these data did not achieve a significance threshold (p=", signif(p.ccl27_ccr10_prot), " for CCL27:CCR10 ratio, n=",dim(samps)[1],") we observe several outliers in the progressive group. ",
    "Analysis of PD-L1 (encoded by CD274) and its receptor PD-1 (encoded by PDCD1) is included due to its relevance in clinical practice; again we do not achieve statistically significant results but do observe three marked outliers with PD-L1 expression >25%, all of which progressed to cancer. ",
    "(g) Illustrative immunohistochemistry staining for TNFSF9. ",
    "All p-values are calculated using linear mixed effects modeling to account for samples from the same patient; ***p < 0.001 **p < 0.01 *p<0.05 #p<0.1. Units for gene expression figures represent normalised microarray intensity values."
  ),
  sep, "# End figure 4 legend", sep, ""
))


###############################################################################################
# End Main Figures
###############################################################################################




##############################################################################
# Figure S1
# Tile figure of what was done to each sample
##############################################################################
message("Plotting Figure S1")

df <- rbind(
  data.frame(SampleID = pheno$SampleID, mod = "WGS", val = pheno$Whole.Genome.Sequencing, patient = pheno$Patient.Number),
  data.frame(SampleID = pheno$SampleID, mod = "Methylation", val = pheno$Methylation, patient = pheno$Patient.Number),
  data.frame(SampleID = pheno$SampleID, mod = "Gene expression", val = pheno$Stroma.GXN, patient = pheno$Patient.Number),
  data.frame(SampleID = pheno$SampleID, mod = "IHC", val = pheno$IHC_K, patient = pheno$Patient.Number),
  data.frame(SampleID = pheno$SampleID, mod = "Image analysis", val = pheno$HandE, patient = pheno$Patient.Number)
)
df$outcome <- as.character(pheno$Outcome[match(df$SampleID, pheno$SampleID)])
df$val[which(df$val == FALSE)] <- NA
df$val[which(!is.na(df$val))] <- df$outcome[which(!is.na(df$val))]
df$val <- factor(df$val, levels=c('Progression', 'Regression', 'Control'))
df$nmods <- unlist(lapply(df$SampleID, function(x) {
  length(which(df$SampleID == x & !is.na(df$val)))
}))
df <- df[order(df$outcome, df$patient, df$nmods, decreasing = F),]
# df$outcome <- factor(df$outcome)
df$SampleID <- paste0(" Pt", df$patient, "-", as.character(df$SampleID))
df$SampleID <- factor(as.character(df$SampleID), levels=unique(as.character(df$SampleID)))

fig <- ggplot(df[which(df$outcome %in% c("Progression", "Regression")),], aes(x=mod, y=SampleID, fill=val)) +
  geom_tile(aes(width = 0.9, height = 0.8), size=2) +
  scale_fill_manual(values=plotcols[c(1,2)]) +
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 6), axis.ticks.x = element_blank(), legend.title = element_blank()) +
  ylab("Sample")
fig

save_plot(paste0(figdir, "figS1.tiff"), fig, base_width = fig.width, base_height = fig.maxheight)


##############################################################################
# Figure S2
# Pro-anti inflammatory cytokines
##############################################################################
message("Plotting Figure S2")

# Let's look at the geometric mean of proinflammatory vs anti-inflammatory cytokines
geomean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

x <- gpheno.pair
x$pro.mean.t <- apply(gdata.pair.t[intersect(rownames(gdata.pair.t), gene.lists$cytokines.pro),], 2, geomean)
x$anti.mean.t <- apply(gdata.pair.t[intersect(rownames(gdata.pair.t), gene.lists$cytokines.anti),], 2, geomean)

x$pro.mean.s <- apply(gdata.pair.s[intersect(rownames(gdata.pair.s), gene.lists$cytokines.pro),], 2, geomean)
x$anti.mean.s <- apply(gdata.pair.s[intersect(rownames(gdata.pair.s), gene.lists$cytokines.anti),], 2, geomean)

x$Outcome <- gsub('ression', '.', x$Outcome)
p <- compare.fn(pro.mean.t ~ Outcome + (1 | Patient.Number), data = x)
f1 <- ggplot(x, aes(x=Outcome, y=pro.mean.t)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  ylab("Pro-inflammatory cytokines (epithelium)") +
  theme(axis.title.x = element_blank())

p <- compare.fn(anti.mean.t ~ Outcome + (1 | Patient.Number), data = x)
f2 <- ggplot(x, aes(x=Outcome, y=anti.mean.t)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  ylab("Anti-inflammatory cytokines (epithelium)") +
  theme(axis.title.x = element_blank())

p <- compare.fn(pro.mean.t / anti.mean.t ~ Outcome + (1 | Patient.Number), data = x)
f3 <- ggplot(x, aes(x=Outcome, y=pro.mean.t / anti.mean.t)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  ylab("Pro:Anti ratio (epithelium)") +
  theme(axis.title.x = element_blank())
# Evidence of proinflammatory phenotype in regression

# Demonstrate it's not the case in stroma:
p <- compare.fn(pro.mean.s ~ Outcome + (1 | Patient.Number), data = x)
f4 <- ggplot(x, aes(x=Outcome, y=pro.mean.s)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  ylab("Pro-inflammatory cytokines (stroma)") +
  theme(axis.title.x = element_blank())

p <- compare.fn(anti.mean.s ~ Outcome + (1 | Patient.Number), data = x)
f5 <- ggplot(x, aes(x=Outcome, y=anti.mean.s)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  ylab("Anti-inflammatory cytokines (stroma)") +
  theme(axis.title.x = element_blank())
p <- compare.fn(pro.mean.s / anti.mean.s ~ Outcome + (1 | Patient.Number), data = x)
f6 <- ggplot(x, aes(x=Outcome, y=pro.mean.s / anti.mean.s)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  ylab("Pro:Anti ratio (stroma)") +
  theme(axis.title.x = element_blank())

fig <- plot_grid(f1, f2, f3, f4, f5, f6, ncol=3, labels=c('a','b','c','d','e','f'))

save_plot(paste0(figdir, "figS2.tiff"), fig, base_width = fig.width, base_height = fig.width)

# Print also the p-value for CXCL8:
x$cxcl8 <- as.numeric(gdata.pair.t["CXCL8",])
p <- compare.fn(cxcl8 ~ Outcome + (1 | Patient.Number), data = x)
write(paste0("p-value for CXCL8 (uncorrected): ", p$p), file = opfile, append = T)

##############################################################################
# Figure S3
# Pro-anti inflammatory cytokines (individual plots)
##############################################################################
message("Plotting Figure S3")

genes <- intersect(rownames(gdata.pair.t), c(gene.lists$cytokines.anti, gene.lists$cytokines.pro))
uvv <- limmaCompare(gdata.pair.t, gpheno.pair, fdr_limit = 1)
ps <- lapply(genes, function(x) {
  df <- gpheno.pair
  df$val = as.numeric(gdata.pair.t[x,])
  df$Outcome <- gsub('ression', '.', df$Outcome)
  
  p <- compare.fn(val ~ Outcome + (1 | Patient.Number), data = df)
  return(p)
})
names(ps) <- genes
ps.adj <- p.adjust(sapply(ps, function(x) {x$p}))
for(i in 1:length(ps)) {
  ps[[i]]$p <- signif(ps.adj[i], 2)
}

plots <- lapply(genes, function(x) {
  df <- gpheno.pair
  df$val = as.numeric(gdata.pair.t[x,])
  df$Outcome <- gsub('ression', '.', df$Outcome)
  
  p <- ps[[x]]
  
  fig <- ggplot(df, aes(x=Outcome, y=val)) +
    geom_boxplot(aes(fill=Outcome)) +
    geom_point() +
    guides(fill=F) +
    stat_pvalue_manual(p, label='FDR={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3, vjust = 1.7, label.size = 2.5) +
    ylab(x) +
    theme(axis.title.x = element_blank())
  return(fig)
})
names(plots) <- genes
# Plot pro- and anti-inflammatory cytokines separately:
fig.a <- plot_grid(
  plots[["IFNG"]], plots[["TNF"]], plots[["IL1B"]], plots[["IL2"]], plots[["IL7"]], plots[["CXCL8"]], plots[["IL12A"]], plots[["IL17A"]], plots[["IL23A"]]
)
fig.b <- plot_grid(
  plots[["TGFB1"]], plots[["IL10"]], plots[["IL1RAP"]], plots[["IL4"]], plots[["IL6"]], plots[["IL11"]], plots[["IL13"]]
)
fig <- plot_grid(fig.a, fig.b, labels = c('a','b'), nrow = 1)

save_plot(paste0(figdir, "figS3.tiff"), fig, base_width = fig.width*1.2, base_height = fig.width)


##############################################################################
# Figure S4 
# Smoking status and its effect on progression
##############################################################################

message("Plotting Figure S4")

library(DescTools)

pheno$quit_group <- pheno$SmokingStatus
pheno$quit_group[which(pheno$quit_group == 'Former' & is.na(pheno$YearsSinceQuitting))] <- 'Unknown'
pheno$quit_group[which(pheno$YearsSinceQuitting <= 1)] <- '<1'
pheno$quit_group[which(pheno$YearsSinceQuitting > 1 & pheno$YearsSinceQuitting <= 5)] <- '1-5'
pheno$quit_group[which(pheno$YearsSinceQuitting > 5 & pheno$YearsSinceQuitting <= 10)] <- '5-10'
pheno$quit_group[which(pheno$YearsSinceQuitting > 10)] <- '>10'

pheno$quit_group <- factor(pheno$quit_group, levels=c('Current', '<1', '1-5', '5-10', '>10', 'Never', 'Unknown'))

pheno$smoking_dose <- as.numeric(pheno$quit_group)

# First use all samples to show an increasing trend of progression from never -> former -> current
p <- pheno[which(!is.na(pheno$SmokingStatus) & pheno$Outcome != 'Control' & pheno$SmokingStatus != 'Unknown'),]
# Check infiltration status:
p$Infiltrated <- p$lymphocytes_perArea > median(p$lymphocytes_perArea, na.rm=T)

# We will use a boostrap method to account for multiple samples per patient
# Test multiple combinations of samples from each patient
# Collect the Fisher p-value for each and summarise
# The number of combinations is too large to be calculated, so sample multiple times
set.seed(42)
ps <- lapply(1:1000, function(i) {
  # Pick a sample at random from each patient:
  samps <- sapply(unique(p$Patient.Number), function(pid) {
    sample(p$SampleID[which(p$Patient.Number == pid)], 1)
  })
  # Run three comparions:
  # Smoking status vs outcome
  # Time since quitting vs outcome (in former smoker group)
  # Smoking status vs high/low lymphocytes
  p.tmp <- p[match(samps, p$SampleID),]
  t1 <- table(p.tmp$Outcome, p.tmp$SmokingStatus)
  sel.former <- which(p.tmp$SmokingStatus == 'Former' & !is.na(p.tmp$YearsSinceQuitting))
  t2 <- table(p.tmp$Outcome[sel.former], p.tmp$quit_group[sel.former])
  p.tmp$has_lym <- p.tmp$lymphocytes_perArea > median(p.tmp$lymphocytes_perArea, na.rm=T)
  t3 <- table(p.tmp$has_lym, p.tmp$SmokingStatus)
  return(data.frame(
    run=i,
    p1.chisq = prop.trend.test(as.numeric(t1['Progression',]), colSums(t1))$p.value,
    p1.CA2sided = CochranArmitageTest(t1, alternative = 'two.sided')$p.value, # Identical to chisq trend
    p1.CAincreasing = CochranArmitageTest(t1, alternative = 'increasing')$p.value, # Measures increasing regression,
    
    p2.chisq = prop.trend.test(as.numeric(t2['Progression',]), colSums(t2))$p.value,
    p2.CA2sided = CochranArmitageTest(t2, alternative = 'two.sided')$p.value, # Identical to chisq trend
    p2.CAincreasing = CochranArmitageTest(t2, alternative = 'increasing')$p.value, # Measures increasing regression
    
    p3.chisq = prop.trend.test(as.numeric(t3['TRUE',]), colSums(t3))$p.value,
    p3.CA2sided = CochranArmitageTest(t3, alternative = 'two.sided')$p.value, # Identical to chisq trend
    p3.CAincreasing = CochranArmitageTest(t3, alternative = 'increasing')$p.value # Measures increasing regression
  ))
})
ps <- do.call('rbind', ps)


# CA tests:
t <- table(p$Outcome, p$SmokingStatus)
p.ca <- CochranArmitageTest(t, alternative = 'increasing')$p.value
t <- table(p$Infiltrated, p$SmokingStatus)
p.ca2 <- CochranArmitageTest(t, alternative = 'increasing')$p.value
print(paste('CA test for all samples for progression Current > Former > Never p=', signif(p.ca,2)))
print(paste('CA test for all samples for infiltration Current > Former > Never p=', signif(p.ca2,2)))
print(paste('Boostrapped CA test for increasing progression in Current > Former > Never smokers: mean p=', signif(mean(ps$p1.CAincreasing),2), 'CI', signif(mean(ps$p1.CAincreasing) - 2*sd(ps$p1.CAincreasing),2), '-', signif(mean(ps$p1.CAincreasing) + 2*sd(ps$p1.CAincreasing), 2)))
print(paste('Boostrapped CA test for increasing progression in former smokers with years since quitting: mean p=', signif(mean(ps$p2.CAincreasing),2), 'CI', signif(mean(ps$p2.CAincreasing) - 2*sd(ps$p2.CAincreasing),2), '-', signif(mean(ps$p2.CAincreasing) + 2*sd(ps$p2.CAincreasing), 2)))
print(paste('Boostrapped CA test for increasing lymphocyte infiltration in Never > Former > Current smokers: mean p=', signif(mean(ps$p3.CAincreasing),2), 'CI', signif(mean(ps$p3.CAincreasing) - 2*sd(ps$p3.CAincreasing),2), '-', signif(mean(ps$p3.CAincreasing) + 2*sd(ps$p3.CAincreasing), 2)))




# Bar plot of progression with smoking status:
mytheme <- theme(axis.title = element_text(size = 9))
# Absolute:
a <- ggplot(p, aes(x=SmokingStatus, fill=Outcome)) + 
  geom_bar() + 
  mytheme +
  theme(axis.title.x = element_blank(), legend.direction = 'horizontal', legend.justification = 'center') 
leg <- get_legend(a)
a <- a + guides(fill=F)
# Relative:
p.rel <- p %>% 
  group_by(SmokingStatus, Outcome) %>% 
  summarise(count=n()) %>% 
  mutate(proportion=count/sum(count))
b <- ggplot(p.rel, aes(x = factor(SmokingStatus), y = proportion, fill = factor(Outcome))) +
  geom_bar(stat="identity", width = 0.7) +
  mytheme +
  theme(axis.title.x = element_blank()) +
  guides(fill=F)

# Repeat with lymphocytes:
sel <- which(!is.na(p$Infiltrated))
c <- ggplot(p[sel,], aes(x=SmokingStatus, fill=Infiltrated)) + 
  geom_bar() +
  mytheme +
  theme(axis.title.x = element_blank(), legend.direction = 'horizontal', legend.justification = 'center') +
  scale_fill_manual(values = c("#1C63A3", "#D70017"))
leg2 <- get_legend(c)
c <- c + guides(fill=F)
# Relative:
p.rel2 <- p[sel,] %>% 
  group_by(SmokingStatus, Infiltrated) %>% 
  summarise(count=n()) %>% 
  mutate(proportion=count/sum(count))
d <- ggplot(p.rel2, aes(x = factor(SmokingStatus), y = proportion, fill = factor(Infiltrated))) +
  geom_bar(stat="identity", width = 0.7) +
  mytheme +
  theme(axis.title.x = element_blank()) +
  guides(fill=F) +
  scale_fill_manual(values = c("#1C63A3", "#D70017"))

# Plot a histogram of p-values
e <- ggplot(ps, aes(x=p1.CAincreasing)) +
  geom_histogram(aes(y=..density..), bins = 15, alpha = 0.5) +
  geom_density(bw = sd(ps$p1.CAincreasing)/2) +
  mytheme +
  xlab('p-value for outcome')

f <- ggplot(ps, aes(x=p3.CAincreasing)) +
  geom_histogram(aes(y=..density..), bins = 15, alpha = 0.5) +
  geom_density(bw = sd(ps$p3.CAincreasing)/2) + 
  mytheme +
  xlab('p-value for infiltration')

fig.ab <- plot_grid(a,b,ncol=1, labels = c('a','b'))
fig.cd <- plot_grid(c,d,ncol=1, labels = c('d','e'))

fig <- plot_grid(
  fig.ab, e,
  fig.cd, f,
  ncol=2, nrow=2,
  rel_widths = c(2,1),
  labels = c('', 'c','', 'f')
)

# Add breakdown of ex-smokers
g <- ggplot(p, aes(x=quit_group, fill=Outcome)) + 
  geom_bar() +
  mytheme +
  theme(axis.title.x = element_blank()) + guides(fill=F)
p.rel <- p %>% 
  group_by(quit_group, Outcome) %>% 
  summarise(count=n()) %>% 
  mutate(proportion=count/sum(count))
h <- ggplot(p.rel, aes(x = factor(quit_group), y = proportion, fill = factor(Outcome))) +
  geom_bar(stat="identity", width = 0.7) +
  mytheme +
  theme(axis.title.x = element_blank()) +
  guides(fill=F)

i <- ggplot(p[which(!is.na(p$Infiltrated)),], aes(x=quit_group, fill=Infiltrated)) + 
  geom_bar() +
  mytheme +
  theme(axis.title.x = element_blank()) + guides(fill=F) +
  scale_fill_manual(values = c("#1C63A3", "#D70017"))
p.rel <- p[which(!is.na(p$Infiltrated)),] %>% 
  group_by(quit_group, Infiltrated) %>% 
  summarise(count=n()) %>% 
  mutate(proportion=count/sum(count))
j <- ggplot(p.rel, aes(x = factor(quit_group), y = proportion, fill = factor(Infiltrated))) +
  geom_bar(stat="identity", width = 0.7) +
  mytheme +
  theme(axis.title.x = element_blank()) +
  guides(fill=F) +
  scale_fill_manual(values = c("#1C63A3", "#D70017"))


fig <- plot_grid(
  fig,
  g, h, i,j,
  plot_grid(leg, leg2, nrow = 1),
  ncol = 1,
  rel_heights = c(4,1,1,1,1, 0.3),
  labels = c('', 'g','h', 'i', 'j')
)
fig

save_plot(paste0(figdir, "figS4.tiff"), fig, base_width = fig.width, base_height = fig.maxheight)

do.log(c(
  sep, "# Figure S4 legend", sep,
  paste0("Smoking status was available for ",length(which(pheno$SmokingStatus %in% c("Current", "Former", "Never"))),
         " CIS lesions from ",length(unique(pheno$Patient.Number[which(pheno$SmokingStatus %in% c("Current", "Former", "Never"))])),
         " patients (",length(which(pheno$SmokingStatus %in% c("Current"))),
         " lesions from ",length(unique(pheno$Patient.Number[which(pheno$SmokingStatus %in% c("Current"))])),
         " current smokers; ",length(which(pheno$SmokingStatus %in% c("Former"))),
         " from ",length(unique(pheno$Patient.Number[which(pheno$SmokingStatus %in% c("Former"))])),
         " former smokers; ",length(which(pheno$SmokingStatus %in% c("Never"))),
         " from ",length(unique(pheno$Patient.Number[which(pheno$SmokingStatus %in% c("Never"))])),
         " never smokers). Here we show the absolute (a) and relative (b) numbers of lesions in each group which progressed to cancer. Using a Cochrane-Armitage test to look for a trend from Current to Former to Never smokers, we found a trend towards lower chance of regression (p=",signif(p.ca, 1),
         "). To account for samples from the same patient, we randomly selected one sample per patient and repeated this Cochrane-Armitage test, obtaining a p-value < 0.1 irrespective of sample choice (c). Similarly, we show the absolute (d) and relative (e) numbers of lesions which are “infiltrated”, defined as having above the median number of infiltrating lymphocytes per unit area. Again, we observe a trend towards higher infiltration from Current to Former to Never smokers, which is reasonably robust to sample selection (f)."),
  sep, "# End figure S4 legend", sep, ""
))


##############################################################################
# Figure S5
# Neoantigen plots:
#  a) Predicted neoantigens vs mutational burden
#  b) Clonal neoantigens P vs R
#  c) Proportion clonal P vs R
#  d) Proportion clonal corrected for purity
#  e) Affinity P vs R
#  f) Rank P vs R
#  g) Depletion stats P vs R
##############################################################################
message("Plotting Figure S5")

c <- cor.test(pheno$neoants.strong, pheno$burden)
r2 <- paste0("r2=", signif(c$estimate, 2), "; p=", signif(c$p.value,2))
fig.a <- ggplot(pheno[which(pheno$Whole.Genome.Sequencing),], aes(x=burden, y=neoants.strong)) +
  geom_point(aes(color=Outcome)) +
  xlab("Mutational Burden") + ylab('Strong Neoantigens') +
  guides(color=F) +
  stat_cor()
# geom_text(x=20000,y=400,label=paste0(r2), color='black')

# Calculate any additional measures
pheno$neoants.strong.clonal <- unlist(lapply(pheno$Sample.Number..WGS., function(x) {
  length(which(muts.all$is.neo.strong & muts.all$is.clonal & muts.all$patient == x))
}))
source('private/neoant_depletion_analysis.R')
pheno$depletion.stat <- depletion_stat[pheno$Sample.Number..WGS.]

# Format for plots
plotdata <- pheno[which(pheno$Whole.Genome.Sequencing),]
plotdata$Outcome <- gsub("ression", ".", plotdata$Outcome)
plotdata$patient <- factor(plotdata$Patient.Number)

p <- compare.fn(neoants.strong ~ Outcome + (1 | patient), data = plotdata)
fig.b <- ggplot(plotdata, aes(x=Outcome, y=neoants.strong)) +
  geom_boxplot(aes(fill=Outcome)) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  # geom_signif(comparisons = list(c("Progression", "Regression"))) +
  geom_point() +
  guides(fill=F) +
  ylab("Strong Neoantigens") +
  theme(axis.title.x = element_blank())

p <- compare.fn(neoants.strong.clonal ~ Outcome + (1 | patient), data = plotdata)
fig.c <- ggplot(plotdata, aes(x=Outcome, y=neoants.strong.clonal)) +
  geom_boxplot(aes(fill=Outcome)) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  geom_point() +
  guides(fill=F) +
  ylab("Clonal Strong Neoantigens") +
  theme(axis.title.x = element_blank())

p <- compare.fn(neoants.strong.clonal/neoants.strong ~ Outcome + (1 | patient), data = plotdata)
fig.d <- ggplot(plotdata, aes(x=Outcome, y=(neoants.strong.clonal / neoants.strong))) +
  geom_boxplot(aes(fill=Outcome)) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  geom_point() +
  guides(fill=F) +
  ylab("Proportion Clonal Strong Neoantigens") +
  theme(axis.title.x = element_blank())

# Comparing neoantigen properties P vs R
neos.all <- lapply(1:length(neoantigen.data), function(i) {
  df <- neoantigen.data[[i]]
  df$sample <- names(neoantigen.data)[i]
  df <- df[which(df$is_strong),c("sample", "Peptide", "best.aff", "best.rank", "best.DAI")]
})
neos.all <- do.call('rbind', neos.all)
neos.all$Outcome <- plotdata$Outcome[match(neos.all$sample, plotdata$Sample.Number..WGS.)]
neos.all$patient <- factor(pheno$Patient.Number[match(neos.all$sample, pheno$Sample.Number..WGS.)])

p <- compare.fn(log(best.aff) ~ Outcome + (1 | patient), data = neos.all)
fig.e <- ggplot(neos.all, aes(x=Outcome, y=log(best.aff))) +
  geom_boxplot(aes(fill=Outcome)) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  # geom_signif(comparisons = list(c("Progression", "Regression"))) +
  geom_point() +
  guides(fill=F) +
  ylab("Binding affinity (log)") +
  theme(axis.title.x = element_blank())

p <- compare.fn(best.rank ~ Outcome + (1 | patient), data = neos.all)
fig.f <- ggplot(neos.all, aes(x=Outcome, y=best.rank)) +
  geom_boxplot(aes(fill=Outcome)) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  geom_point() +
  guides(fill=F) +
  ylab("Rank") +
  theme(axis.title.x = element_blank())

p <- compare.fn(best.DAI ~ Outcome + (1 | patient), data = neos.all)
fig.g <- ggplot(neos.all, aes(x=Outcome, y=best.DAI)) +
  geom_boxplot(aes(fill=Outcome)) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  geom_point() +
  guides(fill=F) +
  ylab("DAI") +
  theme(axis.title.x = element_blank())

# Plot depletion stats:
p <- compare.fn(depletion.stat ~ Outcome + (1 | patient), data = plotdata)
fig.h <- ggplot(plotdata, aes(x=Outcome, y=depletion.stat)) +
  geom_boxplot(aes(fill=Outcome)) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  geom_point() +
  guides(fill=F) +
  ylab("Depletion Score") +
  theme(axis.title.x = element_blank())


fig <- plot_grid(fig.a, fig.b, fig.c, fig.d, fig.e, fig.f, fig.g, fig.h, nrow=3, ncol=3, 
                 labels = c('a', 'b','c','d','e','f','g','h'))

save_plot(paste0(figdir, "figS5.tiff"), fig, ncol=3, nrow=3, base_width = fig.width/3, base_height = fig.width/2.5)

##############################################################################
# Figure S6
# Methylation DMRs across the genome
##############################################################################
message("Plotting Figure S6")

# Plot DMRs as a circos plot.
# Skip this plot if the user does not have circos installed.
circos.cmd <- "/Users/adam/Software/circos-0.69-6/bin/circos"
if(file.exists(circos.cmd)){
  
  circos.dir <- "results/circos_tmp/"
  dir.create(circos.dir, recursive = T, showWarnings = F)
  
  conf <- paste(circos.dir, "circos.conf", sep="")
  meth.circos <- paste(circos.dir, "meth.circosdata.txt", sep="")
  meth.circos.tcga <- paste(circos.dir, "meth.circosdata.tcga.txt", sep="")
  circos.labels <- paste(circos.dir, "circos.labels.txt", sep="")
  
  # Label relevant genes. Use biomaRt to find locations
  ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org")
  # Get all genes
  bm <- getBM(
    attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
    filters=c('hgnc_symbol'),
    values=gene.lists$hla.assoc, 
    mart = ensembl
  )
  hla.track <- bm[which(bm$chromosome_name %in% 1:22),]
  hla.track <- data.frame(
    chr=paste0('hs', hla.track$chromosome_name),
    start=hla.track$start_position, end=hla.track$end_position,
    gene=hla.track$hgnc_symbol
  )
  write.table(hla.track, sep="\t", quote=F, col.names=F, row.names = F, file=circos.labels)
  
  
  meth.track <- data.frame(
    chr=gsub("chr", "hs", dmrs$seqnames),
    start=dmrs$start,
    end=dmrs$end,
    value=dmrs$betaAv_Progressive - dmrs$betaAv_Regressive
  )
  write.table(meth.track, sep="\t", quote=F, col.names=F, row.names = F, file=meth.circos)
  
  meth.track.tcga <- data.frame(
    chr=gsub("chr", "hs", dmrs.tcga$seqnames),
    start=dmrs.tcga$start,
    end=dmrs.tcga$end,
    value=dmrs.tcga$betaAv_TCGA.SqCC - dmrs.tcga$betaAv_TCGA.Control
  )
  write.table(meth.track.tcga, sep="\t", quote=F, col.names=F, row.names = F, file=meth.circos.tcga)
  
  # Create the circos config file
  text <- paste("
                karyotype = data/karyotype/karyotype.human.txt
                
                <ideogram>
                
                <spacing>
                default = 0.005r
                </spacing>
                
                radius    = 0.95r
                thickness = 2p
                fill      = no
                color     = dgrey
                
                show_label       = yes
                # see etc/fonts.conf for list of font names
                label_font       = default 
                label_radius     = 0.35r
                label_size       = 30
                label_parallel   = yes
                
                </ideogram>
                
                
                <plots>
                # Labels for drivers
                <plot>
                type = text
                color            = black
                file             = ",circos.labels,"
                r0 = 0.85r
                r1 = 1.5r
                
                show_links     = yes
                link_dims      = 4p,4p,8p,4p,4p
                link_thickness = 2p
                link_color     = red
                
                label_size   = 36p
                label_font   = condensed
                
                padding  = 0p
                rpadding = 0p
                
                label_snuggle = yes
                
                </plot>
                
                
                # CIS methylation DMRs
                <plot>
                type = scatter
                file = ",meth.circos,"
                r1   = 0.85r
                r0   = 0.65r
                min = -0.5
                max = 0.5
                
                <axes>
                show=yes
                <axis>
                color=vdgrey
                position=0.5r
                </axis>
                </axes>
                <backgrounds>
                <background>
                color=vvlyellow
                y0=0.5r
                </background>
                <background>
                color=vvlblue
                y1=0.5r
                </background>
                </backgrounds>
                </plot>
                
                # TCGA methylation DMRs
                <plot>
                type = scatter
                file = ",meth.circos.tcga,"
                r1   = 0.6r
                r0   = 0.4r
                min = -0.5
                max = 0.5
                fill_color = yellow
                thickness = 5
                extend_bin=no
                <axes>
                show=yes
                <axis>
                color=vdgrey
                position=0.5r
                </axis>
                </axes>
                <backgrounds>
                <background>
                color=vvlyellow
                y0=0.5r
                </background>
                <background>
                color=vvlblue
                y1=0.5r
                </background>
                </backgrounds>
                </plot>
                
                </plots>
                ################################################################
                # The remaining content is standard and required. It is imported 
                # from default files in the Circos distribution.
                #
                # These should be present in every Circos configuration file and
                # overridden as required. To see the content of these files, 
                # look in etc/ in the Circos distribution.
                
                <image>
                # Included from Circos distribution.
                <<include etc/image.conf>>
                </image>
                
                # RGB/HSV color definitions, color lists, location of fonts, fill patterns.
                # Included from Circos distribution.
                <<include etc/colors_fonts_patterns.conf>>
                
                # Debugging, I/O an dother system parameters
                # Included from Circos distribution.
                # This file need to be edited to set max_points_per_track to 250000 (for the TCGA CNA track)
                <<include etc/housekeeping.conf>>
                ",sep="")
  
  write(text, file=conf)
  
  wd <- getwd() 
  setwd(circos.dir)
  system(circos.cmd)
  setwd(wd)
  # file.copy(paste0(circos.dir, 'circos.png'), paste0(figdir, "figS4a.png"), overwrite = T)
  file.copy(paste0(circos.dir, 'circos.svg'), paste0(figdir, "circos.tmp.svg"), overwrite = T)
  
  fig.a <- ggdraw() + draw_image(paste0(figdir, "circos.tmp.svg"))
  
  # Plot significant HLA-A probes in a facet_wrap:
  hla.a.probes <- c("cg09803951", "cg05157171", "cg23489273")
  # For each patient, need each of 3 probes, Sample_Group, ID
  plotdata <- mpheno
  plotdata <- cbind(plotdata, t(mdata[hla.a.probes,make.names(plotdata$sampleID)]))
  plotdata$Sample_Group <- factor(plotdata$Sample_Group, levels=c('Progressive', 'Regressive', 'Control'))
  plotdata <- plotdata[order(-as.numeric(plotdata$Sample_Group)),]
  plotdata$index <- 1:dim(plotdata)[1]
  
  
  plotdata2 <- rbind(
    cbind(plotdata, data.frame(probe='cg09803951', val=plotdata$cg09803951)),
    cbind(plotdata, data.frame(probe='cg05157171', val=plotdata$cg05157171)),
    cbind(plotdata, data.frame(probe='cg23489273', val=plotdata$cg23489273))
  )
  
  fig.b <- ggplot(plotdata2, aes(x=index, y=val, color=Sample_Group)) +
    geom_point() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'bottom', legend.justification = 'center', legend.title = element_blank()) +
    facet_wrap(~probe) +
    ylab("Methylation beta value") +
    scale_color_manual(values = plotcols)
  
  fig <- plot_grid(fig.a, fig.b, labels=c('a','b'), ncol=1, nrow=2, rel_heights = c(2,1))
  save_plot(paste0(figdir, "figS6.tiff"), fig, nrow = 2, ncol=1, base_height = 5, base_width=fig.width)
  file.remove(paste0(figdir, "circos.tmp.svg"))
} else {
  message("Warning: circos installation not found. Figure S6 will not be plotted.")
}


##############################################################################
# Figure S7
# HLA silencing in CIS and TCGA (plots of GXN vs methylation)
##############################################################################
message("Plotting Figure S7")

gmcors.tcga <- data.frame(
  gene = gene.lists$all.immune,
  r2 = NA,
  p = NA
)
gmcors.tcga$gene <- as.character(gmcors.tcga$gene)
for(i in 1:dim(gmcors.tcga)[1]) {
  gene <- gmcors.tcga$gene[i]
  if(!(gene %in% rownames(gm.tcga.gdata)) | !(gene %in% rownames(gm.tcga.mdata.genes))) {
    next
  }
  c <- cor.test(as.numeric(log(gm.tcga.gdata[gene, ])), as.numeric(gm.tcga.mdata.genes[gene, ]))
  gmcors.tcga$r2[i] <- c$estimate
  gmcors.tcga$p[i] <- c$p.value
}
gmcors.tcga$fdr <- p.adjust(gmcors.tcga$p)
gmcors.tcga <- gmcors.tcga[order(gmcors.tcga$p),]




tcga.plots <- list()
cis.plots <- list()
ol <- pheno[which(pheno$Methylation & pheno$Gene.expression),]
genes.to.plot <- c("HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "B2M", "WNT5A", "TNFSF9", "CIITA")
for(gene in genes.to.plot) {
  plotdata <- data.frame(
    gene.expression = log(as.numeric(gm.tcga.gdata[gene,]) + 0.5),
    methylation = as.numeric(gm.tcga.mdata.genes[gene,])
  )
  tcga.plots[[gene]] <- ggplot(plotdata, aes(x=methylation, y=gene.expression)) +
    geom_point(size = 0.1, color='orange') +
    geom_smooth(method='lm') +
    stat_cor() +
    xlab("Methylation beta value") +
    ylab("Gene expression (log counts)") +
    ggtitle(paste("TCGA", gene)) +
    theme(axis.title = element_text(size = 8))
  
  
  plotdata <- data.frame(
    gene.expression = as.numeric(gdata[gene,ol$SampleID]),
    methylation = as.numeric(mdata.genes.nosnp[gene, ol$SampleID]),
    lym_grad = ol$lym_grad
  )
  cis.plots[[gene]] <- ggplot(plotdata, aes(x=methylation, y=gene.expression)) +
    geom_point(size = 0.4, color='darkorange') +
    geom_smooth(method='lm') +
    stat_cor() +
    xlab("Methylation beta value") +
    ylab("Gene expression") +
    ggtitle(paste("CIS", gene)) +
    theme(axis.title = element_text(size = 8))
}
fig <- plot_grid(
  tcga.plots[["HLA-A"]], tcga.plots[["HLA-B"]], tcga.plots[["HLA-C"]],
  tcga.plots[["TAP1"]], tcga.plots[["TAP2"]], tcga.plots[["B2M"]], 
   tcga.plots[["TNFSF9"]], tcga.plots[["CIITA"]],
  cis.plots[["HLA-A"]], cis.plots[["HLA-B"]], cis.plots[["HLA-C"]],
  cis.plots[["TAP1"]], cis.plots[["TAP2"]], cis.plots[["B2M"]],
   # cis.plots[["TNFSF9"]], 
  cis.plots[["CIITA"]],
  ncol=4,
  labels = c("a", "", "", "", "", "", "", "", "b")
  # labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
)
save_plot(paste0(figdir, "figS7.tiff"), fig, nrow = 4, ncol=4, base_width = fig.width/2.75, base_height = fig.width/2.75)

##############################################################################
# Figure S8
# Methylation patterns over above genes
##############################################################################
message("Plotting Figure S8")

fig.main <- plot_grid(
  plotMethylForGene("HLA-A") + guides(color=F) + scale_color_manual(values = plotcols), 
  plotMethylForGene("HLA-B") + guides(color=F) + scale_color_manual(values = plotcols), 
  plotMethylForGene("HLA-C") + guides(color=F) + scale_color_manual(values = plotcols), 
  plotMethylForGene("TAP1") + guides(color=F) + scale_color_manual(values = plotcols), 
  plotMethylForGene("B2M") + guides(color=F) + scale_color_manual(values = plotcols), 
  plotMethylForGene("TNFSF9") + guides(color=F) + scale_color_manual(values = plotcols), 
  nrow=3
)
# Add a legend:
x <- ggplot(mpheno, aes(x=Sample_Group, y=Sample_Group, color=Sample_Group)) +
  geom_point() +
  theme(legend.position = 'bottom', legend.justification = 'center', legend.title = element_blank())
x <- cowplot::get_legend(x)

fig <- plot_grid(fig.main, x, ncol=1, rel_heights = c(9.5,0.5))
save_plot(paste0(figdir, "figS8.tiff"), fig, nrow = 3, ncol=2, base_width = fig.width/2)

##############################################################################
# Figure S9
# Boxplots of immune checkpoints with CIS tissue
##############################################################################
message("Plotting Figure S9")

l <- lapply(gene.lists$checkpoints, function(gene) {
  if(gene %in% rownames(gdata.v)) {
    return(data.frame(gene=gene, val=as.numeric(gdata.v[gene,]), Outcome=gpheno.v$Outcone, index=1:dim(gdata.v)[2]))
  }
})
plotdata <- do.call('rbind', l)

fig <- ggplot(plotdata, aes(x=index, y=val, color=Outcome)) + 
  geom_point() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = 'bottom', legend.justification = 'center') +
  facet_wrap(~gene, scales = 'free') +
  ylab('Gene expression')

save_plot(paste0(figdir, "figS9.tiff"), fig, base_width = fig.width, base_height = fig.width*0.85)


##############################################################################
# Figure S10
# Cytokine:receptor analysis
##############################################################################
message("Plotting Figure S10")

ck.ratio <- gdata.pair.t[gene.lists$chemokines$name,] / gdata.pair.t[gene.lists$chemokines$receptor,]
rownames(ck.ratio) <- paste(gene.lists$chemokines$name, gene.lists$chemokines$receptor, sep=".")
ck.ratio.diff <- limmaCompare(ck.ratio, gpheno.pair, fdr_limit = 0.05)

write(paste0("PvR Limma analysis for chemokine:receptor pairs: sig pairs ", paste(paste(rownames(ck.ratio.diff), "FC:", ck.ratio.diff$fc, "FDR:", ck.ratio.diff$fdr, sep=','), sep = '; ')), file = opfile, append = T)

# We see an interesting signal in CCL27:CCR10 - increased signalling relevant in melanoma through AKT/PI3K signalling
# (Ref https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3926818/)
plotdata <- gpheno.pair
plotdata$ccl27.ccr10 <- as.numeric(ck.ratio["CCL27.CCR10",])
plotdata$ccl27 <- as.numeric(gdata.pair.t["CCL27",])
plotdata$ccr10 <- as.numeric(gdata.pair.t["CCR10",])
plotdata$Outcome <- gsub('ression', '.', plotdata$Outcome)
p <- compare.fn(ccl27.ccr10 ~ Outcome + (1 | Patient.Number), data = plotdata)
f1 <- ggplot(plotdata, aes(x=Outcome, y=ccl27.ccr10)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  ylab("CCL27:CCR10 ratio") +
  guides(fill=F) +
  theme(axis.title.x = element_blank())

p <- compare.fn(ccl27 ~ Outcome + (1 | Patient.Number), data = plotdata)
f2 <- ggplot(plotdata, aes(x=Outcome, y=ccl27)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  guides(fill=F) +
  ylab('CCL27 expression') +
  theme(axis.title.x = element_blank())

p <- compare.fn(ccr10 ~ Outcome + (1 | Patient.Number), data = plotdata)
f3 <- ggplot(plotdata, aes(x=Outcome, y=ccr10)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  guides(fill=F) +
  ylab('CCR10 expression') +
  theme(axis.title.x = element_blank())

# This is supposed to act through PIK3CA/AKT - any correlation?
plotdata$pik3ca <- as.numeric(gdata.pair.t["PIK3CA",])
f4 <- ggplot(plotdata, aes(x=ccl27, y=pik3ca)) +
  geom_point() + geom_smooth(method='lm') + stat_cor() +
  ylab('PIK3CA expression') +
  xlab('CCL27 expression')
plotdata$akt1 <- as.numeric(gdata.pair.t["AKT1",])
f5 <- ggplot(plotdata, aes(x=ccl27, y=akt1)) +
  geom_point() + geom_smooth(method='lm') + stat_cor() +
  ylab('AKT1 expression')+
  xlab('CCL27 expression')

# fig <- plot_grid(f1,f2,f3,f4,f5, ncol = 3, labels = c('a','b','c','d','e'))
fig <- plot_grid(f4,f5, ncol = 2, labels = c('a','b'))
save_plot(paste0(figdir, "figS10.tiff"), fig, base_width = fig.width, base_height = fig.width/2)


##############################################################################
# Sup. Data 1
# Table of all lesions and profiling modalities
##############################################################################
message("Plotting Table S1")

pheno.public <- pheno[,c(1,2,3,5,8,10,12,13,17,20,25,26,29,30,31,32,33,34, 6,9,11,14)]
pheno.public$Sample.Number..Meth.[which(pheno.public$Sample.Number..Meth. == 'X')] <- ''
pheno.public$late.progression <- pheno$exclude.reg
# Sort properly by sample number:
pheno.public$tmp <- str_pad(str_extract(pheno.public$SampleID, "[0-9]+"), width = 3, pad = '0')
pheno.public <- pheno.public[order(pheno.public$tmp),]
pheno.public$tmp <- NULL
pheno.public$Previous.lung.cancer <- as.character(pheno.public$Previous.lung.cancer)
WriteXLS('pheno.public', ExcelFileName = paste0(supdir, "TableS1.xlsx"))

##############################################################################
# Sup. Data 2
# Markers used for mIHC
# Created manually
##############################################################################

##############################################################################
# Sup. Data 3
# Danaher and MethylCIBERSORT decomposition data, with p-values for P vs R
##############################################################################
message("Plotting Table S3")

# Remove columns that aren't supported by IHC:
dan.to.rm <- c(2,11,15, 6,13)

df1 <- gdata.danaher.t[,-dan.to.rm]
df1$patient <- factor(gpheno.pair$Patient.Number)
df1$Outcome <- gpheno.pair$Outcome
df1 <- df1[,which(!apply(df1, 2, function(x) {all(is.na(x))}))]
colnames(df1) <- make.names(colnames(df1))
cols <- colnames(df1)
p.PvR <- sapply(cols, function(x){
  if(x %in% c('patient', 'Outcome')) {return(NA)}
  f <- as.formula(paste0(make.names(x), " ~ Outcome + (1 | patient)"))
  return(compare.fn(f, data = df1)$p)
})
FDR.PvR <- p.adjust(p.PvR)
names(p.PvR) <- cols
mean.p <- t(data.frame(apply(df1[which(df1$Outcome == 'Progression'),], 2, function(x) {mean(as.numeric(x))})))
rownames(mean.p) <- "mean.p"
df1 <- rbind(df1, mean.p)
mean.r <- t(data.frame(apply(df1[which(df1$Outcome == 'Regression'),], 2, function(x) {mean(as.numeric(x))})))
rownames(mean.r) <- "mean.r"
df1 <- rbind(df1, mean.r)
df1 <- rbind(df1, t(data.frame(p.PvR)))
df1 <- rbind(df1, t(data.frame(FDR.PvR)))

df1s <- gdata.danaher.s[,-dan.to.rm]
df1s$patient <- factor(gpheno.pair$Patient.Number)
df1s$Outcome <- gpheno.pair$Outcome
df1s <- df1s[,which(!apply(df1s, 2, function(x) {all(is.na(x))}))]
colnames(df1s) <- make.names(colnames(df1s))
cols <- colnames(df1s)
p.PvR <- sapply(cols, function(x){
  if(x %in% c('patient', 'Outcome')) {return(NA)}
  f <- as.formula(paste0(make.names(x), " ~ Outcome + (1 | patient)"))
  return(compare.fn(f, data = df1s)$p)
})
FDR.PvR <- p.adjust(p.PvR)
names(p.PvR) <- cols
mean.p <- t(data.frame(apply(df1s[which(df1s$Outcome == 'Progression'),], 2, function(x) {mean(as.numeric(x))})))
rownames(mean.p) <- "mean.p"
df1s <- rbind(df1s, mean.p)
mean.r <- t(data.frame(apply(df1s[which(df1s$Outcome == 'Regression'),], 2, function(x) {mean(as.numeric(x))})))
rownames(mean.r) <- "mean.r"
df1s <- rbind(df1s, mean.r)
df1s <- rbind(df1s, t(data.frame(p.PvR)))
df1s <- rbind(df1s, t(data.frame(FDR.PvR)))
df1s <- cbind(data.frame(sampleID = rownames(df1s)), df1s)

df2 <- data.frame(methCS[,1:11])
df2$PC.immune <- rowSums(methCS[,methCS.immune.cols])
df2$patient <- factor(mpheno$Patient_ID)
df2$Outcome <- factor(mpheno$Sample_Group)
cols <- colnames(df2)
p.PvR <- sapply(cols, function(x){
  if(x %in% c('patient', 'Outcome')) {return(NA)}
  f <- as.formula(paste0(make.names(x), " ~ Outcome + (1 | patient)"))
  # Exclude control samples here:
  return(compare.fn(f, data = df2[which(df2$Outcome %in% c("Progressive", "Regressive")),])$p)
})
FDR.PvR <- p.adjust(p.PvR)
names(p.PvR) <- cols
mean.p <- t(data.frame(apply(df2[which(df2$Outcome == 'Progressive'),], 2, function(x) {mean(as.numeric(x))})))
rownames(mean.p) <- "mean.p"
df2 <- rbind(df2, mean.p)
mean.r <- t(data.frame(apply(df2[which(df2$Outcome == 'Regressive'),], 2, function(x) {mean(as.numeric(x))})))
rownames(mean.r) <- "mean.r"
df2 <- rbind(df2, mean.r)
df2 <- rbind(df2, t(data.frame(p.PvR)))
df2 <- rbind(df2, t(data.frame(FDR.PvR)))

df1 <- cbind(data.frame(SampleID = rownames(df1)), df1)
df2 <- cbind(data.frame(SampleID = rownames(df2)), df2)


###
# Add mIHC data
df3 <- mihc.cis.norm2
df3['prop.cd8.grzb',] <- as.numeric(mihc.cis.data['CD45+CD3+CD8+Grzb+',]) / as.numeric(mihc.cis.data['CD45+CD3+CD8+',])
df3['prop.cd8.eomes',] <- as.numeric(mihc.cis.data['CD45+CD3+CD8+EOMES+PD1+',]) + as.numeric(mihc.cis.data['CD45+CD3+CD8+EOMES+PD1-',]) / as.numeric(mihc.cis.data['CD45+CD3+CD8+',])

df3$p <- sapply(rownames(df3), function(x) {
  p <- pheno.mihc
  p$val <- as.numeric(df3[x,])
  p2 <- compare.fn(val ~ Outcome + (1 | Patient), data = p)
  p2$p
})
df3 <- t(data.frame(df3))
df3 <- cbind(data.frame(sampleID = rownames(df3)), df3)

# Repeat for stroma
df3.s <- mihc.stroma.norm2
df3.s['prop.cd8.grzb',] <- as.numeric(mihc.stroma.data['CD45+CD3+CD8+Grzb+',]) / as.numeric(mihc.stroma.data['CD45+CD3+CD8+',])
df3.s['prop.cd8.eomes',] <- as.numeric(mihc.stroma.data['CD45+CD3+CD8+EOMES+PD1+',]) + as.numeric(mihc.stroma.data['CD45+CD3+CD8+EOMES+PD1-',]) / as.numeric(mihc.stroma.data['CD45+CD3+CD8+',])
df3.s$p <- sapply(rownames(df3.s), function(x) {
  p <- pheno.mihc
  p$val <- as.numeric(df3.s[x,])
  p2 <- compare.fn(val ~ Outcome + (1 | Patient), data = p)
  p2$p
})
df3.s <- t(data.frame(df3.s))
df3.s <- cbind(data.frame(sampleID = rownames(df3.s)), df3.s)

p <- pheno.mihc
# p$cd8.pc <- as.numeric(mihc.cis.norm2['CD8 T cells',])
p$macro.pc <- as.numeric(mihc.cis.norm2['Macrophages',])
p$macro <- as.numeric(mihc.cis.data['CD45+CD3-CD20-CD11b+CD66b-CD163+CD68+',]) + as.numeric(mihc.cis.data['CD45+CD3-CD20-CD11b+CD66b-CD163-CD68+',])
ggplot(p, aes(x=Outcome, y=macro)) + geom_boxplot() + geom_point()
ggplot(p[which(p$macro > 5),], aes(x=Outcome, y=macro.pc)) + geom_boxplot() + geom_point()
# What if we require a reasonable number of counts to compare?
# Now we see an almost-significant difference with more in prog...
ggplot(p, aes(x=Outcome, y=macro.pc)) + geom_boxplot() + geom_point()
compare.fn(macro.pc ~ Outcome + (1 | Patient), data = p[which(p$macro > 5),])

WriteXLS(list(df1, df2, df3, df1s, df3.s), row.names = F,
         SheetNames = c("CIS Danaher GXN", "CIS MethylCIBERSORT", 'CIS mIHC', 'Stroma Danaher GXN', 'Stroma mIHC'), ExcelFileName = paste0(supdir, "TableS3.xlsx"))

##############################################################################
# Sup. Data 4
# Lists of genes used in this analysis
# Include general immune genes (used for stroma PvsR comparison) and immunomodulators
##############################################################################
message("Plotting Table S4")

df1 <- data.frame(Gene.Name=sort(gene.lists$all.immune))
df2 <- gene.lists$somatic.immune.genes[,c('GENES', 'REGULATION')]
df3 <- data.frame(Gene.Name=sort(gene.lists$checkpoints))
df4 <- rbind(
  data.frame(Gene.Name=sort(gene.lists$cytokines.pro), Type='Pro-inflammatory'),
  data.frame(Gene.Name=sort(gene.lists$cytokines.anti), Type='Anti-inflammatory')
)
df5 <- gene.lists$chemokines
WriteXLS(
  list(df1, df2, df3, df4, df5), 
  ExcelFileName = paste0(supdir, "TableS4.xlsx"), 
  SheetNames = c('All immune genes', 'Antigen presentation genes', 'Immunomodulators', 'Inflammatory cytokines', 'Chemokine-receptor pairs'), col.names = F
)


###########################################################################
# Sup. Table 5
# Illumina vs Affymetrix analysis
# Prog/Reg pathways previously published were based on Illumina microarrays. Immune signals were not highly significant.
# We ask whether this was influenced by microarray design
###########################################################################
message("Plotting Table S5")
# We find two key differences between Illumina and Affymetrix microarrays:
#  1. Affymetrix has many more probes, covering a wider range of transcripts
#  2. Affymetrix often covers multiple transcripts per gene, where Illumina has only one
#
# We propose the following approach to address the impact of these differences on immune profiling:
#   * Reduce the Affymetrix dataset to genes covered by Illumina probes
#   * Discard genes covered by multiple transcripts, which are ambiguous
#   * Also discard these genes from the Illumina data
#   * We hypothesise that the two present similar pathways when analysed P vs R, and these do not include immune pathways

# Reduce affy data to Illumina probes:
gdata.v.r <- gdata.v[which(rownames(gdata.v) %in% rownames(gdata.d)),]

# Find ambiguous genes (multiple transcripts in Affy) and discard from _both_ datasets
ambig.genes <- rownames(gdata.v)[which(duplicated(rownames(gdata.v)))]
gdata.v.r2 <- gdata.v.r[which(!(rownames(gdata.v.r) %in% ambig.genes)),]
gdata.d.r <- gdata.d[which(!(rownames(gdata.d) %in% ambig.genes)),]

# Illumina data does include some late progressors, so exclude these
lateprog <- pheno$SampleID[which(pheno$exclude.reg == 'TRUE')]
gpheno.d2 <- gpheno.d
gpheno.d2$Outcone[which(gpheno.d2$sampleID %in% lateprog)] <- 'Progression'

# For an unbiased analysis of the validation set we need to collapse ambiguous transcripts
gdata.v.collapsed <- aggregate(gdata.v, by=list(rownames(gdata.v)), FUN=mean)
rownames(gdata.v.collapsed) <- gdata.v.collapsed$Group.1
gdata.v.collapsed$Group.1 <- NULL

# Function to compare pathways from a given dataset and pheno array
pathwayComp <- function(gd, gp, fdr_limit=0.05) {
  # gd <- normalizeBetweenArrays(gd)
  uvv <- limmaCompare(gd, gp, fdr_limit = 1)
  ranks <- uvv$t
  names(ranks) <- rownames(uvv)
  gsea <- fgsea(gene.lists$c2.kegg.gsets, stats = ranks, nperm = 1000, maxSize = 500)
  gsea <- gsea[which(gsea$padj < fdr_limit),]
  gsea <- gsea[order(abs(gsea$NES), decreasing = T),]
  gsea$leadingEdge <- sapply(gsea$leadingEdge, function(a) {paste(a, collapse=",")})
  
  g <- gage(gd, gsets = gene.lists$c2.kegg.gsets, samp = which(gp$Outcone == 'Progression'), ref = which(gp$Outcone == 'Regression'), compare = 'unpaired')
  greater = data.frame(g$greater)
  less = data.frame(g$less)
  x <- rbind(
    greater[which(greater$q.val < fdr_limit), 1:4],
    less[which(less$q.val < fdr_limit), 1:4]
  )
  x <- x[order(abs(x$stat.mean), decreasing = T),]
  
  return(list(data.frame(gsea), data.frame(x)))
}

# Original unreduced data:
x <- pathwayComp(gdata.d, gpheno.d)
gsea.d <- x[[1]]
gage.d <- x[[2]]
x <- pathwayComp(gdata.v.collapsed, gpheno.v)
gsea.v <- x[[1]]
gage.v <- x[[2]]

# Break down in to shared, discover only and validation only
gsea.shared <- intersect(gsea.d$pathway, gsea.v$pathway)
gsea.donly <- gsea.d$pathway[which(!(gsea.d$pathway %in% gsea.v$pathway))]
gsea.vonly <- gsea.v$pathway[which(!(gsea.v$pathway %in% gsea.d$pathway))]
gage.shared <- intersect(rownames(gage.d), rownames(gage.v))
gage.donly <- rownames(gage.d)[which(!(rownames(gage.d) %in% rownames(gage.v)))]
gage.vonly <- rownames(gage.v)[which(!(rownames(gage.v) %in% rownames(gage.d)))]

x <- pathwayComp(gdata.d.r, gpheno.d)
gsea.d.r <- x[[1]]
gage.d.r <- x[[2]]
x <- pathwayComp(gdata.v.r2, gpheno.v)
gsea.v.r2 <- x[[1]]
gage.v.r2 <- x[[2]]

gsea.r.shared <- intersect(gsea.d.r$pathway, gsea.v.r2$pathway)
gsea.r.donly <- gsea.d.r$pathway[which(!(gsea.d.r$pathway %in% gsea.v.r2$pathway))]
gsea.r.vonly <- gsea.v.r2$pathway[which(!(gsea.v.r2$pathway %in% gsea.d.r$pathway))]
gage.r.shared <- intersect(rownames(gage.d.r), rownames(gage.v.r2))
gage.r.donly <- rownames(gage.d.r)[which(!(rownames(gage.d.r) %in% rownames(gage.v.r2)))]
gage.r.vonly <- rownames(gage.v.r2)[which(!(rownames(gage.v.r2) %in% rownames(gage.d.r)))]



# Write output to Excel, as a supplementary table S5:
# Include a write-up sheet
x <- data.frame(
  x = c(
    "Supplementary Table 5: Pathway Analysis and Gene Expression Platform Comparison.",
    "",
    "We present pathway analysis between progressive and regressive samples using two methods: Gene Set Enrichment Analysis (GSEA) and GAGE.",
    "In each case we present our previously published Discovery set, profiled on Illumina microarrays, and our Validation set, profiled on Affymetrix microarrays.",
    "To compare the platforms, we include analysis of 'reduced' versions of each dataset, including only genes that are present on both platforms and are unambiguous (one transcript per gene).",
    "We find immune pathways present in the Affymetrix validation set which are not seen in the Illumina discovery set, with weaker signals seen using only the reduced dataset."
  )
)

WriteXLS(
  c("x", "gsea.d", "gsea.v", "gsea.d.r", "gsea.v.r2", "gage.d", "gage.v", "gage.d.r", "gage.v.r2"), 
  SheetNames = c("Description", "GSEA-Illumina", "GSEA-Affymetrix", "GSEA-Illumina-reduced", "GSEA-Affymetrix-reduced", "GAGE-Illumina", "GAGE-Affymetrix", "GAGE-Illumina-reduced", "GAGE-Affymetrix-reduced"),
  ExcelFileName = paste0(supdir, "TableS5.xlsx"), AdjWidth = T, row.names = T
)



##############################################################################
# Additional Calculations
# These are quoted in the main text, though not used in figures.
##############################################################################
message("Performing additional calculations")

# Differential analysis between Prog/Reg limited to checkpoints
uvv <- limmaCompare(gdata.v[intersect(gene.lists$checkpoints, rownames(gdata.v)),], gpheno.v, fdr_limit = 1)
sel <- which(uvv$fdr < 0.05)
write(paste0("PvR Limma analysis limited to immunomodulators: sig genes ", paste(paste(rownames(uvv)[sel], uvv$fdr[sel], sep=','), sep = '; ')), file = opfile, append = T)
write(paste0("Corresponding p-value for receptor TNFSRSF9: ", uvv["TNFRSF9",]$fdr), file = opfile, append = T)

# Comparison of HLA types between prog and reg
# Find all known HLA types - limit to 4 point resolution
all.hlas <- unique(unlist(lapply(pheno$hla.type, function(x) {
  unlist(strsplit(x, "|", fixed = T))
})))
all.hlas <- sort(all.hlas[which(!is.na(all.hlas))])
df <- data.frame(
  hla = all.hlas,
  prog = unlist(lapply(all.hlas, function(hla) {length(grep(hla, pheno.nodups[which(pheno.nodups$Outcome == 'Progression'),]$hla.type))})),
  reg = unlist(lapply(all.hlas, function(hla) {length(grep(hla, pheno.nodups[which(pheno.nodups$Outcome == 'Regression'),]$hla.type))}))
)
df$n.prog <- length(which(pheno.nodups$Outcome == 'Progression' & pheno.nodups$Whole.Genome.Sequencing))
df$n.reg <- length(which(pheno.nodups$Outcome == 'Regression' & pheno.nodups$Whole.Genome.Sequencing))
df$fisher.p <- unlist(lapply(1:dim(df)[1], function(i) {
  m = matrix(c(df$prog[i], df$n.prog[i] - df$prog[i], df$reg[i], df$n.reg[i] - df$reg[i]), ncol = 2)
  return(fisher.test(m)$p.value)
}))
# Demonstrate no HLAs are significantly different between prog and reg:
length(which(df$fisher.p < 0.05))

# Try with 2-digit resolution:
all.hlas <- unique(substr(all.hlas, 1,7))
df <- data.frame(
  hla = all.hlas,
  prog = unlist(lapply(all.hlas, function(hla) {length(grep(hla, pheno.nodups[which(pheno.nodups$Outcome == 'Progression'),]$hla.type))})),
  reg = unlist(lapply(all.hlas, function(hla) {length(grep(hla, pheno.nodups[which(pheno.nodups$Outcome == 'Regression'),]$hla.type))}))
)
df$n.prog <- length(which(pheno.nodups$Outcome == 'Progression' & pheno.nodups$Whole.Genome.Sequencing))
df$n.reg <- length(which(pheno.nodups$Outcome == 'Regression' & pheno.nodups$Whole.Genome.Sequencing))
df$fisher.p <- unlist(lapply(1:dim(df)[1], function(i) {
  m = matrix(c(df$prog[i], df$n.prog[i] - df$prog[i], df$reg[i], df$n.reg[i] - df$reg[i]), ncol = 2)
  return(fisher.test(m)$p.value)
}))
# Demonstrate no HLAs are significantly different between prog and reg:
length(which(df$fisher.p < 0.05))

# What proportion of _patients_ have HLA loss according to LOHHLA? 
# For patients with multiple samples, include that patient if any samples show HLA LOH
pts <- unique(pheno$Patient.Number[which(pheno$Whole.Genome.Sequencing)])
pts.lohhla <- unlist(lapply(pts, function(x) {
  any(pheno[which(pheno$Patient.Number == x & pheno$Whole.Genome.Sequencing),]$hla.loss > 0)
}))
names(pts.lohhla) <- pts
pts.outcome <- unlist(lapply(pts, function(x) {
  pheno$Outcome[which(pheno$Patient.Number == x)][1]
}))
m <- table(pts.lohhla, pts.outcome)
p <- fisher.test(m)$p.value

write(paste0("Patients identified with HLA LOH (progF, progT, regF, regT):", paste(m, collapse=",")), file = opfile, append = T)
write(paste0(" -> ", 100*length(which(pts.lohhla))/length(pts.lohhla), "% of all patients with CIS have HLA LOH"), file = opfile, append=T)
write(paste0("Fisher p-value for prog vs reg: ", p), file=opfile, append=T)

message("Plots generated successfully")








