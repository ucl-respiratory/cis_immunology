# Create all figures and supplementary data for this paper.

version <- "1.1"

# Output directory
figdir <- paste0("~/Dropbox/CIS_Immunology/cis_immunology2/results/v", version, "/figures/")
supdir <- paste0("~/Dropbox/CIS_Immunology/cis_immunology2/results/v", version, "/supdata/")
dir.create(figdir, showWarnings = F, recursive = T)
dir.create(supdir, showWarnings = F, recursive = T)

# Create a text file for other results quoted in the text:
opfile <- paste0(supdir, "results.in.text.txt")
write(c("# Statistical results quoted in the main text:", ""), file = opfile)

# Required libraries:
libs <- c("tidyverse", "ggplot2", "ggsignif", "ggrepel", "ggpubr", "cowplot", "magick", "WriteXLS", "foreach", "ChAMP", "pathview", "pheatmap", "grid", "ggplotify", "biomaRt",
          "httr", "jsonlite", "xml2", "sva", "gdata", "BSgenome.Hsapiens.UCSC.hg19", "Homo.sapiens", "limma", "lme4", "car", "tibble", "dndscv")

for(lib in libs) {
  library(lib, character.only = T)
}

# Export the current versions of required libraries to the results folder
vfile <- paste0(supdir, "package.versions.csv")
version.data <- installed.packages()[libs,]
write.csv(version.data, file = vfile)


# Figure width - assume double column = 183mm (see https://www.nature.com/nature/for-authors/formatting-guide)
# -> 7.204" (needed for save_plot base_width method)
fig.width <- 7.204
fig.maxheight <- 9.72

# Format figures in 12 point Times New Roman (suggested by Nature formatting guide)
theme_set(theme_cowplot(font_size=11, font_family = "Times"))

# Set default colour palette - for consistency, redefine the ggplot function as follows:
# Use Dark2 from Rcolorbrewer, which is colourblind safe
pal <- 'Set2'
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette=pal, direction = -1) + scale_fill_brewer(palette = pal, direction = -1)
# Define manual colours for progressive, regressive, control in that order (used in a couple of plots)
library(RColorBrewer)
plotcols <- brewer.pal(3, pal)[c(2,1,3)]

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

# Minor data fixes
gm.tcga.pheno$PatientBarcode <- substr(gm.tcga.pheno$submitter_id2, 1, 12)

# Include a pheno data frame with no duplicated patients (e.g. for HLA analysis).
# Prioritise multi-omic samples.
pheno.nodups <- pheno[order(pheno$Methylation, pheno$Whole.Genome.Sequencing, pheno$Stroma.GXN, pheno$Gene.expression, pheno$HandE, pheno$Nanostring, pheno$IHC_K, decreasing = T),]
pheno.nodups <- pheno.nodups[which(!duplicated(pheno.nodups$Patient.Number)),]

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

gdata.davoli.t <- do.davoli(gdata.pair.t)
gdata.davoli.s <- do.davoli(gdata.pair.s)

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
# IHC quantification data
##############################################################################
load("data/ihc_quant.RData")
# Recalculate percentages (note that here hems represents cells identified on haematoxylin stain that are negative for other markers):
preinvSum$total.cells <- (preinvSum$hems_IHC + preinvSum$cd4s + preinvSum$cd8s + preinvSum$foxp3s)
preinvSum$cd4_per <- preinvSum$cd4s / preinvSum$total.cells
preinvSum$cd8_per <- preinvSum$cd8s / preinvSum$total.cells
preinvSum$foxp3_per <- preinvSum$foxp3s / preinvSum$total.cells
preinvSum2$total.cells <- (preinvSum2$hems_IHC + preinvSum2$cd4s + preinvSum2$cd8s + preinvSum2$foxp3s)
preinvSum2$cd4_per <- preinvSum2$cd4s / preinvSum2$total.cells
preinvSum2$cd8_per <- preinvSum2$cd8s / preinvSum2$total.cells
preinvSum2$foxp3_per <- preinvSum2$foxp3s / preinvSum2$total.cells

# Match to pheno
preinvSum$SampleID <- pheno$SampleID[match(preinvSum$FileName, paste0(gsub(" ", "_", str_trim(pheno$Sample.Number..IHC_K.)), "IHC.ndpi"))]
preinvSum$patient <- pheno$Patient.Number[match(preinvSum$FileName, paste0(gsub(" ", "_", str_trim(pheno$Sample.Number..IHC_K.)), "IHC.ndpi"))]

# Annotate the pheno data frame for downstream analysis
preinvSum.cis <- preinvSum[which(preinvSum$site == 'CIS'),]
preinvSum.stroma <- preinvSum[which(preinvSum$site == 'Stroma'),]

pheno <- merge(pheno, preinvSum.cis, by='SampleID', all.x=T)
colnames(preinvSum.stroma) <- paste0('stroma_', colnames(preinvSum.stroma))
pheno <- merge(pheno, preinvSum.stroma, by.x='SampleID', by.y='stroma_SampleID', all.x=T)

pheno$cd4_grad <- pheno$cd4_per - pheno$stroma_cd4_per
pheno$cd8_grad <- pheno$cd8_per - pheno$stroma_cd8_per
pheno$stroma_foxp3_per <- pheno$stroma_foxp3_per - pheno$stroma_foxp3_per

# Correct some variables after merging
pheno$Outcome <- pheno$Outcome.x

##############################################################################
# Imaging data
##############################################################################
# Load image data, annotate pheno
# Note that this overwrites the preinvSum and preinvSum2 variables above
load('data/image_data_cache.RData')
# Add patient ID to imaging data:
preinvSum$patient <- pheno$Patient.Number[match(preinvSum$UniqueID, pheno$HandE.SampleID)]
# Sometimes we need to match AltSamples:
sel <- which(is.na(preinvSum$patient))
alts <- lapply(pheno$HandE.AltSamples, function(x) {unlist(strsplit(x, ", "))})
pts <- unlist(lapply(sel, function(i){
  uid <- preinvSum$UniqueID[i]
  samp <- which(unlist(lapply(alts, function(x) {uid %in% x})))
  return(pheno$Patient.Number[samp])
}))
preinvSum$patient[sel] <- pts

preinvSum.cis <- preinvSum[which(preinvSum$site == 'CIS'),]
preinvSum.stroma <- preinvSum[which(preinvSum$site == 'Stroma'),]
pheno$lymphocytes_per <- preinvSum.cis$lymphocyte_per[match(pheno$HandE.SampleID, preinvSum.cis$UniqueID)]
# Also look at lymphocyte gradient tumour - stroma
pheno$lym_grad <- unlist(lapply(pheno$HandE.SampleID, function(sid) {
  preinvSum.cis$lymphocyte_per[match(sid, preinvSum.cis$UniqueID)] - preinvSum.stroma$lymphocyte_per[match(sid, preinvSum.cis$UniqueID)]
}))


##############################################################################
# Figure 1: IHC and Danaher data
# TODO: adjust for patient ID
##############################################################################
message("Plotting Figure 1")
# a/b - prog/reg IHC images
# c - IHC quantified CD8 cells
# d - Danaher TIL scores
# e - MethylCS pc scores
# f - imaging data

fig.a <- ggdraw() + draw_image("resources/fig1_progressive.jpg", scale=0.95) +
  draw_label('Progressive CIS lesion', y=0.13, size = 11)
fig.b <- ggdraw() + draw_image("resources/fig1_regressive.jpg", scale=0.95) +
  draw_label('Regressive CIS lesion', y=0.13, size = 11)
fig.ab <- plot_grid(fig.a, fig.b, labels=c('a','b'))


# Plot IHC data
plotdata <- pheno[which(pheno$IHC_K),]
# plotdata <- ihc
plotdata$Outcome <- gsub("ression", ".", plotdata$Outcome)
p.cd4 <- compare.fn(cd4_perArea ~ Outcome + (1 | patient), data = plotdata)
p.cd8 <- compare.fn(cd8_perArea ~ Outcome + (1 | patient), data = plotdata)
p.foxp3 <- compare.fn(foxp3_perArea ~ Outcome + (1 | patient), data = plotdata)
# Assess _net_ cytotoxic activity (CD8/FOXP3), excluding samples with no FOXP3 cells
p.cd8foxp3 <- compare.fn(cd8_perArea/foxp3_perArea ~ Outcome + (1 | patient), data = plotdata[which(plotdata$foxp3_perArea > 0),])

# Plot only CD8 here (only this is significant)
fig.c <- ggplot(plotdata, aes(x=Outcome, y=cd8_perArea)) + 
  geom_boxplot(aes(fill = Outcome)) +
  # geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.1) +
  geom_point() +
  stat_pvalue_manual(p.cd8, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  # geom_signif(comparisons = list(c("Progression","Regression")), annotations = 'tmp') +
  guides(fill=FALSE) +
  theme(axis.title.x = element_blank()) +
  # facet_wrap(~cell, strip.position = 'right') +
  ylab('CD8 cell count / area')

# Save other comparisons in a text file
write(paste0("CD4 P vs R p=", p.cd4$p), file = opfile, append = T)
write(paste0("CD8 P vs R p=", p.cd8$p), file = opfile, append = T)
write(paste0("FOXP3 P vs R p=", p.foxp3$p), file = opfile, append = T)
write(paste0("CD8:FOXP3 P vs R p=", p.cd8foxp3$p), file = opfile, append = T)

# Include stromal comparisons here
p.cd4s <- compare.fn(stroma_cd4_perArea ~ Outcome + (1 | patient), data = plotdata)
p.cd8s <- compare.fn(stroma_cd8_perArea ~ Outcome + (1 | patient), data = plotdata)
p.foxp3s <- compare.fn(stroma_foxp3_perArea ~ Outcome + (1 | patient), data = plotdata)
p.cd8foxp3s <- compare.fn(stroma_cd8_perArea/stroma_foxp3_perArea ~ Outcome + (1 | patient), data = plotdata[which(plotdata$stroma_foxp3_perArea > 0),])
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


# Danaher TIL scores
plotdata <- data.frame(til.score = gdata.danaher.t$til.score, Outcome = gpheno.pair$Outcome, patient = factor(gpheno.pair$Patient.Number))
plotdata$Outcome <- gsub("ression", ".", plotdata$Outcome)
p <- compare.fn(til.score ~ Outcome + (1 | patient), data = plotdata)
fig.d <- ggplot(plotdata, aes(x=Outcome, y=til.score)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  guides(fill=FALSE) +
  theme(axis.title.x = element_blank()) +
  ylab('Danaher TIL score')


# methylCS scores:
sel <- which(mpheno$Sample_Group %in% c("Progressive", "Regressive"))
plotdata <- data.frame(pc.immune = rowSums(methCS[sel,c(2,3,4,5,6,8,10,11)]), Outcome = mpheno$Sample_Group[sel], patient = factor(mpheno$Patient_ID[sel]))
plotdata$Outcome <- gsub('ive$', 'ion', plotdata$Outcome)
plotdata$Outcome <- gsub("ression", ".", plotdata$Outcome)
p <- compare.fn(pc.immune ~ Outcome + (1 | patient), data = plotdata)
fig.e <- ggplot(plotdata, aes(x=Outcome, y=pc.immune)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  guides(fill=FALSE) +
  theme(axis.title.x = element_blank()) +
  ylab('MethylCIBERSORT % immune cells') 


# Imaging data:
plotdata <- preinvSum.cis
plotdata$Outcome <- gsub("ive", "ion", plotdata$Outcome)
plotdata$Outcome <- gsub("ression", ".", plotdata$Outcome)
p <- compare.fn(lymphocyte_per ~ Outcome + (1 | patient), data = plotdata)
fig.f <- ggplot(plotdata, aes(x=Outcome, y=lymphocyte_per)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  guides(fill=FALSE) +
  theme(axis.title.x = element_blank()) +
  ylab('% lymphocytes on H&E')

fig.cdef <- plot_grid(fig.c, fig.d, fig.e, fig.f, labels=c('c','d','e','f'), ncol=4, nrow=1)

fig <- plot_grid(fig.ab, fig.cdef, ncol=1, nrow=2, rel_heights = c(1,1))
save_plot(paste0(figdir, "fig1.pdf"), fig, nrow = 2, ncol=1, base_width = fig.width, base_height = 2.9)


##############################################################################
# Figure 2: Clustering analyses
# Cluster based on Danaher data and methylCIBERSORT
##############################################################################
message("Plotting Figure 2")

annot.colors <- list(outcome = c('Progression' = plotcols[1], 'Regression' = plotcols[2], 'Progressive' = plotcols[1], 'Regressive' = plotcols[2]))
# Clustering on affy gxn immune genes:
annot.gxn <- data.frame(outcome = factor(gpheno.v$Outcone, levels=c('Regression', 'Progression')), row.names = gpheno.v$sampleID)
# pheatmap(gdata.v[intersect(gene.lists$all.immune, rownames(gdata.v)),], annotation_col = annot.gxn, scale='row', show_rownames = F, show_colnames = F, cutree_cols = 2)

# Similar using deconvoluted Danaher values:
fig.a <- pheatmap(t(gdata.danaher.t[,-c(2,15,16)]), annotation_col = annot.gxn, scale='row', show_colnames = F, legend = F, annotation_legend = F, cutree_cols = 2, treeheight_row = 0, annotation_colors = annot.colors)

# Similar using methylCIBERSORT
# Exclude Controls - these won't have TILs as they are brushings
annot.meth <- data.frame(outcome = factor(mpheno$Sample_Group, levels=c('Regressive', 'Progressive')), row.names = mpheno$sampleID)
sel <- which(annot.meth$outcome %in% c("Progressive", "Regressive"))
plotdata <- t(methCS[sel, 1:11])
# Rename some columns
rownames(plotdata)[which(rownames(plotdata) == 'Eos')] <- 'Eosinophils'
rownames(plotdata)[which(rownames(plotdata) == 'CD4_Eff')] <- 'CD4'
rownames(plotdata)[which(rownames(plotdata) == 'Neu')] <- 'Neutrophils'
fig.b <- pheatmap(plotdata, annotation_col = annot.meth, show_colnames = F, legend = F, cutree_cols = 2, annotation_legend = F, treeheight_row = 0, annotation_colors = annot.colors)

# This shows a cluster with lots of cancer i.e. low immune infiltration

# Show the groups TIL score / percent lymphocytes
groups <- cutree(fig.a$tree_col, k = 2)
plotdata <- data.frame(
  group = factor(groups, levels = c(2,1), labels = c('cold', 'hot')),
  til.score = gdata.danaher.t[names(groups), 'til.score'],
  Outcome = gpheno.pair$Outcome[match(names(groups), gpheno.pair$SampleID)],
  patient = gpheno.pair$Patient.Number[match(names(groups), gpheno.pair$SampleID)]
)
p <- compare.fn(til.score ~ group + (1 | patient), data = plotdata, comparison=list('hot','cold'))
fig.c <- ggplot(plotdata, aes(x=group, y=til.score)) +
  geom_boxplot() +
  geom_point(aes(color=Outcome)) +
  xlab("Cluster") +
  ylab("TIL score") +
  stat_pvalue_manual(p, xmin = "group1", xmax="group2", label = "p={p}", tip.length = 0.01, size=3) +
  theme(axis.title.x = element_blank()) +
  guides(color=F)
  # geom_signif(comparisons = list(c('1','2')))
  # guides(color=F)

groups <- cutree(fig.b$tree_col, k = 2)
plotdata <- data.frame(
  group = factor(groups, levels = c(1,2), labels=c('cold', 'hot')),
  pc.lym = rowSums(methCS[sel,c(2,3,4,5,6,8,10,11)]),
  Outcome = mpheno$Sample_Group[match(names(groups), mpheno$sampleID)],
  patient = mpheno$Patient_ID[match(names(groups), mpheno$sampleID)]
)
p <- compare.fn(pc.lym ~ group + (1 | patient), data = plotdata, comparison = list('hot', 'cold'))
fig.d <- ggplot(plotdata, aes(x=group, y=pc.lym)) +
  geom_boxplot() +
  geom_point(aes(color=Outcome)) +
  xlab("Cluster") +
  ylab("% immune cells") +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  theme(axis.title.x = element_blank()) +
  guides(color=F)

fig.main <- plot_grid(as.grob(fig.a), as.grob(fig.b), fig.c, fig.d, nrow = 2, labels=c('a','b','c','d'), rel_heights = c(5,3))

# Make a separate legend plot under the others
plotdata <- data.frame(
  x=c(0.25,1.1),
  y=c(1,1),
  label = c('Progressive', 'Regressive')
)
leg <- ggplot(plotdata, aes(x=x, y=y, label=label)) +
  geom_point(aes(color=label)) +
  xlim(-1,2.5) +
  geom_label(nudge_x = 0.33, label.size = 0) +
  guides(color=F) +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  annotate('text', x=-0.25, y=1, label="Outcome:")

fig <- plot_grid(fig.main, leg, ncol = 1, rel_heights = c(15,1))

save_plot(paste0(figdir, "fig2.pdf"), fig, nrow = 2, ncol=2, base_width = fig.width/2)

# Clustering of immune MVPs:
# pheatmap(mdata.nosnp[intersect(rownames(mdata.nosnp), mvps.immune), rownames(annot.meth)], annotation_col = annot.meth, scale='row', show_rownames = F, show_colnames = F)

##############################################################################
# Figure 3: Summary of genomic, epigenetic and transcriptomic changes
##############################################################################
message("Plotting Figure 3")

# Find muts in these genes only and check VEP status
mhc.muts <- muts.all[which(muts.all$gene %in% gene.lists$hla.assoc),]
#mhc.muts <- muts.all[which(muts.all$gene %in% tgfb.genes),]

p <- pheno[which((pheno$Whole.Genome.Sequencing | pheno$Methylation | pheno$Stroma.GXN) & pheno$Outcome != 'Control'),]
p <- p[order(p$Whole.Genome.Sequencing, p$Methylation, p$Stroma.GXN),]
hla.changes <- muts.for.plotting(gene.lists$hla.assoc, pheno = p)

# Define a colour scheme (uses help from RColorBrewer):
cols <- c(
  'Overexpression'='#fc8d59',
  'Underexpression'='#91bfdb',
  
  'Hypermethylation'='#af8dc3',
  'Hypomethylation'='#ffffbf',
  
  'Amplification'='#f0027f',
  'CN deletion'='#386cb0',
  
  'LOH'='#386cb0',
  
  'Missense'='#7fc97f',
  'Nonsense'='#beaed4',
  'Frameshift'='#ffff99',
  'Rearrangement'='#fdc086',
  'Start Lost'='#bf5b17',
  
  'Normal'='#eeeeee'
)

plotdata <- hla.changes
plotdata$gene <- factor(plotdata$gene, levels = gene.lists$hla.assoc)
# Remove rows with no changes
plotdata <- plotdata[-which(plotdata$gene %in% c('CNX', 'HSPA', 'HSPC')),]

samps.ordered <- unique(p$SampleID[order((p$Whole.Genome.Sequencing + p$Methylation + p$Stroma.GXN), p$Whole.Genome.Sequencing, p$Methylation, p$Stroma.GXN, decreasing = T)])

fig.main <- ggplot(plotdata, aes(
    # x = factor(sample, levels=unique(plotdata$sample[order(plotdata$n.modalities, plotdata$n.changes, decreasing = T)])), 
    x = factor(sample, levels=samps.ordered),
    y = factor(mod, levels=c('gxn', 'meth', 'genomic'))
  )) +
  geom_tile(aes(fill=factor(type, levels=names(cols))), colour='white', size=0.25) +
  scale_fill_manual(values=cols, na.value='#ffffff') +
  # scale_fill_manual(values=c('#DDDDDD','#56B4E9', '#FF0000')) +
  facet_grid(gene ~ factor(outcome), scales = 'free_x', space='free_x', switch = 'y', ) +
  theme(
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    strip.text.y = element_text(size=8), strip.text.x = element_blank(), strip.background = element_blank(),
    legend.title = element_blank(),
    axis.line = element_blank()
  )
  # guides(fill=F)

# Make a bar chart along the top showing number of changes/modalities:
# plotdata2 <- data.frame(
#   sample = samps.ordered,
#   n.modalities = unlist(lapply(samps.ordered, function(x) {
#     y <- pheno[which(pheno$SampleID == x),]
#     return(sum(y$Whole.Genome.Sequencing, y$Methylation, y$Stroma.GXN))
#   })),
#   n.changes = unlist(lapply(samps.ordered, function(x) {
#     y <- hla.changes[which(hla.changes$sample == x),]
#     return(length(which(y$result == 1)))
#   }))
# )
# plotdata2$n.possible.changes <- plotdata2$n.modalities * length(unique(plotdata$gene))

# Alternate:
plotdata3 <- plotdata[which(plotdata$result >= 0),]
top.bar <- ggplot(plotdata3, aes(factor(sample, levels=samps.ordered))) +
  geom_bar(aes(fill=factor(result))) + 
  facet_grid(~ factor(outcome), scales = 'free_x', space='free_x', switch = 'y') +
  guides(fill=F) +
  theme(
    axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.text.y = element_text(size=8), strip.background = element_blank(),
    legend.title = element_blank(),
    axis.line.x = element_blank(), 
    # plot.margin = margin(0,109, 0,13, "pt") # Use if using no axis
    plot.margin = margin(0,109, 0,1, "pt")
  ) +
  scale_fill_manual(values=c("#EEEEEE", "#ef8a62"))


# Plot together:
fig <- plot_grid(top.bar, fig.main, nrow = 2, rel_heights = c(1,10))


# Alternate with genomic changes only:
# hla.changes.genomic <- hla.changes[which(hla.changes$mod == 'genomic' & !is.na(hla.changes$type)),]
# hla.changes.genomic$sample <- as.character(hla.changes.genomic$sample)
# fig2 <- ggplot(hla.changes.genomic, aes(
#   x = factor(sample), 
#   y = factor(gene, levels = rev(gene.lists$hla.assoc))
# )) +
#   geom_tile(aes(fill=factor(type, levels=names(cols))), colour='white', size=0.25) +
#   scale_fill_manual(values=cols, na.value='#ffffff') +
#   facet_grid( ~ factor(outcome), scales = 'free_x', space = 'free_x') +
#   # scale_fill_manual(values=c('#DDDDDD','#56B4E9', '#FF0000')) +
#   theme(
#     axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
#     strip.text.y = element_text(size=7), strip.background = element_blank(),
#     legend.title = element_blank()
#   )

save_plot(paste0(figdir, "fig3.pdf"), fig, base_width = fig.width, base_height = 9)

# Check which samples have at least one change, and look for P vs R differences
hla.changes.by.sample <- aggregate(hla.changes$result, by=list(hla.changes$sample), FUN=function(x) {any(x == 1)})
rownames(hla.changes.by.sample) <- hla.changes.by.sample$Group.1
hla.changes.by.sample$Group.1 <- NULL
hla.changes.by.sample$Outcome <- pheno$Outcome[match(rownames(hla.changes.by.sample), pheno$SampleID)]
m <- table(hla.changes.by.sample$x, hla.changes.by.sample$Outcome)
chisq.test(m)

# Check if this is true even after correcting for mutational burden
# Here we count the number of changes
hla.changes.by.sample$count <- unlist(lapply(rownames(hla.changes.by.sample), function(x) {
  length(which(hla.changes$sample == x & hla.changes$result == 1))
}))
hla.changes.by.sample$burden <- pheno$burden[match(rownames(hla.changes.by.sample), pheno$SampleID)]
hla.changes.by.sample$count.corr <- hla.changes.by.sample$count / hla.changes.by.sample$burden

# Also check number of changes using appropriate mixed model strategy
hla.changes.by.sample$patient <- factor(pheno$Patient.Number[match(rownames(hla.changes.by.sample), pheno$SampleID)])
p <- compare.fn(count ~ Outcome + (1 | patient), data = hla.changes.by.sample)
write(paste0("Aberrations in HLA genes were more common in progressive samples (p=", p$p, ")"), file = opfile, append = T)

# Corrected for burden (not in the paper as it is not really correct to include meth/gxn in burden correction):
p <- compare.fn(count/burden ~ Outcome + (1 | patient), data = hla.changes.by.sample)

# Repeat above using only genomic changes:
sel <- which(hla.changes$mod == 'genomic')
hla.genomic.by.sample <- aggregate(hla.changes$result[sel], by=list(hla.changes$sample[sel]), FUN=function(x) {any(x == 1)})
rownames(hla.genomic.by.sample) <- hla.genomic.by.sample$Group.1
hla.genomic.by.sample$Group.1 <- NULL
hla.genomic.by.sample$Outcome <- pheno$Outcome[match(rownames(hla.genomic.by.sample), pheno$SampleID)]
hla.genomic.by.sample$count <- unlist(lapply(rownames(hla.genomic.by.sample), function(x) {
  length(which(hla.changes$sample == x & hla.changes$result == 1))
}))
hla.genomic.by.sample$burden <- pheno$burden[match(rownames(hla.genomic.by.sample), pheno$SampleID)]
# hla.genomic.by.sample$burden2 <- unlist(lapply(rownames(hla.genomic.by.sample), function(x) {
#   length(which(muts.all$sampleID == x))
# }))
hla.genomic.by.sample <- hla.genomic.by.sample[which(rownames(hla.genomic.by.sample) %in% pheno$SampleID[which(pheno$Whole.Genome.Sequencing)]),]

hla.genomic.by.sample$patient <- factor(pheno$Patient.Number[match(rownames(hla.genomic.by.sample), pheno$SampleID)])
p <- compare.fn(count ~ Outcome + (1 | patient), data = hla.genomic.by.sample)
write(paste0("*Genomic* Aberrations in HLA genes were more common in progressive samples (p=", p$p, ")"), file = opfile, append = T)
p <- compare.fn(count/burden ~ Outcome + (1 | patient), data = hla.genomic.by.sample)
write(paste0("Genomic Aberrations in HLA genes were more common in progressive samples even after correcting for overall mutational burden (p=", p$p, ")"), file = opfile, append = T)

# Add some additional comparisons:
# HLA-A GXN changes:
plotdata <- data.frame(
  hla.a = as.numeric(gdata.v["HLA-A",]),
  hla.b = as.numeric(gdata.v["HLA-B",]),
  hla.c = as.numeric(gdata.v["HLA-C",]),
  Outcome = gsub('ression', '.', gpheno.v$Outcone),
  patient = factor(gpheno.v$Patient)
)
p <- compare.fn(hla.a ~ Outcome + (1 | patient), data = plotdata)
write(paste0("Expression of HLA-A was correspondingly reduced in progressive compared to regressive samples (p=", p$p, ")"), file = opfile, append = T)
p <- compare.fn(hla.b ~ Outcome + (1 | patient), data = plotdata)
write(paste0("Expression of HLA-B was correspondingly reduced in progressive compared to regressive samples (p=", p$p, ")"), file = opfile, append = T)
p <- compare.fn(hla.c ~ Outcome + (1 | patient), data = plotdata)
write(paste0("Expression of HLA-C was correspondingly reduced in progressive compared to regressive samples (p=", p$p, ")"), file = opfile, append = T)

##############################################################################
# Figure 4: TIL plots PvsR; TIL gradient vs FTGFB; IHC + hotspot analysis
# Awaiting data for second part
##############################################################################
message("Plotting Figure 4")

# a) Stromal volcano plot (showing no differences)
uvv <- limmaCompare(gdata.pair.s, gpheno.pair, fdr_limit = 1)
fig.a <- ggplot(uvv, aes(x=logratio, y=-log(fdr))) +
  geom_point() +
  ylim(0, -log(0.05)) +
  geom_hline(yintercept=-log(0.05), color='grey') +
  annotate("text", x=1, y=2.9, label="FDR < 0.05") +
  xlab('Log fold change') +
  ylab('-log(False Discovery Rate)')

# b-d) TGFB story 
samples <- pheno[which(pheno$Stroma.GXN),]
pdan.data <- rbind(
  cbind(gdata.danaher.t[samples$SampleID,], data.frame(type='tissue', Outcome=samples$Outcome, sampleID=samples$SampleID, patient=samples$Patient.Number)),
  cbind(gdata.danaher.s[samples$SampleID,], data.frame(type='stroma', Outcome=samples$Outcome, sampleID=samples$SampleID, patient=samples$Patient.Number))
)

# Add FTGFB and TIL gradient calculations to samples data frame
samples$til.gradient <- as.numeric(gdata.danaher.t[samples$SampleID,]$til.score) - as.numeric(gdata.danaher.s[samples$SampleID,]$til.score)
samples$ftgfb <- do.ftgfb(gdata.pair.s[,samples$SampleID])

# Paired comparison of regressive patients, TIL score tissue vs stroma
# Calculate a p-value for all samples, p and r:
p <- compare.fn(til.score ~ type + (1 | sampleID) + (1 | patient), data=pdan.data)
write(paste0("Paired comparison of TIL score tissue vs stroma in ALL samples, p= ", signif(p$p, 3)), file = opfile, append = T)
p <- compare.fn(til.score ~ type + (1 | sampleID) + (1 | patient), data=pdan.data[which(pdan.data$Outcome == 'Regression'),])
write(paste0("Paired comparison of TIL score tissue vs stroma in REGRESSIVE samples, p= ", signif(p$p, 3)), file = opfile, append = T)
p <- compare.fn(til.score ~ type + (1 | sampleID) + (1 | patient), data=pdan.data[which(pdan.data$Outcome == 'Progression'),])
write(paste0("Paired comparison of TIL score tissue vs stroma in PROGRESSIVE samples, p= ", signif(p$p, 3)), file = opfile, append = T)

plotdata <- pdan.data[pdan.data$Outcome == 'Regression',]
plotdata$type <- str_to_title(plotdata$type)
p <- compare.fn(til.score ~ type + (1 | sampleID) + (1 | patient), data = plotdata, comparison = list("Tissue", "Stroma"))
fig.b <- ggplot(plotdata, aes(x=type, y=til.score)) + 
  geom_boxplot(aes(fill = type)) +
  geom_line(aes(x=type, y=til.score, group=sampleID), colour='grey') +
  geom_point() +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01) +
  guides(fill=FALSE) +
  theme(axis.title.x = element_blank()) +
  ylab("TIL score")

plotdata <- pdan.data[pdan.data$Outcome == 'Progression',]
plotdata$type <- str_to_title(plotdata$type)
p <- compare.fn(til.score ~ type + (1 | sampleID) + (1 | patient), data = plotdata, comparison = list("Tissue", "Stroma"))
fig.c <- ggplot(plotdata, aes(x=type, y=til.score)) + 
  geom_boxplot(aes(fill = type)) +
  geom_line(aes(x=type, y=til.score, group=sampleID), colour='grey') +
  geom_point() +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01) +
  guides(fill=FALSE) +
  theme(axis.title.x = element_blank()) +
  ylab("TIL score")

fig.d <- ggplot(samples, aes(x=ftgfb, y=til.gradient)) +
  geom_point(aes(color=Outcome)) +
  xlab("FTGFB signature") + ylab("TIL gradient") +
  guides(color=FALSE) +
  geom_smooth(method = 'lm') +
  stat_cor()

# e-f) TNFSF9 boxplots (CIS and TCGA)
plotdata <- data.frame(
  tnfsf9 = as.numeric(gdata.pair.t["TNFSF9",]),
  tnfrsf9 = as.numeric(gdata.pair.t["TNFRSF9",]),
  Outcome = gsub("ression", ".", gpheno.pair$Outcome),
  patient = factor(gpheno.pair$Patient.Number)
)
p <- compare.fn(tnfsf9 ~ Outcome + (1 | patient), data = plotdata)
fig.e <- ggplot(plotdata, aes(x=Outcome, y=tnfsf9)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01) +
  ylab('CIS TNFSF9 expression') +
  theme(axis.title.x = element_blank())

# Also check the correlation of TNFSF9 with predicted NK cells
samps <- pheno$SampleID[which(pheno$Stroma.GXN)]
ctypes = colnames(gdata.danaher.t)[-c(2,15)]
df <- data.frame(
  ctype = ctypes,
  r2 = unlist(lapply(ctypes, function(x) {
    cor.test(gdata.danaher.t[samps, x], as.numeric(gdata.pair.t["TNFSF9",samps]))$estimate
  })),
  p = unlist(lapply(ctypes, function(x) {
    cor.test(gdata.danaher.t[samps, x], as.numeric(gdata.pair.t["TNFSF9",samps]))$p.value
  }))
)
df$p.adj <- p.adjust(df$p)
rownames(df) <- df$ctype

write(paste0("Correlation of Cytotoxic cells with TNFSF9 expression, r2 = ", df["Cytotoxic cells",'r2']," p= ", signif(df["Cytotoxic cells",'p'], 3)), file = opfile, append = T)
write(paste0("Correlation of NK cells with TNFSF9 expression, r2 = ", df["NK cells",'r2']," p= ", signif(df["NK cells",'p'], 3)), file = opfile, append = T)


load('data/gdata.tcga.lusc.RData')
plotdata <- data.frame(
  tnfsf9 = log(as.numeric(gdata.tcga.lusc["TNFSF9",])),
  tnfrsf9 = log(as.numeric(gdata.tcga.lusc["TNFRSF9",])),
  outcome = gpheno.tcga.lusc$sample_type, 
  stringsAsFactors = F
)
plotdata$outcome[which(plotdata$outcome == 'Solid Tissue Normal')] <- 'Control'
plotdata$outcome <- factor(plotdata$outcome, levels=c('Primary Tumor', 'Control'))

# Manually create a variable for stat_pvalue_manual using a Wilcox test (no samples here from theh same patient)
p <- tibble(
  .y. = 'tnfsf9',
  group1 = 'Primary Tumor',
  group2 = 'Control',
  p=signif(wilcox.test(plotdata$tnfsf9 ~ plotdata$outcome)$p.value, 2),
  y.position = max(plotdata$tnfsf9) + max(plotdata$tnfsf9)/10
)
fig.f <- ggplot(plotdata, aes(x=outcome, y=tnfsf9)) +
  geom_boxplot(aes(fill=outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01) +
  ylab('TCGA TNFSF9 expression') +
  theme(axis.title.x = element_blank())


fig <- plot_grid(fig.a, fig.b, fig.c, fig.d, fig.e, fig.f, 
          labels=c('a','b','c','d','e','f'),
          ncol = 3, nrow=2)

save_plot(paste0(figdir, "fig4.pdf"), fig, nrow = 2, ncol=3, base_width = fig.width/3)

##############################################################################
# Figure S1
# Tile figure of what was done to each sample
##############################################################################
message("Plotting Figure S1")

df <- rbind(
  data.frame(SampleID = pheno$SampleID, mod = "WGS", val = pheno$Whole.Genome.Sequencing),
  data.frame(SampleID = pheno$SampleID, mod = "Gene expression", val = pheno$Stroma.GXN),
  data.frame(SampleID = pheno$SampleID, mod = "Methylation", val = pheno$Methylation),
  data.frame(SampleID = pheno$SampleID, mod = "IHC", val = pheno$IHC_K),
  data.frame(SampleID = pheno$SampleID, mod = "Image analysis", val = pheno$HandE)
)
df$outcome <- as.character(pheno$Outcome[match(df$SampleID, pheno$SampleID)])
df$val[which(df$val == FALSE)] <- NA
df$val[which(!is.na(df$val))] <- df$outcome[which(!is.na(df$val))]
df$val <- factor(df$val, levels=c('Progression', 'Regression', 'Control'))
df$nmods <- unlist(lapply(df$SampleID, function(x) {
  length(which(df$SampleID == x & df$val == TRUE))
}))
df <- df[order(df$outcome, df$nmods, decreasing = T),]
# df$outcome <- factor(df$outcome)
df$SampleID <- factor(as.character(df$SampleID), levels=unique(as.character(df$SampleID)))

fig <- ggplot(df[which(df$outcome %in% c("Progression", "Regression")),], aes(x=mod, y=SampleID, fill=val)) +
  geom_tile(aes(width = 0.8, height = 0.8), size=2) +
  scale_fill_manual(values=plotcols[c(1,2)]) +
  theme(axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.title = element_blank()) +
  ylab("Sample")

save_plot(paste0(figdir, "figS1.pdf"), fig, base_width = fig.width)


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

library(ggplot2)
library(ggsignif)
library(cowplot)
x$Outcome <- gsub('ression', '.', x$Outcome)
p <- compare.fn(pro.mean.t ~ Outcome + (1 | Patient.Number), data = x)
f1 <- ggplot(x, aes(x=Outcome, y=pro.mean.t)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  ylab("Pro-inflammatory cytokines (tissue)") +
  theme(axis.title.x = element_blank())

p <- compare.fn(anti.mean.t ~ Outcome + (1 | Patient.Number), data = x)
f2 <- ggplot(x, aes(x=Outcome, y=anti.mean.t)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  ylab("Anti-inflammatory cytokines (tissue)") +
  theme(axis.title.x = element_blank())

p <- compare.fn(pro.mean.t / anti.mean.t ~ Outcome + (1 | Patient.Number), data = x)
f3 <- ggplot(x, aes(x=Outcome, y=pro.mean.t / anti.mean.t)) +
  geom_boxplot(aes(fill=Outcome)) +
  geom_point() +
  guides(fill=F) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  ylab("Pro:Anti ratio (CIS)") +
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
  ylab("Pro:Anti ratio (Stroma)") +
  theme(axis.title.x = element_blank())

fig <- plot_grid(f1, f2, f3, f4, f5, f6, ncol=3, labels=c('a','b','c','d','e','f'))

save_plot(paste0(figdir, "figS2.pdf"), fig, base_width = fig.width, base_height = fig.width)

# Print also the p-value for CXCL8:
x$cxcl8 <- as.numeric(gdata.pair.t["CXCL8",])
p <- compare.fn(cxcl8 ~ Outcome + (1 | Patient.Number), data = x)
write(paste0("p-value for CXCL8 (uncorrected): ", p$p), file = opfile, append = T)

##############################################################################
# Figure S3
# Pro-anti inflammatory cytokines (individual plots)
##############################################################################
message("Plotting Figure S3")

# Plot the individual cytokines as separate boxplots:
plotdata <- lapply(intersect(rownames(gdata.pair.t), c(gene.lists$cytokines.anti, gene.lists$cytokines.pro)), function(x) {
  data.frame(
    gene=x,
    is.pro = x %in% gene.lists$cytokines.pro,
    val = as.numeric(gdata.pair.t[x,]),
    sampleID = gpheno.pair$SampleID,
    Outcome = gpheno.pair$Outcome
  )
})
plotdata <- do.call('rbind', plotdata)
plotdata$Outcome <- gsub('ression', '.', plotdata$Outcome)
scaleFUN <- function(x) sprintf("%.1f", x)
pro <- ggplot(plotdata[which(plotdata$is.pro),], aes(x=Outcome, y=val)) +
  geom_boxplot(aes(fill=Outcome)) + geom_point()  +
  geom_signif(comparisons = list(c('Prog.', 'Reg.')), map_signif_level = T, vjust = 1.5) +
  guides(fill=F) +
  facet_wrap(~gene, scales = 'free_y') +
  theme(axis.title.x = element_blank()) +
  ylab("CIS gene expression values for pro-inflammatory cytokines") +
  scale_y_continuous(labels=scaleFUN)
anti <- ggplot(plotdata[which(!plotdata$is.pro),], aes(x=Outcome, y=val)) +
  geom_boxplot(aes(fill=Outcome)) + geom_point()  +
  geom_signif(comparisons = list(c('Prog.', 'Reg.')), map_signif_level = T, vjust = 1.5) +
  facet_wrap(~gene, scales = 'free_y') + guides(fill=F) +
  theme(axis.title.x = element_blank()) +
  ylab("CIS gene expression values for anti-inflammatory cytokines") +
  scale_y_continuous(labels=scaleFUN)
fig <- plot_grid(pro, anti, ncol = 2, labels=c('a','b'))

save_plot(paste0(figdir, "figS3.pdf"), fig, base_width = fig.width*1.2, base_height = fig.width)

##############################################################################
# Figure S4
# Neoantigen plots:
#  a) Predicted neoantigens vs mutational burden
#  b) Clonal neoantigens P vs R
#  c) Proportion clonal P vs R
#  d) Proportion clonal corrected for purity
#  e) Affinity P vs R
#  f) Rank P vs R
#  g) Depletion stats P vs R
##############################################################################
message("Plotting Figure S4")

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
  # geom_signif(comparisons = list(c("Progression", "Regression"))) +
  geom_point() +
  guides(fill=F) +
  ylab("Clonal Strong Neoantigens") +
  theme(axis.title.x = element_blank())

p <- compare.fn(neoants.strong.clonal/neoants.strong ~ Outcome + (1 | patient), data = plotdata)
fig.d <- ggplot(plotdata, aes(x=Outcome, y=(neoants.strong.clonal / neoants.strong))) +
  geom_boxplot(aes(fill=Outcome)) +
  stat_pvalue_manual(p, label='p={p}', xmin = "group1", xmax="group2", tip.length = 0.01, size=3) +
  # geom_signif(comparisons = list(c("Progression", "Regression"))) +
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
  # geom_signif(comparisons = list(c("Progression", "Regression"))) +
  geom_point() +
  guides(fill=F) +
  ylab("Rank") +
  theme(axis.title.x = element_blank())

p <- compare.fn(best.DAI ~ Outcome + (1 | patient), data = neos.all)
fig.g <- ggplot(neos.all, aes(x=Outcome, y=best.DAI)) +
  geom_boxplot(aes(fill=Outcome)) +
  # geom_signif(comparisons = list(c("Progression", "Regression"))) +
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
  # geom_signif(comparisons = list(c("Progression", "Regression"))) +
  geom_point() +
  guides(fill=F) +
  ylab("Depletion Score") +
  theme(axis.title.x = element_blank())


fig <- plot_grid(fig.a, fig.b, fig.c, fig.d, fig.e, fig.f, fig.g, fig.h, nrow=3, ncol=3, 
                 labels = c('a', 'b','c','d','e','f','g','h'))

save_plot(paste0(figdir, "figS4.pdf"), fig, ncol=3, nrow=3, base_width = fig.width/3, base_height = fig.width/2.5)

##############################################################################
# Figure S5
# Methylation DMRs across the genome
##############################################################################
message("Plotting Figure S5")

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
  # library(GenomicRanges)
  # library(Homo.sapiens)
  # library(BSgenome.Hsapiens.UCSC.hg19)
  
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
    guides(color=FALSE) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    facet_wrap(~probe) +
    ylab("Methylation beta value") +
    scale_color_manual(values = plotcols)
  
  fig <- plot_grid(fig.a, fig.b, labels=c('a','b'), ncol=1, nrow=2, rel_heights = c(2,1))
  save_plot(paste0(figdir, "figS5.pdf"), fig, nrow = 2, ncol=1, base_height = 5, base_width=fig.width)
  file.remove(paste0(figdir, "circos.tmp.svg"))
} else {
  message("Warning: circos installation not found. Figure S4 will not be plotted.")
}


##############################################################################
# Figure S6
# HLA silencing in CIS and TCGA (plots of GXN vs methylation)
##############################################################################
message("Plotting Figure S6")

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
    ggtitle(gene)
  
  
  plotdata <- data.frame(
    gene.expression = as.numeric(gdata[gene,ol$SampleID]),
    methylation = as.numeric(mdata.genes.nosnp[gene, ol$SampleID]),
    lym_grad = ol$lym_grad
  )
  cis.plots[[gene]] <- ggplot(plotdata, aes(x=methylation, y=gene.expression)) +
    geom_point() +
    geom_smooth(method='lm') +
    stat_cor() +
    xlab("Methylation beta value") +
    ylab("Gene expression (log counts)") +
    ggtitle(gene)
}
fig <- plot_grid(
  tcga.plots[["HLA-A"]], tcga.plots[["HLA-B"]], tcga.plots[["HLA-C"]],
  tcga.plots[["TAP1"]], tcga.plots[["TAP2"]], tcga.plots[["B2M"]], 
  tcga.plots[["WNT5A"]], tcga.plots[["TNFSF9"]], tcga.plots[["CIITA"]],
  # cis.plots[["HLA-A"]], cis.plots[["HLA-B"]], cis.plots[["HLA-C"]],
  # cis.plots[["TAP1"]], cis.plots[["TAP2"]], cis.plots[["B2M"]], 
  # cis.plots[["WNT5A"]], cis.plots[["TNFSF9"]], cis.plots[["CIITA"]],
  nrow=3,
  labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i") #, "j", "k", "l")
)
save_plot(paste0(figdir, "figS6.pdf"), fig, nrow = 3, ncol=2, base_width = fig.width/2, base_height = fig.width/2.5)

# fig <- plot_grid(
#   cis.plots[["HLA-A"]], cis.plots[["HLA-B"]], cis.plots[["HLA-C"]],
#   cis.plots[["TAP1"]], cis.plots[["TAP2"]], cis.plots[["B2M"]],
#   nrow=6,
#   labels = c("a", "b", "c", "d", "e", "f")
# )
# save_plot(paste0(figdir, "figXb.pdf"), fig, nrow = 3, ncol=2, base_width = 9)

##############################################################################
# Figure S7
# Methylation patterns over above genes
##############################################################################
message("Plotting Figure S7")

fig <- plot_grid(
  plotMethylForGene("HLA-A") + guides(color=F) + scale_color_manual(values = plotcols), 
  plotMethylForGene("HLA-B") + guides(color=F) + scale_color_manual(values = plotcols), 
  plotMethylForGene("HLA-C") + guides(color=F) + scale_color_manual(values = plotcols), 
  plotMethylForGene("TAP1") + guides(color=F) + scale_color_manual(values = plotcols), 
  plotMethylForGene("B2M") + guides(color=F) + scale_color_manual(values = plotcols), 
  plotMethylForGene("TNFSF9") + guides(color=F) + scale_color_manual(values = plotcols), 
  nrow=3
)
save_plot(paste0(figdir, "figS7.pdf"), fig, nrow = 3, ncol=2, base_width = fig.width/2)


##############################################################################
# Figure S8
# Cytokine:receptor analysis
##############################################################################
message("Plotting Figure S8")

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

fig <- plot_grid(f1,f2,f3,f4,f5, ncol = 3, labels = c('a','b','c','d','e'))
save_plot(paste0(figdir, "figS8.pdf"), fig, base_width = fig.width, base_height = fig.width)

##############################################################################
# Figure S9
# Boxplots of immune checkpoints with CIS tissue
##############################################################################
message("Plotting Figure S7")

l <- lapply(gene.lists$checkpoints, function(gene) {
  # if(gene %in% rownames(gdata)) {
  #   return(data.frame(gene=gene, val=as.numeric(gdata[gene,]), Outcome=gpheno$Outcone, index=1:dim(gdata)[2]))
  # } else {
  if(gene %in% rownames(gdata.v)) {
    return(data.frame(gene=gene, val=as.numeric(gdata.v[gene,]), Outcome=gpheno.v$Outcone, index=1:dim(gdata.v)[2]))
  }
  # }
})
plotdata <- do.call('rbind', l)

fig.main <- ggplot(plotdata, aes(x=index, y=val, color=Outcome)) + 
  geom_point() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~gene, scales = 'free') +
  ylab('Gene expression') +
  guides(color=FALSE) 

# Re-use legend defined in figure 2:
fig <- plot_grid(fig.main, leg, ncol = 1, rel_heights = c(15,1))

save_plot(paste0(figdir, "figS9.pdf"), fig, base_width = fig.width, base_height = fig.width*0.85)

##############################################################################
# Sup. Data 1
# Table of all lesions and profiling modalities
##############################################################################
message("Plotting Table S1")

pheno.public <- pheno[,c(1,2,3,5,8,10,12,13,17,20,25,28,29,30)]
# Sort properly by sample number:
pheno.public$tmp <- str_pad(str_extract(pheno.public$SampleID, "[0-9]+"), width = 3, pad = '0')
pheno.public <- pheno.public[order(pheno.public$tmp),]
pheno.public$tmp <- NULL
WriteXLS('pheno.public', ExcelFileName = paste0(supdir, "TableS1.xlsx"))

##############################################################################
# Sup. Data 2
# Danaher and MethylCIBERSORT decomposition data, with p-values for P vs R
##############################################################################
message("Plotting Table S2")

df1 <- gdata.danaher.t
df1$patient <- factor(gpheno.pair$Patient.Number)
df1$Outcome <- gpheno.pair$Outcome
df1 <- df1[,which(!apply(df1, 2, function(x) {all(is.na(x))}))]
colnames(df1) <- make.names(colnames(df1))
cols <- colnames(df1)
p.PvR <- lapply(cols, function(x){
  if(x %in% c('patient', 'Outcome')) {return(NA)}
  f <- as.formula(paste0(make.names(x), " ~ Outcome + (1 | patient)"))
  return(compare.fn(f, data = df1)$p)
})
p.PvR <- unlist(p.PvR)
names(p.PvR) <- cols
mean.p <- t(data.frame(apply(df1[which(df1$Outcome == 'Progression'),], 2, function(x) {mean(as.numeric(x))})))
rownames(mean.p) <- "mean.p"
df1 <- rbind(df1, mean.p)
mean.r <- t(data.frame(apply(df1[which(df1$Outcome == 'Regression'),], 2, function(x) {mean(as.numeric(x))})))
rownames(mean.r) <- "mean.r"
df1 <- rbind(df1, mean.r)
df1 <- rbind(df1, t(data.frame(p.PvR)))


df2 <- data.frame(methCS[,1:11])
df2$patient <- factor(mpheno$Patient_ID)
df2$Outcome <- factor(mpheno$Sample_Group)
cols <- colnames(df2)
p.PvR <- lapply(cols, function(x){
  if(x %in% c('patient', 'Outcome')) {return(NA)}
  f <- as.formula(paste0(make.names(x), " ~ Outcome + (1 | patient)"))
  # Exclude control samples here:
  return(compare.fn(f, data = df2[which(df2$Outcome %in% c("Progressive", "Regressive")),])$p)
})
p.PvR <- unlist(p.PvR)
names(p.PvR) <- cols
mean.p <- t(data.frame(apply(df2[which(df2$Outcome == 'Progressive'),], 2, function(x) {mean(as.numeric(x))})))
rownames(mean.p) <- "mean.p"
df2 <- rbind(df2, mean.p)
mean.r <- t(data.frame(apply(df2[which(df2$Outcome == 'Regressive'),], 2, function(x) {mean(as.numeric(x))})))
rownames(mean.r) <- "mean.r"
df2 <- rbind(df2, mean.r)
df2 <- rbind(df2, t(data.frame(p.PvR)))

WriteXLS(list(df1, df2), row.names = T,
         SheetNames = c("Danaher GXN deconvolution", "MethylCIBERSORT deconvolution"), ExcelFileName = paste0(supdir, "TableS2.xlsx"))

##############################################################################
# Sup. Data 3
# Lists of genes used in this analysis
# Include general immune genes (used for stroma PvsR comparison) and immunomodulators
##############################################################################
message("Plotting Table S3")

df1 <- data.frame(Gene.Name=sort(gene.lists$all.immune))
df2 <- data.frame(Gene.Name=sort(gene.lists$hla.assoc))
df3 <- data.frame(Gene.Name=sort(gene.lists$checkpoints))
df4 <- rbind(
  data.frame(Gene.Name=sort(gene.lists$cytokines.pro), Type='Pro-inflammatory'),
  data.frame(Gene.Name=sort(gene.lists$cytokines.anti), Type='Anti-inflammatory')
)
df5 <- gene.lists$chemokines
WriteXLS(
  list(df1, df2, df3, df4, df5), 
  ExcelFileName = paste0(supdir, "TableS3.xlsx"), 
  SheetNames = c('All immune genes', 'Antigen presentation genes', 'Immunomodulators', 'Inflammatory cytokines', 'Chemokine-receptor pairs'), col.names = F
)

##############################################################################
# Sup. Data 4
# Table of all genomic/gxn/methylation changes in immune genes
##############################################################################
message("Plotting Table S4")

# Summary of GXN and meth:
gs <- data.frame(t(gdata.zs.v[gene.lists$hla.assoc[which(gene.lists$hla.assoc %in% rownames(gdata.zs.v))], match(pheno$SampleID, colnames(gdata.zs.v))]))
rownames(gs) <- pheno$SampleID
ms <- data.frame(t(mdata.zs[gene.lists$hla.assoc[which(gene.lists$hla.assoc %in% rownames(mdata.zs))], match(pheno$SampleID, colnames(mdata.zs))]))
rownames(ms) <- pheno$SampleID
# Analagous table of genomic changes
geno.sum <- matrix(NA, nrow = dim(pheno)[1], ncol=length(gene.lists$hla.assoc))
for(i in 1:dim(pheno)[1]) {
  if(!pheno$Whole.Genome.Sequencing[i]){ next }
  for(j in 1:length(gene.lists$hla.assoc)) {
    sel <- which(as.character(mhc.muts$gene) == gene.lists$hla.assoc[j] & as.character(mhc.muts$sampleID) == pheno$SampleID[i])
    if(length(sel) > 0){
      geno.sum[i,j] <- paste(mhc.muts$type[sel], collapse="|")
    } else {
      geno.sum[i,j] <- 0
    }
  }
}
rownames(geno.sum) <- pheno$sampleID
colnames(geno.sum) <- gene.lists$hla.assoc
geno.sum <- data.frame(geno.sum)
mhc.muts <- mhc.muts[order(mhc.muts$patient, mhc.muts$gene, mhc.muts$class, mhc.muts$type),]

WriteXLS(list(pheno.public, gs, ms, geno.sum, mhc.muts), row.names = T,
         SheetNames = c('Summary', 'GXN', 'Meth', 'Mutations', "MutationDetail"), ExcelFileName = paste0(supdir, "TableS4.xlsx"))



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
