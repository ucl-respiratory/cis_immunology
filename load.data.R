# Load necessary datasets.
#
# This script does assume significant pre-processing has already been done. It is shared for information and is not straightforward to run.
# The user must perform the following:
#   Download and process raw gene expression, methylation and whole genome data as described in this manuscript: https://doi.org/10.1038/s41591-018-0323-0
#   Save gene expression data in data/gdata.RData
#   Save methylation data in data/mdata.RData
#   Download matched stromal data from GEO as described in the accompanying manuscript
#   Call mutations from WGS data as described previously (https://doi.org/10.1038/s41591-018-0323-0), save results as a data frame muts.all in data/wgsdata.Rdata
#   Run optitype on germline blood samples to establish patient HLA types: https://github.com/FRED-2/OptiType
#   Run LOHHLA on WGS samples - see https://bitbucket.org/mcgranahanlab/lohhla/src/master/
#   Run netMHC to call neoantigens from WGS samples following methods described here: https://doi.org/10.1038/s41586-019-1032-7
#
# Note that the above steps are not straightforward and do require significant bioinformatic expertise. 
# However, all methods above are published and well-established.
#
# Data is stored in data frames as follows:
#   pheno - details of each sample, and which analyses were performed. Includes some calculated metrics as below.
#   gdata - gene expression dataset with rows as genes, columns as sample names (n=51). Combines illumina and affymetrix data using ComBat.
#   gdata.d, gdata.v - separate illumina/affymetrix datasets used to study genes not present in one or the other.
#   gdata.s - stromal gene expression dataset
#   mdata - methylation 450k dataset with rows as probes, columns as sample names (n=87)
#   muts.all - details of all mutations passing filters. Includes CN deletions derived from LOHHLA analysis.
#   ihc - manually quantified IHC data of CD4+, CD8+, FOXP3+ cells
#   gdata.danaher.v, gdata.danaher.t, gdata.danaher.s - Danaher analysis of validation set / Tissue / Stroma
#   gm.tcga.pheno / gm.tcga.gdata / gm.tcga.mdata / gm.tcga.tcga.mdata.genes - TCGA gxn/methylation data 
#   methCS, methCS.tcga - methylCIBERSORT calls for CIS data / TCGA data
#   
# pheno additionally includes some calculated metrics:
#   hla.a.loss (and b,c) - whether HLA LOH was called using LOHHLA
#   neoantigen counts, separated into weak and strong
#
# The pheno table is included in the resources folder of this repository. It is included as supplementary data with the paper.
# gdata and mdata can be downloaded from GEO

####################################################################################
# Combine datasets
# Use existing data structures and rearrange to match master pheno described above
####################################################################################
message("Combining datasets")

library(gdata)
library(stringr)
library(ChAMP)
pheno <- read.xls('data/pheno.xls', stringsAsFactors=F)

# Fix date formats
for(i in 1:dim(pheno)[1]) {
  if(grepl("[0-9]+/[0-9]+/[0-9]+", pheno$Biopsy.Date[i])) {
    pheno$Biopsy.Date[i] <- as.character(as.Date(pheno$Biopsy.Date[i], format = "%d/%m/%Y"))
  }
}

# Convert 'YES/NO' column to booleans for ease of referencing
pheno$Gene.expression <- pheno$Gene.expression == 'YES'
pheno$Methylation <- pheno$Methylation == 'YES'
pheno$Whole.Genome.Sequencing <- pheno$Whole.Genome.Sequencing == 'YES'
pheno$Stroma.GXN <- pheno$Stroma.GXN == 'YES'
pheno$IHC <- pheno$IHC == 'YES'
pheno$HandE <- pheno$HandE == 'YES'
pheno$IHC_K <- pheno$IHC_K == 'YES'
pheno$Nanostring <- pheno$Nanostring == 'YES'

pheno$Sample.Number..Meth. <- make.names(pheno$Sample.Number..Meth.)
pheno$progression <- ifelse(pheno$Outcome == 'Progression', 1, 0)

# Load gene expression and methylation data
load('data/gdata.RData')
load('data/mdata.RData')

# Merge pheno data
pheno$Gender <- NA
pheno$COPD <- NA
pheno$Pack.years <- NA
pheno$Age.at.bronchoscopy <- NA
pheno$Previous.lung.cancer <- NA

pheno$Gender[match(gpheno$name, pheno$Sample.Number..GXN.)] <- as.character(gpheno$Gender)
pheno$Gender[match(mpheno$Sample_Name, pheno$Sample.Number..Meth.)] <- as.character(mpheno$Gender)
pheno$COPD[match(gpheno$name, pheno$Sample.Number..GXN.)] <- as.character(gpheno$COPD)
pheno$COPD[match(mpheno$Sample_Name, pheno$Sample.Number..Meth.)] <- as.character(mpheno$COPD)
pheno$Pack.years[match(gpheno$name, pheno$Sample.Number..GXN.)] <- as.numeric(gpheno$Pack.years)
pheno$Pack.years[match(mpheno$Sample_Name, pheno$Sample.Number..Meth.)] <- as.numeric(mpheno$Pack.years)
pheno$Age.at.bronchoscopy[match(gpheno$name, pheno$Sample.Number..GXN.)] <- as.numeric(gpheno$Age.at.specimen.profiled)
pheno$Age.at.bronchoscopy[match(mpheno$Sample_Name, pheno$Sample.Number..Meth.)] <- as.numeric(mpheno$Age.at.specimen.collected)
pheno$Previous.lung.cancer[match(gpheno$name, pheno$Sample.Number..GXN.)] <- as.character(gpheno$Prev.History.of.LC)
pheno$Previous.lung.cancer[match(mpheno$Sample_Name, pheno$Sample.Number..Meth.)] <- as.character(mpheno$Previous.Lung.CA.........No..0..Yes..1.)

pheno$COPD[which(pheno$COPD == "Y")] <- "YES"
pheno$COPD[which(pheno$COPD == "N")] <- "NO"

pheno$Previous.lung.cancer <- pheno$Previous.lung.cancer == '1'

# Match additional data for patients not previously published
clinical <- read.xls('resources/clinical_summaries.xls', stringsAsFactors=F)
clinical$ID <- as.numeric(clinical$ID)
clinical$Sex <- substr(clinical$Sex, 1, 1)
clinical$uuid <- str_pad(clinical$ID, 3,pad = '0')

x <- pheno
x$uuid <- str_pad(x$Patient.Number, 3, pad='0')
x <- merge(x, clinical, by='uuid', all.x=T)

# Fill in gaps
sel <- which(is.na(x$Gender) & !is.na(x$Sex))
x$Gender[sel] <- x$Sex[sel]
sel <- which(is.na(x$Pack.years) & !is.na(x$PackYears))
x$Pack.years[sel] <- x$PackYears[sel]
x$age2 <- as.numeric(floor((as.Date(x$Biopsy.Date) - as.Date(x$DOB))/365))
sel <- which(is.na(x$Age.at.bronchoscopy) & !is.na(x$age2))
x$Age.at.bronchoscopy[sel] <- x$age2[sel]
x$COPD.y <- ifelse(x$COPD.y, 'YES', 'NO')
sel <- which(is.na(x$COPD.x) & !is.na(x$COPD.y))
x$COPD.x[sel] <- x$COPD.y[sel]

# Find outstanding issues - some without smoking history:
# View(x[which((is.na(x$Gender) | is.na(x$Pack.years) | is.na(x$Age.at.bronchoscopy) | is.na(x$COPD.x)) & x$Outcome != 'Control'),])

# Copy back to pheno
sel <- match(pheno$SampleID, x$SampleID)
pheno$Gender <- x$Gender[sel]
pheno$Pack.years <- x$Pack.years[sel]
pheno$Age.at.bronchoscopy <- x$Age.at.bronchoscopy[sel]
pheno$COPD <- x$COPD.x[sel]

pheno$QuitDate[which(pheno$QuitDate == '')] <- NA
pheno$QuitDate <- as.Date(pheno$QuitDate)
pheno$YearsSinceQuitting <- as.numeric(floor((as.Date(pheno$Biopsy.Date) - pheno$QuitDate) / 365))
sel <- which(pheno$YearsSinceQuitting < 0)
if(length(sel) > 0) {
  pheno$YearsSinceQuitting[sel] <- NA
  pheno$SmokingStatus[sel] <- 'Current'
}

# Manual corrections (see private/pheno_consistency_check.R):
# Patient 6
pheno$COPD[which(pheno$Patient.Number == 6)] <- 'NO'
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 6)] <- F
# Patient 11
pheno$COPD[which(pheno$Patient.Number == 11)] <- 'YES'
pheno$Pack.years[which(pheno$Patient.Number == 11)] <- 40
# Patient 12 
pheno$COPD[which(pheno$Patient.Number == 12)] <- 'NO'
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 12)] <- T
# Patient 26
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 26)] <- T
# Patient 76 
pheno$COPD[which(pheno$Patient.Number == 76)] <- 'NO'
pheno$Pack.years[which(pheno$Patient.Number == 76)] <- 60
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 76)] <- T
# Patient 89
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 89)] <- F
# Patient 95
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 95)] <- F
# Patient 102
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 102)] <- F
# Patient 103 
pheno$COPD[which(pheno$Patient.Number == 103)] <- 'NO'
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 103)] <- F
# Patient 104
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 104)] <- T
# Patient 105 
pheno$Pack.years[which(pheno$Patient.Number == 105)] <- 9
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 105)] <- F
# Patient 108
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 108)] <- T
# Patient 109
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 109)] <- F
# Patient 110
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 110)] <- F
# Patient 111 
pheno$COPD[which(pheno$Patient.Number == 111)] <- 'YES'
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 111)] <- F
# Patient 114
pheno$Pack.years[which(pheno$Patient.Number == 114)] <- 90
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 114)] <- T
# Patient 119
pheno$COPD[which(pheno$Patient.Number == 119)] <- 'NO'
pheno$Pack.years[which(pheno$Patient.Number == 119)] <- 36
# Patient 123
pheno$Pack.years[which(pheno$Patient.Number == 123)] <- 60
# Patient 124 
pheno$COPD[which(pheno$Patient.Number == 124)] <- 'NO'
pheno$Previous.lung.cancer[which(pheno$Patient.Number == 124)] <- F
# Others needing manual input:
pheno$Previous.lung.cancer[which(pheno$Patient.Number %in% c(25, 131, 134, 135, 136, 138))] <- T
pheno$Previous.lung.cancer[which(pheno$Patient.Number %in% c(88, 125, 132, 139, 140, 141, 142, 144))] <- F


# Re-label data to match pheno
colnames(gdata) <- pheno$SampleID[match(colnames(gdata), pheno$Sample.Number..GXN.)]
gpheno$sampleID <- pheno$SampleID[match(gpheno$name, pheno$Sample.Number..GXN.)]
colnames(gdata.d) <- pheno$SampleID[match(colnames(gdata.d), pheno$Sample.Number..GXN.)]
gpheno.d$sampleID <- pheno$SampleID[match(gpheno.d$name, pheno$Sample.Number..GXN.)]
colnames(gdata.v) <- pheno$SampleID[match(colnames(gdata.v), pheno$Sample.Number..GXN.)]
gpheno.v$sampleID <- pheno$SampleID[match(gpheno.v$name, pheno$Sample.Number..GXN.)]
colnames(mdata) <- pheno$SampleID[match(colnames(mdata), pheno$Sample.Number..Meth.)]
mpheno$sampleID <- pheno$SampleID[match(make.names(mpheno$Sample_Name), pheno$Sample.Number..Meth.)]

# Add corresponding gdata stroma data
# Here we load tissue and stroma CEL files and normalise together 
# Thus values will be different to gdata.v, even though these are the same samples
# Hence we store as gdata.pair.t and gdata.pair.s for tissue and stroma data respectively
load("data/gdata.paired.RData")
gpheno.pair <- pheno[which(pheno$Stroma.GXN),]
gdata.pair.t <- gdata.paired[,paste0('T_', gpheno.pair$Sample.Number..GXN.)]
colnames(gdata.pair.t) <- gpheno.pair$SampleID
gdata.pair.s <- gdata.paired[,paste0('S_', gpheno.pair$Sample.Number..Stroma.)]
colnames(gdata.pair.s) <- gpheno.pair$SampleID

# Copy some variables from pheno to gpheno
for(var in c('SmokingStatus', 'YearsSinceQuitting')) {
  gpheno[,var] <- pheno[match(gpheno$sampleID, pheno$SampleID), var]
  gpheno.d[,var] <- pheno[match(gpheno.d$sampleID, pheno$SampleID), var]
  gpheno.v[,var] <- pheno[match(gpheno.v$sampleID, pheno$SampleID), var]
  gpheno.pair[,var] <- pheno[match(gpheno.pair$SampleID, pheno$SampleID), var]
}

# Load gene lists
load('resources/gene.lists.RData')


####################################################################################
# Load WGS data
####################################################################################
load('data/wgsdata.RData')

# Load annotated pheno including LOHHLA and neoantigen results
# TODO: include LOHHLA annotation code
load('data/wgsPhenoAnnotated.RData')

# Annotate with HLA type:
pheno$hla.type <- NA
sel <- which(pheno$Whole.Genome.Sequencing)
for(i in sel) {
  samp <- pheno$Sample.Number..WGS.[i]
  f <- paste0("~/Scratch/cshpc/cis/data/optitype_output/", samp, "/Scratch/", samp, "/", samp, ".optitype.output.for.netMHC.txt")
  if(file.exists(f)) {
    x <- readLines(f)
    pheno$hla.type[i] <- paste(x, collapse="|")
  } else {
    message(paste("File not found:", f))
  }
}



####################################################################################
# Combine mutation and copy number data
# Add CNAs to muts.all
# Include deletion (CN = 0), gain (CN/ploidy > 1), amplification (CN/ploidy > 2)
# Also LOH, defined as minor CN in tumour == 0 (and !=0 in normal)
####################################################################################
message("Combining mutation and copy number data")

cna.summary.list$ploidy <- as.numeric(wgs.pheno$ploidy[match(cna.summary.list$sample, wgs.pheno$name)])
# cna.summary.list$ploidy <- as.numeric(pheno$ploidy[match(cna.summary.list$sample, pheno$Sample.Number..WGS.)])
cna.summary.list$cn_ploidy_corrected <- cna.summary.list$total.copy.number.inTumour / cna.summary.list$ploidy
cna.summary.list$type <- 'Normal'
cna.summary.list2 <- cna.summary.list
cna.summary.list$type[which(cna.summary.list$cn_ploidy_corrected == 0)] <- "Deletion"
cna.summary.list$type[which(cna.summary.list$cn_ploidy_corrected > 0 & cna.summary.list$cn_ploidy_corrected < 1)] <- "Loss"
cna.summary.list$type[which(cna.summary.list$cn_ploidy_corrected > 1)] <- "Gain"
cna.summary.list$type[which(cna.summary.list$cn_ploidy_corrected > 2)] <- "Amplification"
# Add LOH
# Need to be careful as some samples have LOH and loss together
cna.summary.list2$type[which(cna.summary.list2$minor.copy.number.inTumour == 0 & cna.summary.list2$minor.copy.number.inNormal != 0)] <- "LOH"
cna.summary.list <- rbind(cna.summary.list, cna.summary.list2[which(cna.summary.list2$type == 'LOH'),])

# For each non-normal region, find relevant genes and add to muts.all
# Takes a long time -> cache results
cn.cache.file <- "data/cna.cache.RData"
if(file.exists(cn.cache.file)) {
  load(cn.cache.file)
} else {
  print('Generating CNA cache')
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(Homo.sapiens)
  genome <- BSgenome.Hsapiens.UCSC.hg19
  gene_coord <- genes(Homo.sapiens, columns=c("GENEID", "SYMBOL"))
  source("utility_functions/parallel.setup.R")
  cnas.all <- foreach(i = 1:dim(cna.summary.list)[1], .combine=rbind) %do% {
    if(i %% 100 == 0) { print(paste(i, "/", dim(cna.summary.list)[1])) }
    # if(i %% 10 == 0){print(paste(i, "/", dim(cna.summary.list)[1]))}
    my.na <- data.frame(patient=NA, gene=NA, class=NA, type=NA, ref=NA, alt=NA, chr=NA, start=NA, end=NA, filters=NA, asmd=NA, clpm=NA, vaf=NA,ref.reads=NA,alt.reads=NA,tumour.reads=NA,depth=NA,exonic=NA,protein.change=NA,cds.mut=NA,filters.passed=T,mid=NA,translocation.partner=NA,chr2=NA,start2=NA,end2=NA)
    if(cna.summary.list$type[i] == 'Normal') { return(my.na) }
    genes <- gene_coord[which(
      as.character(seqnames(gene_coord)) == cna.summary.list$Chromosome[i] &
        start(gene_coord) < cna.summary.list$chromEnd[i] & 
        end(gene_coord) > cna.summary.list$chromStart[i]
    ),]
    genes <- as.character(genes$SYMBOL)
    n = length(genes)
    if(n == 0) { return(my.na) }
    df <- data.frame(
      patient=rep(cna.summary.list$sample[i], n),
      gene=genes,
      class=cna.summary.list$type[i],
      type=cna.summary.list$type[i],
      ref=NA,
      alt=NA,
      chr=gsub("chr", "", cna.summary.list$Chromosome[i]),
      start=cna.summary.list$chromStart[i],
      end=cna.summary.list$chromEnd[i],
      filters="PASS",
      asmd=140,
      clpm=0,
      vaf=NA,ref.reads=NA,alt.reads=NA,tumour.reads=NA,depth=NA,exonic=NA,protein.change=NA,cds.mut=NA,filters.passed=T,mid=NA,translocation.partner=NA,chr2=NA,start2=NA,end2=NA
    )
    df <- df[which(!is.na(df$gene)),]
    # muts.all <- rbind(muts.all, df)
    return(df)
  }
  cnas.all <- cnas.all[which(!is.na(cnas.all$gene)),]
  
  
  # Genes which appear twice must have breakpoints within the gene. Relabel these as breaks:
  for(pt in unique(cnas.all$patient)){
    sel.pt <- which(cnas.all$patient == pt)
    dup.genes <- unique(cnas.all[sel.pt,]$gene[which(duplicated(cnas.all[sel.pt,]$gene))])
    # Change the first one to "CN break" type and delete the rest
    for(gene in dup.genes) {
      sel <- which(cnas.all$patient == pt & cnas.all$gene == gene)
      # Skip Loss/LOH combos - we expect these
      if(length(sel) == 2 & 'LOH' %in% cnas.all$type[sel] & 'Loss' %in% cnas.all$type[sel]) {next}
      cnas.all$type[sel[1]] <- "CN break"
      cnas.all$class[sel[1]] <- "CN break"
      cnas.all <- cnas.all[-sel[2:length(sel)],]
    }

  }
  
  save(cnas.all, file=cn.cache.file)
}


muts.all <- rbind(muts.all, cnas.all)

####################################################################################
# LOHHLA
####################################################################################
message("Adding LOHHLA")

# Add LOHHLA results as LOH events to muts.all 
for(i in 1:dim(wgs.pheno)[1]){
  # dels <- cnas.dels[[wgs.pheno$name[i]]]
  dels <- data.frame(
    seqnames=character(0), start=numeric(0), end=numeric(0), width=numeric(0), strand=character(0), GENEID=numeric(0), SYMBOL=character(0), sample=character(0)
  )
  if(wgs.pheno$hla.a.loss[i]) {
    df <- data.frame(seqnames='chr6', start=as.integer(29909037), end=29913661, width=4624, strand='+', GENEID=3105, SYMBOL='HLA-A', row.names = sample(10**10, 1), sample=wgs.pheno$name[i])
    dels <- rbind(dels, df)
  }
  if(wgs.pheno$hla.b.loss[i]) {
    df <- data.frame(seqnames='chr6', start=31321649, end=31324965, width=3316, strand='-', GENEID=3106, SYMBOL='HLA-B', row.names = sample(10**10, 1), sample=wgs.pheno$name[i])
    dels <- rbind(dels, df)
  }
  if(wgs.pheno$hla.c.loss[i]) {
    df <- data.frame(seqnames='chr6', start=31236526, end=31239907, width=3381, strand='-', GENEID=3107, SYMBOL='HLA-C', row.names = sample(10**10, 1), sample=wgs.pheno$name[i])
    dels <- rbind(dels, df)
  }
  # cnas.dels[[wgs.pheno$name[i]]] <- dels
  if(i == 1) {
    lohhla.dels <- dels
  } else {
    lohhla.dels <- rbind(lohhla.dels, dels)
  }
}
# Append LOHHLA results to muts.all as LOH events
df <- data.frame(
  patient=lohhla.dels$sample,
  gene=as.character(lohhla.dels$SYMBOL),
  class='LOH',
  type='LOH',
  ref=NA,
  alt=NA,
  chr=gsub("chr", "", lohhla.dels$seqnames),
  start=lohhla.dels$start,
  end=lohhla.dels$end,
  filters="PASS",
  asmd=140,
  clpm=0,
  vaf=NA,ref.reads=NA,alt.reads=NA,tumour.reads=NA,depth=NA,exonic=NA,protein.change=NA,cds.mut=NA,filters.passed=T,mid=NA,translocation.partner=NA,chr2=NA,start2=NA,end2=NA
)
muts.all <- rbind(muts.all, df)

####################################################################################
# Load Neoantigen data
# Data profiled using netMHCpan4.0
# Store all neoantigen data, and add counts to wgs.pheno
####################################################################################
message("Adding neoantigens")

load('data/neoantigen.data.RData')
wgs.pheno$neoants.weak <- NA
wgs.pheno$neoants.weak[match(names(neoantigen.data), wgs.pheno$name)] <- unlist(lapply(neoantigen.data, function(x) {length(which(x$is_weak))}))
wgs.pheno$neoants.strong <- NA
wgs.pheno$neoants.strong[match(names(neoantigen.data), wgs.pheno$name)] <- unlist(lapply(neoantigen.data, function(x) {length(which(x$is_strong))}))

# Annotate muts.all with a flag identifying strong/weak neoantigens
neos.all <- do.call('rbind', lapply(neoantigen.data, function(x) {
  x[,c("Peptide", "ID", "sample", "full.id", "is_mt", "is_strong", "is_weak", "best.rank", "best.aff", "best.DAI", "best.rank.allele", "best.aff.allele", "best.DAI.allele")]
}))
# Extract genomic position info from neoantigen full_id data
neos.all$mid2 <- unlist(lapply(as.character(neos.all$full.id), function(x) {
  unlist(strsplit(x, '--'))[3]
}))
# Coerce mutation mid into a matching format
muts.all$mid2 <- unlist(lapply(muts.all$mid, function(x) {
  if(is.na(x)) {return(NA)}
  if(!str_detect(x, "[0-9XY]+-[0-9]+-[0-9]+-[A-Z]+-[A-Z]+")) {return(NA)} # Removes Rearrangements which have a different ID format
  s <- unlist(strsplit(x, "-"))
  if(nchar(s[4]) == 1 & nchar(s[5]) == 1) {
    
  } else {
    if(nchar(s[4]) == 1) {
      s[4] <- '-'
    } else {
      s[4] <- str_sub(s[4], start=2)
    }
    if(nchar(s[5]) == 1) {
      s[5] <- '-'
    } else {
      s[5] <- str_sub(s[5], start=2)
    }
    if(nchar(s[4]) > nchar(s[5])) {
      s[2] <- as.numeric(s[2]) + 1
    }
  }
  
  return(paste(s[1], s[2], s[4], s[5], sep=":"))
}))
# All neoantigens should be represented in mid - but note multiple potential neoantigens per mut
strong_mids <- unique(neos.all$mid2[which(neos.all$is_strong)])
weak_mids <- unique(neos.all$mid2[which(neos.all$is_weak)])

if(!all(strong_mids %in% muts.all$mid2)) {message("WARNING: not all strong neoantigen mutations are in the mutation database")}
if(!all(weak_mids %in% muts.all$mid2)) {message("WARNING: not all weak neoantigen mutations are in the mutation database")}

muts.all$is.neo.strong <- 0
muts.all$is.neo.strong[which(muts.all$mid2 %in% strong_mids)] <- 1

muts.all$is.neo.weak <- 0
muts.all$is.neo.weak[which(muts.all$mid2 %in% weak_mids)] <- 1


####################################################################################
# Clonality
# Define a mutation as clonal if CCF is around 1
####################################################################################
message("Adding clonality estimates")

# First, find the local copy number of each mutation (tumour and normal)
muts.all$CNt <- NA
muts.all$CNn <- NA
for(i in 1:dim(cna.summary.list)[1]){
  sel <- which(paste0("chr", as.character(muts.all$chr)) == cna.summary.list$Chromosome[i] & muts.all$start <= cna.summary.list$chromEnd[i] & muts.all$end >= cna.summary.list$chromStart[i] & muts.all$patient == cna.summary.list$sample[i])
  muts.all$CNt[sel] <- cna.summary.list$total.copy.number.inTumour[i]
  muts.all$CNn[sel] <- cna.summary.list$total.copy.number.inNormal[i]
}
# Now apply CCF function to each
# (warning - takes a long time! Cache the CCF data frame)
# ccf.cache.file <- 'data/ccf.RData'
# if(file.exists(ccf.cache.file)){
#   message('Loading CCF cache')
#   load(ccf.cache.file)
# } else {
#   message('Generating CCF file (may take some time)')
#   source('utility_functions/parallel.setup.R')
#   m.tmp <- muts.all[,c("alt.reads", "tumour.reads", "patient", "CNt")]
#   t1 <- Sys.time()
#   ccf <- foreach(i=1:dim(m.tmp)[1], .combine = rbind, .packages = character(0)) %dopar% {
#     if(any(is.na(m.tmp[i, c('alt.reads', 'tumour.reads', 'CNt')]))) { 
#       return(data.frame(ccf.lower.ci=NA,ccf.est=NA,ccf.upper.ci=NA,prob.subclonal=NA,prob.clonal=NA,is.clonal=NA,is.clonal.byprob=NA)) 
#     }
#     source('utility_functions/absolute.cancer.cell.fraction.R')
#     c <- absolute.cancer.cell.fraction(
#       n.alt=m.tmp$alt.reads[i], 
#       depth=m.tmp$tumour.reads[i], 
#       purity=wgs.pheno$purity[match(m.tmp$patient[i], wgs.pheno$name)], 
#       local.copy.number=m.tmp$CNt[i]
#     )
#     return(data.frame(c))
#   }
#   t2 <- Sys.time()
#   print(t2 - t1)
#   # Make sure we can match back to muts.all, even if we change muts.all:
#   ccf$patient <- muts.all$patient
#   ccf$mid <- muts.all$mid
#   save(ccf, file=ccf.cache.file)
#   rm(m.tmp)
# }
# 
# muts.all$ccf <- ccf$ccf.est
# muts.all$is.clonal <- ccf$is.clonal

# New version runs much more efficiently:
source('utility_functions/calcCCF.R')
muts.all <- calcCCF(muts.all, pheno)

# Add sampleID to muts.all
muts.all$sampleID <- pheno$SampleID[match(muts.all$patient, pheno$Sample.Number..WGS.)]

# Add calculated WGS data
for(attr in c('purity', 'ploidy', 'burden', 'neoants.strong', 'neoants.weak', 
              'hla.loss', 'hla.a.loss', 'hla.b.loss', 'hla.c.loss', 'hla.lost.alleles')) {
  pheno[,attr] <- NA
  pheno[match(wgs.pheno$name, pheno$Sample.Number..WGS.), attr] <- wgs.pheno[,attr]
}


####################################################################################
# Add WGII data to pheno data frame
####################################################################################
message("Adding WGII")
# Need to reload copy number data in correct format:
load('data/cna.summary.list.RData')
getWgii <- function(data, ploidy=2, type='all'){
  # data should have columns chr, start, end, cn
  # For each chromosome (except X,Y) find the percentage for which CN !=2
  # Then take the mean across all chromosomes
  chrdata <- c()
  data$chr <- gsub('chr', '', as.character(data$chr))
  for(chr in unique(data$chr)){
    if(!(as.numeric(chr) %in% 1:22)){ next }
    segs <- data[which(data$chr == chr),]
    seglengths <- as.numeric(segs$end) - as.numeric(segs$start)
    if(type == 'all') {
      x <- sum(seglengths[which(segs$cn != ploidy)]) / sum(seglengths)
    }
    if(type == 'loss') {
      x <- sum(seglengths[which(segs$cn < ploidy)]) / sum(seglengths)
    }
    if(type == 'loh') {
      x <- sum(seglengths[which(segs$cn.minor.t == 0 & segs$cn.minor.n != 0)]) / sum(seglengths)
    }
    if(type == 'gain') {
      x <- sum(seglengths[which(segs$cn > ploidy)]) / sum(seglengths)
    }
    chrdata <- c(chrdata, x)
  }
  return(mean(chrdata))
}
pheno$wgii <- sapply(1:dim(pheno)[1], function(i) {
  x <- pheno$Sample.Number..WGS.[i]
  if(is.na(x)) {return(NA)}
  ploidy <- as.numeric(pheno$ploidy[i])
  df <- cna.summary.list[which(cna.summary.list$sample == x),]
  df <- data.frame(
    chr = df$Chromosome, start = df$chromStart, end = df$chromEnd, cn = df$total.copy.number.inTumour, stringsAsFactors = F
  )
  return(getWgii(df, ploidy))
})
# Note this measures LOH rather than any loss
pheno$wgii.loh <- sapply(1:dim(pheno)[1], function(i) {
  x <- pheno$Sample.Number..WGS.[i]
  if(is.na(x)) {return(NA)}
  ploidy <- as.numeric(pheno$ploidy[i])
  df <- cna.summary.list[which(cna.summary.list$sample == x),]
  df <- data.frame(
    chr = df$Chromosome, start = df$chromStart, end = df$chromEnd, cn = df$total.copy.number.inTumour, cn.minor.t = df$minor.copy.number.inTumour, cn.minor.n = df$minor.copy.number.inNormal, stringsAsFactors = F
  )
  return(getWgii(df, ploidy, type = 'loh'))
})
# Here we measure any loss (CN < ploidy)
pheno$wgii.loss <- sapply(1:dim(pheno)[1], function(i) {
  x <- pheno$Sample.Number..WGS.[i]
  if(is.na(x)) {return(NA)}
  ploidy <- as.numeric(pheno$ploidy[i])
  df <- cna.summary.list[which(cna.summary.list$sample == x),]
  df <- data.frame(
    chr = df$Chromosome, start = df$chromStart, end = df$chromEnd, cn = df$total.copy.number.inTumour, stringsAsFactors = F
  )
  return(getWgii(df, ploidy, type = 'loss'))
})
pheno$wgii.gain <- sapply(1:dim(pheno)[1], function(i) {
  x <- pheno$Sample.Number..WGS.[i]
  if(is.na(x)) {return(NA)}
  ploidy <- as.numeric(pheno$ploidy[i])
  df <- cna.summary.list[which(cna.summary.list$sample == x),]
  df <- data.frame(
    chr = df$Chromosome, start = df$chromStart, end = df$chromEnd, cn = df$total.copy.number.inTumour, stringsAsFactors = F
  )
  return(getWgii(df, ploidy, type = 'gain'))
})


####################################################################################
# z-score calculation
####################################################################################
message("Calculating z-scores")
# For gene expression and methylation data, calculate z-scores
# Base this on a reference calculated from TCGA control samples
# where z-score = (expression in sample - mean expression in reference) / standard deviation of expression in reference
source("utility_functions/runComBat.R")
gdata.ref <- tcga.gdata[,which(tcga.gpheno$progression == 0)]
# Batch adjust for a) all gene expression data b) validation set only
# This function also reduces the input dataset to overlapping genes only
c <- runComBat(gdata, gdata.ref)
gdata.1 <- data.frame(c[[1]])
gdata.ref.1 <- data.frame(c[[2]])
c <- runComBat(gdata.v, gdata.ref)
gdata.2 <- data.frame(c[[1]])
gdata.ref.2 <- data.frame(c[[2]])

# sel.ref <- which(gpheno$progression == 0)
ref.mean <- apply(gdata.ref.1, 1, mean)
ref.sd <- apply(gdata.ref.1, 1, sd)
gdata.zs <- apply(gdata.1, 2, function(x){
  (x - ref.mean) / ref.sd
})

# Repeat for the validation set only (Affy data) - use this for genes not expressed in the combined dataset
# sel.ref <- which(gpheno.v$progression == 0)
ref.mean <- apply(gdata.ref.2, 1, mean)
ref.sd <- apply(gdata.ref.2, 1, sd)
gdata.zs.v <- apply(gdata.2, 2, function(x){
  (x - ref.mean) / ref.sd
})

# Try same method for methylation
# Here we will consider average methylation over the gene (crude but reduces noise)
# Use control samples as our reference dataset
data('probe.features')
genes <- as.character(probe.features$gene[match(rownames(mdata), rownames(probe.features))])
sel <- which(genes != '')
mdata.genes <- aggregate(mdata[sel,], by=list(genes[sel]), FUN=mean)
rownames(mdata.genes) <- mdata.genes$Group.1
mdata.genes$Group.1 <- NULL

sel.ref <- which(mpheno$Sample_Group == "Control")
ref.mean <- apply(mdata.genes[,sel.ref], 1, mean)
ref.sd <- apply(mdata.genes[,sel.ref], 1, sd)
mdata.zs <- apply(mdata.genes, 2, function(x){
  (x - ref.mean) / ref.sd
})

####################################################################################
# Methylation SNP removal 
####################################################################################

data(hm450.manifest.hg38)
data(probe.features)
mvps <- rownames(hm450.manifest.hg38)[which(!hm450.manifest.hg38$MASK.general)]
mvps.immune <- rownames(probe.features)[which(probe.features$gene %in% gene.lists$all.immune)]
mvps.hla <- rownames(probe.features)[which(probe.features$gene %in% gene.lists$hla.assoc)]
# Remove SNP-associated
mvps.immune <- mvps.immune[which(mvps.immune %in% mvps)]
genes <- as.character(probe.features[mvps,]$gene)
mdata.nosnp <- mdata[mvps,]
mdata.genes.nosnp <- aggregate(mdata[mvps,], by=list(genes), FUN=function(x) {mean(x, na.rm=T)})
mdata.genes.nosnp <- mdata.genes.nosnp[which(mdata.genes.nosnp$Group.1 != ''),]
rownames(mdata.genes.nosnp) <- mdata.genes.nosnp$Group.1
mdata.genes.nosnp$Group.1 <- NULL
# Remove NA genes (with no probes):
sel.rm <- which(apply(mdata.genes.nosnp, 1, function(x) {any(is.na(x))}))
if(length(sel.rm) > 0){
  mdata.genes.nosnp <- mdata.genes.nosnp[-sel.rm,]
}

####################################################################################
# IHC data
# Load raw IHC data
####################################################################################
load('data/ihc.RData')
ihc$sampleID <- pheno$SampleID[match(ihc$Samples, pheno$Sample.Number..IHC.)]

####################################################################################
# TCGA data
# Data for lung squamous cell carcinoma (LUSC) was downloaded from TCGA using GenomicDataCommons software.
# Samples with overlapping methylation and gene expression data are saved to data/tcga.overlap.RData
# Methylation data is imputed using champ.impute from the ChAMP package to correct for missing probe values
# CNAs are also called from TCGA data using ASCAT, and summarised in a table (tcga.cn.table)
# LOHHLA calls from TCGA are also included.
####################################################################################
load("data/tcga.overlap.RData")
# Rename some variables:
gm.tcga.mdata.genes <- gm.tcga.mdata
gm.tcga.mdata <- gm.tcga.mdata.imp

load('data/lusc.cnas.by.gene.RData')

# TCGA LOHHLA calls (From Rachel)
load('data/tcga.lohhla.calls.RData')
tcga.lohhla <- tmp
gm.tcga.pheno$PatientBarcode <- substr(as.character(gm.tcga.pheno$submitter_id2), 1, 12)
gm.tcga.pheno <- merge(gm.tcga.pheno, tmp, by.x = 'PatientBarcode', by.y='row.names', all.x=T)


####################################################################################
# Danaher data
####################################################################################
source('utility_functions/do.danaher.R')
gdata.danaher.v <- do.danaher(gdata.v)
gdata.danaher.t <- do.danaher(gdata.pair.t)
gdata.danaher.s <- do.danaher(gdata.pair.s)

####################################################################################
# MethylCibersort data
####################################################################################
library(MethylCIBERSORT)
data('V2_Signatures')
sig <- Signatures$lung_NSCLC_squamous_cell_carcinoma_v2_Signature.txt
Prep.CancerType(Beta = mdata, Probes = sig$NAME, fname = "methCS_mdata")
Prep.CancerType(Beta = gm.tcga.mdata, Probes = sig$NAME, fname = "methCS_mdata_tcga")
write.table(sig, file = "methCS_Signature_LUSC", sep = "\t", row.names = FALSE, quote = FALSE)
# CIBERSORT code must be downloaded into the resources folder from here:
# https://cibersort.stanford.edu
source("resources/CIBERSORT.R")
set.seed(42)
methCS <- CIBERSORT('methCS_Signature_LUSC', 'methCS_mdata.txt', perm=100)
methCS.tcga <-CIBERSORT('methCS_Signature_LUSC', 'methCS_mdata_tcga.txt', perm=100) 

####################################################################################
# mIHC data
####################################################################################
mihc <- read.xls('data/mIHC_raw_data.xls', sheet = 1, header=F, stringsAsFactors = F)

pheno.mihc <- data.frame(
  Outcome = as.character(mihc[1,-1]),
  Batch = as.character(mihc[2,-1]),
  Patient = as.character(mihc[3, -1]),
  total.events = as.numeric(as.character(mihc[4,-1])),
  stroma.events = as.numeric(as.character(mihc[5,-1])),
  stroma.16p = as.numeric(as.character(mihc[6,-1])),
  
  cis.events = as.numeric(as.character(mihc[76,-1])),
  cis.16p = as.numeric(as.character(mihc[77,-1]))
)
# pheno.mihc$name <- paste0("S", 1:dim(pheno.mihc)[1])
pheno.mihc$Patient <- gsub("Pt", "Patient ", as.character(pheno.mihc$Patient))
pheno.mihc$name <- pheno$SampleID[match(pheno.mihc$Patient, pheno$Nanostring.ID)]
pheno.mihc$Patient <- pheno$Patient.Number[match(pheno.mihc$name, pheno$SampleID)]

rows.cis <- 78:146
mihc.cis.data <- apply(mihc[rows.cis,-1], c(1,2), as.numeric)
rownames(mihc.cis.data) <- mihc[rows.cis,1]
colnames(mihc.cis.data) <- pheno.mihc$name
# Tidy the rownames (they have Carcinoma1 in them!)
library(stringr)
rownames(mihc.cis.data) <- gsub("(Carcinoma1)", "", rownames(mihc.cis.data), fixed = T)
rownames(mihc.cis.data) <- gsub("Carcinoma1", "", rownames(mihc.cis.data), fixed = T)
rownames(mihc.cis.data) <- str_trim(rownames(mihc.cis.data))

rows.stroma <- 7:75
mihc.stroma.data <- apply(mihc[rows.stroma,-1], c(1,2), as.numeric)
rownames(mihc.stroma.data) <- mihc[rows.stroma,1]
colnames(mihc.stroma.data) <- pheno.mihc$name
rownames(mihc.stroma.data) <- str_trim(rownames(mihc.stroma.data))

# Collapse to common markers:
ctypes <- intersect(rownames(mihc.cis.data), rownames(mihc.stroma.data))
mihc.cis.data <- mihc.cis.data[ctypes,]
mihc.stroma.data <- mihc.stroma.data[ctypes,]

# Label cells we know about:
labs <- list(
  "CD45+" = "All CD45",
  "CD45+CD3+" = "All T cells", 
  "CD45+CD3+CD8+" = "CD8 T cells",
  "CD45+CD3+CD8+EOMES+PD1+" = "CD8 Exhausted T cells",
  "CD45+CD3+CD8+EOMES+PD1-" = "CD8 late effector T cells",
  "CD45+CD3+CD8+Grzb+" = "Cytotoxic T cells",
  "CD45+CD3+CD8-" = "CD4 T cells",
  "CD45+CD3+CD8-FoxP3+" = "Treg",
  "CD45+CD3-CD20+" = "B cells",
  # "CD45+CD3-CD20-CD11b+CD66b+" = "Myeloid cells",
  # "CD45+CD3-CD20-CD11b+CD66b-" = "Myeloid cells",
  "CD45+CD3-CD20-CD11b+CD66b+" = "Neutrophils",
  "CD45+CD3-CD20-CD11b+CD66b-HLADR+DCLamp+" = "Dendritic Cells",
  "CD45+CD3-CD20-CD11b+CD66b-HLADR+DCLamp-" = "Dendritic Cells",
  "CD45+CD3-CD20-CD11b+CD66b-CD163+CD68+" = "Macrophages",
  "CD45+CD3-CD20-CD11b+CD66b-CD163-CD68+" = "Macrophages"
)
# Export these for use in a sup. table:
df.mihc <- data.frame(V1=names(labs), V2=as.character(labs))
WriteXLS::WriteXLS(df.mihc, ExcelFileName = "results/mihc.markers.xlsx", AdjWidth = T)

x <- mihc.cis.data[which(rownames(mihc.cis.data) %in% names(labs)),]
mihc.cis.data2 <- aggregate(x, by=list(unlist(labs[rownames(x)])), FUN=sum)
rownames(mihc.cis.data2) <- mihc.cis.data2$Group.1
mihc.cis.data2$Group.1 <- NULL

x <- mihc.stroma.data[which(rownames(mihc.stroma.data) %in% names(labs)),]
mihc.stroma.data2 <- aggregate(x, by=list(unlist(labs[rownames(x)])), FUN=sum)
rownames(mihc.stroma.data2) <- mihc.stroma.data2$Group.1
mihc.stroma.data2$Group.1 <- NULL

# Now let's normalise everything to be a percent of total >16px events:
mihc.cis.norm <- mihc.cis.data
mihc.stroma.norm <- mihc.stroma.data
mihc.cis.norm2 <- mihc.cis.data2
mihc.stroma.norm2 <- mihc.stroma.data2
for(i in 1:dim(pheno.mihc)[1]) {
  mihc.cis.norm[,i] <- mihc.cis.norm[,i] / pheno.mihc$cis.16p[i]
  mihc.stroma.norm[,i] <- mihc.stroma.norm[,i] / pheno.mihc$stroma.16p[i]
  mihc.cis.norm2[,i] <- mihc.cis.norm2[,i] / pheno.mihc$cis.16p[i]
  mihc.stroma.norm2[,i] <- mihc.stroma.norm2[,i] / pheno.mihc$stroma.16p[i]
}

##########################################################
# VEP annotation of muts.all:
muts.all$mid3 <- paste0(muts.all$patient, "-", muts.all$mid)
vep.input <- data.frame(
  muts.all$chr, 
  muts.all$start,
  muts.all$end,
  paste0(muts.all$ref, '/', muts.all$alt),
  '*',
  mid = muts.all$mid3
)
vep.input <- vep.input[which(!is.na(muts.all$mid)),]
write.table(vep.input, file='data/vep.input.txt', sep='\t', quote = FALSE, col.names = F, row.names = F)

# Now we run Ensembl VEP command line version:
# (Note this can be tricky to install...)
# ~/Software/ensembl-vep/vep -i ~/Dropbox/CIS_Immunology/cis_immunology2/data/vep.input.txt -o ~/Dropbox/CIS_Immunology/cis_immunology2/data/vep.output.txt --cache --offline
# 
vep.output <- read.table('data/vep.output.txt', sep='\t', stringsAsFactors = F)
vep.output$impact <- str_extract(vep.output$V14, "IMPACT=[A-Z,]+")
vep.output$impact <- gsub('IMPACT=', '', vep.output$impact)

# Aggregate the impact scores by V1:
vep.output <- aggregate(vep.output[,c('V7', 'impact')], by=list(vep.output$V1), FUN=function(x) {
  paste(sort(unique(x)), collapse=',')
})

# vep.output$patient <- str_extract(vep.output[,1], "PD[0-9]+[a-z]+")
# vep.output$mid <- gsub('PD[0-9]+[a-z]+-', '', vep.output[,1])

muts.all$vep.impact <- vep.output$impact[match(muts.all$mid3, vep.output$Group.1)]
muts.all$vep.consequence <- vep.output$V7[match(muts.all$mid3, vep.output$Group.1)]

# Which are relevant?
# Include amplifications and deletions here
# Exclude gains and LOH _except_ HLA LOH
muts.all$vep.sig <- grepl('HIGH|MODERATE', muts.all$vep.impact) | muts.all$type %in% c('Amplification', 'Deletion', 'CN break', 'Rearrangement', 'LOH', 'Gain') # | muts.all$type == 'LOH' & grepl('^HLA-', muts.all$gene)

muts.all$mid3 <- NULL



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

preinvSum$tmp <- 100*(preinvSum$cd4s + preinvSum$cd8s + preinvSum$foxp3s + preinvSum$hems_IHC)/(preinvSum$tissuearea_grayImage*(0.227*0.227)*2*2)

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
# Imaging data (lymphocytes quantified from H&E)
##############################################################################
# Load image data, annotate pheno
# Note that this overwrites the preinvSum and preinvSum2 variables above
# load('data/image_data_cache.RData')
# # Add patient ID to imaging data:
# preinvSum$patient <- pheno$Patient.Number[match(preinvSum$UniqueID, pheno$HandE.SampleID)]
# # Sometimes we need to match AltSamples:
# sel <- which(is.na(preinvSum$patient))
# alts <- lapply(pheno$HandE.AltSamples, function(x) {unlist(strsplit(x, ", "))})
# pts <- unlist(lapply(sel, function(i){
#   uid <- preinvSum$UniqueID[i]
#   samp <- which(unlist(lapply(alts, function(x) {uid %in% x})))
#   return(pheno$Patient.Number[samp])
# }))
# preinvSum$patient[sel] <- pts

# PerArea calculation (as per Khalid)
# preinvSum$lymphocytes_perArea <- 100*(preinvSum$lymphocytes)/(preinvSum$tissuearea_grayImage*(0.227*0.227)*2*2)
# 
# preinvSum.cis <- preinvSum[which(preinvSum$site == 'CIS'),]
# preinvSum.stroma <- preinvSum[which(preinvSum$site == 'Stroma'),]
# pheno$lymphocytes_per <- preinvSum.cis$lymphocyte_per[match(pheno$HandE.SampleID, preinvSum.cis$UniqueID)]
# pheno$lymphocytes_perArea <- preinvSum.cis$lymphocytes_perArea[match(pheno$HandE.SampleID, preinvSum.cis$UniqueID)]
# pheno$stroma_lymphocytes_per <- preinvSum.stroma$lymphocyte_per[match(pheno$HandE.SampleID, preinvSum.stroma$UniqueID)]
# pheno$stroma_lymphocytes_perArea <- preinvSum.stroma$lymphocytes_perArea[match(pheno$HandE.SampleID, preinvSum.stroma$UniqueID)]
# # Also look at lymphocyte gradient tumour - stroma
# pheno$lym_grad <- unlist(lapply(pheno$HandE.SampleID, function(sid) {
#   preinvSum.cis$lymphocyte_per[match(sid, preinvSum.cis$UniqueID)] - preinvSum.stroma$lymphocyte_per[match(sid, preinvSum.cis$UniqueID)]
# }))

x <- read.csv("~/Dropbox/HandE_for_Khalid/edited_output/2_sum_patientSite.csv")
x$UniqueID <- str_extract(x$FileName, "^S[0-9]+")

# Percentage calculations
x$total_cells <- x$stromals + x$lymphocytes + x$tumors + x$othercells
x$lymphocyte_per <- x$lymphocytes / x$total_cells
x$stromal_per <- x$stromals / x$total_cells
x$tumor_per <- x$tumors / x$total_cells
x$othercells_per <- x$othercells / x$total_cells

preinvSum <- x

# Copy data to pheno
preinvSum.cis <- preinvSum[which(preinvSum$site == 'CIS'),]
preinvSum.stroma <- preinvSum[which(preinvSum$site == 'Stroma'),]
pheno$lymphocytes_per <- preinvSum.cis$lymphocyte_per[match(pheno$SampleID, preinvSum.cis$UniqueID)]
pheno$lymphocytes_perArea <- preinvSum.cis$lymphocyte_perArea[match(pheno$SampleID, preinvSum.cis$UniqueID)]
pheno$stroma_lymphocytes_per <- preinvSum.stroma$lymphocyte_per[match(pheno$SampleID, preinvSum.stroma$UniqueID)]
pheno$stroma_lymphocytes_perArea <- preinvSum.stroma$lymphocyte_perArea[match(pheno$SampleID, preinvSum.stroma$UniqueID)]
# Also look at lymphocyte gradient tumour - stroma
pheno$lym_grad <- unlist(lapply(pheno$SampleID, function(sid) {
  preinvSum.cis$lymphocyte_per[match(sid, preinvSum.cis$UniqueID)] - preinvSum.stroma$lymphocyte_per[match(sid, preinvSum.cis$UniqueID)]
}))

pheno$HandE <- pheno$SampleID %in% preinvSum$UniqueID


###########

# Minor data fixes
gm.tcga.pheno$PatientBarcode <- substr(gm.tcga.pheno$submitter_id2, 1, 12)

# Add follow up data for regressive cases
pheno.surv <- read.xls("data/pheno.with.survival.xls")
pheno$last.seen <- as.character(pheno.surv$LastContact[match(pheno$SampleID, pheno.surv$SampleID)])
pheno$cancer.date <- as.character(pheno.surv$CancerDate[match(pheno$SampleID, pheno.surv$SampleID)])
pheno$exclude.reg <- as.character(pheno.surv$ExcludeAsRegressive[match(pheno$SampleID, pheno.surv$SampleID)]) == 'TRUE'
pheno$time.from <- unlist(lapply(as.character(pheno$Biopsy.Date), function(x) {
  as.character(as.Date(x, tryFormats = c("%d/%m/%Y", "%Y-%m-%d")))
}))
pheno$fu.time <- difftime(as.Date(pheno$last.seen, format = "%Y-%m-%d"), pheno$time.from)
pheno$fu.time <- as.numeric(pheno$fu.time) / (24*365.25) # Convert to years


# Include a pheno data frame with no duplicated patients (e.g. for HLA analysis).
# Prioritise multi-omic samples.
pheno.nodups <- pheno[order(pheno$Methylation, pheno$Whole.Genome.Sequencing, pheno$Stroma.GXN, pheno$Gene.expression, pheno$HandE, pheno$Nanostring, pheno$IHC_K, decreasing = T),]
pheno.nodups <- pheno.nodups[which(!duplicated(pheno.nodups$Patient.Number)),]

####################################################################################
# Save output
####################################################################################
save(
  pheno, pheno.nodups,
  gdata, gpheno, gdata.d, gpheno.d, gdata.v, gpheno.v,
  gdata.pair.t, gdata.pair.s, gpheno.pair,
  mdata, mpheno,
  mdata.genes, mdata.genes.nosnp, mdata.nosnp, mvps.immune, mvps.hla,
  gdata.zs, gdata.zs.v, mdata.zs, mdata.genes,
  ihc,
  muts.all,
  neoantigen.data,
  gm.tcga.gdata, gm.tcga.mdata, gm.tcga.mdata.genes, gm.tcga.pheno,
  tcga.cn.genes, tcga.cn.table, tcga.cn.pheno,
  gdata.danaher.v, gdata.danaher.t, gdata.danaher.s,
  methCS, methCS.tcga,
  pheno.mihc, mihc.cis.data, mihc.stroma.data, mihc.cis.norm, mihc.cis.norm2, mihc.stroma.norm, mihc.stroma.norm2,
  file="data/cached.RData"
)
