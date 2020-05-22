source("utility_functions/vep.annotate.R")
muts.for.plotting <- function(genes, pheno, do.vep=T, excluded.muttypes = c('upstream', 'downstream', 'intronic', '3prime_UTR_variant', 'silent')){
  mhc.muts <- muts.all[which(muts.all$gene %in% genes),]
  tmp <- mhc.muts
  tmp$Chromosome <- mhc.muts$chr
  tmp$Start.Position <- mhc.muts$start
  tmp$Reference.Allele <- mhc.muts$ref
  tmp$Variant.Allele <- mhc.muts$alt
  
  # Filter out low impact variants
  # This can take ages with lots of genes so allow the user to skip this step
  if(do.vep) {
    # x <- vep.annotate(tmp[which(tmp$class %in% c("SNV", 'DI', 'R')),])
    # mhc.muts$vep.impact <- NA
    # mhc.muts$vep.impact[match(x$mid, mhc.muts$mid)] <- x$vep.impact
    # 
    # # Ignore those which are modifier/low impact and those without a prediction
    # mhc.muts.sig <- mhc.muts
    # sel.rm <- which(
    #   mhc.muts.sig$vep.impact %in% c("MODIFIER", "LOW", "LOW,MODIFIER") |
    #     (is.na(mhc.muts.sig$vep.impact) & as.character(mhc.muts.sig$type) %in% excluded.muttypes)
    # )
    # if(length(sel.rm) > 0) { mhc.muts.sig <- mhc.muts.sig[-sel.rm,] }
    mhc.muts.sig <- mhc.muts[which(mhc.muts$vep.sig),]
  } else {
    mhc.muts.sig <- mhc.muts
    sel.rm <- which(
      (is.na(mhc.muts.sig$vep.impact) & as.character(mhc.muts.sig$type) %in% c('upstream', 'downstream', 'intronic', '3prime_UTR_variant', 'silent'))
    )
    if(length(sel.rm) > 0) { mhc.muts.sig <- mhc.muts.sig[-sel.rm,] }
  }
  
  # Create a data frame showing the presence of changes
  hla.changes <- foreach(i=1:dim(pheno)[1], .combine=rbind, .export = c('pheno', 'gdata.zs', 'gdata.zs.v', 'mdata.zs')) %do% {
    # print(i)
    df <- data.frame(
      gene <- character(0)
    )
    sample <- pheno$SampleID[i]
    for(gene in genes) {
      df <- rbind(df, data.frame(
        gene = rep(gene, 3),
        mod = c('genomic', 'meth', 'gxn'),
        id = paste(rep(gene, 3), c('genomic', 'meth', 'gxn'), sep="."),
        sample = rep(sample, 3),
        outcome = rep(pheno$Outcome[i]),
        result = -1,
        type = NA
      , stringsAsFactors = F))
      # Check if they have genomic data, set to 0 or 1 if they do
      # Include the mutation type
      if(pheno$Whole.Genome.Sequencing[i]) {
        muts.s <- mhc.muts.sig[which(as.character(mhc.muts.sig$sampleID) == as.character(pheno$SampleID[i]) & as.character(mhc.muts.sig$gene) == gene),]
        if(dim(muts.s)[1] > 0) {
          df$result[which(df$gene == gene & df$mod == 'genomic')] <- 1
          df$type[which(df$gene == gene & df$mod == 'genomic')] <- as.character(muts.s[1,]$type) #paste(unique(muts.s$type), collapse = "|")
        } else {
          df$result[which(df$gene == gene & df$mod == 'genomic')] <- 0
          df$type[which(df$gene == gene & df$mod == 'genomic')] <- 'Normal'
        }
      }
      # Check if they have methylation data, set to 0 or 1 if they do
      if(pheno$Methylation[i] & gene %in% rownames(mdata.zs) & as.character(pheno$SampleID[i]) %in% colnames(mdata.zs)) {
        z <- mdata.zs[gene,pheno$SampleID[i]]
        if(abs(z) >= 2){
          df$result[which(df$gene == gene & df$mod == 'meth')] <- 1
          df$type[which(df$gene == gene & df$mod == 'meth')] <- ifelse(z > 0, 'Hypermethylation', 'Hypomethylation')
        } else {
          df$result[which(df$gene == gene & df$mod == 'meth')] <- 0
          df$type[which(df$gene == gene & df$mod == 'meth')] <- 'Normal'
        }
      }
      # Check if they have expression data, set to 0 or 1 if they do
      if(pheno$Stroma.GXN[i]) {
      # if(pheno$Gene.expression[i]) {
        z.use <- NULL
        # if(gene %in% rownames(gdata.zs)) {
        #   z.use <- gdata.zs
        # } else {
          if(gene %in% rownames(gdata.zs.v) & as.character(pheno$SampleID[i]) %in% colnames(gdata.zs.v)) {
            z.use <- gdata.zs.v
          }
        # }
        if(!is.null(z.use)) {
          z <- z.use[gene,as.character(pheno$SampleID[i])]
          if(abs(z) >= 2){
            df$result[which(df$gene == gene & df$mod == 'gxn')] <- 1
            df$type[which(df$gene == gene & df$mod == 'gxn')] <- ifelse(z > 0, 'Overexpression', 'Underexpression')
          } else {
            df$result[which(df$gene == gene & df$mod == 'gxn')] <- 0
            df$type[which(df$gene == gene & df$mod == 'gxn')] <- 'Normal'
          }
        }
      }
      
    }
    return(df)
  }
  
  # Don't plot Control samples
  hla.changes <- hla.changes[which(hla.changes$outcome %in% c('Progression', 'Regression')),]
  
  hla.changes$o <- unlist(lapply(hla.changes$sample, function(x) {sum(hla.changes$result[which(hla.changes$sample == x)])}))
  hla.changes$n.modalities <- unlist(lapply(hla.changes$sample, function(x) {
    y <- pheno[which(pheno$SampleID == x),]
    return(length(which(c(y$Stroma.GXN, y$Methylation, y$Whole.Genome.Sequencing))))
  }))
  hla.changes$n.changes <- unlist(lapply(hla.changes$sample, function(x) {sum(hla.changes$result[which(hla.changes$sample == x & hla.changes$result == 1)])}))
  
  hla.changes <- hla.changes[which(hla.changes$n.modalities > 0),]
  
  # Capitalise:
  hla.changes$type[which(hla.changes$type == 'start_lost')] <- "Start Lost"
  sel <- which(hla.changes$type %in% c('missense', 'nonsense', 'frameshift'))
  hla.changes$type[sel] <- str_to_title(hla.changes$type[sel])
  
  return(hla.changes)
}
