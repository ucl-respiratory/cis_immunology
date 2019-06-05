# Make a simple mutation table from a list of genes
simpleMutTable <- function(genes, p) {
  muts <- muts.all[which(muts.all$gene %in% genes & muts.all$sampleID %in% p$SampleID),]
  # Filter out irrelevant ones
  sel.rm <- which(muts$type %in% c('intronic', 'silent', 'upstream', 'downstream'))
  if(length(sel.rm) > 0) {
    muts <- muts[-sel.rm,]
  }
  
  df <- matrix(NA, ncol=length(p$SampleID), nrow=length(genes))
  for(i in 1:dim(df)[1]) {
    for(j in 1:dim(df)[2]) {
      sel <- which(muts$gene == genes[i] & muts$sampleID == p$SampleID[j])
      if(length(sel) == 0){ next }
      m <- muts[sel,]
      df[i,j] <- paste(m$type, collapse="|")
    }
  }
  df <- data.frame(df)
  rownames(df) <- genes
  colnames(df) <- p$SampleID
  return(df)
}
