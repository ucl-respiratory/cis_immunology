# Here we estimate immune cell populations using the Davoli method
# Reference: https://www.ncbi.nlm.nih.gov/pubmed/28104840 - see table S4
# Load genes used in the Davoli method (direct download of their table S4):
cache.file = 'data/davoli.genes.RData'
if(file.exists(cache.file)) {
  load(cache.file)
} else {
  library(gdata)
  s4 <- read.xls("~/Dropbox/CIS_Immunology/cis_immunology/resources/NIHMS893889-supplement-TableS4.xlsx", sheet = "Table S4d", skip = 2)
  celltypes <- colnames(s4)
  davoli.genes <- lapply(celltypes, function(x) {
    y <- unique(as.character(s4[,x]))
    y <- y[which(y != '')]
    return(y)
  })
  names(davoli.genes) <- celltypes
  save(davoli.genes, file = cache.file)
}

do.davoli <- function(data) {
  df <- data.frame(matrix(nrow = dim(data)[2], ncol=length(davoli.genes)))
  rownames(df) <- colnames(data)
  colnames(df) <- names(davoli.genes)
  
  for(celltype in names(davoli.genes)) {
    genes <- davoli.genes[[celltype]]
    missing <- genes[which(!(genes %in% rownames(data)))]
    if(length(missing) > 0){
      print(paste("Missing genes for", celltype, ":", paste(missing, collapse=", ")))
      genes <- genes[which(genes %in% rownames(data))]
    }
    df[,celltype] <- apply(data[genes,], 2, mean)
  }
  
  return(df)
}
