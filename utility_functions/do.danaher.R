# Method to create a matrix of Danaher TIL expression from input gene expression data
til.genes <- list(
  "B-cells"=c("BLK", "CD19", "FCRL2", "MS4A1", "KIAA0125", "TNFRSF17", "TCL1A", "SPIB", "PNOC"),
  "CD45"=c("PTRPC"),
  "Cytotoxic cells"=c("PRF1", "GZMA", "GZMB", "NKG7", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "GNLY"),
  "DC"=c("CCL13", "CD209", "HSD11B1"),
  "Exhausted CD8"=c("LAG3", "CD244", "EOMES", "PTGER4"),
  "Macrophages"=c("CD68", "CD84", "CD163", "MS4A4A"),
  "Mast cells"=c("TPSB2", "TPSAB1", "CPA3", "MS4A2", "HDC"),
  "Neutrophils"=c("FPR1", "SIGLEC5", "CSF3R", "FCAR", "FCGR3B", "CEACAM3", "S100A12"),
  "NK CD56dim cells"=c("KIR2DL3", "KIR3DL1", "KIR3DL2", "IL21R"),
  "NK cells"=c("XCL1", "XCL2", "NCR1"),
  "T-cells"=c("CD6", "CD3D", "CD3E", "SH2D1A", "TRAT1", "CD3G"),
  "Th1 cells"=c("TBX21"),
  "Treg"=c("FOXP3"),
  "CD8 T cells"=c("CD8A", "CD8B"),
  "CD4 cells"=c()
)
til.score.types <- c("B-cells", "Cytotoxic cells", "Exhausted CD8", "Macrophages", "Neutrophils", "NK CD56dim cells", "NK cells", "T-cells", "Th1 cells", "CD8 T cells")

do.danaher <- function(data){
  
  
  df <- data.frame(matrix(NA, ncol=length(til.genes), nrow=dim(data)[2]))
  rownames(df) <- colnames(data)
  colnames(df) <- names(til.genes)
  library(limma)
  for(i in 1:length(til.genes)){
    genes <- alias2Symbol(til.genes[[i]])
    sel <- which(!(genes %in% rownames(data)))
    if(length(sel) > 0){
      print(paste("Missing genes for",names(til.genes)[i],":", paste(genes[sel], collapse = ",")))
      genes <- genes[-sel]
    }
    if(length(genes) == 1){
      score <- data[genes,]
    }else{
      score <- apply(data[genes,], 2, mean)
    }
    df[,names(til.genes)[i]] <- as.numeric(score)
  }
  
  # TIL score - average of cell scores excluding DC Tregs mast cells
  df$til.score <- apply(df[,til.score.types], 1, mean)
  
  return(df)
}
