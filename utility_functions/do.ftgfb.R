# Define the fibroblast TGFB response signature from https://www.nature.com/articles/nature25501
# This function applies this to a gene expression array
ftgfb.genes <- c("ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CTGF", "CTPS1", "FAM101B", "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A", "TAGLN", "TGFBI", "TNS1", "TPM1")

do.ftgfb <- function(gene.data) {
  genes <- ftgfb.genes[which(ftgfb.genes %in% rownames(gene.data))]
  ftgfb <- apply(gene.data[genes,], 2, mean)
  return(ftgfb)
}
