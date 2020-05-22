# Function to run ESTIMATE
#
# To install ESTIMATE:
# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)

library(estimate)

do.estimate <- function(data, platform = 'affymetrix') {
  f1 <- "data/estimate.gdata.in.txt"
  f2 <- "data/estimate.gdata.out.gct"
  f3 <- "data/estimate.gdata.res.gct"
  write.table(data, file = f1, sep='\t', quote = F)
  filterCommonGenes(input.f = f1, output.f = f2, id = 'GeneSymbol')
  estimateScore(input.ds = f2, output.ds = f3, platform = platform)
  x <- read.table(f3, sep='\t', skip=2, header = T, row.names = 1)
  x <- data.frame(t(x[,-1]))
  return(x)
}
