# Add stars to a plot at position x, based on p-value p
ann_fun <- function(geom, x, p, m=NA) {
  if(is.na(m)) {
    # Extract the variable on the y-axis and find the max value
    y.lab <- as.character(geom$mapping$y)[2]
    m <- max(geom$data[,y.lab], na.rm = T)
  }
  size <- 7
  if(p > 0.1) {return(geom)}
  if(p < 0.001) { lab <- '***' } else {
    if(p < 0.01) { lab <- '**' } else {
      if(p < 0.05) { lab <- '*' } else {
        lab <- '#'
        size <- 4
      }
    }
  }
  
  return( geom + annotate("text", x=x, y=m, label=lab, size=size, color = '#2648FE') )
}
