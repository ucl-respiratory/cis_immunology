plotMethylForGene <- function(gene, do.plot=T) {
  mvps <- probe.features[which(probe.features$gene == gene),]
  for(i in 1:dim(mpheno)[1]) {
    df <- data.frame(
      MAPINFO = mvps$MAPINFO,
      feature = mvps$feature,
      sample = mpheno$sampleID[i],
      Sample_Group = mpheno$Sample_Group[i],
      val = as.numeric(mdata[rownames(mvps),i]),
      label = NA
    )
    if(i == 1) {
      df$label <- df$feature
      df$label[which(duplicated(df$label))] <- NA
      df2 <- df
    } else {
      df2 <- rbind(df2, df)
    }
  }

  p <- ggplot(df2, aes(x=MAPINFO, y=val, label=label, color=Sample_Group))
  
  p <- p + geom_point(cex=0.2) +
    geom_smooth() +
    ylim(-0.2,1) +
    xlab("Position") +
    ylab("Beta value") +
    ggtitle(gene)
  for(feature in as.character(unique(df$feature))) {
    
    p <- p + annotate(
      "segment", 
      x = min(df$MAPINFO[which(df$feature == feature)]), xend = max(df$MAPINFO[which(df$feature == feature)]), 
      y = -0.1, yend = -0.1,
      colour = "#cccccc", size = 1.5
    )
    p <- p + annotate(
      'text', 
      x = mean(df$MAPINFO[which(df$feature == feature)]),
      y = -0.15,
      label = c(feature),
      cex = 3
    )
  }

  if(do.plot) {
    print(p)
  }
  
  return(p)
}