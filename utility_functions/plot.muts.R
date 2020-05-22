# Function to plot deleterious mutations and copy number alterations in a list of genes
#
# hla.genes <- unique(c(somatic.genes.lusc, gene.lists$hla.assoc))
# hla.genes[which(hla.genes == 'HSPA')] <- 'HSPBP1'
# hla.genes[which(hla.genes == 'CNX')] <- 'CANX'
#
# hla.genes <- unique(c(wellenstein.genes, thorsson.drivers, thorsson.pos, thorsson.neg))

plot.muts <- function(
  hla.genes,
  gene.groups = NA,
  reorder = T
) {
  
  if(is.na(gene.groups)) {
    # Group genes by function
    # Use oncogene/TSG definitions from COSMIC
    x <- read.csv('data/Cancer_Gene_Census_Sanger_Jun12_2018.csv', stringsAsFactors = F)
    gene.groups <- list(
      'Onc' = x$Gene.Symbol[grep('oncogene', x$Role.in.Cancer)],
      'TSG' = x$Gene.Symbol[grep('TSG', x$Role.in.Cancer)],
      'AP' = gene.lists$hla.assoc, 
      'TD' = gene.lists$thorsson.drivers,
      'LF+' = gene.lists$thorsson.pos, 
      'LF-' = gene.lists$thorsson.neg,
      'AM' = gene.lists$wellenstein.genes
    )
  }
  gene.groups <- gene.groups[which(!sapply(gene.groups, is.null))]
  
  # Check if groups are unique:
  groups.uniq <- !any(as.numeric(table(unlist(gene.groups))) > 1)
  
  # Plots for heatmaps:
  hmcol.high = '#4575B4'
  hmcol.mid = '#FEFEBF'
  hmcol.low = '#D73027'
  
  # We don't care about those with no changes at all 
  # (safety check for bad gene names)
  hla.genes <- hla.genes[which(hla.genes %in% muts.all$gene)]
  
  # Find muts in these genes only and check VEP status
  mhc.muts <- muts.all[which(muts.all$gene %in% hla.genes),]
  #mhc.muts <- muts.all[which(muts.all$gene %in% tgfb.genes),]
  
  # p <- pheno[which((pheno$Whole.Genome.Sequencing | pheno$Methylation | pheno$Stroma.GXN) & pheno$Outcome != 'Control'),]
  p <- pheno[which((pheno$Whole.Genome.Sequencing) & pheno$Outcome != 'Control'),]
  p <- p[order(p$Whole.Genome.Sequencing, p$Methylation, p$Stroma.GXN),]
  # hla.changes <- muts.for.plotting(hla.genes, pheno = p, excluded.muttypes = excluded.muttypes)
  
  # Create a matrix of changes genes x samples
  hla.changes <- data.frame(
    gene = rep(hla.genes, dim(p)[1]),
    sampleID = unlist(lapply(p$SampleID, function(x) {rep(x, length(hla.genes))})),
    stringsAsFactors = F
  )
  
  mut.types <- c('SNV', 'D', 'I', 'DI', 'R', 'CN break', 'LOH')
  cna.types <- c('Deletion', 'Amplification', 'Loss', 'Gain')
  
  hla.changes$outcome <- p$Outcome[match(hla.changes$sampleID, p$SampleID)]
  hla.changes$n_muts <- sapply(1:dim(hla.changes)[1], function(i) {
    length(which(
      mhc.muts$sampleID == hla.changes$sampleID[i] & 
      mhc.muts$gene == hla.changes$gene[i] & 
      mhc.muts$vep.sig &
      mhc.muts$class %in% mut.types
    ))
  })
  hla.changes$n_cnas <- sapply(1:dim(hla.changes)[1], function(i) {
    length(which(
      mhc.muts$sampleID == hla.changes$sampleID[i] & 
        mhc.muts$gene == hla.changes$gene[i] & 
        mhc.muts$class %in% cna.types
    ))
  }) 
  hla.changes$mutation = 'Normal'
  hla.changes$cn_status = 'Normal'
  for(i in 1:dim(hla.changes)[1]) {
    if(hla.changes$n_muts[i] > 0) {
      mtypes <- mhc.muts$type[which(mhc.muts$sampleID == hla.changes$sampleID[i] & 
                                      mhc.muts$gene == hla.changes$gene[i] & 
                                      mhc.muts$vep.sig &
                                      mhc.muts$class %in% c('SNV', 'D', 'I', 'DI', 'R', 'CN break', 'LOH'))]
      hla.changes$mutation[i] <- paste(unique(mtypes), collapse = '/')
    }
    if(hla.changes$n_cnas[i] > 0) {
      mtypes <- mhc.muts$type[which(mhc.muts$sampleID == hla.changes$sampleID[i] & 
                                      mhc.muts$gene == hla.changes$gene[i] & 
                                      mhc.muts$class %in% c('Deletion', 'Amplification', 'Loss', 'Gain'))]
      hla.changes$cn_status[i] <- paste(unique(mtypes), collapse = '/')
    }
  }
  
  
  # For updated version, use only genomic changes
  # hla.changes <- hla.changes[which(hla.changes$mod == 'genomic'),]
  
  # Tidy up names
  hla.changes$mutation <- gsub('inframe|frameshift', 'Frameshift', hla.changes$mutation)
  hla.changes$mutation <- gsub('missense', 'Missense', hla.changes$mutation)
  hla.changes$mutation <- gsub('nonsense', 'Nonsense', hla.changes$mutation)
  hla.changes$mutation <- gsub('splice_region|ess_splice', 'Splice variant', hla.changes$mutation)
  hla.changes$mutation <- gsub('intronic', 'Intronic', hla.changes$mutation)
  hla.changes$mutation <- gsub('start_lost', 'Start Lost', hla.changes$mutation)
  # hla.changes$mutation[grep('/', hla.changes$mutation)] <- 'Multiple'
  
  library(data.table)
  # dt <- data.table(melt(hla.changes))
  # dt <- dt[, strsplit(as.character(mutation), "/"), by=list(gene, sampleID)]  # this expands "X/Y/Z" into three rows
  dt <- hla.changes
  # sel <- grep('/', dt$mutation)
  # if(length(sel) > 0) {
  #   for(i in sel) {
  #     dt <- rbind(dt, data.table(
  #       gene = dt$gene[i],
  #       sampleID = dt$sampleID[i],
  #       outcome = dt$outcome[i],
  #       n_muts = dt$n_muts[i],
  #       n_cnas = dt$n_cnas[i],
  #       mutation = unlist(strsplit(unlist(dt$mutation[i]), '/')),
  #       cn_status = dt$cn_status[i]
  #     ))
  #   }
  #   
  #   dt <- dt[-sel,]
  # }
  for(i in 1:dim(dt)[1]) {
    if(grepl('/', dt$mutation[i])) {
      mtypes <- unlist(strsplit(unlist(dt$mutation[i]), '/'))
    } else {
      if(dt$mutation[i] == 'Normal') {
        mtypes <- 'Normal'
      } else {
        mtypes <- c(dt$mutation[i], 'Normal')
      }
      
    }
    x <- data.table(
      gene = dt$gene[i],
      sampleID = dt$sampleID[i],
      outcome = dt$outcome[i],
      n_muts = dt$n_muts[i],
      n_cnas = dt$n_cnas[i],
      mutation = mtypes,
      cn_status = dt$cn_status[i]
    )
    if(i == 1) {
      y = x
    } else {
      y = rbind(y, x)
    }
  }
  dt <- y
  
  dt <- data.table(dt)
  
  dt[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(gene, sampleID)]
  dt[, height:=1/.N, by=list(gene, sampleID)]
  # ggplot(dt, aes(sampleID,y=as.numeric(factor(gene)) + shift, fill=mutation, height=height)) + 
  #   geom_tile(color="yellow", size=1) + xlab('Patient') + ylab('Gene')
  hla.changes <- data.frame(dt)
  
  # Define a colour scheme (uses help from RColorBrewer):
  cols <- c(
    # 'Overexpression'='#fc8d59',
    # 'Underexpression'='#91bfdb',
    # 
    # 'Hypermethylation'='#af8dc3',
    # 'Hypomethylation'='#ffffbf',
    
    'Amplification'='#f0027f',
    'Gain'='#fe9acf',
    'Loss'='#386cb0',
    'LOH'='#80b1d3',
    'CN deletion'='#386cb0',
    
    'Missense'='#8dd3c7',
    'Frameshift'='#ffffb3',
    'Nonsense'='#bebada',
    'Start Lost'='#fb8072',
    'Splice variant'='#b15928',
    'Rearrangement'='#fdb462',
    
    'Intronic' = '#b3de69',
    'CN break'='#fccde5',
    
    
    'Multiple' = '#bc80bd',
    
    'Normal'='#eeeeee'
  )
  
  plotdata <- hla.changes
  # plotdata$gene <- factor(plotdata$gene, levels = hla.genes)
  # How many mutations in this entire sample?
  p$n_muts <- sapply(p$SampleID, function(x) {
    length(which(mhc.muts$sampleID == x & mhc.muts$vep.sig & mhc.muts$gene %in% hla.genes & mhc.muts$class %in% mut.types))
  })
  p$n_cnas <- sapply(p$SampleID, function(x) {
    length(which(mhc.muts$sampleID == x & mhc.muts$gene %in% hla.genes & mhc.muts$class %in% cna.types))
  })
  samps.ordered <- unique(p$SampleID[rev(order(p$n_muts + p$n_cnas))])
  plotdata$muts.in.sample <- p$n_muts[match(plotdata$sampleID, p$SampleID)]
  # Order mutations by frequency (include CN here)
  df <- data.frame(
    gene = hla.genes,
    n.muts = sapply(hla.genes, function(x) {length(which(mhc.muts$gene == x & mhc.muts$class %in% mut.types & mhc.muts$vep.sig))}),
    n.cnas = sapply(hla.genes, function(x) {length(which(mhc.muts$gene == x & mhc.muts$class %in% cna.types))}),
    stringsAsFactors = F
  )
  if(groups.uniq) {
    df$group <- sapply(df$gene, function(x) {
      names(which(sapply(gene.groups, function(y) {x %in% y})))
    })
    df$group <- factor(df$group, levels = rev(names(gene.groups)))
    df <- df[order(df$group, df$n.muts + df$n.cnas, decreasing = T),]
  } else {
    df <- df[order(df$n.muts + df$n.cnas, decreasing = T),]
  }
  
  if(reorder) {
    genes.ordered <- df$gene
  } else {
    genes.ordered <- hla.genes
  }
  
  # genes.ordered <- names(sort(table(as.character(mhc.muts$gene[which(mhc.muts$vep.sig)])), decreasing = T))
  plotdata$gene <- factor(as.character(plotdata$gene), levels = genes.ordered)
  
  # p$n.muts <- sapply(p$SampleID, function(x) {length(which(plotdata$sample == x & plotdata$type != 'Normal'))})
  # Remove rows with no changes
  # plotdata <- plotdata[-which(plotdata$gene %in% c('CNX', 'HSPA', 'HSPC')),]
  # plotdata$n.muts.per.gene <- sapply(plotdata$gene, function(x) { length(which(plotdata$gene == x & plotdata$type != 'Normal'))})
  # plotdata <- plotdata[which(plotdata$n.muts.per.gene > 0),]
  
  # samps.ordered <- unique(p$SampleID[order((p$Whole.Genome.Sequencing + p$Methylation + p$Stroma.GXN), p$Whole.Genome.Sequencing, p$Methylation, p$Stroma.GXN, decreasing = T)])
  # samps.ordered <- unique(p$SampleID[rev(order(p$n.muts))])
  
  # plotdata$outcome[which(pheno$exclude.reg[match(plotdata$sampleID, pheno$SampleID)])] <- 'LateProg'
  
  
  # Gene-level metrics
  # For each gene, does its expression correlate with Danaher TIL score?
  # What is the prog/reg z-score?
  # What is the methylation z-score?
  tcga.dan <- do.danaher(gm.tcga.gdata)
  gp <- gpheno.pair
  gp$lympocytes_perArea <- pheno$lymphocytes_perArea[match(gp$SampleID, pheno$SampleID)]
  gene.metrics <- data.frame(
    gene = genes.ordered,
    n.muts = sapply(genes.ordered, function(x) {
      sum(plotdata$n_muts[which(plotdata$gene == x)])
    }),
    n.cnas = sapply(genes.ordered, function(x) {
      sum(plotdata$n_cnas[which(plotdata$gene == x)])
    }),
    til.cor = sapply(genes.ordered, function(x) {
      if(!(x %in% rownames(gdata.pair.t))) {return(NA)}
      cor(as.numeric(gdata.pair.t[x,]), gdata.danaher.t$til.score)
    }),
    til.cor.sig = sapply(genes.ordered, function(x) {
      if(!(x %in% rownames(gdata.pair.t))) {return(NA)}
      cor.test(as.numeric(gdata.pair.t[x,]), gdata.danaher.t$til.score)$p.value < 0.05
    }),
    lym.cor = sapply(genes.ordered, function(x) {
      if(!(x %in% rownames(gdata.pair.t))) {return(NA)}
      cor(as.numeric(gdata.pair.t[x,gp$SampleID]), gp$lympocytes_perArea)
    }),
    lym.cor.sig = sapply(genes.ordered, function(x) {
      if(!(x %in% rownames(gdata.pair.t))) {return(NA)}
      cor.test(as.numeric(gdata.pair.t[x,gp$SampleID]), gp$lympocytes_perArea)$p.value < 0.05
    }),
    til.cor.tcga = sapply(genes.ordered, function(x) {
      if(!(x %in% rownames(gm.tcga.gdata))) {return(NA)}
      cor(as.numeric(gm.tcga.gdata[x,]), tcga.dan$til.score)
    }),
    til.cor.sig.tcga = sapply(genes.ordered, function(x) {
      if(!(x %in% rownames(gm.tcga.gdata))) {return(NA)}
      cor.test(as.numeric(gm.tcga.gdata[x,]), tcga.dan$til.score)$p.value < 0.05
    }),
    # Rough guide of net direction of gains/losses
    cn.direction = sapply(genes.ordered, function(x) {
      length(which(plotdata$gene == x & plotdata$cn_status %in% c('Gain', 'Amplification'))) -
        length(which(plotdata$gene == x & plotdata$cn_status %in% c('LOH', 'Deletion')))
    }),
    stringsAsFactors = F
  )
  # Add cytoband
  x <- data.frame(org.Hs.egSYMBOL)
  gene.metrics$entrez <- x$gene_id[match(gene.metrics$gene, x$symbol)]
  x <- data.frame(org.Hs.egMAP)
  gene.metrics$cytoband <- x$cytogenetic_location[match(gene.metrics$entrez, x$gene_id)]
  # For long band names, use the first part
  gene.metrics$cytoband <- str_extract(gene.metrics$cytoband, '[0-9XY]+[p|q][0-9]+([.][0-9]+)?')
  
  # GXN FC/FDR
  uvv <- limmaCompare(gdata.pair.t, gpheno.pair, fdr_limit = 1)
  gene.metrics$gxn.fc <- uvv[gene.metrics$gene,]$fc
  gene.metrics$gxn.fdr <- uvv[gene.metrics$gene,]$fdr
  gene.metrics$gxn.sig <- gene.metrics$gxn.fdr < 0.01
  # Meth - return FDR and deltaBeta from most significant probe for this gene (max |t|)
  # gene.metrics$meth.fdr <- NA
  # gene.metrics$meth.db <- NA
  # for(i in 1:dim(gene.metrics)[1]) {
  #   gene <- gene.metrics$gene[i]
  #   probes <- dmps[which(dmps$gene == gene),]
  #   if(dim(probes)[1] == 0) {next}
  #   probes <- probes[order(abs(probes$t), decreasing = T),]
  #   gene.metrics$meth.fdr[i] <- probes$adj.P.Val[1]
  #   gene.metrics$meth.db[i] <- probes$deltaBeta[1]
  # }
  # gene.metrics$meth.sig <- gene.metrics$meth.fdr < 0.01 & abs(gene.metrics$meth.db) >= 0.3
  
  # For methylation we want to use ChAMP to compare averaged beta values over this gene set
  # We need to name the genes as probes or ChAMP will fail:
  # x <- mdata.genes[hla.genes,]
  # data("probe.features")
  # y <- probe.features[which(!duplicated(probe.features$gene)),]
  # rownames(x) <- rownames(y[match(rownames(x), y$gene),])
  # x <- x[which(!is.na(rownames(x))),]
  # x.dmp <- champ.DMP(x, pheno = mpheno$Sample_Group)
  # x.dmp <- x.dmp$Progressive_to_Regressive
  # gene.metrics$meth.fdr <- x.dmp$adj.P.Val[match(gene.metrics$gene, x.dmp$gene)]
  # gene.metrics$meth.db <- x.dmp$deltaBeta[match(gene.metrics$gene, x.dmp$gene)]
  # gene.metrics$meth.sig <- gene.metrics$meth.fdr < 0.01 # & abs(gene.metrics$meth.db) >= 0.3
  
  
  ###########
  # Create plots
  
  #### Add gene groups
  gene.metrics$genegroup <- NA
  plotdata$genegroup <- NA
  if(!is.na(gene.groups)) {
    for(i in 1:length(gene.groups)) {
      n = names(gene.groups)[i]
      gene.metrics$genegroup[which(gene.metrics$gene %in% gene.groups[[i]])] <- n
      plotdata$genegroup[which(plotdata$gene %in% gene.groups[[i]])] <- n
    }
  }
  
  #### Gene Ordering
  # Order by mutations + CNAs (default)
  plotdata$gene <- factor(as.character(plotdata$gene), levels = genes.ordered)
  # Optional - order by CN state instead of by #mutations
  # genes.ordered <- gene.metrics$gene[order(gene.metrics$cn.direction)]
  # plotdata$gene <- factor(as.character(plotdata$gene), levels = genes.ordered)
  
  # If there are gene groups defined, order by them:
  # if(!all(is.na(plotdata$genegroup))) {
  #   genes.ordered <- gene.metrics$gene[order(gene.metrics$genegroup, gene.metrics$n.muts + gene.metrics$n.cnas, decreasing = T)]
  #   plotdata$gene <- factor(as.character(plotdata$gene), levels = genes.ordered)
  # }

  # Need to fix multiple mutations in a gene...  
  # x <- plotdata[grep('/', plotdata$mutation),]
  # x1 <- x
  # x1$mutation <- sapply(x1$mutation, function(y) {unlist(strsplit(y, '/'))[1]})
  # x2 <- x
  # x2$mutation <- sapply(x2$mutation, function(y) {unlist(strsplit(y, '/'))[2]})
  # plotdata <- rbind(rbind(plotdata, x1), x2)
  # plotdata <- plotdata[-grep('/', plotdata$mutation),]
  
  
  # Matrix of which genes are in which gene group(s)
  groupmat <- lapply(names(gene.groups), function(group) {
    data.frame(
      gene = genes.ordered, 
      group = group, 
      val = sapply(genes.ordered, function(gene) {gene %in% gene.groups[[group]]}),
      stringsAsFactors = F
    )
  })
  groupmat <- do.call('rbind', groupmat)
  groupmat$gene <- factor(groupmat$gene, levels = rev(levels(plotdata$gene)))
  groupmat$group <- factor(groupmat$group, levels = names(gene.groups))
  
  # Alternative if groups are unique:
  groupmat2 <- lapply(1:length(gene.groups), function(i){
    data.frame(
      group = names(gene.groups)[i],
      gene = gene.groups[[i]],
      stringsAsFactors = F
    )
  })
  groupmat2 <- do.call('rbind', groupmat2)
  groupmat2$gene <- factor(groupmat2$gene, levels = rev(genes.ordered))
  groupmat2$group <- factor(groupmat2$group, levels = names(gene.groups))
  
  
  # Ignore those with no changes
  gene.metrics <- gene.metrics[which(gene.metrics$n.muts > 0 | gene.metrics$n.cnas > 0),]
  plotdata <- plotdata[which(plotdata$gene %in% gene.metrics$gene),]
  # plotdata <- plotdata[which(plotdata$n_muts > 0 | plotdata$n_cnas > 0),]
  groupmat <- groupmat[which(groupmat$gene %in% gene.metrics$gene),]
  groupmat2 <- groupmat2[which(groupmat2$gene %in% gene.metrics$gene),]
  
  fig.main <- ggplot(plotdata, aes(
    # x = factor(sample, levels=unique(plotdata$sample[order(plotdata$n.modalities, plotdata$n.changes, decreasing = T)])), 
    x = factor(sampleID, levels=samps.ordered),
    y = 1
  )) +
    geom_tile(
      aes(
        y = shift,
        # x = 1, y=1+shift,
        # Fill based on mutation type
        fill=factor(mutation, levels=names(cols)), 
        # Border based on CN status
        color = factor(cn_status, levels = names(cols)),
        
        height = height
        
      ),
      width = 0.7, #height = 0.7,
      size = 0.65
      # size=0.25
    ) +
    scale_fill_manual(values=cols, na.value='#ffffff') +
    scale_color_manual(values=cols, na.value='#ffffff') +
    # scale_fill_manual(values=c('#DDDDDD','#56B4E9', '#FF0000')) +
    facet_grid(gene ~ factor(outcome), scales = 'free_x', space='free_x', switch = 'y') +
    theme(
      axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(size = 8, angle = 90),
      strip.text.y = element_text(size=8, angle=-180), strip.text.x = element_text(size=8), strip.background = element_blank(),
      legend.title = element_blank(),
      axis.line = element_blank(),
      legend.direction = 'horizontal',
      legend.text = element_text(size = 8)
    )
  # fig.main
  
  # Extract the legend them remove it
  leg <- get_legend(fig.main)
  # as_ggplot(leg)
  fig.main <- fig.main + theme(legend.position = 'none')
  
  # Highlight which group(s) each gene is in
  fig.group <- ggplot(groupmat, aes(
    x = group,
    y = gene
  )) +
    geom_tile(
      aes(
        # Fill based on mutation type
        fill=factor(val, levels=c('TRUE', 'FALSE'))
      ),
      width = 0.7, height = 0.7,
      size = 0.65
      # size=0.25
    ) +
    theme(
      legend.position = 'none',
      axis.title = element_blank(),
      axis.text.x = element_text(size = 8, angle = 90),
      axis.text.y = element_blank(),
      # axis.text.y = element_text(size = 8), # Useful to debug!
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(size = 8),
      plot.margin = margin(t = 17, r = 0, b = 5, l = 0, unit = "pt")
    ) 
  # fig.group <- ggdraw(switch_axis_position(fig.group, axis = 'y'))
  
  # Alternate version - if the groups are unique:
  
  fig.group2 <- ggplot(groupmat2, aes(
    x = 1,
    y = gene
  )) +
    geom_tile(
      aes(
        # Fill based on mutation type
        fill=group
      ),
      width = 0.7, height = 0.7,
      size = 0.65
      # size=0.25
    ) +
    xlab('group') +
    theme(
      # legend.position = 'none',
      legend.text = element_text(size = 8),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 8),
      # axis.text.x = element_text(size = 8, angle = 90),
      axis.text = element_blank(),
      # axis.text.y = element_text(size = 8), # Useful to debug!
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(size = 8),
      plot.margin = margin(t = 17, r = 0, b = 13, l = 0, unit = "pt")
    ) +
    geom_text(aes(label = group), size = 8 / .pt)
  group.leg <- get_legend(fig.group2)
  fig.group2 <- fig.group2 + theme(legend.position = 'none')
  
  # plot_grid(fig.main, fig.group, nrow = 1, rel_widths = c(5,1))
    # scale_fill_manual(values=c('#DDDDDD','#56B4E9', '#FF0000')) +
    # facet_grid(genegroup ~ ., scales = 'free_x', space='free_x', switch = 'y') +
    # theme(
    #   axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
    #   axis.text.x = element_text(size = 8, angle = 90),
    #   strip.text.y = element_text(size=8, angle=-180), strip.text.x = element_text(size=8), strip.background = element_blank(),
    #   legend.title = element_blank(),
    #   axis.line = element_blank(),
    #   # Hide legend (optional):
    #   legend.position = 'none'
    # )
  
  
  
  mytheme <- theme(
    axis.title.y = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.x = element_text(size = 8),
    strip.text.y = element_text(size=8, angle=-180), strip.text.x = element_text(size=8), strip.background = element_blank(),
    legend.title = element_blank(),
    axis.line = element_blank(),
    legend.position = 'none',
    plot.title = element_text(size = 8),
    plot.margin = margin(b = 13, r = 0, t = 17, l = 0, unit = "pt")
  )
  
  # Illustrative CNA gain/loss plot
  m <- max(abs(gene.metrics$cn.direction))
  fig.cnas <- ggplot(gene.metrics, aes(
    x = 1,
    y = factor(gene, levels = rev(genes.ordered))
  )) +
    geom_tile(
      aes(
        fill=cn.direction
      ),
      width = 0.7, height = 0.7,
      size = 0.95
    ) +
    scale_fill_gradient2(low = '#386cb0', high = '#f0027f', mid = '#ffffff', limits = c(-m,m)) +
    # scale_color_manual(values=c('TRUE' = '#0000ff', 'FALSE' = '#ffffff'), na.value='#ffffff') 
    mytheme +
    geom_text(aes(label = cytoband), size = 8 / .pt)+
    xlab('CNAs')
    # ggtitle('CNAs')
  
  fig.genes <- ggplot(gene.metrics, aes(
    x = 1,
    y = factor(gene, levels = rev(genes.ordered))
  )) +
    geom_tile(
      aes(
        # Fill based on mutation type
        fill=til.cor, 
        color = factor(til.cor.sig)
        # Border based on CN status
        # color = factor(cn_status, levels = names(cols))
      ),
      width = 0.7, height = 0.7,
      size = 0.95
      # size=0.25
    ) +
    # scale_fill_continuous(low = '#ff0000', high = '#00ff00', limits = c(-1,1)) +
    scale_fill_gradient2(low = hmcol.low, mid = hmcol.mid, high = hmcol.high, limits = c(-1,1)) +
    scale_color_manual(values=c('TRUE' = '#0000ff', 'FALSE' = '#ffffff'), na.value='#ffffff') +
    # scale_fill_manual(values=c('#DDDDDD','#56B4E9', '#FF0000')) +
    # facet_grid(gene ~ factor(outcome), scales = 'free_x', space='free_x', switch = 'y') +
    mytheme +
    geom_text(aes(label = round(til.cor, 2)), size = 8 / .pt)+
    xlab('TILcor')
  
  fig.lymcor <- ggplot(gene.metrics, aes(
    x = 1,
    y = factor(gene, levels = rev(genes.ordered))
  )) +
    geom_tile(
      aes(
        # Fill based on mutation type
        fill=lym.cor, 
        color = factor(lym.cor.sig)
        # Border based on CN status
        # color = factor(cn_status, levels = names(cols))
      ),
      width = 0.7, height = 0.7,
      size = 0.95
      # size=0.25
    ) +
    # scale_fill_continuous(low = '#ff0000', high = '#00ff00', limits = c(-1,1)) +
    scale_fill_gradient2(low = hmcol.low, mid = hmcol.mid, high = hmcol.high, limits = c(-1,1)) +
    scale_color_manual(values=c('TRUE' = '#0000ff', 'FALSE' = '#ffffff'), na.value='#ffffff') +
    # scale_fill_manual(values=c('#DDDDDD','#56B4E9', '#FF0000')) +
    # facet_grid(gene ~ factor(outcome), scales = 'free_x', space='free_x', switch = 'y') +
    mytheme +
    geom_text(aes(label = round(lym.cor, 2)), size = 8 / .pt)+
    xlab('Lym.cor')
  
  fig.tcga <- ggplot(gene.metrics, aes(
    x = 1,
    y = factor(gene, levels = rev(genes.ordered))
  )) +
    geom_tile(
      aes(
        # Fill based on mutation type
        fill=til.cor.tcga, 
        color = factor(til.cor.sig.tcga)
        # Border based on CN status
        # color = factor(cn_status, levels = names(cols))
      ),
      width = 0.7, height = 0.7,
      size = 0.95
      # size=0.25
    ) +
    # scale_fill_continuous(low = '#ff0000', high = '#00ff00', limits = c(-1,1)) +
    scale_fill_gradient2(low = hmcol.low, mid = hmcol.mid, high = hmcol.high, limits = c(-1,1)) +
    scale_color_manual(values=c('TRUE' = '#0000ff', 'FALSE' = '#ffffff'), na.value='#ffffff') +
    # scale_fill_manual(values=c('#DDDDDD','#56B4E9', '#FF0000')) +
    # facet_grid(gene ~ factor(outcome), scales = 'free_x', space='free_x', switch = 'y') +
    mytheme +
    geom_text(aes(label = round(til.cor.tcga, 2)), size = 8 / .pt)+
    xlab('TCGA.TILcor')
  
  # GXN
  m <- max(abs(gene.metrics$gxn.fc))
  fig.gxn <- ggplot(gene.metrics, aes(
    x = 1,
    y = factor(gene, levels = rev(genes.ordered))
  )) +
    geom_tile(
      aes(
        fill=gxn.fc, 
        color = factor(gxn.sig)
      ),
      width = 0.7, height = 0.7,
      size = 0.95
      # size=0.25
    ) +
    # scale_fill_continuous(low = '#ff0000', high = '#00ff00', limits = c(-m,m)) +
    scale_fill_gradient2(low = hmcol.low, mid = hmcol.mid, high = hmcol.high, limits = c(-m,m)) +
    scale_color_manual(values=c('TRUE' = '#0000ff', 'FALSE' = '#ffffff'), na.value='#ffffff') +
    # scale_fill_manual(values=c('#DDDDDD','#56B4E9', '#FF0000')) +
    # facet_grid(gene ~ factor(outcome), scales = 'free_x', space='free_x', switch = 'y') +
    mytheme + 
    geom_text(aes(label = round(gxn.fc, 2)), size = 8 / .pt)+
    xlab('GXN PvR')
  
  
  # m <- max(abs(gene.metrics$meth.db))
  # fig.meth <- ggplot(gene.metrics, aes(
  #   x = 1,
  #   y = factor(gene, levels = rev(genes.ordered))
  # )) +
  #   geom_tile(
  #     aes(
  #       fill=meth.db, 
  #       color = factor(meth.sig)
  #     ),
  #     width = 0.7, height = 0.7,
  #     size = 0.95
  #     # size=0.25
  #   ) +
  #   # scale_fill_continuous(low = '#ff0000', high = '#00ff00', limits = c(-m,m)) +
  #   scale_fill_gradient2(low = hmcol.low, mid = hmcol.mid, high = hmcol.high, limits = c(-m,m)) +
  #   scale_color_manual(values=c('TRUE' = '#0000ff', 'FALSE' = '#ffffff'), na.value='#ffffff') +
  #   # scale_fill_manual(values=c('#DDDDDD','#56B4E9', '#FF0000')) +
  #   # facet_grid(gene ~ factor(outcome), scales = 'free_x', space='free_x', switch = 'y') +
  #   mytheme +
  #   geom_text(aes(label = round(meth.db, 2)), size = 8 / .pt) +
  #   xlab('Meth PvR')
  
  # fig <- plot_grid(
  #   fig.main, 
  #   fig.cnas, fig.genes, fig.gxn, 
  #   fig.group,
  #   nrow = 1, 
  #   rel_widths = c(14,1,1,1,3)
  # )
  # fig
  # save_plot('results/mutation.plot.demo.pdf', fig, base_width = fig.width*2, base_height = fig.width*1.5)
  
  # Compare progressive and regressive
  # Use a mixed effects model
  # Account for patient, mutational burden and wgii score
  # Scale variables:
  # p.scaled <- p
  # p.scaled$n_muts <- scale(p.scaled$n_muts)
  # p.scaled$n_cnas <- scale(p.scaled$n_cnas)
  # p.scaled$burden <- scale(p.scaled$burden)
  # p.scaled$wgii <- scale(p.scaled$wgii)
  # # maxburden <- max(pheno$burden, na.rm=T)
  # compare.fn(n_muts ~ Outcome + (1 | Patient.Number), data = p.scaled)
  # compare.fn(n_muts/burden ~ Outcome + (1 | Patient.Number), data = p.scaled)
  # compare.fn(n_cnas ~ Outcome + (1 | Patient.Number), data = p.scaled)
  # compare.fn(n_cnas/wgii ~ Outcome + (1 | Patient.Number), data = p.scaled)
  # lmm <- lmer(n_muts ~ Outcome + burden + wgii + (1 | Patient.Number), data = p.scaled, REML = FALSE)
  # a <- car::Anova(lmm)
  # a
  # lmm <- lmer(n_cnas ~ Outcome + burden + wgii + (1 | Patient.Number), data = p.scaled, REML = FALSE)
  # a <- car::Anova(lmm)
  # a
  # lmm <- lmer(n_muts + n_cnas ~ Outcome + burden + wgii + (1 | Patient.Number), data = p.scaled, REML = FALSE)
  # a <- car::Anova(lmm)
  # a
  
  
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
  
  # Plot the proportions with mutations and CNAs, coloured by type
  line.col = '#aaaaaa'

  theme.topbar <- theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = 'none',
    strip.text = element_blank(), strip.background = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 6), 
    plot.margin = margin(l = 37, r = 7),
    axis.line = element_line(colour = line.col),
    axis.ticks.y = element_line(color = line.col)
  )
  # lim <- aggregate(plotdata[,c(4,5)], by=list(plotdata$sampleID), FUN=sum)
  # lim = max(lim[,c(2,3)])
  # Limit = max number of mutations or CNAs per sample
  sel.mut <- which(plotdata$mutation != 'Normal')
  sel.cna <- which(plotdata$cn_status != 'Normal')
  lim = max(
    aggregate(plotdata[sel.mut,4], by=list(plotdata$sampleID[sel.mut]), FUN=sum)$x,
    aggregate(plotdata[sel.cna,5], by=list(plotdata$sampleID[sel.cna]), FUN=sum)$x
  )
  
  # Fix to make sure all samples are displayed:
  # Add NAs to samples not in the list (i.e. those with no mutations)
  x.muts <- plotdata[sel.mut,]
  s <- samps.ordered[which(!(samps.ordered %in% x.muts$sampleID))]
  if(length(s) > 0) {
    newx <- plotdata[as.numeric(sapply(s, function(sid) {which(plotdata$sampleID == sid)[1]})),]
    newx$sampleID <- s
    newx$mutation <- NA
    x.muts <- rbind(x.muts, newx)
  }
  topbar.muts <- ggplot(x.muts, aes(x = factor(sampleID, levels = samps.ordered), y=n_muts, fill=mutation)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values=cols, na.value='#ffffff') +
    facet_grid(~ factor(outcome), scales = 'free_x', space='free_x', drop = F) +
    ylab('mutations') +
    ylim(0,lim) +
    theme.topbar 
  
  # Repeat for CNAs
  x.cna <- plotdata[sel.cna,]
  s <- samps.ordered[which(!(samps.ordered %in% x.cna$sampleID))]
  if(length(s) > 0) {
    newx <- plotdata[as.numeric(sapply(s, function(sid) {which(plotdata$sampleID == sid)[1]})),]
    newx$sampleID <- s
    newx$cn_status <- NA
    x.cna <- rbind(x.cna, newx)
  }
  topbar.cnas <- ggplot(x.cna, aes(x = factor(sampleID, levels = samps.ordered), y=n_cnas, fill=cn_status)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values=cols, na.value='#ffffff') +
    facet_grid(~ factor(outcome), scales = 'free_x', space='free_x') +
    ylab('CNAs') +
    ylim(0,lim) +
    theme.topbar 
  
  
  # Make equivalent sidebars by gene rather than sample
  theme.sidebar <- theme(
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = 'none',
    strip.text = element_blank(), strip.background = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    plot.margin = margin(t=15, b=3),
    axis.line = element_line(colour = line.col),
    axis.ticks.x = element_line(color = line.col)
  )
  # lim <- aggregate(plotdata[,c(4,5)], by=list(plotdata$gene), FUN=sum)
  # lim = max(lim[,c(2,3)])
  lim = max(
    aggregate(plotdata[sel.mut,4], by=list(plotdata$gene[sel.mut]), FUN=sum)$x,
    aggregate(plotdata[sel.cna,5], by=list(plotdata$gene[sel.cna]), FUN=sum)$x
  )
  
  # Here we need to make sure the data frame has an entry for every gene to be plotted:
  x.muts <- plotdata[sel.mut,]
  genes <- as.character(unique(plotdata$gene))
  g <- genes[which(!(genes %in% x.muts$gene))]
  if(length(g) > 0) {
    # newx <- plotdata[as.numeric(sapply(g, function(sid) {which(plotdata$gene == sid)[1]})),]
    newx <- plotdata[1:length(g),]
    newx$gene <- g
    newx$mutation <- NA
    x.muts <- rbind(x.muts, newx)
  }
  # x.muts <- x.muts[which(!is.na(x.muts$gene)),]
  x.muts$gene <- factor(as.character(x.muts$gene), levels = rev(genes.ordered))
  sidebar.muts <- ggplot(x.muts, aes(x = factor(gene, levels = rev(genes.ordered)), y=n_muts, fill=mutation)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values=cols, na.value='#ffffff') +
    ylab('mutations') +
    ylim(0,lim) + 
    coord_flip() +
    theme.sidebar
  
  x.cnas <- plotdata[sel.cna,]
  genes <- as.character(unique(plotdata$gene))
  g <- genes[which(!(genes %in% x.cnas$gene))]
  if(length(g) > 0) {
    # newx <- plotdata[as.numeric(sapply(g, function(sid) {which(plotdata$gene == sid)[1]})),]
    newx <- plotdata[1:length(g),]
    newx$gene <- g
    newx$cn_status <- NA
    x.cnas <- rbind(x.cnas, newx)
  }
  x.cnas <- x.cnas[which(!is.na(x.cnas$gene)),]
  x.cnas$gene <- factor(as.character(x.cnas$gene), levels = rev(genes.ordered))
  sidebar.cnas <- ggplot(x.cnas, aes(x = factor(gene, levels = rev(genes.ordered)), y=n_cnas, fill=cn_status)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values=cols, na.value='#ffffff') +
    ylab('CNAs') +
    ylim(0,lim) +
    coord_flip() +
    theme.sidebar
  

  # plot_grid(
  #   plot_grid(topbar.muts, topbar.cnas, fig.main, ncol = 1, rel_heights = c(1,1,8)),
  #   plot_grid(NULL, plot_grid(
  #     sidebar.muts, sidebar.cnas
  #   ), rel_heights = c(2,8), ncol = 1),
  #   ncol=2,
  #   rel_widths = c(8,2)
  # )
  
  
  
  # Now we can try to paste these all together...
  if(groups.uniq) {
    group.plot <- fig.group2
    gpsize = 0.75
  } else {
    group.plot <- fig.group
    gpsize = 2
  }
  # fig <- 
  fig <- plot_grid(
    plot_grid(
      plot_grid(topbar.muts, topbar.cnas, ncol = 1),
      as_ggplot(leg),
      ncol=2,
      rel_widths = c(12, 4.25 + gpsize)
      # rel_widths = c(12, 6.25)
    ),
    plot_grid(
      fig.main, 
      sidebar.muts, sidebar.cnas,
      fig.cnas, fig.genes, fig.gxn, 
      group.plot,
      nrow = 1, 
      rel_widths = c(12,1,1,0.75,0.75,0.75,gpsize)
    ),
    ncol = 1,
    rel_heights = c(1,8)
  )
  return(fig)
  # fig
  
  
  # leg <- ggplot(plotdata, aes(x=x, y=y, label=label)) +
  #   geom_point(color=c(hotcol, coldcol)) +
  #   xlim(-1,2.5) +
  #   geom_label(nudge_x = 0.33, label.size = 0) +
  #   guides(color=F) +
  #   theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  #   annotate('text', x=-0.25, y=1, label="Cluster:")
  
  
  # # Alternate:
  # plotdata3 <- plotdata[which(plotdata$result >= 0),]
  # top.bar <- ggplot(plotdata3, aes(factor(sample, levels=samps.ordered))) +
  #   geom_bar(aes(fill=factor(result))) + 
  #   facet_grid(~ factor(outcome), scales = 'free_x', space='free_x', switch = 'y') +
  #   guides(fill=F) +
  #   theme(
  #     axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  #     strip.text.y = element_text(size=8), strip.background = element_blank(),
  #     legend.title = element_blank(),
  #     axis.line.x = element_blank(), 
  #     # plot.margin = margin(0,109, 0,13, "pt") # Use if using no axis
  #     plot.margin = margin(0,90, 0,37, "pt")
  #   ) +
  #   scale_fill_manual(values=c("#EEEEEE", "#777777"))
  # # scale_fill_manual(values=c("#EEEEEE", "#ef8a62"))
  # 
  # 
  # # Do some stats:
  # # Are changes in these genes more common in progressive?
  # # How about if we ignore CNAs and correct for mutational burden?
  # x <- hla.changes[which(hla.changes$type != 'Normal'),]
  # p <- pheno[which(pheno$Whole.Genome.Sequencing),]
  # p$n_changes <- table(x$sample)[p$SampleID]
  # p$n_changes[which(is.na(p$n_changes))] <- 0
  # p1 <- compare.fn(n_changes ~ Outcome + (1 | Patient.Number), data = p, comparison = c('Progression', 'Regression'))
  # ggplot(p, aes(x=Outcome, y=n_changes)) + geom_boxplot() + geom_point() + stat_pvalue_manual(p1)
  # p2 <- compare.fn(n_changes/burden ~ Outcome + (1 | Patient.Number), data = p, comparison = c('Progression', 'Regression'))
  # ggplot(p, aes(x=Outcome, y=n_changes/burden)) + geom_boxplot() + geom_point() + stat_pvalue_manual(p2)
  # print(paste("P-value PvR for all changes:", p1$p))
  # print(paste("P-value PvR for all changes correcting for burden:", p2$p))
  # 
  # # Repeat without LOH/Gain
  # # Do include LOH specifically for HLA genes as polymorphism is important here
  # x <- hla.changes[which(!(hla.changes$type %in% c('Normal', 'Gain'))),]
  # sel.rm <- which(x$type == 'LOH' & !grepl('^HLA-', x$gene))
  # if(length(sel.rm) > 0) {
  #   x <- x[-sel.rm,]
  # }
  # p <- pheno[which(pheno$Whole.Genome.Sequencing),]
  # p$n_changes <- table(x$sample)[p$SampleID]
  # p$n_changes[which(is.na(p$n_changes))] <- 0
  # p1 <- compare.fn(n_changes ~ Outcome + (1 | Patient.Number), data = p, comparison = c('Progression', 'Regression'))
  # ggplot(p, aes(x=Outcome, y=n_changes)) + geom_boxplot() + geom_point() + stat_pvalue_manual(p1)
  # p2 <- compare.fn(1e4*n_changes/burden ~ Outcome + (1 | Patient.Number), data = p, comparison = c('Progression', 'Regression'))
  # ggplot(p, aes(x=Outcome, y=1e4*n_changes/burden)) + geom_boxplot() + geom_point() + stat_pvalue_manual(p2)
  # print(paste("P-value PvR without LOH/Gain:", p1$p))
  # print(paste("P-value PvR without LOH/Gain correcting for burden:", p2$p))
  # sort(p$n_changes)
  # t <- table(p$n_changes > 0, p$Outcome)
  # if(all(dim(t) == c(2,2))) {
  #   f <- fisher.test(t)
  #   print(paste("P-value PvR without LOH/Gain at least 1 change:", f$p.value))
  # }
  # 
  # 
  # # Plot together:
  # fig <- plot_grid(top.bar, fig.main, nrow = 2, rel_heights = c(1,10))
  # fig
}

