# Function to compare mutational burden and CNAs between Prog and Reg for a set of genes
# 
compareBurden <- function(genes, do.dnds=T) {
  
  mhc.muts <- muts.all[which(muts.all$gene %in% genes),]
  p <- pheno[which((pheno$Whole.Genome.Sequencing) & pheno$Outcome != 'Control'),]
  
  mut.types <- c('SNV', 'D', 'I', 'DI', 'R', 'CN break')
  cna.types <- c('Deletion', 'Amplification', 'LOH', 'Loss', 'Gain')
  
  p$n_muts <- sapply(p$SampleID, function(x) {
    length(which(mhc.muts$sampleID == x & mhc.muts$vep.sig & mhc.muts$gene %in% genes & mhc.muts$class %in% mut.types))
  })
  p$n_cnas <- sapply(p$SampleID, function(x) {
    length(which(mhc.muts$sampleID == x & mhc.muts$gene %in% genes & mhc.muts$class %in% cna.types))
  })
  p$n_loss <- sapply(p$SampleID, function(x) {
    length(which(mhc.muts$sampleID == x & mhc.muts$gene %in% genes & mhc.muts$class %in% c('Loss', 'Deletion')))
  })
  p$n_loh <- sapply(p$SampleID, function(x) {
    length(which(mhc.muts$sampleID == x & mhc.muts$gene %in% genes & mhc.muts$class %in% c('LOH')))
  })
  p$n_gain <- sapply(p$SampleID, function(x) {
    length(which(mhc.muts$sampleID == x & mhc.muts$gene %in% genes & mhc.muts$class %in% c('Gain', 'Amplification')))
  })
  
  # Compare progressive and regressive
  # Use a mixed effects model
  # Account for patient, mutational burden and wgii score
  # Scale variables:
  p.scaled <- p
  p.scaled$n_muts <- scale(p.scaled$n_muts)
  p.scaled$n_cnas <- scale(p.scaled$n_cnas)
  p.scaled$n_loss <- scale(p.scaled$n_loss)
  p.scaled$n_loh <- scale(p.scaled$n_loh)
  p.scaled$n_gain <- scale(p.scaled$n_gain)
  p.scaled$burden <- scale(p.scaled$burden)
  p.scaled$wgii <- scale(p.scaled$wgii)
  p.scaled$wgii.loss <- scale(p.scaled$wgii.loss)
  p.scaled$wgii.gain <- scale(p.scaled$wgii.gain)
  
  # Uncorrected comparisons (correcting for patient only)
  # compare.fn(n_muts ~ Outcome + (1 | Patient.Number), data = p.scaled)
  # compare.fn(n_muts/burden ~ Outcome + (1 | Patient.Number), data = p.scaled)
  # compare.fn(n_cnas ~ Outcome + (1 | Patient.Number), data = p.scaled)
  # compare.fn(n_cnas/wgii ~ Outcome + (1 | Patient.Number), data = p.scaled)
  
  # Comparisons of mutations, CNAs, and the sum of both ("changes")
  if(sum(p$n_muts) > 0 & sum(p$n_cnas) > 0) {
    lmm <- lmer(n_muts + n_cnas ~ Outcome + burden + wgii + (1 | Patient.Number), data = p.scaled, REML = FALSE)
    a.both <- anova(lmm)
    # a.both <- car::Anova(lmm)
  } else {
    a.both <- NA
  }
  
  # Repeat but uncorrected for burden/wgii
  if(sum(p$n_muts) > 0 & sum(p$n_cnas) > 0) {
    lmm <- lmer(n_muts + n_cnas ~ Outcome + (1 | Patient.Number), data = p.scaled, REML = FALSE)
    a.uncor <- anova(lmm)
    # a.both <- car::Anova(lmm)
  } else {
    a.uncor <- NA
  }
  
  # Mutations corrected for burden
  if(sum(p$n_muts) > 0) {
    lmm <- lmer(n_muts ~ Outcome + burden + (1 | Patient.Number), data = p.scaled, REML = FALSE)
    a.muts <- anova(lmm)
    
    lmm <- lmer(n_muts ~ Outcome + (1 | Patient.Number), data = p.scaled, REML = FALSE)
    a.muts.uncor <- anova(lmm)
  } else {
    a.muts <- NA
    a.muts.uncor <- NA
  }

  # CNAs corrected for WGII
  if(sum(p$n_cnas) > 0) {
    lmm <- lmer(n_cnas ~ Outcome + wgii + (1 | Patient.Number), data = p.scaled, REML = FALSE)
    a.cnas <- anova(lmm)
    
    lmm <- lmer(n_cnas ~ Outcome + (1 | Patient.Number), data = p.scaled, REML = FALSE)
    a.cnas.uncor <- anova(lmm)
  } else {
    a.cnas <- NA
    a.cnas.uncor <- NA
  }
  
  # Losses and gains independently
  if(sum(p$n_loh) > 0) {
    lmm <- lmer(n_loh ~ Outcome + wgii.loh + (1 | Patient.Number), data = p.scaled, REML = FALSE)
    a.loh <- anova(lmm)
  } else {
    a.loh <- NA
  }
  if(sum(p$n_loss) > 0) {
    lmm <- lmer(n_loss ~ Outcome + wgii.loss + (1 | Patient.Number), data = p.scaled, REML = FALSE)
    a.loss <- anova(lmm)
  } else {
    a.loss <- NA
  }
  if(sum(p$n_gain) > 0) {
    lmm <- lmer(n_gain ~ Outcome + wgii.gain + (1 | Patient.Number), data = p.scaled, REML = FALSE)
    a.gain <- anova(lmm)
  } else {
    a.gain <- NA
  }
  
  
  # Run dndscv on the mutations (for these genes only)
  # Limit to those in RefCDS
  # Use a try-catch block
  dnds <- NA
  dnds.global <- NA
  if(do.dnds) {
    tryCatch({
      data('refcds_hg19')
      valid.genes <- sapply(RefCDS, function(x) {
        x[[1]]
      })
      dndscv.input <- mhc.muts[which(mhc.muts$gene %in% valid.genes & !(mhc.muts$class %in% c('Deletion', 'LOH', 'Gain', 'Amplification', 'R', 'CN break'))), c("patient", "chr", "start", "ref", "alt")]
      # Remove duplicates from samples from the same patient
      colnames(dndscv.input) <- c("sampleID", "chr", "pos", "ref", "mut")
      dndscv.input$uuid <- paste(substr(dndscv.input$sampleID, 1, 7), dndscv.input$chr, dndscv.input$pos, dndscv.input$ref, dndscv.input$mut, sep="-")
      dndscv.input <- dndscv.input[-which(duplicated(dndscv.input$uuid)),]
      dndscv.input$uuid <- NULL
      
      dndsout = dndscv(dndscv.input, gene_list = intersect(genes, valid.genes), outmats = T)
      ci <- geneci(dndsout = dndsout, gene_list = intersect(genes, valid.genes))
      
      # See if any appear significant and return those from sel_cv
      sel_cv <- dndsout$sel_cv
      cols <- colnames(sel_cv)[grep('^q', colnames(sel_cv))]
      sel <- which(apply(sel_cv[,cols], 1, function(x) {any(x < 0.1)}))
      if(length(sel) > 0) {
        dnds <- sel_cv[sel,]
      } else {
        dnds <- NA
      }
      dnds.global <- dndsout$globaldnds
      # return(dnds)
    }, error = function(e) {
      dnds <- NA
      dnds.global <- NA
    })
  }
  
  
  # Check for a significant subset in expression data
  # Look for genes with significant gene expression differences, correcting only for this set
  uvv <- limmaCompare(gdata.pair.t[intersect(rownames(gdata.pair.t), genes),], gpheno.pair, fdr_limit = 0.05)
  
  
  return(list(
    'both' = a.both,
    'uncorrected' = a.uncor,
    'muts.uncorrected' = a.muts.uncor,
    'cnas.uncorrected' = a.cnas.uncor,
    'muts' = a.muts,
    'cnas' = a.cnas,
    'loh' = a.loh,
    'loss' = a.loss,
    'gain' = a.gain,
    'dnds' = dnds,
    'dnds.global' = dnds.global,
    'sig.expression' = uvv
  ))
}
