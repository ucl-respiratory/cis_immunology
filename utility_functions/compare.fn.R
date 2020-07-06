# Library to compare two means using a mixed effects model.
# Outputs a data frame in the format expected by stat_pvalue_manual
# Limited - will only allow one outcome variable

library(lme4)
library(lmerTest)
library(tibble)
compare.fn <- function(formula, data, comparison=list('Prog.', 'Reg.'), y.position=NULL) {
  lmm <- lmer(formula, data = data, REML = FALSE)
  a <- anova(lmm)
  # Return results in a format expected by stat_pvalue_manual:
  t <- table(lmm@frame[,1], lmm@frame[,2])
  if(is.null(y.position)) {
    y.position = max(lmm@frame[,1]) + (max(lmm@frame[,1]) - min(lmm@frame[,1]))/10
  }
  df <- tibble(
    .y. = colnames(lmm@frame)[1],
    group1 = comparison[[1]],
    group2 = comparison[[2]],
    # p=signif(a$`Pr(>Chisq)`, 2),
    p=signif(a$`Pr(>F)`, 2),
    y.position = y.position
  )
  return(df)
}
