# Library to compare two means using a mixed effects model.
# Outputs a data frame in the format expected by stat_pvalue_manual
# Limited - will only allow one outcome variable

library(lme4)
library(lmerTest)
# library(car)
library(tibble)
compare.fn <- function(formula, data, comparison=list('Prog.', 'Reg.'), y.position=NULL) {
  lmm <- lmer(formula, data = data, REML = FALSE)
  # a <- car::Anova(lmm)
  a <- anova(lmm)
  # Return results in a format expected by stat_pvalue_manual:
  t <- table(lmm@frame[,1], lmm@frame[,2])
  if(is.null(y.position)) {
    y.position = max(lmm@frame[,1]) + (max(lmm@frame[,1]) - min(lmm@frame[,1]))/10
  }
  df <- tibble(
    # group1=list(rownames(t)[which(t[,1] == 1)]),
    # group2=list(rownames(t)[which(t[,2] == 1)]),
    .y. = colnames(lmm@frame)[1],
    group1 = comparison[[1]],
    group2 = comparison[[2]],
    # p=signif(a$`Pr(>Chisq)`, 2),
    p=signif(a$`Pr(>F)`, 2),
    y.position = y.position
  )
  return(df)
}
# library(lme4)
# lmm <- lmer(til.score ~ outcome + (1 | patient), data = plotdata, REML = FALSE)
# library(car)
# a <- Anova(lmm)
# a$`Pr(>Chisq)`

# aov.ponly <- function(formula, data) {
#   a <- aov(formula, data)
#   return(summary(a)[[1]][["Pr(>F)"]])[1]
# }