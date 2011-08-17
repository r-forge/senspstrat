calc.pvalue <- function(x, var) {
  a <- pnorm(x/sqrt(var))
  2 * ifelse(a > 0.5, 1 - a, a)
}
