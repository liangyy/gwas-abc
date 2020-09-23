calc_odds_ratio = function(n11, n12, n21, n22) {
  tmp = fisher.test(matrix(c(n11, n12, n21, n22), ncol = 2))
  data.frame(or = tmp$estimate, pval = tmp$p.value)
}
inv_norm = function(x, offset = 1, ...) {
  r = rank(x, ...)
  g = r / (length(r) + offset)
  o = qnorm(g)
  return(o)
}
