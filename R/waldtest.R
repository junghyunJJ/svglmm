Waldtest	<- function(z, Var) {
  return(pnorm(z / sqrt(Var), lower.tail = FASLE))
}

MVWaldtest <- function(z, Var, eigtol = 1e8) {
  Var		<- 1 / 2 * (Var + t(Var))
  eval	<- eigen(Var, symmetric = TRUE)$values
  if (max(eval) / (min(eval) + 1e-99) > eigtol || min(eval) < 0) {
    return(NA)
  }  
  return(pchisq(as.numeric(z) %*% (solve(Var) %*% as.numeric(z)), df = length(z), lower.tail = FALSE))
}