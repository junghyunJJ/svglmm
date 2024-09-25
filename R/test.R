LRtest <- function(ll1, ll0, df) {
  pchisq(2 * (ll1 - ll0), df = df, lower.tail = FALSE)
}

Waldtest	<- function(z, var) {
  return(pnorm(z / sqrt(var), lower.tail = FASLE))
}

MVWaldtest <- function(z, var, eigtol = 1e8) {
  var		<- 1 / 2 * (var + t(var))
  eval	<- eigen(var, symmetric = TRUE)$values
  if (max(eval) / (min(eval) + 1e-99) > eigtol || min(eval) < 0) {
    return(NA)
  }  
  return(pchisq(as.numeric(z) %*% (solve(var) %*% as.numeric(z)), df = length(z), lower.tail = FALSE))
}