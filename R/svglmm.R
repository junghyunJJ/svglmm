# 1) 'spatial_model' is equivalent to 'spatialDE' and models only the mean effects of the environment.
# 2) 'II model' allows environment-specific spatial variance and noise but assumes that th ese values are constant across every environment.
# 3) 'III allows the spatial and noise levels to freely vary between environments.

# scale_cov <- function(K) {
#   w <- mean(diag(K))
#   K <- K / w
#   return(K = K, w = w)
# }

check_K <- function(K, X0) {
  Xa <- X0 %*% solve(t(X0) %*% X0)
  KXXa <- (K %*% X0) %*% t(Xa)
  K <- K - KXXa - t(KXXa) + X0 %*% (t(Xa) %*% (KXXa))
  w	<- mean(diag(K))
  K <- K / w
 
  keep <- which(lower.tri(K, diag = TRUE), arr.ind = TRUE)
  non_mis <- rep(1, length(keep[, 1]))
  keep <- cbind(keep, non_mis, K[keep])
  keep <- keep[order(keep[, 1], keep[, 2]), ]
  browser()
  #!!!!!!!!!
  return(list(K = keep, w = w))
}

fitlmm <- function(y, X, K, df, stype, etype, rescaling) {

  X <- cbind(1, X)

  cat("[", format(Sys.time()), "]", " - Rescaling 'spatial' kernel\n", sep = "")
  K <- check_K(K, X)
  ws <- mean(diag(K))
  K <- K / ws
  K_list <- list(K)

  # stype
  if (stype == "iid") {
    cat("[", format(Sys.time()), "]", " - Rescaling 'stype = iid' kernel\n", sep = "")
    K2 <- K * tcrossprod(Z) # K * (Z Z^T)
    res_check_K2 <- check_K(K2, X)
    
    K_list <- append(K_list, list(res_check_K2$K))
    ws <- c(ws, res_check_K2$w)
  } else if (stype == "free") {
    cat("[", format(Sys.time()), "]", " - Rescaling 'stype = free' kernel\n", sep = "")
    for (k in seq_len(df)) {
      sel_Z <- Z[, k, drop = FALSE]
      K2 <- K * tcrossprod(sel_Z) # K * (Zk Zk^T)
      res_check_K2 <- check_K(K2, X)
    
      K_list <- append(K_list, list(res_check_K2$K))
      ws <- c(ws, res_check_K2$w)
    }
  }

  # etype
  if (etype == "free") {
    cat("[", format(Sys.time()), "]", " - Rescaling 'etype = free' kernel\n", sep = "")
    for (k in seq_len(df - 1)) {
      K3 <- Z[, k]^2 # I * (Zk Zk^T)
      res_check_K3 <- check_K(K3, X)
    
      K_list <- append(K_list, list(res_check_K3$K))
      ws <- c(ws, res_check_K3$w)
    }
  }
  browser()

  fit <- qgg::greml(y = y, X = X, GRM = K_list)
  vc <- fit$theta
  asd <- fit$asd
  browser()

  # rescaling
  # TO DO: we need to check the rescaling
  ws <- c(ws, 1 - ncol(X) / nrow(X))
  if (rescaling) {
    vc <- vc / ws
    # sig2ses <- sig2ses / w
    asd <- diag(1 / ws) %*% asd %*% diag(1 / ws)
  }

  # # names of output
  # if (model == "s") {
  #   stype <- "hom"
  #   etype <- "hom"
  #   names(vc) <- c("v_s", "v_e")
  # } else if (model == "sc") {
  #   stype <- "iid"
  #   etype <- "hom"
  #   names(vc) <- c("v_s", "v_sc", "v_e")
  # } else if (model == "sc") {
  #   stype <- "iid"
  #   etype <- "iid"
  #   # names(vc) <- c("v_s", "v_sc", "v_e")
  # } else {
  #   warning("Please set 'modle = s' or 'modle = sc'")
  # }
  return(list(vc = vc, asd = asd, ll = fit$ll, beta = fit$b, v_beta = fit$vb, stype = stype, etype = etype))
}

cal_cov <- function(res_fit, df) {
  sig2out <- sig2map(res_fit$vc, res_fit$stype, res_fit$etype, df)
  h2Covmat <- msm::deltamethod(sig2out$form, res_fit$vc, res_fit$asd, ses = FALSE)

  if (length(sig2out$h2) == 1) {
    h2 <- sig2out$h2
  } else {
    h2 <- c(sig2out$h2, tot = sum(sig2out$h2))
  }
  return(list(h2 = h2, h2Covmat = h2Covmat))
}

cal_p <- function(res_h2, model) {
  if (model == "s") {
    p <- GxEMM::Waldtest(res_h2$h2, res_h2$h2Covmat[1, 1])
  } else if (model == "sc") {
    p <- sapply(seq_len(length(res_h2$h2) - 1), function(i) {
      as.numeric(GxEMM::Waldtest(res_h2$h2[i], res_h2$h2Covmat[i, i]))
    })
    p_tot <- as.numeric(GxEMM::Waldtest(res_h2$h2['tot'], sum(diag(res_h2$h2Covmat))))
    p <- c(p, p_tot)
    names(p) <- c('s', 'sc', 'tot')
  } else {
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  }
  names(p) <- paste0('p_', names(p))
  
  return(p)
}


svglmm <- function(y, X, K, stype = c('hom', 'iid', 'free')[1], etype = c('hom', 'free')[1], rescaling = TRUE, resfull = TRUE) {
  y <- as.matrix(scale(y))
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  K <- as.matrix(K)
  df <- ncol(Z)

  stopifnot(all.equal(nrow(y), nrow(K)))
  stopifnot(all.equal(nrow(y), nrow(Z)))

  # 1. fitting
  cat("[", format(Sys.time()), "]", " - Fitting liner mixed model\n", sep = "")
  res_fit <- fitlmm(y, X, K, df, stype, etype, rescaling)
  vc <- res_fit$vc
  asd <- res_fit$asd
  
  # 2. cal cov matrix using Delta method
  cat("[", format(Sys.time()), "]", " - Calulating covariance matrix\n", sep = "")
  res_h2 <- cal_cov(res_fit, df)
  h2 <- res_h2$h2
  
  # 3. cal pvalue using Wald
  cat("[", format(Sys.time()), "]", " - Calulating p-value\n", sep = "")
  p <- cal_p(res_h2, model)
  
  if (resfull) {
    final_res <- list(
      vc = vc, ll = as.numeric(res_fit$ll), df = df,
      h2 = h2, h2Covmat = res_h2$h2Covmat, h2_p = p, 
      beta = as.numeric(res_fit$beta), v_beta = as.matrix(res_fit$v_beta)
    ) 
  } else {
    final_res <- list(
      vc = vc, ll = as.numeric(res_fit$ll), df = df,
      h2 = h2, h2_p = p, 
      beta = as.numeric(res_fit$beta), v_beta = as.matrix(res_fit$v_beta)
    )
  }
  return(final_res)
}
