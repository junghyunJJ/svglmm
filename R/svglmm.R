# 1) 'hom' is equivalent to 'spatialDE' and models only the mean effects of the environment.
# 2) 'iid' allows cancer-specific (i.g., output from METI) spatial variance and noise but assumes that th ese values are constant across every environment.
# 3) 'free' allows the spatial and noise levels to freely vary between cancer.


# check_K <- function(K, X0) {
#   Xa <- X0 %*% solve(t(X0) %*% X0)
#   KXXa <- (K %*% X0) %*% t(Xa)
#   K <- K - KXXa - t(KXXa) + X0 %*% (t(Xa) %*% (KXXa))
#   return(K)
# }

summary_res <- function(vc, asd, df, Z, stype, etype) {
  
  names(vc) <- NULL
  names(vc)[1] <- "sig2s0"
  names(vc)[length(vc)] <- "sig2e"
  
  if (stype == "hom") {
    sig2s <- vc[1]
    sig2e <- vc[2]
  } else if (stype == "iid") {
    names(vc)[2] <- "sig2s"
    sig2s <- vc[1:2]
    sig2e <- vc[3]
  } else if (stype == "free") {
    if (etype == "hom") {
      idx_v <- seq.int(from = 2, length.out = df)
      v <- vc[idx_v]
      names(v) <- paste0("sig2s_", colnames(Z))
      sig2s <- c(vc[1], v)
      sig2e <- vc[length(vc)]

      names(vc)[c(idx_v)] <- paste0("sig2s_", seq_len(length(c(idx_v))))

    } else if (etype == "free") {
      idx_v <- seq.int(from = 2, length.out = df)
      idx_w <- seq.int(from = (max(idx_v) + 1), to = length(vc) - 1)

      v <- vc[idx_v]
      w <- c(vc[idx_w], 0) + vc[length(vc)]

      names(v) <- paste0("sig2s_", colnames(Z))
      names(w) <- paste0("sig2e_", colnames(Z))
      sig2s <- c(vc[1], v)
      sig2e <- w

      names(vc)[c(idx_v, idx_w)] <- paste0("sig2s_", seq_len(length(c(idx_v, idx_w))))
    }
  }
  
  return(list(vc = vc, sig2s = sig2s, sig2e = sig2e))
}
  

fitlmm <- function(y, X, K, Z, n_Z, stype, etype, rescaling, maxit) {

  X <- cbind(1, X)
  # K <- check_K(K, X)
  ws <- mean(diag(K))
  K <- K / ws
  K_list <- list(K)

  # stype
  if (stype == "iid") {
    cat("[", format(Sys.time()), "]", " - Rescaling 'stype = iid' kernel\n", sep = "")
    K2 <- K * tcrossprod(Z) # K * (Z Z^T)
    # K2 <- check_K(K2, X)
    w	<- mean(diag(K2))
    K2 <- K2 / w
    K_list <- append(K_list, list(K2))
    ws <- c(ws, w)
  } else if (stype == "free") {
    cat("[", format(Sys.time()), "]", " - Rescaling 'stype = free' kernel\n", sep = "")
    for (k in seq_len(n_Z)) {
      sel_Z <- Z[, k, drop = FALSE]
      K2 <- K * tcrossprod(sel_Z) # K * (Zk Zk^T)
      # K2 <- check_K(K2, X)
      w	<- mean(diag(K2))
      K2 <- K2 / w
      K_list <- append(K_list, list(K2))
      ws <- c(ws, w)
    }
  }
  
  # etype
  if (etype == "free") {
    cat("[", format(Sys.time()), "]", " - Rescaling 'etype = free' kernel\n", sep = "")
    for (k in seq_len(n_Z - 1)) {
      K3 <- diag(Z[, k]^2) # I * (Zk Zk^T)
      # K3 <- check_K(K3, X)
      w	<- mean(diag(K3))
      K3 <- K3 / w
      K_list <- append(K_list, list(K3))
      ws <- c(ws, w)
    }
  }

  # fit lmm
  fit <- qgg::greml(y = y, X = X, GRM = K_list, maxit = maxit)
  vc <- fit$theta
  asd <- fit$asd
  df <- length(ws)
  
  # rescaling
  # TO DO: we need to check the rescaling
  ws <- c(ws, 1 - ncol(X) / nrow(X))
  if (rescaling) {
    vc <- vc / ws
    asd <- diag(1 / ws) %*% asd %*% diag(1 / ws)
  }
  
  # summary res
  summary_sig2 <- summary_res(vc, asd, n_Z, Z, stype, etype)
  
  return(list(
    sig2s = summary_sig2$sig2s, sig2e = summary_sig2$sig2e, vc = summary_sig2$vc,
    asd = asd, ll = fit$ll, df = df, beta = fit$b, v_beta = fit$vb, 
    stype = stype, etype = etype
  ))
}


cal_cov <- function(res_fit, n_Z, Z, stype, etype) {
  # browser()
  sig2out <- sig2map(res_fit$vc, res_fit$stype, res_fit$etype, n_Z)
  h2Covmat <- msm::deltamethod(sig2out$form, res_fit$vc, res_fit$asd, ses = FALSE)
  h2 <- sig2out$h2

  if (stype == "hom") {
    names(h2) <- "h2"
  } else if (stype == "iid") {
    names(h2) <- c("h2_0", "h2")
  } else if (stype == "free") {
    # if (etype == "hom") {
    # } else if (etype == "free") {
    # }
    names(h2) <- paste0("h2_", colnames(Z))
  }

  if (length(sig2out$h2) != 1) {
    h2 <- c(h2, tot = sum(sig2out$h2))
  }
  return(list(h2 = h2, h2Covmat = h2Covmat))
}


cal_h2_p <- function(res_h2, stype) {
  h2 <- res_h2$h2
  h2Covmat <- res_h2$h2Covmat

  if (stype == "hom") {
    p <- GxEMM::Waldtest(h2, h2Covmat[1, 1])
  } else if (stype == "iid") {
    p <- sapply(seq_len(length(h2) - 1), function(i) {
      as.numeric(GxEMM::Waldtest(h2[i], h2Covmat[i, i]))
    })
    p_tot <- as.numeric(GxEMM::Waldtest(res_h2$h2['tot'], sum(diag(res_h2$h2Covmat))))
    p <- c(p, p_tot)
  } else if (stype == "free") {
    p <- sapply(seq_len(length(h2) - 1), function(i) {
      as.numeric(GxEMM::Waldtest(h2[i], h2Covmat[i, i]))
    })
    p_tot <- as.numeric(GxEMM::Waldtest(h2['tot'], sum(diag(h2Covmat))))
    p <- c(p, p_tot)
  }
  names(p) <- paste0('p_', names(h2))
  return(p)
}

svglmm <- function(y, X, K, Z, stype = c('hom', 'iid', 'free')[1], etype = c('hom', 'free')[1], rescaling = TRUE, resfull = TRUE, maxit = 100) {
  cat("[", format(Sys.time()), "] - START (stype : ", stype, " / etype : ", etype, ")\n", sep = "")
  y <- as.matrix(scale(y))
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  K <- as.matrix(K)
  n_Z <- ncol(Z)

  stopifnot(all.equal(nrow(y), nrow(K)))
  stopifnot(all.equal(nrow(y), nrow(Z)))

  # fitting
  cat("[", format(Sys.time()), "]", " - Fitting liner mixed model\n", sep = "")
  res_fit <- fitlmm(y, X, K, Z, n_Z, stype, etype, rescaling, maxit)
  sig2s <- res_fit$sig2s
  vc <- res_fit$vc
  asd <- res_fit$asd
  df <- res_fit$df

  # cal p value
  cat("[", format(Sys.time()), "]", " - Calulating p-value\n", sep = "")
  p_sig2s <- sapply(seq_len(length(sig2s)), function(ii) {
    GxEMM::Waldtest(sig2s[ii], asd[ii, ii])
  })
  names(p_sig2s) <- paste0("p_", names(p_sig2s))
  
  # # tests for spatial heterogeneity using Free model
  # sel_idx <- seq(from = 2, length.out = length(res_fit$sig2s[2:3]))
  # GxEMM::MVWaldtest(res_fit$sig2s[sel_idx], res_fit$asd[sel_idx, sel_idx]) 

  # summary resutls
  if (resfull) {
    cat("[", format(Sys.time()), "]", " - Calulating h2 matrix\n", sep = "")
    res_h2 <- cal_cov(res_fit, n_Z, Z, stype, etype)
    p_h2 <- cal_h2_p(res_h2, stype)
    
    final_res <- list(
      sig2s = sig2s, p_sig2s = p_sig2s, sig2e = res_fit$sig2e,  
      vc = vc,
      asd = asd, ll = as.numeric(res_fit$ll), df = df,
      h2 = res_h2$h2, h2Covmat = res_h2$h2Covmat, p_h2 = p_h2, 
      beta = as.numeric(res_fit$beta), v_beta = diag(as.matrix(res_fit$v_beta))
    )
  } else {
    final_res <- list(
      sig2s = sig2s, p_sig2s = p_sig2s, sig2e = res_fit$sig2e,  
      vc = vc,
      asd = asd, ll = as.numeric(res_fit$ll), df = df
    )
  }
  cat("[", format(Sys.time()), "]", " - END\n", sep = "")
  return(final_res)
}