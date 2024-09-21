rm(list = ls())
library(tidyverse)
library(data.table)
library(jjutil)

# library(GxEMM)

set.seed(1234)
N <- 1e3 # sample size
S <- 1e2 # number SNPs

Z1  <- rbinom(N, 1, .5)
Z   <- cbind(Z1, 1 - Z1) ### two discrete environments

snps  <- scale(matrix(rbinom(N * S, 2, .1), N, S))
K <- snps %*% t(snps) / S
K %>% h

X <- Z[, -1] # fixed effect covariates. Must include Z! Column is dropped here so that cbi
X %>% ll

#genetic variances--assumed heterogeneous in this simulation
sig2hom <- 0
sig2het <- c(.1, .4)

# noise--assumed homogeneous in this simulation
epsilon <- sqrt(1 - sig2hom - sum(colMeans(Z^2) * sig2het)) * rnorm(N)

# heterogeneous SNP effects--details of this expression are not so important
betas <- sapply(sig2het, function(sig) rnorm(S, sd = sqrt(sig / S)))
uhet <- sapply(seq_len(nrow(Z)), function(i) snps[i, ] %*% (betas %*% Z[i, ]))

y  <- as.numeric(Z %*% c(.1, -.5) + uhet + epsilon)

#wget http://dougspeed.com/wp-content/uploads/ldak5.mac_.zip
ldak_loc  <- "gxemm-master/ldak5.linux "

######################################################################
### hom ##############################################################
######################################################################

out_hom <- GxEMM::GxEMM(y, X, K, Z, gtype = 'hom', ldak_loc = ldak_loc)
# $h2
#       hom
# 0.1118206

# $sig2g
# [1] 0.1034026

# $sig2e
# [1] 0.8213166

# $df
# [1] 1

# $sig2s
# [1] 0.1034026 0.8213166

# $ll
# [1] -1363.064

# $h2Covmat
#              [,1]
# [1,] 0.0008360986

# $sig2Var
#               [,1]          [,2]
# [1,]  0.0007149648 -0.0007149107
# [2,] -0.0007149107  0.0007148566

# $betas
# [1]  0.265808 -0.551468

# GxEMM::Waldtest(out_hom$h2, out_hom$h2Covmat[1, 1])
# 5.505362e-05

# 1) 'spatial_model' is equivalent to 'spatialDE' and models only the mean effects of the environment.
# 2) 'II model' allows environment-specific spatial variance and noise but assumes that th ese values are constant across every environment.
# 3) 'III allows the spatial and noise levels to freely vary between environments.

source("R/sig2map.R")
svglmm <- function(y, X, K, model = c('s', 'sc')[1], resfull = FALSE) {
  y <- as.matrix(y)
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  K <- as.matrix(K)
  K0 <- ncol(Z)

  stopifnot(all.equal(nrow(y), nrow(K)))
  stopifnot(all.equal(nrow(y), nrow(Z)))

  # rescaling
  # TO DO: we need to check the rescaling
  w <- mean(diag(K))
  K <- K / w
  ws <- c(w, 1 - ncol(X) / nrow(X))

  if (model == "s") {
    # spatial_model
    fit <- qgg::greml(y = y, X = cbind(1, X), GRM = list(K))
  
  } else if (model == "sc") {
    # spatial celltype model
    K2 <- K * (Z %*% t(Z))
    w <- mean(diag(K2))
    K2 <- K2 / w
    ws <- c(ws, 1 - ncol(X) / nrow(X))
    fit <- qgg::greml(y = y, X = cbind(1, X), GRM = list(K, K2))
  }
  vc <- fit$theta
  asd <- fit$asd 
  
  browser()

  # rescaling
  # TO DO: we need to check the rescaling
  vc <- vc / ws
  # sig2ses <- sig2ses / w
  asd <- diag(1 / ws) %*% asd %*% diag(1 / ws)

  if (model == "s") {
    gtype <- etype <- "hom"
    names(vc) <- c("vs", "ve")
  } else if (model == "sc") {
    gtype <- "iid"
    etype <- "hom"
    names(vc) <- c("vs", "vsc", "ve")
  } else {
    warning("Please set 'modle = s' or 'modle = sc'")
  }

  sig2out <- sig2map(vc, gtype, etype, K0)
  h2Covmat <- msm::deltamethod(sig2out$form, vc, asd, ses = FALSE)

  if (model == "s") {
    p <- GxEMM::Waldtest(sig2out$h2, out_hom$h2Covmat[1, 1])
  } else if (model == "sc") {
    # TODO!!!!
    p <- GxEMM::Waldtest(sig2out$h2, out_hom$h2Covmat[1, 1])
  }
  
  if (resfull == TRUE) {
    final_res <- list(h2 = sig2out$h2, vc = vc, p = p)
  } else {
    final_res <- list(h2 = sig2out$h2, vc = vc, ll = fit$ll, p = p)
  }

  return(final_res)
}

# svglmm(y, X, K, model = "s")
svglmm(y, X, K, model = "sc")



######################################################################
### iid ##############################################################
######################################################################


out_iid	<- GxEMM(y, X, K, Z, gtype = 'iid', ldak_loc = ldak_loc) ### need to add etype='iid' for non-discrete environments
# $h2
#        hom        het 
# 0.02423972 0.20236330 

# $sig2g
# [1] 0.02271815 0.18966065

# $sig2e
# [1] 0.7248497

# $df
# [1] 2

# $sig2s
# [1] 0.02271815 0.18966065 0.72484970

# $ll
# [1] -1347.182

# $h2Covmat
#              [,1]         [,2]
# [1,]  0.001795932 -0.001689620
# [2,] -0.001689620  0.002966507

# $sig2Var
#               [,1]         [,2]          [,3]
# [1,]  1.577542e-03 -0.001484154 -9.338019e-05
# [2,] -1.484154e-03  0.002605806 -1.121566e-03
# [3,] -9.338019e-05 -0.001121566  1.214855e-03

# $betas
# [1]  0.265808 -0.551468
Waldtest(out_iid$h2[1], out_iid$h2Covmat[1, 1])
Waldtest(out_iid$h2[2], out_iid$h2Covmat[2, 2])
Waldtest(out_iid$h2, out_iid$h2Covmat)


iid_grm_1 <- genio::read_grm("gxemm_tmp_iid/K.1")
iid_grm_2 <- genio::read_grm("gxemm_tmp_iid/K.2")

microbenchmark::microbenchmark(
  qgg::greml(y = y, X = cbind(1, X), GRM = list(iid_grm_1$kinship, iid_grm_2$kinship)),
  qgg::greml(y = y, X = cbind(1, X), GRM = list(iid_grm_1$kinship, iid_grm_2$kinship), ncores = 2),
  times = 5
)

fitG <- qgg::greml(y = y, X = cbind(1, X), GRM = list(iid_grm_1$kinship, iid_grm_2$kinship))
fitG$theta
#         G1         G2          E 
# 0.02507015 0.20930025 0.79997022 

K2 <- K * (Z %*% t(Z))

# qgg::greml(y = y, X = cbind(1, X), GRM = list(K, K2))$theta
# #         G1         G2          E 
# # 0.02512229 0.20973556 0.79997022 


K0	<- ncol(Z)
binary <- FALSE
gtype <- 'hom'
etype <- 'hom'


sig2out		<- sig2map(fitG$theta, gtype, etype, K0, binary = binary)
h2Covmat	<- msm::deltamethod(sig2out$form, fitG$theta, out_hom$sig2Var, ses = FALSE)

list(
  h2 = sig2out$h2, sig2g = sig2out$sig2g, sig2e = sig2out$sig2e, df = r,
  sig2s = gout$sig2s, ll = gout$ll, h2Covmat = h2Covmat, sig2Var = gout$sig2Var, betas = gout$betas
)



######################################################################
### free #############################################################
######################################################################

out_free <- GxEMM(y, X, K, Z, gtype = "free", etype = "free", ldak_loc = ldak_loc)
# $h2
#       h2_1       h2_2 
# 0.04465943 0.36073450 

# $sig2g
#   sig2g_hom     sig2g_1     sig2g_2 
# 0.030651619 0.001945418 0.390421777 

# $sig2e
# [1] 0.6973056 0.7461934

# $df
# [1] 4

# $sig2s
# [1]  0.030651619  0.001945418  0.390421777 -0.048887814  0.746193387

# $ll
# [1] -1330.965

# $h2Covmat
#              [,1]         [,2]
# [1,] 1.252950e-03 1.562555e-05
# [2,] 1.562555e-05 3.155530e-03

# $sig2Var
#               [,1]          [,2]          [,3]          [,4]          [,5]
# [1,]  0.0011477595 -1.077593e-03 -0.0009516815  0.0000329504 -1.486050e-04
# [2,] -0.0010775932  1.685847e-03  0.0007636436 -0.0004802894  9.011915e-05
# [3,] -0.0009516815  7.636436e-04  0.0065201949  0.0009781194 -3.124657e-03
# [4,]  0.0000329504 -4.802894e-04  0.0009781194  0.0058625682 -3.300605e-03
# [5,] -0.0001486050  9.011915e-05 -0.0031246571 -0.0033006051  3.335328e-03

# $betas
# [1]  0.265808 -0.551468
MVWaldtest(out_free$sig2s[2:3], out_free$sig2Var[2:3, 2:3]) 


free_grm_1 <- genio::read_grm("gxemm_tmp_free/K.1")
free_grm_2 <- genio::read_grm("gxemm_tmp_free/K.2")
free_grm_3 <- genio::read_grm("gxemm_tmp_free/K.3")
free_grm_4 <- genio::read_grm("gxemm_tmp_free/K.4")

fitG <- qgg::greml(y = y, X = cbind(1, X), GRM = list(free_grm_1$kinship, free_grm_2$kinship, free_grm_3$kinship, free_grm_4$kinship))
fitG$theta
#          G1          G2          G3          G4           E 
# 0.033957296 0.001075617 0.209995439 0.000000001 0.823509985 



K_free <- lapply(seq_len(ncol(Z)), function(k) {
  K * (Z[, k, drop = FALSE] %*% t(Z[, k, drop = FALSE]))
})
K_free %>% str


E_free <- diag(Z[, 1]^2)


fitG <- qgg::greml(y = y, X = X, GRM = list(K, f_free[[1]], f_free[[2]], E_free))
fitG$theta
