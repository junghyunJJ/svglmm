bw_scott <- function(dd) {
  return(stats::bw.SJ(dd, method = "dpi"))
}

bw_silverman <- function(dd) {
  return(stats::bw.nrd0(dd))
}

bw_mlcv <- function(dd) {
  # https://cran.r-project.org/web/packages/kedd/vignettes/kedd.pdf
  return(kedd::h.mlcv(dd, kernel = "gaussian")$h)
}

bandwidth_select <- function(expr, method = method, nthread = nthread) {

  res_bw <- safe_mclapply(seq_len(nrow(expr)), function(i) {
    if (method == "Scott") {
      bw_scott(expr[i, ])
    } else if (method == "Silverman") {
      bw_silverman(expr[i, ])
    } else if (method == "MLCV") {
      bw_mlcv(expr[i, ])
    } else {
      warning("Please select 'Scott', 'Siverman', 'MLCV'")
    }
  }, mc.cores = nthread, stop.on.error = FALSE)

  # return(unlist(res_bw))
  return(median(na.omit(unlist(res_bw))))
}

kernel_Gaussian <- function(coord, bandwidth) {
  pairwise_sq_dists <- as.matrix(dist(coord)^2)
  K <- exp(-1 * pairwise_sq_dists / bandwidth)
  return(K)
}

cal_spatial_kernel <- function(exp, coord, bandwidthtype = "MLCV", nthread = 1) {

  # Please check the deim for exp and coord
  stopifnot(ncol(exp) == nrow(coord))
  
  # 1. standadization
  expr <- t(scale(t(exp)))

  # 2. cal bandwidth
  bandwidth <- bandwidth_select(expr, method = bandwidthtype, nthread = nthread)
  cat(paste0("The bandwidth is: ", round(bandwidth, 3), " (", bandwidthtype, ")\n"))

  # 3. Calculate the kernel matrix using the bandwidth
  coord_normalized <- scale(coord)
  kernelmat <- kernel_Gaussian(coord = coord_normalized, bandwidth = bandwidth)

  return(kernelmat)
}
