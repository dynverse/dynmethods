cell <- colnames(scent_data$expMC)[[1]]


scent_data$expMC[, cell]







exp.v = scent_data$expMC[, cell]
adj.m = as.matrix(scent_data$adjMC)
maxSR = max_entropy


scent::CompSRana(exp.v, adj.m, maxSR=maxSR)

test <- function(exp.v, adj.m, maxSR) {
  start.time <- Sys.time()

  sumexp.v <- as.vector(adj.m %*% matrix(exp.v, ncol = 1)) # (Ax)_i
  p.m <- t(t(adj.m) * exp.v)/sumexp.v # p_ij

  S.v <- apply(p.m, 1, CompS) # Sum of (p_ij * log(p_ij))

  invP.v <- exp.v * sumexp.v
  nf <- sum(invP.v)
  invP.v <- invP.v/nf # pi_i

  SR <- sum(invP.v * S.v) # Sr(x)

  if (is.null(maxSR) == FALSE) {
    SR <- SR/maxSR
  }
  if (local) {
    NS.v <- apply(p.m, 1, CompNS)
  } else {
    NS.v <- NULL
  }
  list(sr = SR, inv = invP.v, s = S.v, ns = NS.v, time=Sys.time()-start.time)
}


CompS <- function (p.v) {
  tmp.idx <- which(p.v > 0)
  S <- -sum(p.v[tmp.idx] * log(p.v[tmp.idx]))
  return(S)
}

