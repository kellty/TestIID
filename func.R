sup_test <- function(X,h) {
  X <- as.matrix(X)
  n <- nrow(X)
  hxx <- sapply(1:n, function(i) sapply(1:n, function(j) h(X[i,],X[j,])))
  diag(hxx) <- 0; hxx <- hxx - sum(hxx)/(n*(n-1));  diag(hxx) <- 0
  dist_idx <- abs(row(hxx) - col(hxx))
  T_0 <- max(abs(sapply(1:n, function(d) sum(hxx[dist_idx<d]))))#/(n*(n-1))
  rho <- 2*sapply(1:n, function(k) sapply(1:(n-1), function(d) min(d,k-1)+min(d,n-k)))
  rho_0 <- rho - (2/n) * (1:(n-1))*((2*n-2):n)
  hx1 <- rowSums(hxx)/(n-1)
  T_b <- replicate(1e3, max(abs(rho_0 %*% (hx1 * rnorm(n)))))
  return(pval = 1 - mean(T_0 >= T_b))
}


pca_test <- function(X,q) {
  X <- as.matrix(X)
  n <- nrow(X); n_half <- n/2
  idx <- 2*(1:n_half)
  pca <- eigen(t(X[idx,])%*%X[idx,])$vectors[,1:q]
  T_0 <- sum((X[-idx,] - X[-idx,]%*%pca%*%t(pca))^2)
  T_perm <- replicate(1e2, {
    idx <- sample(n, n_half)
    pca <- eigen(t(X[idx,])%*%X[idx,])$vectors[,1:q]
    return(sum((X[-idx,] - X[-idx,]%*%pca%*%t(pca))^2))
  })
  return(pval = 1 - mean(T_0 >= T_perm))
}


cp_test <- function(X,h) {
  X <- as.matrix(X)
  n <- nrow(X)
  hxx <- sapply(1:n, function(i) sapply(1:n, function(j) h(X[i,],X[j,])), simplify='array')
  hx_s <- sapply(1:(n-1), function(i) rowSums(as.matrix(hxx[,i,(i+1):n])))
  T_0 <- max(abs(rowSums(as.matrix(hx_s))))
  T_b <- replicate(1e3, max(abs(hx_s %*% rnorm(n-1))))
  return(pval = 1 - mean(T_0 >= T_b))
}


adc_test <- function(X,lambda, is.corr=TRUE) {
  X <- as.matrix(X)
  n <- nrow(X)
  maxlag <- floor(3*n^lambda)
  # z <- (1:(n-1))/maxlag
  # kern <- ifelse( abs(z) <= 1, 1 - abs(z), 0 )  #Bartlett
  crossdist <- sapply(1:n, function(i) sapply(1:n, function(j) abs(X[i,]-X[j,])), simplify='array')
  if (is.corr) {
    A_mean1 <- rowMeans(crossdist, dims=2)
    A_mean2 <- rowMeans(A_mean1)
    D2 <- rowMeans(sapply(1:n, function(i) sapply(1:n, function(j)
      (crossdist[,i,j] - A_mean1[,i] - A_mean1[,j] + A_mean2)^2), simplify='array'), dims=1)
    for (k in 1:ncol(X)) {
      crossdist[k,,] <- crossdist[k,,] / D2[k]^(1/2)
  }}
  AB_list <- NULL
  for (l in 1:maxlag) {
    AB <- matrix(nrow=n-l, ncol=n-l)
    A_mean1 <- rowMeans(crossdist[,1:(n-l),1:(n-l)], dims=2)
    A_mean2 <- rowMeans(A_mean1)
    B_mean1 <- rowMeans(crossdist[,(l+1):n,(l+1):n], dims=2)
    B_mean2 <- rowMeans(B_mean1)
    for(i in 1:(n-l)) for(j in 1:(n-l)) {
      Aij <- crossdist[,i,j] - A_mean1[,i] - A_mean1[,j] + A_mean2
      Bij <- crossdist[,i+l,j+l] - B_mean1[,i] - B_mean1[,j] + B_mean2
      AB[i,j] <- sum(Aij) * sum(Bij)
    }
    AB_list[[l]] <- AB
  }
  adcor2sum <- lapply(AB_list, mean)
  T_0 <- sum(sapply(1:maxlag, function(l) (n-l) * (1-l/maxlag)^2 * adcor2sum[[l]]))
  T_b <- replicate(1e2, {
    adcor2sum_b <- lapply(AB_list, function(AB){
      w <- rnorm(nrow(AB))
      mean(diag(w)%*%AB%*%diag(w))
    })
    gc()
    return(sum(sapply(1:maxlag, function(l) (n-l) * (1-l/maxlag)^2 * adcor2sum_b[[l]])))
  })
  return(pval = 1 - mean(T_0 >= T_b))
}

