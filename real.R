source('func.R')
iid_pval <- function(x) {
  pval <- c(sup_test(x,function(x1,x2){ 1 / (sum((x1 - x2)^2)^2 + 1) }),
            cp_test(x,function(x1,x2){ x1 - x2 }),
            cp_test(x,function(x1,x2){ x1^2 - x2^2 }),
            adc_test(x,0.2, FALSE))
  names(pval) <- c('odsup', 'cp1', 'cp2', 'adcv')
  return(pval)
}


mydata <- NULL

# library(dslabs)
# mnist <- read_mnist()
load("mnist.RData")
mydata[['mnist']] <- mnist$train$images[1:800,]
mnist_v <- matrix(nrow=800, ncol=28)
for (i in 1:800) {
  mnist_v[i,] <- svd(matrix(mnist$train$images[i,],nrow=28))$u[,1]
}
mydata[['mnist_v']] <- mnist_v

library(R.matlab)
# https://github.com/matteoiacopini/ZIL-T-MS
email <- readMat('Data_email.mat')
financial <- readMat('Data_financial.mat')
mydata[['financial']] <- t(apply(financial$Xt, 4, c))
mydata[['email']] <- t(apply(email$Xt, 4, c))

set.seed(1234)
print(sapply(mydata, iid_pval))


library(foreach)
library(doParallel)

for (x in mydata) {
  n <- nrow(x)
  rej <- 0
  for (M in 0:49) {
    gc()
    cl <- makeCluster(20)
    registerDoParallel(cl)
    pvals <- foreach(mc=M*20+(1:20), .combine=rbind) %dopar% {
      set.seed((618+mc)^2)
      iid_pval(x[sample(n,replace=TRUE),])
    }
    stopCluster(cl)
    rej <- rej + colSums(pvals<0.05)
  }
  print(rej)
}

