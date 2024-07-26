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

air <- as.matrix(read.csv("PRSA_Data_Aotizhongxin_20130301-20170228.csv"))
air <- air[rowSums(is.na(air))==0,]
mydata[['air']] <- apply(air[1:400,6:11], 2, as.numeric)

# library(dslabs)
# mnist <- read_mnist()
# mnist <- mnist$train$images[1:800,]
load("mnist.RData")
mydata[['mnist']] <- mnist[1:800,]
mnist_v <- matrix(nrow=800, ncol=28)
for (i in 1:800) {
  mnist_v[i,] <- svd(matrix(mnist[i,],nrow=28))$u[,1]
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

