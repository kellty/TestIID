MC <- 1e3
p <- 5
mu_drift_list <- seq(0,15e-4,15e-5) #11e-4
sigma_change_list <- seq(1,1.5,.05) #1.4
ar_coef_list <- seq(0,1,.1) #0.9
ma_coef_list <- seq(0,10,1) #8

n_list <- c(400,800); len_n <- length(n_list)
m_list <- c('odsup','pca','cp1','cp2','adcr','adcv')
source('func.R')
iid_pval <- function(x,method) {
  if(method=='odsup'){ sup_test(x,function(x1,x2){ exp(-sqrt(sum((x1-x2)^2))) }) }
  else if(method=='pca'){ pca_test(x,2) }
  else if(method=='knn'){ knn_test(x,5) }
  else if(method=='cp1'){ cp_test(x,function(x1,x2){ x1 - x2 }) }
  else if(method=='cp2'){ cp_test(x,function(x1,x2){ x1^2 - x2^2 }) }
  else if(method=='adcr'){ adc_test(x,0.2,TRUE) }
  else if(method=='adcv'){ adc_test(x,0.2,FALSE) }
}
iid_result <- function(model,method,pow,t0) {
  print(paste(model,'|', method,'| n:',n,'| time(secs):',ceiling(proc.time()[3]-t0)))
  print(100*pow)
}


library(foreach)
library(doParallel)
cl <- makeCluster(100) # detectCores()
registerDoParallel(cl)

power_exch <- NULL
x_pop <- t(sapply(1:1000, function(i) rep(i/200,p)))
for(n in n_list) for(method in m_list) {
  t0 <- proc.time()[3]
  pval <- foreach(mc=1:MC, .combine=c) %dopar% {
    set.seed(999+(mc+55)^2)
    x <- x_pop[sample(1000,n,replace=FALSE),]
    iid_pval(x,method)
  }
  pow <- mean(pval<0.05)
  power_exch[[method]] <- c(power_exch[[method]], pow)
  iid_result("exchangeable", method, pow, t0)
}
save(power_exch, file="power_exch.RData")

power_mean <- NULL
for(n in n_list) for(method in m_list) {
t0 <- proc.time()[3]
pval <- foreach(mu_drift=mu_drift_list, .combine=cbind) %:% foreach(mc=1:MC, .combine=c) %dopar% {
  set.seed(999+(mc+55)^2)
  x <- t(sapply(1:n, function(i) rep(i*mu_drift,p) + rnorm(p)))
  iid_pval(x,method)
}
pow <- colMeans(pval<0.05)
names(pow) <- mu_drift_list * 1e4
power_mean[[method]] <- rbind(power_mean[[method]], pow)
iid_result("mean drift", method, pow, t0)
}
save(power_mean, file="power_mean.RData")

power_var <- NULL
for(n in n_list) for(method in m_list) {
t0 <- proc.time()[3]
pval <- foreach(sigma_change=sigma_change_list, .combine=cbind) %:% foreach(mc=1:MC, .combine=c) %dopar% {
  set.seed(999+(mc+55)^2)
  x <- t(sapply(1:n, function(i) ifelse(i>n/2,sigma_change,1) * rnorm(p)))
  iid_pval(x,method)
}
pow <- colMeans(pval<0.05)
names(pow) <- sigma_change_list
power_var[[method]] <- rbind(power_var[[method]], pow)
iid_result("variance changepoint", method, pow, t0)
}
save(power_var, file="power_var.RData")

power_ar <- NULL
for(n in n_list) for(method in m_list) {
t0 <- proc.time()[3]
pval <- foreach(ar_coef=ar_coef_list, .combine=cbind) %:% foreach(mc=1:MC, .combine=c) %dopar% {
  set.seed(999+(mc+55)^2)
  x <- matrix(nrow=n, ncol=p)
  x[1,] <- rnorm(p)
  x[2,] <- rnorm(p)
  for (i in 3:n) {
    x[i,] <- ar_coef * (x[i-1,]-x[i-2,]) + rnorm(p)
  }
  iid_pval(x,method)
}
pow <- colMeans(pval<0.05)
names(pow) <- ar_coef_list
power_ar[[method]] <- rbind(power_ar[[method]], pow)
iid_result("autoregressive model", method, pow, t0)
}
save(power_ar, file="power_ar.RData")

power_ma <- NULL
  for(n in n_list) for(method in m_list) {
t0 <- proc.time()[3]
pval <- foreach(ma_coef=ma_coef_list, .combine=cbind) %:% foreach(mc=1:MC, .combine=c) %dopar% {
  set.seed(999+(mc+55)^2)
  x <- matrix(rnorm(n*p), nrow=n, ncol=p, byrow=TRUE)
  x[-1,] <- x[-1,] + ma_coef * x[-n,]
  iid_pval(x,method)
}
pow <- colMeans(pval<0.05)
names(pow) <- ma_coef_list
power_ma[[method]] <- rbind(power_ma[[method]], pow)
iid_result("moving average model", method, pow, t0)
}
save(power_ma, file="power_ma.RData")


power_ma_mean <- NULL
for(ma_coef in c(7,10)) for(n in n_list) for(method in m_list) {
t0 <- proc.time()[3]
pval <- foreach(mu_drift=10*mu_drift_list, .combine=cbind) %:% foreach(mc=1:MC, .combine=c) %dopar% {
  set.seed(999+(mc+55)^2)
  x <- matrix(rnorm(n*p), nrow=n, ncol=p, byrow=TRUE)
  x[-1,] <- x[-1,] + ma_coef * x[-n,]
  x <- x + t(sapply(1:n, function(i) rep(i*mu_drift,p)))
  iid_pval(x,method)
}
pow <- colMeans(pval<0.05)
names(pow) <- mu_drift_list * 1e4
power_ma_mean[[as.character(ma_coef)]][[method]] <- rbind(power_ma_mean[[as.character(ma_coef)]][[method]], pow)
iid_result("moving average model with mean drift", method, pow, t0)
}
save(power_ma_mean, file="power_ma_mean.RData")

power_ma_var <- NULL
for(ma_coef in c(7,10)) for(n in n_list) for(method in m_list) {
t0 <- proc.time()[3]
pval <- foreach(sigma_change=sigma_change_list, .combine=cbind) %:% foreach(mc=1:MC, .combine=c) %dopar% {
  set.seed(999+(mc+55)^2)
  x <- matrix(rnorm(n*p), nrow=n, ncol=p, byrow=TRUE)
  x[-1,] <- x[-1,] + ma_coef * x[-n,]
  x[-(1:(n/2))] <- sigma_change * x[-(1:(n/2))]
  iid_pval(x,method)
}
pow <- colMeans(pval<0.05)
names(pow) <- sigma_change_list
power_ma_var[[as.character(ma_coef)]][[method]] <- rbind(power_ma_var[[as.character(ma_coef)]][[method]], pow)
iid_result("moving average model with variance changepoint", method, pow, t0)
}
save(power_ma_var, file="power_ma_var.RData")


stopCluster(cl)


model_list <- c('mean','var','ar','ma','ma7_mean','ma10_mean','ma7_var','ma10_var')
power_plot <- function(mdl, len_m=6, is.eps=FALSE) {
  if(mdl==1){param<-mu_drift_list; power<-power_mean; xlab<-expression(paste("mean drift ",mu))}
  if(mdl==2){param<-sigma_change_list; power<-power_var; xlab<-expression(paste("changed std. dev. ",sigma))}
  if(mdl==3){param<-ar_coef_list; power<-power_ar; xlab<-"autoregressive coefficient a"}
  if(mdl==4){param<-ma_coef_list; power<-power_ma; xlab<-"moving average coefficient b"}
  if(mdl==5){param<-10*mu_drift_list; power<-power_ma_mean[['7']]; xlab<-expression(paste("mean drift ",mu))}
  if(mdl==6){param<-10*mu_drift_list; power<-power_ma_mean[['10']]; xlab<-expression(paste("mean drift ",mu))}
  if(mdl==7){param<-sigma_change_list; power<-power_ma_var[['7']]; xlab<-expression(paste("changed std. dev. ",sigma))}
  if(mdl==8){param<-sigma_change_list; power<-power_ma_var[['10']]; xlab<-expression(paste("changed std. dev. ",sigma))}
  if (is.eps) {
    setEPS()
    postscript(paste0('power_',model_list[mdl],'.eps'), width=8*len_n, height=8)
  } else {
    png(paste0('power_',model_list[mdl],'.png'), width=600*len_n, height=600)
  }
  par(mfrow=c(1,len_n), mar=c(5,5,1,1))
  lt <- sapply(1:len_m, function(idx_m) idx_m+(idx_m>2)-(idx_m==6))
  clr <- c('black','#1C9A35','#1B4E8C','#925E9F','#EDB85D','#BB4D0E')
  for (idx_n in 1:len_n) {
    plot(param, power[[m_list[1]]][idx_n,], type='l', lwd=2+is.eps, ylim=c(0,1), xlab=xlab, ylab=bquote(power~(n == .(n_list[idx_n]))))
    for (idx_m in 2:len_m) {
      lines(param, power[[m_list[idx_m]]][idx_n,], lwd=2+is.eps, lty=lt[idx_m], col=clr[idx_m])
    }
    abline(h=0.05, col='lightgray')
    legend('bottomright', m_list[1:len_m], lwd=2+is.eps, lty=lt, col=clr[1:len_m])
  }
  dev.off()
}

sapply(1:length(model_list), function(mdl) power_plot(mdl,length(m_list)))
sapply(1:length(model_list), function(mdl) power_plot(mdl,6,TRUE))

