## This package replicates the programming code for implementing the test to detect changes in time series tails sequentially.
#  Reference: Hoga, Y., Wied, D. (2017). Sequential monitoring of the tail behavior of dependent data. Journal of Statistical Planning and Inference, 182, 29-49.

# install.packages('sde')
## The VW_limit function simulates the asymptotic limits of the detectors V_t and W_t, given in Theorem 1.
# input: reps - the number of replications to obtain the limiting distribution.
#        n - the size of training sample
#        T - the ratio between the whole sample size and the training sample size
#        kn - tunning parameter, Hoga and Wied (2017) suggests to choose kn = 0.2.
# output: emdis - empirical distribution of the V and W statistics, given the first column for the limit of V_t and the second column for the limit of W_t.
VW_limit<-function(reps,n,T,kn){
  int_approx <- function(x){
        temp_n = length(x)
        return((1/(temp_n)) * sum(x))}
  t0p=kn
  
  vdis=matrix(NA,reps,1)
  wdis=matrix(NA,reps,1)
  for(i in 1:reps){
     bx=BM(x=0, t0=0, T, N=n)
     vnumerator=max((bx[(n/T+n/T*t0p):n]-(((n/T+n/T*t0p):n)/(n/T))*bx[n/T])^2)
     vdenominator=int_approx((bx[(n/T*t0p):(n/T)]-(((n/T*t0p):(n/T))/(n/T))*bx[n/T])^2)
     vdis[i]=vnumerator/vdenominator
     
     wnumerator=max((bx[(n/T+n/T*t0p):n]-bx[(n/T):(n-n/T*t0p)]-t0p*bx[n/T])^2)
     wdenominator=int_approx((bx[(n/T*t0p):(n/T)]-bx[1:(n/T-n/T*t0p)]-t0p*bx[n/T])^2)
     wdis[i]=wnumerator/wdenominator
  }
  emdis=cbind(vdis,wdis)
  return(emdis)
}


## The tail_stat function monitors the stability of the tail behaviour of a time-series sequence of interest.
# input: xtrain - training sample, split from the whole sample
#        xtest - test sample, i.e., whole sample - traning sample
#        kn - tunning parameter, Hoga and Wied (2017) suggests to choose kn = 0.2.
#        p - significance level, e.g., p = 0.05.
# output: retable - test statistics and detected changes for V and W tests.
#         wvec - w statistics
#         emw - the (limiting) threshold for the w test.
# , given a close-end sample x, composed by training sample and testing sample
tail_stat<-function(xtrain,xtest,kn,p){
  n=length(xtrain)
  mn=length(xtest)-n*kn
  x=c(xtrain,xtest)
  T=(length(xtrain)+length(xtest))/length(xtrain)

  emdis=limits(reps=5000,n,T,kn)
  emv=quantile(emdis[,1],0.95)
  emw=quantile(emdis[,2],0.95)

  vvec=matrix(NA,mn,1)
  wvec=matrix(NA,mn,1)
  for(i in 1:mn){
     vvec[i]=Vxdet(x,xtrain,kn,p,i)
     wvec[i]=Wxdet(x,xtrain,kn,p,i)
  }
  
  vind=which(vvec-emv > 0)[1]
  wind=which(wvec-emw > 0)[1]

  vchange=floor(vind+n*kn)
  wchange=floor(wind+n*kn)
  retable=matrix(c(vvec[vind],vchange,wvec[wind],wchange),2,2)

  plot(wvec,main = "",xlab="Date",ylab="Detector",col="blue",lwd=1.5, type="l")
  abline(h=emw, col="goldenrod3",lty=2)

  return(list(retable,wvec,emw))
}
