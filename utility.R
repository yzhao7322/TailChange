## This file provides the functions for the feasibility of the main functions in cv_stat.R
# Sequential monitoring of the tail behavior of dependent data
# Hoga and Wied (2016), JSPI.
#library('sde')

### detectors statistics
tailindex <- function(x,kn){
  n=length(x)
  k=round(n*kn)
  orderx=sort(x,decreasing = TRUE)
  taili=0

  for(i in 1:k){
  	taili=taili+log(orderx[i]/orderx[k+1])
  }
  tai=taili/k
  return(tai)
}

# for equation 7
Vdet <- function(x,xtrain,kn,tind){
  int_approx <- function(x){
        temp_n = length(x)
        return((1/(temp_n)) * sum(x))}

  n=length(xtrain)
  t0=kn
  xtest=x[(n+1):((n+n*kn)+tind)]
  
  numerator=(((n+n*kn+tind)/n-1)*(tailindex(xtest,kn)-tailindex(xtrain,kn)))^2
  numx=matrix(NA,(n-n*kn),1)
  for(s in floor(n*t0+1):n){
      sf=s/n
      xs=x[1:s]
      numx[s-floor(n*t0)]=(sf*(tailindex(xs,kn)-tailindex(xtrain,kn)))^2
  }
  denominator=int_approx(numx)
  return(numerator/denominator)
}


Wdet <- function(x,xtrain,kn,tind){
  int_approx <- function(x){
        temp_n = length(x)
        return((1/(temp_n)) * sum(x))}

  n=length(xtrain)
  t0=kn
  xtest=x[((n+n*kn)+tind-n*kn):((n+n*kn)+tind)]
  
  numerator=(t0*(tailindex(xtest,kn)-tailindex(xtrain,kn)))^2
  numx=matrix(NA,(n-n*t0),1)
  for(s in floor(n*kn+1):n){
      sf=s/n
      xs=x[(s-floor(n*kn)):s]
      numx[s-floor(n*kn)]=(t0*(tailindex(xs,kn)-tailindex(xtrain,kn)))^2
  }
  denominator=int_approx(numx)
  return(numerator/denominator)
}


# for equation 8
Vxdet <-function(x,xtrain,kn,p,tind){
  int_approx <- function(x){
    temp_n = length(x)
    return((1/(temp_n)) * sum(x))}

  n=length(xtrain)
  t0=kn

  xp<-function(data,p,n,k,s,t){
     xk=sort(data[floor(s*n+1):floor(t*n)],decreasing = TRUE)[floor(k*(t-s)*1)+1]
     return(xk*(n*p/k)^(-tailindex(data[floor(s*n+1):floor(t*n)],k/n)))
  }

  numerator=(((n+n*kn+tind)/n-1)*log(xp(x,p,n,n*kn,1,(n+n*kn+tind)/n)/xp(x,p,n,n*kn,0,1)))^2
  numx=matrix(NA,(n-n*kn),1)
  for(s in floor(n*t0+1):n){
      sf=s/n
      xs=x[1:s]
      numx[s-floor(n*t0)]=(sf*log(xp(x,p,n,n*kn,0,sf)/xp(x,p,n,n*kn,0,1)))^2
  }
  denominator=int_approx(numx)
  return(numerator/denominator)
}

Wxdet <- function(x,xtrain,kn,p,tind){
  int_approx <- function(x){
        #temp_n = length(!is.na(x))
        #return((1/(temp_n)) * sum(x,na.rm = TRUE))}
        temp_n = length(x)
        return((1/(temp_n)) * sum(x))}
  n=length(xtrain)
  t0=kn
  
  xp<-function(data,p,n,k,s,t){
     xk=sort(data[floor(s*n+1):floor(t*n)],decreasing = TRUE)[floor(k*(t-s)*1)+1]
     
     #xorder=sort(data[floor(s*n+1):floor(t*n)],decreasing = TRUE)
     #pk=which(xorder <= 0)[1]-1
     #if (xk<=0){xk=xorder[pk]}
     
     return(xk*(n*p/k)^(-tailindex(data[(floor(s*n)+1):floor(t*n)],k/n)))
  }

  numerator=(t0*log(xp(x,p,n,n*kn,(n+tind)/n,(n+n*kn+tind)/n)/xp(x,p,n,n*kn,0,1)))^2
  numx=matrix(NA,(n-n*t0),1)

  for(s in floor(n*kn+1):n){
      sf=s/n
      xs=x[(s-floor(n*kn)):s]
      numx[s-floor(n*kn)]=(t0*log(xp(x,p,n,floor(n*kn),sf-t0,sf)/xp(x,p,n,n*kn,0,1)))^2
  }

  denominator=int_approx(numx)
  return(numerator/denominator)
}


