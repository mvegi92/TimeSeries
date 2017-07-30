############################################################################################
##                                         ACF                                            ##
############################################################################################

## Function to calcule simple autocorrelation (acf)
## pk = gamma(k)/gamma(0) = Corr(Yt,Yt-k) = sum_{t=k+1}^{n}(Y_t-Y)(Y_(t-k)-Y)/sum_{t=1}^{n}(Y_t-Y)^2

acf.m<-function(s,p,lwd){
  
  n=length(s)
  mu<-mean(s)
  j=n*p
  num_pk1<-numeric(j)
  num=numeric(p)
  
  ## Calculation of the numerator of the correlation coefficient of order k
  g=1           
  for(k in 0:(p-1)){
    for(t in 1:(n-k)){
      num_pk1[g]<-((s[t]-mu)*(s[t+k]-mu))
      g=g+1
    }
  }
  
  a=1
  b=n
  i=1
  for(j in 0:(p-1)){
    num[i]<-sum(num_pk1[a:b])
    a=a+n-j
    b=a+n-j-2
    i=i+1
  }
  
  ## We check that the length of our vector is equal to the order k
  if(length(num)==p){
    TRUE
  }else{
    FALSE
  }
  
  ## Calculation of the denominator of the correlation coefficient of order k
  den_pk1=numeric(n)
  for(t in 1:n){
    den_pk1[t]<-((s[t]-mu)^2)
  }
  den=sum(den_pk1)
  
  ## Calculation of the correlation coefficient of order k
  Rk=numeric(p)
  for(k in 1:p){
    Rk[k]<-num[k]/den
  }
  
  ## Calculation the x variable to show the acf of the time series
  x=numeric(p)
  for(i in 1:p){
    x[i]<-(i-1)*(sqrt(n)/n)
  }
  
  ## Boundary confiance
  lim1=(2/sqrt(n))
  lim2=(-2/sqrt(n))
  
  ## Ranges
  if(min(Rk)<lim2){
    range1=min(Rk)
  }else{
    range1=lim2
  }
  range2=max(Rk)
  
  ## Show the acf graphic
  plot(x,Rk,type="h",lwd=lwd,xlab="Lag",ylab="ACF",main="Simple Autocorrelation",ylim=c(range1,range2))
  abline(h=lim1,col="blue",lty=2)
  abline(h=lim2,col="blue",lty=2)
  abline(h=0,col="black",lty=1)
  
}

############################################################################################
