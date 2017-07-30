############################################################################################
##                                   AUTOCOVARIANCES                                      ##
############################################################################################

## Function to calculate the autocovariances
## pk = gamma(0) = (1/n)*sum_{t=1}^{n}(Y_t-Y)^2

autocovarianzas<-function(s,p,lwd){
  
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
  
  ## Calculation the autocovariance
  num2=numeric(p)
  num2=num/n
  
  ## Calculation the x variable to show the autocovariance of the time series
  x=numeric(p)
  for(i in 1:p){
    x[i]<-(i-1)*(sqrt(n)/n)
  }
  
  ## Ranges
  range1=min(num2)
  range2=max(num2)
  
  ## Show the autocovariances graphic
  plot(x,num2,type="h",lwd=lwd,xlab="Lag",ylab="Covariance",main="Autocovariances",ylim=c(range1,range2))
}

############################################################################################