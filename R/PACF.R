############################################################################################
##                                        PACF                                            ##
############################################################################################

## Function to calculate the partial autocorrelation (PACF)
## phi_kk = Corr(Y_t,Y_(t-k)|Y_(t-1),Y_(t-2),...,Y_(t-k+1))

pacf.m<-function(s,p,lwd){

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
  
  ## ## Calculation of the denominator of the correlation coefficient of order k 
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
  
  ## Create a matrix where we add te coefficients
  phi<-matrix(0, nrow=p, ncol=p, byrow=TRUE)
  vn=numeric(p)
  a=0
  for(k in 1:(p-1)){
    a=a+k
    grado<-a
  }
  phi2=numeric(grado)
  
  ## Initial values known
  phi[1,1]=(num[2]/den)
  vn[1]=den
  
  ## Calculation the coefficients
  l=0
  k=0
  h=1
  for(i in 2:p){
    for(j in 1:(i-1)){
      phi2[h]<-phi[i-1,j]*num[i+1-j]
      h=h+1
    }
    vn[i]<-vn[i-1]*(1-((phi[i-1,i-1])^2))
    phi[i,i]<-(num[i+1]-(sum(phi2[(i-1+l):(i-1+k)])))/vn[i]
    if(i!=j){
      for(j in 1:(i-1)){
        phi[i,j]<-phi[i-1,j]-(phi[i,i]*phi[i-1,i-j])
      }
    }
    l=l+i-2
    k=k+i-1
  }
  
  ## Create a list where add diagonal matrix which have the coefficients
  phikk=numeric(p)
  phikk=diag(phi)
  
  ## If there are NA values, we allocate a zero
  phikk[is.na(phikk)] <- 0
  
  ## Calculation the x variable to show the pacf of the time series
  x=numeric(p)
  for(i in 1:p)
    x[i]<-(i)*(sqrt(n)/n)
  }
  
  ## Boundary confiance
  lim1=(2/sqrt(n))
  lim2=(-2/sqrt(n))
  
  ## Ranges
  if(min(phikk)<lim2){
    range1=min(phikk)
  }else{
    range1=lim2
  }
  range2=max(phikk)
  
  ## Show the pacf graphic
  plot(x,phikk,type="h",lwd=1,xlab="Lag",ylab="Partial ACF",main="Partial Autocorrelation",ylim=c(range1,range2))
  abline(h=lim1,col="blue",lty=2)
  abline(h=lim2,col="blue",lty=2)
  abline(h=0,col="black",lty=1)
}

############################################################################################