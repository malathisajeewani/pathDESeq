pgbEst <-
  function(data,state,m,n){
    
    gamma.est_MM <- function(x) {
      mu <- mean(x); v <- var(x) 
      return(c(mu^2/v,mu/v))
    }
    
    #estimate Gamma parameters considering only EE genes
    data.0<-subset(data,state==0)
    lambda<-apply(data.0,1,mean,na.rm = TRUE)
    gamma_est<-gamma.est_MM (lambda)
    
    #estimate Beta parameters considering DR genes
    data.dr<-subset(data,state==-1)
    ctl<-data.dr[,1:m]
    trt<-data.dr[,(m+1):(m+n)]
    trt.mean<-apply(trt,1,mean,na.rm = TRUE)
    ctl.mean<-apply(ctl,1,mean,na.rm = TRUE)
    
    #Estimate the Beta random variable "delta" as a fraction of group means 
    delta<-ifelse(trt.mean/ctl.mean==Inf,NA,trt.mean/ctl.mean)
    max.int<-exp(700)
    delta<-ifelse(delta>max.int,max.int,delta)
    mu1<-mean(delta,na.rm=TRUE)
    v1<-var(delta,na.rm=TRUE)
    a1<- (mu1*(1-mu1)-v1)*mu1/v1
    b1<-a1/mu1-a1
    a1;b1
    
    #estimate Beta parameters considering UR genes
    data.ur<-subset(data,state==1)
    ctl<-data.ur[,1:m]
    trt<-data.ur[,(m+1):(m+n)]
    
    trt.mean<-apply(trt,1,mean,na.rm = TRUE)
    ctl.mean<-apply(ctl,1,mean,na.rm = TRUE)
    
    #Estimate the Beta random variable "delta" as a fraction of group means 
    delta<-ifelse(ctl.mean/trt.mean==Inf,NA,ctl.mean/trt.mean)
    delta<-ifelse(delta>max.int,max.int,delta)
    mu2<-mean(delta,na.rm=TRUE)
    v2<-var(delta,na.rm=TRUE)
    a2<- (mu2*(1-mu2)-v2)*mu2/v2
    b2<-a2/mu2-a2
    a2;b2
    
    #take the average
    beta_est<-c(mean(c(a1,a2)),mean(c(b1,b2)))
    
    #create an object with parameter estimates
    #est<-c(Gamma_Shape=gamma_est[1],Gamma_Rate=gamma_est[2],Beta_Shape1=beta_est[1],Beta_Shape2=beta_est[2])
    
    return(c(gamma_est,beta_est))
    
  }
