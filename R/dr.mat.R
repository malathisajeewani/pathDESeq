dr.mat <-function(lambda,delta,y1,y2){
        
    logdpois1 <- apply(as.matrix(lambda),1,calc.ldpois,y=y1)
    logdpois.1<-apply(logdpois1,2,sum)
    
    logdpois2 <- apply(as.matrix(lambda*delta$nodes),1,calc.ldpois,y=y2)
    logdpois.2<-apply(logdpois2,2,sum)
    dr.mat<-exp(logdpois.1)*exp(logdpois.2)*(delta$weights)
    dr.mat
  }
