estDE <-function(data,m,n,PGB.parameters,MRF.parameters,neib.matrix,state,k=40){

  #data=gene expression data corrresponding to 2 groups
  #-first m columns reperesent data for control group and
  #-the next n columns for treatment group
  #m = number of subjects in control group
  #n = number of subjects in treatment group
  #PGB.parameters= vector of parameters for the Poisson-Gamma-Beta model
  #- alpha,kappa,a and b.
  #MRF.parameters=vector of parameters for the MRF model;gamma.i where i={u,d}, beta.1 and beta.2
  #neib.matrix =matrix with 0's and 1's represents the  neighbourhood
  #- information of the pathway genes, i.e neighbours=1, non-neighbours=0
  #state= vector of DE states of genes


  #rate=Gamma rate parameter to calculate negative log likelihood for genes in ctrl
  #-group
  #shape=Gamma shape parameter to calculate the negative loglikelihood for all genes

  #gamma.i where i={u,d} and beta are MRF parameters
  #gamma.i = an arbitraty parameter corresponding to the DE state  i
  #beta.1 = an arbitrary parameter which encourages neighbouring genes to have same
  #-DE states
  #beta.2= an arbitrary parameter which discourages neighbouring genes to have DE different states
  #beta.1 and beta.2 parameters should be always greater than 0


  #calculating log likelihood functions
  #loglik1=log likelihood for genes at state=0
  #loglik2=log likelihood for genes at state=1
  #loglik3=log likelihood for genes at state=-1
  #require("statmod")

  gauss.quad.prob<-statmod::gauss.quad.prob

  #initiating parameter values
  shape<-PGB.parameters[1]
  rate<-PGB.parameters[2]
  a<-PGB.parameters[3]
  b<-PGB.parameters[4]

  gamma.u <- MRF.parameters[1]
  gamma.d <- MRF.parameters[2]
  beta1.u <- MRF.parameters[3]; beta2.u = MRF.parameters[4]
  beta1.d <- MRF.parameters[5]; beta2.d = MRF.parameters[6]
  G=nrow(data)

  #calculating summations over groups
  #yi.m= sumation gene expression for ith gene over control group
  #yi.n= sumation gene expression for ith gene over treatment group
  #yi..= sumation gene expression for ith gene over both groups
  #sum.log.fac.y = summation of log factorial yij for each gene over m+n replicates

  yi.m<-apply(data[,1:m],1,sum,na.rm=T)
  yi.n<-apply(data[,(m+1):(m+n)],1,sum,na.rm=T)
  yi..<-apply(data[,1:(m+n)],1,sum,na.rm=T)
  sum.log.fac.y <- apply(data,1,function(y){
    sum(lfactorial(y),na.rm=T)})
  max.int=exp(700)
  min.int=-exp(700)

  estimated.x<-array(0,dim=c(G,1))


  for( i in 1:G){


    #calculating neihbourhood information for an individual gene
    #neib.ee=neigbours with current state==0
    #neib.ur=neigbours with current state==1
    #neib.dr=neigbours with current state==-1

    neib.ee <- sum(as.matrix(neib.matrix)[i,-i]==1 & as.matrix(state)[-i]==0)
    neib.ur <- sum(as.matrix(neib.matrix)[i,-i]==1 & as.matrix(state)[-i]==1)
    neib.dr <- sum(as.matrix(neib.matrix)[i,-i]==1 & as.matrix(state)[-i]==-1)


    #neib.u= proportion of neibhours with current state==1
    #neib.e= proportion of neibhours with current state==0
    #neib.d= proportion of neibhours with current state==-1

    neib.u<-ifelse(neib.ur!=0,neib.ur/(sum(neib.matrix[i,])-1),0)
    neib.e<-ifelse(neib.ee!=0,neib.ee/(sum(neib.matrix[i,])-1),0)
    neib.d<-ifelse(neib.dr!=0,neib.dr/(sum(neib.matrix[i,])-1),0)

    # calculating conditional probabilities for MRF model compared to the EE state


    A = 0
    B = gamma.u + beta1.u * neib.u - beta2.u * neib.e
    C = gamma.d + beta1.d * neib.d - beta2.d * neib.e



    #yvec1=data for ith gene in control group
    #yvec2=data for ith gene in treatment group
    y<-as.matrix(data[i,])
    yvec1<-as.matrix(data[i,1:m])
    yvec2<-as.matrix(data[i,(m+1):(m+n)])

    #generating gaussian quadrature points using 'statmod' package
    #quad.pt1= k number of gaussian quadrature points from Gamma distribution
    #quad.pt2= k number of gaussian quadrature points from Beta distribution
    #lambda=k number of gaussian quadrature nodes from Gamma distribution
    #delta=k number of gaussian quadrature nodes from Beta distribution

    quad.pt1 <- gauss.quad.prob(k,dist='gamma',alpha=shape,beta=1/rate)
    quad.pt2 <- gauss.quad.prob(k,dist='beta',alpha=a,beta=b)

    #lambda is a vector of guassian nodes from Gamma distribution
    lambda<-as.matrix(quad.pt1$nodes)

    #delta is a vector of guassian nodes and weights from Beta distribution
    delta=quad.pt2

    ################################################
    #calculating loglikelihood for EE genes
    ################################################
    # using closed form formula for Poisson-Gamma joint density
    loglike1<-(log(rate^shape)-lfactorial(shape-1)+lfactorial(yi..[i]+shape-1)-
                 (yi..[i]+shape)*log(m+n+rate)+A)-sum.log.fac.y[i]


    ################################################
    #calculating loglikelihood for UR genes
    ################################################
    matrix2<-apply(as.matrix(lambda),1,ur.mat,delta=delta,y1=yvec1,y2=yvec2)

    matrix2.col.sum<-apply(matrix2,2,sum)
    loglike2<-log(sum(matrix2.col.sum*(quad.pt1$weights)))
    loglike2<-ifelse(loglike2<min.int,-max.int,loglike2)
    loglike2<-loglike2+B

    ################################################
    #calculating loglikelihood for UR genes
    ################################################
    matrix3<-apply(as.matrix(lambda),1,dr.mat,delta=delta,y1=yvec1,y2=yvec2)

    matrix3.col.sum<-apply(matrix3,2,sum)
    loglike3<-log(sum(matrix3.col.sum*(quad.pt1$weights)))
    loglike3<-ifelse(loglike3<min.int,-max.int,loglike3)
    loglike3<-loglike3+C

    #estimate the DE state based on loglikelihood values and group mean
    x.up<-ifelse(loglike2>loglike1 & loglike2>loglike3 & (yi.m[i]/m)< (yi.n[i]/n) ,1,0)
    x.down<-ifelse(loglike3>loglike1 & loglike3>loglike2 & (yi.m[i]/m)> (yi.n[i]/n),1,0)


    #assign state i as optimal DE state for ith gene
    state[i]<-(1*x.up-1*x.down)
  }
  #the vector "state" represent the optimal gene DE states for all genes
  state
}
