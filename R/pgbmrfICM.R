pgbmrfICM<-function(genes,data,interactions,m,n,sig=0.05,k=40,pgb.start=c(log(10),log(0.2),log(2),log(3)),iterations=12){

  t.results<-ttest(data=data,m=m,n=n,sig=sig)
  tval<-t.results[,1]
  pval<-t.results[,2]


  #state1= 1 if state==1, 0 otherwise
  #state2= 1 if state==-1, 0 otherwise
  state1<-ifelse(tval<=(qt(sig/2,(m+n-2))),1,0)
  state2<-ifelse(tval>=(-qt(sig/2,(m+n-2))),-1,0)

  #vector "t.state" represents the initial DE states of genes as three states(-1,0,1)
  t.state<-state1+state2
  #obtain a summary of t test DE states
  #table(t.state)

  #update final data set by including initial DE states
  data.final<-data.frame(genes,data,t.state)
  #write.table(t.state,"t state.txt")

  #sort data frame by gene names
  #gene expression data and the neib.matrix (in step 3) should be in the same gene order.
  data.final <- data.final[order(data.final[,1]),]
  write.table(data.final,"selected dataset.txt")
  #dim(data.final)
  #store selected gene names for the analysis
  selected.Reactome.genes<-data.final[,1]
  #store selected gene expression data for the analysis
  data<-data.final[,2:(m+n+1)]
  #store t states
  t.state<-data.final[,(m+n+2)]
  #table(t.state)

  #adjust pvalue using FDR
  pval.adj <- p.adjust(pval, method='fdr')
  # how many genes significant at FDR=0.05
  #sum(pval.adj<sig, na.rm=T)

  #The vector "t.state.fdr" represents the DE states of genes as
  #`three states(0,1,-1) using Benjamini & Horchberg method.
  t.state.fdr1<-ifelse(pval.adj<sig & tval<=(qt(sig/2,(m+n-2))),1,0)
  t.state.fdr2<-ifelse(pval.adj<sig & tval>=(-qt(sig/2,(m+n-2))),-1,0)
  t.state.fdr<-t.state.fdr1+t.state.fdr2
  #obtain a summary of t test DE states corrected for FDR
  #table(t.state.fdr)

  ttest.result<-data.frame(selected.Reactome.genes,tval,pval,t.state,pval.adj,t.state.fdr)
  #store t test summary table for further use
  write.table(ttest.result,"ttest result.txt")

  #After filtering the data we can select a list of unique genes for futher analysis
  #selected.Reactome.genes=list of unique gene names mapped with Reactome pathway genes
  length(selected.Reactome.genes)
  cat(paste(length(selected.Reactome.genes)," genes were selected for the PGBMRF analysis."),"\n")
  #save as a character vector
  selected.Reactome.genes<-as.character(selected.Reactome.genes)
  cat(paste("Creating the neigbourhood matrix...."),"\n")
  neib.matrix<-neibMat(pathway.genes=selected.Reactome.genes,interactions=interactions)

  #assing t.states as initial DE states for genes
  state<-t.state
  #create a data frame to store the summary result
  summary<-data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0)
  colnames(summary)<-c("kappa","alpha","a","b","gamma.u","gamma.d",
                       "beta1.u","beta2.u","beta1.d","beta2.d","EE.genes", "UR.Genes","DR.Genes")

  #creating a dataframe to store results with initial parameters
  result<-data.frame(0,0,0,0,0,0,0,0,0,0,sum(t.state==0),sum(t.state==1), sum(t.state==-1))
  colnames(result)<-c("kappa","alpha","a","b","gamma.u","gamma.d",
                      "beta1.u","beta2.u","beta1.d","beta2.d","EE.genes", "UR.Genes","DR.Genes")

  #Let's perform the parameter estimation for the PGBMRF model using ICM algorithm.
  #Follow the iterative steps of optimization untill the convergence of estimated DE states


  #1. estimate Poisson-Gamma-Beta parameters using Method of Moment(MoM)
  #2. estimate Markov Random Feild model parameters using MLE's
  #3. estimate optimal DE states for given PGBMRF model parameters
  #4. follow step 1-3 until the convergence of estimated DE states
  old.state<-c(rep(0,length(state)))
  conv1 <- FALSE
  conv2 <- FALSE

  #ICM algorithm
  while(conv1==FALSE & conv2==FALSE){
    old.state1<-old.state
    old.state<-state
    #calculate PGB MOM estimates using "pgbEst" function
    #- a and b for given DE states
    #parameter_estimates1[1] = MOM estimate for kappa
    #parameter_estimates1[2] = MOM estimatefor alpha
    #parameter_estimates1[3] = MOM estimate for a
    #parameter_estimates1[4] = MOM estimate for b

    cat(paste('Running ICM iteration',(nrow(result)),"...",sep=" "),"\n")
    #step 1: parameter estimation for PGB using MOM
    parameter_estimates1<-pgbEst(data=data,state=state,m=m,n=n)
    #print(parameter_estimates1)

    #step 2: parameter estimation for MRF model using Maximum Likelighood method
    #optimize the "mrf_nloglik" function and obtaining MLE's for Markov Random Field model.
    #parameter_estimates2[1] = MLE's for gamma.u
    #parameter_estimates2[2] = MLE's for gamma.d
    #parameter_estimates2[3] = MLE's for beta.1
    #parameter_estimates2[4] = MLE's for beta.2

    #NEW method for estimating MRF pars
    neib.ur = neib.matrix %*% as.matrix(state == 1)
    neib.ur = neib.ur - as.numeric(state == 1)
    neib.ee = neib.matrix %*% as.matrix(state == 0)
    neib.ee = neib.ee - as.numeric(state == 0)
    neib.dr = neib.matrix %*% as.matrix(state == -1)
    neib.dr = neib.dr - as.numeric(state == -1)
    #prop1= proportion of DR genes in the neigborhood
    #prop2=proportion of EE genes in the neigborhood
    #prop3=proportion of UR genes in the neigborhood


    prop1 = ifelse(neib.dr != 0, neib.dr/(apply(neib.matrix, 1, sum) - 1), 0)
    prop2= ifelse(neib.ee != 0, neib.ee/(apply(neib.matrix, 1, sum) - 1), 0)
    prop3 = ifelse(neib.ur != 0, neib.ur/(apply(neib.matrix, 1, sum) - 1), 0)
    MRF.data <- data.frame(state,prop1,prop2,prop3)


    model1 <- glm(abs(state)~prop1+prop2,family='binomial',subset=state<1,data=MRF.data)
    model2 <- glm(abs(state)~prop3+prop2,family='binomial',subset=state>=0,data=MRF.data)
    #according to the model 1
    gamma.d <- coef(model1)[1]
    beta1.d <-  coef(model1)[2]
    beta2.d <- -coef(model1)[3]

    #according to the model 2
    gamma.u <- coef(model2)[1]
    beta1.u <-  coef(model2)[2]
    beta2.u <- -coef(model2)[3]

    parameter_estimates2<-c(gamma.u,gamma.d,beta1.u,beta2.u,beta1.d,beta2.d)


    #print(parameter_estimates2)


    #step 3:estimating DE states for PGBMRF model using PGB and MRF estimates
    state<-estDE(data=data,m=m,n=n,PGB.parameters=parameter_estimates1,
                 MRF.parameters=parameter_estimates2,neib.matrix=neib.matrix,state=state,k=k)

    #store the results of current iteration as current.results dataframe
    current.result<-data.frame(t(parameter_estimates1),t(parameter_estimates2),sum(state==0),sum(state==1),sum(state==-1))
    colnames(current.result)<-c("kappa","alpha","a","b","gamma.u","gamma.d",
                                "beta1.u","beta2.u","beta1.d","beta2.d","EE.genes", "UR.Genes","DR.Genes")
    #print(current.result)
    #combine the results of current iteration with previous results
    result<-rbind(result,current.result)


    # check for convergence of DE states
    crit <- sum(old.state!=state)/length(state)
    crit1<-sum(old.state1!=state)/length(state)

    # if < 0.1% of genes change state, it converges
    conv1 <- ifelse((crit< 0.001 |crit1< 0.001),TRUE,FALSE)
    conv2 <- ifelse(nrow(result)>iterations,TRUE,FALSE)

  }
  cat(paste("ICM algorithm converged after",(nrow(result)-1),"iterations.See the following result table."),"\n")

  print(result[-1,])
  #write result table as a text file
  write.table(result[-1,],"PGBMRF results.txt")
  #write final DE states as a text file
  write.table(state,"PGBMRF states.txt")
  UR.genes<-selected.Reactome.genes[state==1]
  DR.genes<-selected.Reactome.genes[state==-1]
  write.table(UR.genes,"PGBMRF identified UR genes.txt")
  write.table(DR.genes,"PGBMRF identified DR genes.txt")
  cat(paste('See ',getwd()," for details.",sep=" "),"\n")

}
