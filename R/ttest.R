#This function will perform the independent two sample t test and will store the t statistics and p values.
#m= number of replicates for control group
#n= number of replicates for treatment group
#sig= level of significance

#t.test fails when both groups have constant counts. This function would still work in those situations.
ttest<-function(data,m,n,sig=0.05){
  G=nrow(data)
  tval<-array(0,dim=c(G,1))
  pval<-array(0,dim=c(G,1))

  constant.count=0
  tval<-apply(data,1,function(x){
    tryCatch(t.test(x[1:m],x[(m+1):(m+n)],"two.sided",conf.level = 1-sig)$statistic,
             error = function(e) sign(mean(x[1:m]) - mean(x[(m+1):(m+n)]))*Inf)

  })
  constant.count <- constant.count + sum(tval == -Inf | tval == Inf)
  pval<-apply(data,1,function(x){
    tryCatch(t.test(x[1:m],x[(m+1):(m+n)],"two.sided",conf.level = 1-sig)$p.value,
             error = function(e) as.numeric((mean(x[1:m]) == mean(x[(m+1):(m+n)]))*1))

  })
  return(data.frame(tval,pval))
}
