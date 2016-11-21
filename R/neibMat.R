neibMat<-function(pathway.genes,interactions){
  interactions$Gene.1=toupper(interactions$Gene.1)
  interactions$Gene.2=toupper(interactions$Gene.2)
  subTable <- interactions[(interactions[,1] %in% pathway.genes) & (interactions[,2] %in% pathway.genes),]
  single.genes <- data.frame(Gene.1=pathway.genes,Gene.2=pathway.genes)
  subTable = rbind(subTable,single.genes)
  colnames(subTable) <- c('Source','Target')

  subTable.sort <- t(apply(subTable,1,sort))
  subTable.sort=as.data.frame(subTable.sort)
  colnames(subTable.sort)=colnames(subTable)
  #neib.matrix conntain the neib matrix
  neib.matrix=adj2mat(subTable.sort)
  #store neib.matrix as a text file
  write.table(neib.matrix,"neib_matrix.txt")
  neib.matrix
}
