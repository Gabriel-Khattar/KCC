#' Impact_List_function

#' @param x A numeric matrix or data.frame (sites × species).
#' @param ranking Logical; whether to rank impacts before SES calculation.
#' @param ncores Integer; number of cores for parallel computation.
#' @return Impact list.#Note that the calculations performed in the function Impact_List_function should not be considered as estimates of community keystoneness. See main text.
#'
#'@export

##x= Species by site matrix
#Ranking is a Boolean to decide whether to rank communities in the impact list or not
#Ranking is advisable because improves statistical tractability. See main text.
#ncores: defines the number of cores used to estimate keystoneness.
#Depending on the size of the sp x sites matrix, it can be useful
Impact_List_function<-function(x,ranking=T,ncores){
  x<-x[which(rowSums(x)!=0),which(colSums(x)!=0)]
  Intact<-Second_smallest_eign_Weighted_Laplacian(x)
  if(ncores==1){
    Sencond_eign<-lapply(1:nrow(x),function(i)Second_smallest_eign_Weighted_Laplacian(x[-i,]))
  }else{
    Sencond_eign<-parallel::mclapply(1:nrow(x),function(i)Second_smallest_eign_Weighted_Laplacian(x[-i,]),mc.cores = ncores)

  }

  names(Sencond_eign)<-rownames(x)

  Impact_matrix<-t(sapply(Sencond_eign,function(x)(Intact-x),USE.NAMES = T))
  colnames(Impact_matrix)<-rownames(Intact)

  if(ranking==T){
    Impact_matrix<-apply(Impact_matrix,2,rank)
    Results<-list(Impact_matrix,Intact)
  }else{
    Results<-list(Impact_matrix,Intact)
  }
  return(Results)
}
