#' Keystoneness Calculation
#'
#' Runs the Keystoneness framework described in main paper: builds null distributions of impacts,
#' computes mean & SD of null impacts, and returns observed impact,
#' expected impact, SD, and standardized effect size (keystoneness).
#'
#' @param x A numeric matrix or data.frame (sites × species).
#' @param rep Integer; number of null‐model replicates keeping fixed local community size and species regional abundance distribution
#' @param ranking Logical; whether to rank impacts before SES calculation.
#' @param ncores Integer; number of cores for parallel computation.
#' @return A list with components:
#'   - `ObS`: observed impact matrix
#'   - `Expected`: mean null impacts (sites × indices)
#'   - `SD`: null impact standard deviations
#'   - `SES`: (Observed – Expected) / SD
#'
#'@export

Keystoneness_calculation<-function(x,rep,ranking,ncores){
  x<-x[which(rowSums(x)>0),which(colSums(x)>0)]
  names_Index<-c("chao","hel","canberra","chord","kulczynski")
  MM<-Mod_matrix(x)
  List_of_null_impacts<-list()
  Intact<-matrix(0,length(names_Index),rep)
  for(i in 1:rep){
    samp=sample(MM$Elev,length(MM$Elev),replace=F)#
    matrix_null=data.frame("Elevs"=samp,"sp"=MM$SpR)
    matrix_null_beta=reshape2::dcast(matrix_null,Elevs~sp,fun.aggregate =length,value.var="sp")
    rownames(matrix_null_beta)<-matrix_null_beta[,1]
    matrix_null_beta<-matrix_null_beta[,-1]
    Null_Impact_matrix<-Impact_List_function(matrix_null_beta,ranking = ranking,ncores=ncores)
    Null_Impact_matrix[[1]]<-Null_Impact_matrix[[1]][rownames(x),]
    List_of_null_impacts[[i]]<-Null_Impact_matrix[[1]]
    Intact[,i]<-Null_Impact_matrix[[2]]
    print(i)
  }


  #Separating Null impacts per index

  List_of_null_impacts<-lapply(1:length(names_Index),function(cc)sapply(List_of_null_impacts,function(zz)zz[,cc]))
  names(List_of_null_impacts)<-names_Index


  # Ranking and calculating Mean ranked position and SD of ranked position
  Mean_Null_impact<-sapply(List_of_null_impacts,function(xx)rowMeans(xx))
  SD_Null_impact<-sapply(List_of_null_impacts,function(xx)apply(xx,1,sd))


  ObS_Impact<-Impact_List_function(x,ranking = ranking,ncores=ncores)[[1]]



  Keystoneness_per_index<-(ObS_Impact-Mean_Null_impact)/SD_Null_impact

  Results<-list(ObS=ObS_Impact,Expected=Mean_Null_impact,SD=SD_Null_impact,SES=Keystoneness_per_index)
  return(Results)
}
