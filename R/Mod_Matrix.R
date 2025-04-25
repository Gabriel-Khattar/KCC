
#' Modify Species x Site Matrix for Null Models
#'
#' Prepares a species-by-site matrix for randomization by expanding abundance data.
#'
#' @param x A species-by-site matrix (sites as rows, species as columns)
#' @return A data frame with expanded rows by species occurrences
#' @export



Mod_matrix<-function(x){#Function created to modify spXsite in order to facilitate the following randomization.
  x<-x[which(rowSums(x)!=0),which(colSums(x)!=0)]
  Local=row.names(x)
  local_abund<-apply(x,1,sum)
  total_abund<-sum(local_abund)
  names_sp<-colnames(x)
  Elev<-rep(Local,local_abund)
  local<-list()
  for(i in 1:nrow(x)){
    n<-rep(names_sp,as.numeric(x[i,]))
    local[[i]]<-n
  }
  SpR<-unlist(local,use.names=F)
  matrix_mod<-as.data.frame(cbind(Elev,SpR))
}
