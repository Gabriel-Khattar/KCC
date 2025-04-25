
#' Second Smallest Eigenvalue of the Weighted Laplacian for Beta Diversity
#'
#' Calculates the second‐smallest eigenvalue of the weighted graph Laplacian
#' derived from multiple beta‐diversity distance indices.
#'
#' @param x A numeric matrix or data.frame (sites × species).
#' @return A one‐column matrix with row names for each index
#'   (`"chao"`, `"hel"`, `"canberra"`, `"chord"`, `"kulczynski"`)
#'   and the second smallest eigenvalue of the Laplacian.
#' @details
#'   - Uses `vegdist()` (vegan) and `chaodist()` for distance calculations.
#'   - Returns a message and stops if the second eigenvalue is exactly zero
#'    @export



Second_smallest_eign_Weighted_Laplacian<-function(x){
  x<-x[which(rowSums(x)!=0),which(colSums(x)!=0)]
  Index<-c("chao","hel","canberra","chord","kulczynski")
  Matrix_results<-matrix(NA,nrow= length(Index),ncol=1,
                         dimnames=list(c(Index),"Second Smallest Eign"))
  names(Matrix_results)<-Index

  for(i in 1:length(Index)){
    index_used<-Index[i]

    if(index_used=="hel"){
      matrix_dist<-as.matrix(sqrt(2)-vegan::vegdist(decostand(x,"hel"),"euclidean"))
    }else{
      if(index_used=="chord"){
        matrix_dist<-as.matrix(sqrt(2)-vegan::vegdist(x,"chord"))
      }else{
        if(index_used=="Chao_sor"){
          matrix_dist<-as.matrix(1-vegan::chaodist(x, method = "1 - 2*U*V/(U+V)"))
        }else{
          matrix_dist<-as.matrix(1-vegan::vegdist(x,method=index_used))
        }
      }
    }
    Diagonal_laplacian<-rowSums(matrix_dist)
    Lapplacian<--1*matrix_dist
    diag(Lapplacian)<-Diagonal_laplacian

    ### MAtrix spectrum
    Eigenvalues_lap<-eigen(Lapplacian)$values
    Second_small_eignv<-sort(Eigenvalues_lap)[2]
    Matrix_results[index_used,]<-Second_small_eignv
  }

  if(Second_small_eignv==0){
    print("Stopping calculations: Second_small_eignv = 0")
  }

  return(Matrix_results)
}
