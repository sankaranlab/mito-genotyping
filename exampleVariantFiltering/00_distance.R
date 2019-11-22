distMito <- function(af3, cov3, minCov = 100){
  ncells <- dim(af3)[2]
  stopifnot(dim(cov3)[2] == ncells)
  
  covboo <- cov3 > minCov
  ncomp_mat <- dist_mat <- matrix(0, nrow = ncells, ncol = ncells)
  rownames(ncomp_mat) <- colnames(ncomp_mat) <- rownames(dist_mat) <- colnames(dist_mat) <- colnames(af3)
  for(i in 1:ncells){
    for(j in 1:ncells){
      ncomp <- sum(covboo[,i] * covboo[,j])
      ncomp_mat[i,j] <- ncomp
      dist_mat[i,j] <- sum(abs(sqrt(af3[,i]*covboo[,i]) - sqrt(af3[,j]*covboo[,j])))/ncomp
    }
  }
  return(list(dist_mat, ncomp_mat))
}