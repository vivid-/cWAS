#!/usr/bin/R
#'
#' functions to calculate intermediate variables
#' @title calc_A
#' @param signature signature matrix
#' @param tol tolerancce of the inverse matrix calculation
#'
#' @export
#'

calc_A <- function(signature,tol){
  library(MASS)
  s <- as.matrix(signature[,-1])
  tmp <- s %*% t(s)
  h <- try(chol2inv(chol(tmp)))
  if(inherits(h,'try-error')){
    h <- ginv(tmp, tol)
  }
  return(h)
}


#'
#' function to get the variance of cell fraction for each
#' cell type in reference database
#' @title calc_Mc
#' @param betaS
#' @param A
#'
#' @export
#'
calc_Mc <- function(betaS,A){
  result_df <- data.frame()
  A <- as.matrix(A)
  betaS <- as.matrix(betaS)
  for(i in 1:ncol(A)){
    #a <- expr_beta[,-c(1:3)]
    #b <- signature[,-c(1)]
    tmp <- betaS %*% A[,i]
    result_df <- rbind(result_df,t(tmp))
  }

  return(result_df)
}



