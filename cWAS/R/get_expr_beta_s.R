#!/usr/bin/R

#' function to match imputation weights and
#'  get the products of weights and signature
#'  expression
#'
#' @title get_expr_beta_s
#'
#' @param expr_beta imputation weights after matching
#' @param signature signature matrix after matching
#'
#' @export
#'

get_expr_beta_s <- function(expr_beta,signature){
  # match gene in expr_beta and signature
  inter_gene <- intersect(colnames(expr_beta)[-c(1:3)],colnames(signature)[-1])

  # match intersected genes
  expr_beta_1 <- expr_beta[,c(1:3,match(inter_gene,colnames(expr_beta)))]
  signature_1 <- signature[,c(1,match(inter_gene,colnames(signature)))]
  #\hat\beta %*% S
  S <- t(signature_1[,-1])
  beta <- expr_beta_1[,-c(1:3)]
  result <- data.frame(as.matrix(beta) %*% as.matrix(S))
  result$rsID <- expr_beta_1[,1]
  result$ref <- expr_beta_1[,2]
  result$alt <- expr_beta_1[,3]
  colnames(result) <- c(as.character(signature_1[,1]),"rsID","ref","alt")

  return(result)
}
