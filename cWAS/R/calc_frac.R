#!/usr/bin/R
#library("genio")

#' function to calculate the reference frequency
#'
#' @title calc_frac
#'
#' @param profile_f directory location of profile files
#' @param S_all Signature matrix
#'
#' @export
#'

calc_frac <- function(profile_f,S_all){
  # read imputed expression files
  result <- data.frame()
  files <- list.files(profile_f,pattern="\\.profile")
  for(f in files){
    gene <- gsub("\\.profile","",f)
    tmp_f <- paste0(profile_f,"/",f)

    tmp <- read.table(tmp_f,header=T,stringsAsFactors=F)
    if(length(result)==0){
      result <- tmp[,c("IID","SCORESUM")]
      colnames(result) <- c("IID",gene)
    }else{
      names <- colnames(result)
      result <- cbind(result,tmp[,"SCORESUM"])
      colnames(result) <- c(names,gene)
    }
  }

  pred_expr <- result
  genes <- colnames(S_all)[-1]
  signa <- S_all[,-1]
  sig_genes <- genes
  sig_cells = S_all[,1]

  signa <- t(signa)
  matched_ind <- match(sig_genes,colnames(pred_expr))
  input_extract  <- data.frame(pred_expr[,c(1,na.omit(matched_ind))])
  sig_extrac <- data.frame(signa[match(colnames(input_extract)[-1],sig_genes),])
  cell_fracs <- data.frame()

  for(i in 1:dim(input_extract)[1]){
    test <- cbind(t(input_extract[i,-1]),sig_extrac)
    colnames(test)[1] <- "y"
    a <- lm(y~.,data=test)
    cell_fracs <- rbind(cell_fracs,coef(a)[-1])
    colnames(cell_fracs) <- names(coef(a)[-1])
  }
  colnames(cell_fracs) = sig_cells
  cell_fracs$FID <- pred_expr$IID
  return(cell_fracs)
}
