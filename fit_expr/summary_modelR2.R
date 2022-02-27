#!/usr/bin/R
result_dir <- "/gpfs/ysm/pi/zhao-data/wl382/cWAS/data/fit_expr"

all_df <- data.frame()
for(chr in 1:22){

  tmp_dir <- paste0(result_dir,"/chr",chr,"/")
  files <- list.files(tmp_dir)
  
  for(f in files){
    tmp_f_dir <- paste0(tmp_dir,f,"/")
    
    txt_files <- list.files(tmp_f_dir,pattern="expr_r2_bothM0M1_M0M1.txt")

    for(txt_f in txt_files){
       tmp_f <- paste0(tmp_f_dir,txt_f) 
       if(!file.exists(tmp_f)){
          next()
       }


       tmp <- read.table(tmp_f,header=T,stringsAsFactors=F,sep="\t")
       tmp$CHR = paste0("chr",chr)
       all_df <- rbind(all_df,tmp)
    }
  }
}

write.table(all_df,file="fit_expr_r2_summary.txt",col.names=T,row.names=F,sep="\t",quote=F)


