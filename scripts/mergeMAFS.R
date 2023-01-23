library(maftools)
library(R.utils)
my_path <- "/media/rafael/WORK/WES_pipeline/Somatic/results"
setwd(my_path)
group <- list.dirs(my_path, full.names = FALSE, recursive = FALSE)
for (i in 1:length(group)) {
  directory <- paste(my_path, "/", group[i], "/Merge/maf", sep = "")
  setwd(directory)
  nam=paste("all_VEP_custom",i,sep="")
  MAF <- assign(paste0('all_VEP_custom',i,sep=""),read.maf(maf='all_VEP_custom.maf'))
} 

merge_mafs(maf= MAF, verbose = TRUE)

                 