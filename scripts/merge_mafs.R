### Merge all maf files in one maf cohort
### 
# real all all_VEP_custom.maf
files <- list.files("/media/rafael/WORK/WES_pipeline/Somatic/results/",recursive = TRUE, pattern = "all_VEP_custom.maf")

library(maftools)
list_mafs <- list()
for (i in 1:length(files)) {
  list_mafs[[i]] <- read.maf(paste0("/media/rafael/WORK/WES_pipeline/Somatic/results/",files[i] ))
}

cohort <- merge_mafs(list_mafs)

