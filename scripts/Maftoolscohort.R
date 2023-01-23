#####Previous steps VCF2MAF snakefile##########
##########################
library(maftools)
library(R.utils)
library(argparse)
#####
parser <- ArgumentParser(description="Directory")
parser$add_argument("--directory",help="put here your input directory")
parser$print_help()
args <- parser$parse_args()
directory <- args$directory #for this script the directory will be in results
######
getwd()
setwd(directory)
dir.create("ResultsMaftoolsCohort")
### Merge all maf files in one maf cohort
# real all all_VEP_custom.maf
files <- list.files("/media/rafael/WORK/WES_pipeline/Somatic/results/",recursive = TRUE, pattern = "all_VEP_custom.maf")
list_mafs <- list()
for (i in 1:length(files)) {
  list_mafs[[i]] <- read.maf(paste0("/media/rafael/WORK/WES_pipeline/Somatic/results/",files[i] ))
}
cohort <- merge_mafs(list_mafs)
getwd()
setwd(paste0(directory,"/ResultsMaftoolsCohort/"))
#Create list of 10 most mutated genes
gene <- getGeneSummary(cohort)
head(gene)
geneorder <- gene[order(-gene$total)]
list <- head(geneorder,10)
genelist <- list[,1]
genelist
gene.vec <- as.vector(genelist$Hugo_Symbol)
length(gene.vec)
#We don't have the TCGA barcode so we use this function to extract it
for (i in seq_along(gene.vec)) {
  print(genesToBarcodes(maf = cohort, genes = gene.vec[i]))
}
getSampleSummary(cohort)
##Shows gene summary.
getGeneSummary(cohort)
#shows clinical data associated with samples
getClinicalData(cohort)
cohort
#Shows all fields in MAF
#getFields(x = cohort) #clinicalinformation
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = cohort, basename = 'mafall')
cohort
########VISUALIZATION##########
png(paste0(directory,"/ResultsMaftoolsCohort/plotmafSummary"))
plotmafSummary(maf = cohort, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
#For all the cohort
#Oncoplot
#oncoplot for top ten mutated genes.
cohort
png(paste0(directory,"/ResultsMaftoolsCohort/oncoplot"))
suppressWarnings(oncoplot(maf = cohort, top = 10, showTumorSampleBarcodes = TRUE))
dev.off()
#Transition and Transversions
MAF.titv = titv(maf = cohort, plot = FALSE, useSyn = TRUE)
#plot titv summary
png(paste0(directory,"/ResultsMaftoolsCohort/plottitv"))
plotTiTv(res = MAF.titv)
dev.off()
#Lollipop plots for amino acid changes for 10 most mutated genes in the MAF file
for (i in seq_along(gene.vec)){  
  png(paste0(directory,"/ResultsMaftoolsCohort/lollipopplot", gene.vec[i]))
  tryCatch(expr = {
    lollipopPlot(
      maf = cohort,
      gene = gene.vec[i],
      showMutationRate = TRUE)
  },
  error = function(e){
    print("The protein structure for this gene is not available")},
  warning = function(w){        # Specifying warning message
    message("There was a warning message.")
  },
  finally = {                   # Specifying final message
    message("Lollipopplot is finished.")
  }
  )
  dev.off()
}

#refseqID for 10 most mutated genes in the MAF file
for (i in seq_along(gene.vec)){  
  png(paste0(directory,"/ResultsMaftoolsCohort/plotprotein", gene.vec[i]))
  tryCatch(expr = {
    plotProtein(gene = gene.vec[i])
  },
  error = function(e){
    print("The protein structure for this gene is not available")},
  warning = function(w){        # Specifying warning message
    message("There was a warning message.")
  },
  finally = {                   # Specifying final message
    message("PlotProtein is finished.")
  }
  )
  dev.off()
}

#Rainfall plots
png(paste0(directory,"/ResultsMaftoolsCohort/rainfallplot"))
rainfallPlot(maf = cohort , detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38", width = 10, height = 5)
dev.off()

#Compare mutation load against TCGA cohorts
png(paste0(directory,"/ResultsMaftoolsCohort/TCGAcohorts"))
MAF.mutload = tcgaCompare(maf = cohort, cohortName = 'Cohort', logscale = TRUE)
dev.off()

#This function will return warning messages if it didn't find genes above the threshold (fdrCutOff)

tryCatch(expr = {
   MAF.sig = oncodrive(maf = cohort, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore') #minMut minimum number of mutations required for a gene to be included in analysis
   head(MAF.sig)
   png(paste0(directory,"/ResultsMaftools/oncodrive"))
   plotOncodrive(res = MAF.sig , fdrCutOff = 0.1, useFraction = TRUE, labelSize = 1) #fdr cutoff to call a gene as a driver.
   dev.off()
  },
  error = function(e){
    print("There's no enough mutations")},
  warning = function(w){        # Specifying warning message
    message("There's no enough mutations")
  },
  finally = {                   # Specifying final message
    message("Oncodrive is finished.")
  }
  )

#Adding and summarizing pfam domains
MAF.pfam = pfamDomains(maf = cohort, AACol = 'HGVSp', top = 10)
#Protein summary
MAF.pfam$proteinSummary[,1:7, with = FALSE]
#Domain summary
MAF.pfam$domainSummary[,1:3, with = FALSE]
#Clinical enrichment analysis
#fab.ce = clinicalEnrichment(maf = cohort, clinicalFeature = 'Tumor_Sample_Barcode')
#getClinicalData(cohort)
#fab.ce$groupwise_comparision[p_value < 0.05]
#plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)
###Drug-Gene Interactions
 tryCatch(expr = {
   dgi = drugInteractions(maf = cohort,fontSize = 0.75)
   dnmt3a.dgi = drugInteractions(genes = gene.vec, drugs = TRUE)
   dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]
     },
  error = function(e){
    print("There's no enough mutations")},
  warning = function(w){        # Specifying warning message
    message("There's no enough mutations")
  },
  finally = {                   # Specifying final message
    message("Drug-Gene Interactions is finished.")
  }
  )

#Oncogenic Signaling Pathways
 tryCatch(expr = {
   OncogenicPathways(maf = cohort, pathways = NULL)
   png(paste0(directory,"/ResultsMaftools/oncogenicpathways"))
   PlotOncogenicPathways(maf = cohort)
     },
  error = function(e){
    print("There's no enough mutations")},
  warning = function(w){        # Specifying warning message
    message("There's no enough mutations")
  },
  finally = {                   # Specifying final message
    message("OncogenicPathways is finished.")
  }
  )
  
#Tumor heterogeneity and MATH scores
#Heterogeneity in sample 
library("mclust")

 tryCatch(expr = {
  heterogeneity = inferHeterogeneity(maf = cohort, top = 5, vafCol = 'gnomAD_AF')
  print(heterogeneity$clusterMeans)
  png(paste0(directory,"/ResultsMaftoolsCohort/inferHeterogeneity"))
  plotClusters(clusters = heterogeneity)
  dev.off()
     },
  error = function(e){
    print("There's no enough mutations")},
  warning = function(w){        # Specifying warning message
    message("There's no enough mutations")
  },
  finally = {                   # Specifying final message
    message("Tumor heterogeneity is finished.")
  }
  )
#################GISTIC
all.lesions <- read.table(paste0(directory,"/Gistic2/all_lesions.conf_99.txt"), header = TRUE, sep = "\t")
amp.genes <- read.table(paste0(directory,"/Gistic2/amp_genes.conf_99.txt"), header = TRUE, sep = "\t")
del.genes <- read.table(paste0(directory,"/Gistic2/del_genes.conf_99.txt"), header = TRUE, sep = "\t")
scores.gis <- read.table(paste0(directory,"/Gistic2/scores.gistic"), header = TRUE, sep = "\t")
laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
laml.gistic
genomeplot
png(paste0(directory,"ResultsMaftools/gisticgenomeplot"))
gisticChromPlot(gistic = laml.gistic, markBands = "all")
dev.off()
bubbleplot
png(paste0(directory,"ResultsMaftools/gisticbubbleplot"))
gisticBubblePlot(gistic = laml.gistic)
dev.off()
png(paste0(directory,"ResultsMaftools/gisticoncoplot"))
gisticOncoPlot(gistic = laml.gistic, clinicalData = getClinicalData(x = laml), clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, top = 10)
dev.off()
print("Maftoolscohort has finished")