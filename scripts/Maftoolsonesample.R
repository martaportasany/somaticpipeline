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
directory <- args$directory
######
getwd()
setwd(directory)
getwd()
#path to MAF file
setwd(paste0(directory,"/Merge/maf/"))
getwd()
MAF = read.maf(maf = "all_VEP_custom.maf")
setwd(directory)
dir.create("ResultsMaftools")
setwd(paste0(directory,"/ResultsMaftools/"))
#Create list of 10 most mutated genes
gene <- getGeneSummary(MAF)
head(gene)
geneorder <- gene[order(-gene$total)]
list <- head(geneorder,10)
genelist <- list[,1]
genelist
gene.vec <- as.vector(genelist$Hugo_Symbol)
is.vector(gene.vec)
length(gene.vec)
#We don't have the TCGA barcode so we use this function to extract it
for (i in seq_along(gene.vec)) {
  print(genesToBarcodes(maf = MAF, genes = gene.vec[i]))
}
getSampleSummary(MAF)
##Shows gene summary.
getGeneSummary(MAF)
#shows clinical data associated with samples
getClinicalData(MAF)
MAF
#Shows all fields in MAF
#getFields(x = MAF) #de momento no funciona
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = MAF, basename = 'maf')
########VISUALIZATION##########
png(paste0(directory,"/ResultsMaftools/plotmafSummary"))
plotmafSummary(maf = MAF, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
#Esto ya requiere más de una muestra
#Oncoplot
#oncoplot for top ten mutated genes.
#oncoplot(maf = MAF, top = 10) #ESTO TIENE QUE SER PARA VARIAS MUESTRAS
MAF
#Transition and Transversions
MAF.titv = titv(maf = MAF, plot = TRUE , useSyn = TRUE)
#plot titv summary
png(paste0(directory,"/ResultsMaftools/plottitv"))
plotTiTv(res = MAF.titv)
dev.off()
#Lollipop plots for amino acid changes for 10 most mutated genes in the MAF file
for (i in seq_along(gene.vec)){  
  png(paste0(directory,"/ResultsMaftools/lollipopplot", gene.vec[i]))
  tryCatch(expr = {
    lollipopPlot(
      maf = MAF,
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
  png(paste0(directory,"/ResultsMaftools/plotprotein", gene.vec[i]))
  tryCatch(expr = {
    plotProtein(gene = gene.vec[i])
  },
  error = function(e){
    print("The protein structure for this gene is not available")},
  warning = function(w){        # Specifying warning message
    message("There was a warning message.")
  },
  finally = {                   # Specifying final message
    message("Plotprotein is finished.")
  }
  )
  dev.off()
}
#Rainfall plots for the most mutated sample
png(paste0(directory,"/ResultsMaftools/rainfallplot"))
rainfallPlot(maf = MAF, detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")
dev.off()
#Compare mutation load against TCGA cohorts
png(paste0(directory,"/ResultsMaftools/TCGAcohorts"))
MAF.mutload = tcgaCompare(maf = MAF, cohortName = 'MAF', logscale = TRUE)
dev.off()

#Plotting VAF
png(paste0(directory,"/ResultsMaftools/plotVaf"))
plotVaf(maf = MAF, vafCol = 'gnomAD_AF', keepGeneOrder = TRUE, showN = TRUE) #the upper index shows the number of mutated genes
dev.off()
#Analysis
#Somatic Interactions laml is needed
#somaticInteractions(maf = laml, top = 25, pvalue = c(0.05,0.1))
#laml.pfam = pfamDomains(maf = MAF, top = 10) #protein change information needed
#This function will return warning messages if it didn't find genes above the threshold (fdrCutOff)
 tryCatch(expr = {
   MAF.sig = oncodrive(maf = MAF, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
   head(MAF.sig)
   png(paste0(directory,"/ResultsMaftools/oncodrive"))
   plotOncodrive(res = MAF.sig , fdrCutOff = 0.1, useFraction = TRUE, labelSize = 1) #fdr cutoff to call a gene as a driver.
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
png(paste0(directory,"/ResultsMaftools/pfamdomains"))
MAF.pfam = pfamDomains(maf = MAF, AACol = 'HGVSp_Short', top = 10)
dev.off()
#Protein summary
MAF.pfam$proteinSummary[,1:7, with = FALSE]
#Domain summary
MAF.pfam$domainSummary[,1:3, with = FALSE]
#Clinical enrichment analysis
#fab.ce = clinicalEnrichment(maf = MAF, clinicalFeature = 'Tumor_Sample_Barcode') 
#getClinicalData(MAF)
#fab.ce$groupwise_comparision[p_value < 0.05]
#plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)
####Drug-Gene Interactions
dgi = drugInteractions(maf = MAF, fontSize = 0.75)
dnmt3a.dgi = drugInteractions(genes = gene.vec, drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]


#Oncogenic Signaling Pathways
 tryCatch(expr = {
   sink("Oncogenicsignalingpathways")
   OncogenicPathways(maf = MAF)
   sink()
   png(paste0(directory,"/ResultsMaftools/oncogenicpathways"))
   PlotOncogenicPathways(maf = MAF)
   dev.off()
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
heterogeneity = inferHeterogeneity(maf = MAF, top = 5, vafCol = 'gnomAD_AF')
print(heterogeneity$clusterMeans)
png(paste0(directory,"/ResultsMaftools/inferHeterogeneity"))
plotClusters(clusters = heterogeneity)
dev.off()
#Ignoring variants in copy number altered regions
####
##Mutational Signatures
#library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
#laml.tnm = trinucleotideMatrix(maf = MAF, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
#plotApobecDiff(tnm = laml.tnm, maf = MAF, pVal = 0.2)
#library('NMF')
#laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6)
#plotCophenetic(res = laml.sign)
#laml.sig = extractSignatures(mat = laml.tnm, n = 3)
#Compate against original 30 signatures 
#laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")
#Compate against updated version3 60 signatures 
#laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")
#library('pheatmap')
#pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
#maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "SBS")
print("Maftoolsonesample has finished")