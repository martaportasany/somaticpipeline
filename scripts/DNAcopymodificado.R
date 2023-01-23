suppressMessages(library(DNAcopy))
suppressMessages(library(argparse))
parser <- ArgumentParser(description= 'This program perfors segmentation')

parser$add_argument('--input', '-i', help= 'input sample')
parser$add_argument('--output', '-o', help= 'input sample')
xargs <- parser$parse_args()
input<- xargs$input
output <- xargs$output
#setwd("/media/rafael/WORK/WES_pipeline/Somatic/results/")
#Manual
#cn <- read.delim("/media/rafael/WORK/WES_pipeline/Somatic/results/641487/output.copynumber", header=T)
#For try in R studio
#input <- paste0("/media/rafael/WORK/WES_pipeline/Somatic/results/641487/output.copynumber")
cn <- read.delim(input, header=T)
CNA.object <-CNA( genomdat = cn[,7], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
seg.pvalue <- segments.p(segs, ngrid=100, tol=1e-6, alpha=0.05, search.range=100, nperm=1000)
Sample <-sub("/media/rafael/WORK/WES_pipeline/Somatic/results/", "", input)
Sample<- sub("/output.copynumber", "", Sample)
seg.pvalue$ID <- Sample
write.table (seg.pvalue, file=output, row.names=F, col.names=F, quote=F, sep="\t")

