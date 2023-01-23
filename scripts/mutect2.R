#!/usr/bin/env Rscript
suppressMessages(library(argparse))
setwd("/media/rafael/WORK/WES_pipeline/Somatic/")
parser <- ArgumentParser(description= 'This progrom run Mutec2')

parser$add_argument('--input_normal', '-in', help= 'normal input file')
parser$add_argument('--input_tumor', '-it', help= 'tumor input file')
parser$add_argument('--name_normal', '-nn', help= 'normal name')
parser$add_argument('--name_tumor', '-nt', help= 'tumor name')
parser$add_argument('--pon', '-p', help= 'panel of normals')
parser$add_argument('--intervals', '-i', help= 'intervals list')
parser$add_argument('--ref', '-r', help= 'Reference Genome')
parser$add_argument('--output', '-o', help= 'output file')
parser$add_argument('--fir2', '-f', help = 'fir file')
xargs<- parser$parse_args()


intervals <-xargs$intervals
output <- xargs$output
pon <- xargs$pon
input_normal <- xargs$input_normal
genome <-xargs$ref
name_normal <- xargs$name_normal
input_tumor <-xargs$input_tumor
name_tumor <-xargs$name_tumor
fir <- xargs$fir2

log <- xargs$log
####################################################################
## Prueba manual
# intervals <-"Regions_GATK.interval_list"
# output <- "results/574409/Mutec2/unfiltered.vcf.gz"
# pon <- "pon/pon.vcf.gz"
# input_normal <- "574409/recal/574409-H.bam"
# genome <- "genome.fasta"
# name_normal <- "574409-H"
# input_tumor <- "574409/recal/574409-D.bam"
# name_tumor <- "574409-D"
#################### Constantes que se deben cambiar en otros equipos    
GATK.path<-"/home/rafael/Aplicaciones/GATK/gatk-4.2.5.0/gatk Mutect2"
#bamsdir<- "/media/rafael/Elements/BAMs_Jorge/"
#refGenomeDir<- "/media/rafael/DATA/Genome_GRCh38_100/resources/"  
outdir <- substr(output,1,nchar(output)-17)
dir.create(outdir)
#####################################################################

#intervals<- paste0(refGenomeDir,intervals)
#genome <-paste0(refGenomeDir,genome)

comand <- paste0(GATK.path, " -L ", intervals, " -R ",genome, " -I ", input_normal, " -normal ", name_normal, 
                 " -I ", input_tumor, " -tumor ", name_tumor, " -pon ", pon, " -f1r2-tar-gz ",fir,  " --output ", output )

system(comand)

