import pandas as pd
from snakemake.utils import validate, min_version

import glob,os
import numpy as np
configfile: "config/config.yaml"

callers=["Mutect2", "strelka", "VarScan"]
callers_2 =["strelka", "VarScan"]
vartypes=["snvs", "indels"]
sample_table    = pd.read_table(config['samples'], sep='\t', lineterminator='\n', dtype=str).set_index(["sample", "family"], drop=False)
sample_table.index = sample_table.index.set_levels([i.astype(str) for i in sample_table.index.levels])  # enforce str in index
sample_types = ['normal', 'map']
intervals =  ["{:04d}".format(n) for n in range(40)]

wildcard_constraints:
    sample = "|".join(sample_table['sample']),
    group = "|".join(sample_table['family']),
    vartype = "|".join(vartypes),
    sample_type="|".join(sample_types)

filter_flags: "--min-reads-per-strand 1"
############# Functions ##################

# Genomic Regions
regions_bed = config['genomic_regions']
regions_gatk = os.path.basename(regions_bed).replace('.bed', '.interval_list')
regions_gatk = os.path.join('interval-files', regions_gatk)

num_workers= 40

def get_intervals():
    ints = []
    for i in range(num_workers):
        num_zeros = 4 - len(str(i))
        interval = '0' * num_zeros + str(i)
        ints.append(interval)
    return ints

def get_interval_files():
    ints = get_intervals()
    files = [i + '-scattered.interval_list' for i in ints]
    files = [os.path.join("interval-files", f) for f in files]
    return files

interval_files = get_interval_files()



def get_bams(wildcards):
    import re

    list_bams = config['BAMdir'] + "results/{group}/recal/{sample}.bam",
    return  [a for a in list_bams if re.match("|".join(['-D.bam', '-H.bam']), a)]


def get_mutect2_input(wildcards):
    files = {}
    files['map'] = f"bams/{wildcards.patient}.map.bam"
    if use_pon:
        files['pon'] = pon_vcf
    if not tumor_only:
        files['normal'] = f"bams/{wildcards.patient}.normal.bam"
    return files

#########################################################################################################################################
############################## DICT PARA MERGED #########################
familias = {k: list(v) for k, v in sample_table.groupby('Group')['sample']}
familias_Comb = dict((el,callers) for el in sample_table['family'].unique())
familias_Combine = {
    key: ["results/" + key + "/" + el + "/filtered/all.norm.filtered.vcf.gz" for el in value]
    for key, value in familias_Comb.items()
} 
familias_w = {
    key: [config['bamsdir'] + key + "/recal/" + e + ".bam" for e in value]
    for key, value in familias.items()
}

familias_w_bai = {
    key: [config['bamsdir']  + key + "/recal/" + e + ".bam.bai" for e in value]
    for key, value in familias.items()
}
familias_Comb = dict((el,callers) for el in sample_table['family'].unique())
familias_Combine = {
    key: ["results/" + key + "/" + el + "/filtered/all.norm.filtered_2.vcf.gz" for el in value]
    for key, value in familias_Comb.items()
} 

familias_Combine_idx = {
    key: ["results/" + key + "/" + el + "/filtered/all.norm.filtered_2.vcf.gz.csi" for el in value]
    for key, value in familias_Comb.items()
} 

familias_Comb_2 = dict((el,callers_2) for el in sample_table['family'].unique())
familias_Combine_2 = {
    key: ["results/" + key + "/" + el + "/filtered/all.norm.filtered_2.vcf.gz" for el in value]
    for key, value in familias_Comb_2.items()
} 

familias_Combine_idx_2 = {
    key: ["results/" + key + "/" + el + "/filtered/all.norm.filtered_2.vcf.gz.csi" for el in value]
    for key, value in familias_Comb_2.items()
} 

#######################################################################################################################################################################################################

rule all:
    input:
        config['refGenomeDir'] + "Regions_GATK.interval_list",
        config['refGenomeDir'] + "somatic-hg38_af-only-gnomad.hg38.vcf.gz",
        interval_files,
        expand("pon/{group}.pon.vcf.gz",  group=sample_table['family']),
        "pon/pon_db",
        "pon/pon.vcf.gz",
        expand("BAMs/{group}/normal/{group}.normal.bam" , group=sample_table['family']),
        expand("BAMs/{group}/map/{group}.map.bam" , group=sample_table['family']),
        #expand("results/{group}/Mutect2/unfiltered.vcf.gz",  group=sample_table['family'].unique()),
        #expand("vcfs/{group}.annotated.vcf.gz",  group=sample_table['family'].unique()), #este probablemente se pueda eliminar
        expand("results/{group}/tumorgetpilesummaries.table", group=sample_table['family'].unique()),
        expand("results/{group}/normalgetpilesummaries.table", group=sample_table['family'].unique()),
        expand("results/{group}/tumorcalculatecontamination.table", group=sample_table['family'].unique()),
        expand("results/{group}/normalcalculatecontamination.table", group=sample_table['family'].unique()),
        expand("results/{group}/matchedcontamination.table", group=sample_table['family'].unique()),
        #artifacts prior
        #expand("results/{group}/Mutect2/filtered/all.filtered.vcf.gz", group=sample_table['family'].unique()),
      # expand("results/{group}/Mutect2/{group}.totalmutationrecords.txt",  group=sample_table['family'].unique()),
      # expand("results/{group}/Mutect2/{group}.missenserecords.txt",  group=sample_table['family'].unique()),
       # expand("results/{group}/Manta/filtered_diploidSV.vcf.gz",  group=sample_table['family'].unique()),
       # expand("results/{group}/Strelka/strelka/results/variants/somatic.snvs.vcf.gz" , group=sample_table['family'].unique()),
       # expand("results/{group}/Strelka/strelka/results/variants/somatic.indels.vcf.gz", group=sample_table['family'].unique()),
       # expand("results/{group}/strelka/filtered/all.filtered.vcf.gz", group=sample_table['family'].unique()),
     #	expand("allelecountSNVs/{group}/.somatic.snv.vcf.gz" , group=sample_table['family'].unique()),
     # 	expand("allelecountindels/{group}/.somatic.indels.vcf.gz" , group=sample_table['family'].unique()),
        expand("results/{group}/VarScan/filtered/all.filtered.vcf.gz", group=sample_table['family'].unique()),
       #expand("results/{group}/{caller}/filtered/all.norm.filtered.vcf.gz.csi",caller=callers, group=sample_table['family'].unique()),
       #expand("results/{group}/Filter_tables/", group=sample_table['family'].unique()),
       #expand("results/{group}/matchedcontamination.table"	, group=sample_table['family'].unique()),
       #expand("results/{group}/Merge/calls.tsv.gz", group=sample_table['family'].unique()),
       #expand("results/{group}/Merge/all_VEP_custom.tsv.gz", group=sample_table['family'].unique()),
        expand("results/{group}/Merge/all_VEP_custom.vcf.gz", group=sample_table['family'].unique()),
       #expand("results/{group}/Combination_2/0000.vcf.gz", group=sample_table['family'].unique()),
       #expand("results/{group}/Cohort/Maftools/MaftoolsRall.txt", group=sample_table['family'].unique()),
       #expand("results/{group}/Merge/csq_phased_all.vcf.gz" , group=sample_table['family'].unique()),
       #expand("results/{group}/Merge/all_VEP_custom.vcf", group=sample_table['family'].unique()),
       #expand("results/{group}/Merge/all_VEP_custom.tsv", group=sample_table['family'].unique()),
       # expand("results/{group}/Merge/maf/all_VEP_custom.maf", group=sample_table['family'].unique()),
       #expand("results/{group}/Maftools/maftoolsR.txt, group=sample_table['family'].unique()),
       #"results/Cohort/Maftools",
       #expand("results/Cohort/Maftools/MaftoolsRall.txt", group=sample_table['family'].unique()),
       # expand("/media/rafael/WORK/WES_pipeline/Somatic/results/{group}/output.copynumber", group=sample_table['family'].unique()),
       # expand("results/{group}/Gistic2/varScan.copynumber.caller", group=sample_table['family'].unique()),
       #expand("/media/rafael/WORK/WES_pipeline/Somatic/results/{group}/cbs_P1_pvalue", group=sample_table['family'].unique()),
       #expand("results/{group}/Msisensor/report", group=sample_table['family'].unique()),
       # expand("results/{group}/ResultsMaftools/plotmafSummary", group=sample_table['family'].unique()),
      #  "results/ResultsMaftoolsCohort/plotmafSummary",
        #"results/Gistic2/amp_qplot.pdf"
        #"results/RTG/SDF"
        		

        
        
 
       
       

		
			
     

############################################################## Pasos previos #########################################################################################################

rule make_gatk_regions:
	input:
		bed= config['refGenomeDir'] + "Regions.bed",
		d= config['refGenomeDir'] + "genome.dict"
	output:
		intlist=config['refGenomeDir'] + "Regions_GATK.interval_list"
	conda:
		"envs/gatk.yaml"
	log:
		"logs/gatk/make_regions.log"
	shell:
		"""
		gatk BedToIntervalList \
			-I {input.bed} \
			-SD {input.d} \
			-O {output} &> {log}
		"""
		 
rule split_intervals:
    input:
        ref=config['refGenomeDir'] + "genome.fasta",
        intervals=config['refGenomeDir'] + "Regions_GATK.interval_list",
    output:
        interval_files
    log:
    	"logs/split_intervals.log"
    params:
        N=num_workers,
        d="interval-files"
    
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk SplitIntervals -R {input.ref} -L {input.intervals} \
            --scatter-count {params.N} -O {params.d} \
            --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION
            &> {log}
        """
############################################################### PON #############################################################################################################

rule mutect2_pon:
    input:
        bam=config['bamsdir'] + "{group}/recal/{group}-H.bam",
        intervals= config['refGenomeDir'] + "Regions_GATK.interval_list",
        ref=config['refGenomeDir'] + "genome.fasta",
    output:
        vcf="pon/{group}.pon.vcf.gz",
        idx="pon/{group}.pon.vcf.gz.tbi",
        stats="pon/{group}.pon.vcf.gz.stats"
    conda:
        "envs/gatk.yaml"
    threads:
        4
    log:
    	"logs/{group}/mutect2_pon.log"
    shell:
        """
        gatk Mutect2 \
            -I {input.bam} \
            -R {input.ref} \
            -O {output.vcf} \
            -L {input.intervals} \
            -ip 100 \
            --max-mnp-distance 0 &> {log} #max-np-distance tiene que ser 0 pq sino causa interferencias con el GenomicsDBImport
        """


rule gather_variants:
    input:
        vcfs=expand("pon/{group}.pon.vcf.gz", group=sample_table['family'].unique()),
        intervals=config['refGenomeDir'] + "Regions_GATK.interval_list",
        ref=config['refGenomeDir'] + "genome.fasta",
    output:
        directory("pon/pon_db")
    params:
        vcfs=lambda wildcards, input:  [f" -V {v}" for v in input["vcfs"]]
    conda:
        "envs/gatk.yaml"
    threads:
        32
    #log:
   		#"logs/{group}/gather_variants.log"
    shell:
        """
          gatk GenomicsDBImport {params.vcfs} --genomicsdb-workspace-path {output} -L {input.intervals} -ip 100 --merge-input-intervals true 
        """

rule create_pon:
    input:
        var="pon/pon_db",
        ref=config['refGenomeDir'] + "genome.fasta",
    output:
        vcf="pon/pon.vcf.gz",
        idx="pon/pon.vcf.gz.tbi"
    conda:
        "envs/gatk.yaml"
    log:
    	"logs/create_pon.log"
    shell:
        """
        gatk CreateSomaticPanelOfNormals -R {input.ref} -V gendb://{input.var} -O {output.vcf} &> {log}
        """
####################################################################################################

rule samtools_index_normal:
	input:
		normal= config['bamsdir'] + "{group}/recal/{group}-H.bam"
	output:
		normal= "BAMs/{group}/normal/{group}.normal.bam"
	
	shell:
		"""
		samtools index {input.normal} > {output.normal}
		"""
rule samtools_index_tumor:
	input:
		map= config['bamsdir'] + "{group}/recal/{group}-D.bam"
	output:
		map= "BAMs/{group}/map/{group}.map.bam"
	shell:
		"""
		samtools index {input.map} > {output.map}
		
		"""        


rule mutect2_Rscript:
	input:
		normal = config['bamsdir'] + "{group}/recal/{group}-H.bam",
		tumor = config['bamsdir'] + "{group}/recal/{group}-D.bam",
		pon = "pon/pon.vcf.gz",
		intervals = config['refGenomeDir'] + "Regions_GATK.interval_list",
		genome = config['refGenomeDir'] + "genome.fasta",
	params:
		normal = "{group}-H",
		tumor = "{group}-D"
	output:
		vcf = "results/{group}/Mutect2/unfiltered.vcf.gz",
		fir = "results/{group}/Mutect2/f1r2.tar.gz"
	threads:
		32
	log:
		"logs/{group}/mutect2/mutect2.log"
	shell:
		"scripts/mutect2.R --input_normal {input.normal} --name_normal {params.normal} --input_tumor {input.tumor} --name_tumor {params.tumor} --pon {input.pon} --intervals {input.intervals} --ref {input.genome} --fir2 {output.fir} --output {output.vcf} 2> {log}"

  	
    	

rule GetPileupSummaries_tumorsamples:
	input:
		tumor= config['bamsdir'] + "{group}/recal/{group}-D.bam",
		vcf = "config/af-only-gnomad.hg38.common_biallelic.chr1-22XY.vcf",
		intervals = config['refGenomeDir'] + "Regions_GATK.interval_list",
	output:
		"results/{group}/tumorgetpilesummaries.table"
	log:
		"logs/{group}/mutect2/getpileupsummariestumor.log"
	conda:
		"envs/gatk.yaml"
	shell:
		"gatk GetPileupSummaries -I {input.tumor} -V {input.vcf} -L {input.intervals} -O {output} &> {log}"

rule GetPileupSummaries_normalsamples:
	input:
		tumor= config['bamsdir'] + "{group}/recal/{group}-H.bam",
		vcf = "config/af-only-gnomad.hg38.common_biallelic.chr1-22XY.vcf",
		intervals = config['refGenomeDir'] + "Regions_GATK.interval_list",
	output:
		"results/{group}/normalgetpilesummaries.table"
	log:
		"logs/{group}/mutect2/getpileupsummariesnormal.log"
	conda:
		"envs/gatk.yaml"
	shell:
		"gatk GetPileupSummaries -I {input.tumor} -V {input.vcf} -L {input.intervals} -O {output} &> {log}"

rule CalculateContamination_tumor:
	input:
		"results/{group}/tumorgetpilesummaries.table"
	output:
		"results/{group}/tumorcalculatecontamination.table"
	conda:
		"envs/gatk.yaml"
	log:
		"logs/{group}/mutect2/calculatecontamination_tumor.log"
	shell:
		"gatk CalculateContamination -I {input} -O {output} &> {log}"

rule CalculateContamination_normal:
	input:
		"results/{group}/normalgetpilesummaries.table"
	output:
		"results/{group}/normalcalculatecontamination.table"
	conda:
		"envs/gatk.yaml"
	log:
		"logs/{group}/mutect2/calculatecontamination_normal.log"
	shell:
		"gatk CalculateContamination -I {input} -O {output} &> {log}"

rule CalculateContamination_matched:
	input:
		tumor= "results/{group}/tumorgetpilesummaries.table",
		normal= "results/{group}/normalgetpilesummaries.table"
	output:
		"results/{group}/matchedcontamination.table"
	log:
		"logs/{group}/mutect2/calculatecontamination_matched.log"
	conda:
		"envs/gatk.yaml"
	shell:
		"gatk CalculateContamination -I {input.tumor} -matched {input.normal} -O {output} &> {log}"

#rule CollectF1R2Counts:
#	input:
#		ref=config['refGenomeDir'] + "genome.fasta",
#		tumor= config['bamsdir'] + "{group}/recal/{group}-D.bam"	
#	output:
#		"f1r2.tar.gz"
#	log:
#		"logs/mutect2/CollectF1R2Counts.log"
#	##aquí necesito un group yo creo
#	shell:
#		"gatk CollectF1R2Counts -R {input.ref} -I {input.tumor} -O f1r2.tar.gz &> {log}"

rule learnorientationmodel:
	input:
		f1r2="results/{group}/Mutect2/f1r2.tar.gz"
	output:
		artifacts="results/{group}/Mutect2/artifacts_prior.tar.gz"
	log:
		"logs/{group}/mutect2/learnorientationmodel.log"
	conda:
		"envs/gatk.yaml"
	shell:
		"gatk LearnReadOrientationModel -I {input.f1r2} -O {output.artifacts} &> {log}"
		
rule filtermutectcalls:
	input:
		vcf="results/{group}/Mutect2/unfiltered.vcf.gz",
		artifacts="results/{group}/Mutect2/artifacts_prior.tar.gz",
		table = "results/{group}/matchedcontamination.table",
		ref = config['refGenomeDir'] + "genome.fasta"
	output:
		"results/{group}/Mutect2/filtered/all.filtered.vcf.gz"
	log:
		"logs/{group}/mutect2/filtermutectcalls.log"
	conda:
		"envs/gatk.yaml"
	shell:
		"gatk FilterMutectCalls -V {input.vcf} -R {input.ref} --contamination-table {input.table} --orientation-bias-artifact-priors {input.artifacts} -O {output} &> {log} "

#rule totalmutationrecords:
#	input:
#		vcf="results/{group}/Mutect2/filtered/all.filtered.vcf.gz"
#	output:
#		"results/{group}/Mutect2/{group}.totalmutationrecords.txt"
#	log:
#		"logs/{group}/mutect2/totalmutationrecords.log"
#	shell:
#		"zgrep -v "^#" {input.vcf} | wc > {output} &> {log}"
	
#rule numberofmissenserecords:
#	input:
#		vcf="results/{group}/Mutect2/filtered/all.filtered.vcf.gz"
#	output:
#		"results/{group}/Mutect2/{group}.missenserecords.txt"
#	log:
#		"logs/{group}/mutect2/numberofmissenserecods.log"
#	shell:
#		"zgrep -v "^#" {input.vcf} | grep "|MISSENSE|" | grep PASS | wc > {output} &> {log}"
		

################################################################################ STRELKA ####################################################################################################
##### MANTA #####
rule manta:  
    input:
    	nbam=config['bamsdir'] + "{group}/recal/{group}-H.bam",
	tbam=config['bamsdir'] + "{group}/recal/{group}-D.bam",
	ref=config['refGenomeDir'] + "genome.fasta"
	
    params:
        manta = "results/{group}/Manta/manta/mi",
    
    output:
        "results/{group}/Manta/manta/results/variants/diploidSV.vcf.gz",
        "results/{group}/Manta/manta/results/variants/candidateSmallIndels.vcf.gz",
        
    
    conda:
        "envs/manta.yaml"
    log:
    	"logs/{group}/Manta/manta.log"
    
    threads:
        32
    resources:
        mem_mb = 4096,
          
    shell:
        """
        set -xe
        OUTDIR="$(dirname "{params.manta}")"
              
            configManta.py \
            	--normalBam {input.nbam}\
            	--tumorBam {input.tbam}\
            	--referenceFasta {input.ref}\
            	--runDir "${{OUTDIR}}"\
            	--exome 
            
            /media/rafael/WORK/WES_pipeline/Somatic/results/{wildcards.group}/Manta/manta/runWorkflow.py \
                --quiet \
                -m local \
                -j {threads}  
                     
        """
        
rule filter_manta:
    input:
        "results/{group}/Manta/manta/results/variants/diploidSV.vcf.gz"
    output:
        "results/{group}/Manta/filtered_diploidSV.vcf.gz"
    log:
    	"logs/{group}/Manta/filtermanta.log"
    conda:
        "envs/bcftools.yaml"
    shell:
        """
          bcftools view -i'FILTER="PASS"' {input} -Oz -o {output} &> {log}
        """ 
        
############## STRELKA2 #######################

rule strelka: 
	input:
		nbam=config['bamsdir'] + "{group}/recal/{group}-H.bam",
		tbam=config['bamsdir'] + "{group}/recal/{group}-D.bam",
		ref=config['refGenomeDir'] + "genome.fasta",
		manta="results/{group}/Manta/manta/results/variants/candidateSmallIndels.vcf.gz"
	output:
	    "results/{group}/Strelka/strelka/results/variants/somatic.snvs.vcf.gz",
		"results/{group}/Strelka/strelka/results/variants/somatic.indels.vcf.gz"
		
	
	params:
		strelka = "results/{group}/Strelka/strelka/mi",
	log:
		"logs/{group}/strelka/strelka_general.log"
	threads:
		32			
	shell:
		"""
        set -xe
        OUTDIR="$(dirname "{params.strelka}")"
              
            /home/rafael/Aplicaciones/Strelka2/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
            	--normalBam {input.nbam}\
            	--tumorBam {input.tbam}\
            	--referenceFasta {input.ref}\
            	--indelCandidates {input.manta}\
            	--runDir "${{OUTDIR}}"\
            	--exome 
            
            /media/rafael/WORK/WES_pipeline/Somatic/results/{wildcards.group}/Strelka/strelka/runWorkflow.py \
                --quiet \
                -m local \
                -j {threads}  
                &> {log}     
        """



##################Strelka2

rule merge_strelka:
	input:
		snv="results/{group}/Strelka/strelka/results/variants/somatic.snvs.vcf.gz",
		indel="results/{group}/Strelka/strelka/results/variants/somatic.indels.vcf.gz"
		
	output:
		"results/{group}/strelka/filtered/all.filtered.vcf.gz"
	log:
		"logs/{group}/strelka/merge_strelka.log"
	shell:
		"bcftools merge {input.snv} {input.indel} --force-samples > {output}"


rule allelecount_strelka_SNVs:
	input:
		"results/{group}/Strelka/strelka/results/variants/somatic.snvs.vcf.gz"
	output:
		"results/{group}/Strelka/strelka/results/allelecountsSNVs/{group}.somatic.snv.vcf.gz"
	log:
		"logs/{group}/strelka/allelecount_strelka_SNVs.log"
	script:
		"scripts/allele_counts_from_strelka.py"
		
rule allelecount_strelka_indels:
	input:
		"results/{group}/Strelka/strelka/results/variants/somatic.indels.vcf.gz"
	output:
		"results/{group}/Strelka/strelka/results/allelecountsindels/{group}.somatic.indels.vcf.gz"
	log:
		"logs/{group}/strelka/allelecount_strelka_indels.log"
	script:
		"scripts/allele_counts_from_strelka.py"

###################################################################VarScan#################################################################

###########################SamTools#####################

rule normal_samtools_mpileup:
	input:
		nbam=config['bamsdir'] + "{group}/recal/{group}-H.bam",
		ref=config['refGenomeDir'] + "genome.fasta",
		interval = config['refGenomeDir'] + "Regions_GATK.interval_list"
	output:
		nmpileup="results/{group}.normal.mpileup"
	log:
		"logs/{group}/Samtools/normal_samtools_mpileup.log"
	conda:
		"envs/samtools.yaml"	
	shell:
		"samtools mpileup {input.nbam} -f {input.ref} -o {output.nmpileup} &> {log}"

rule tumor_samtools_mpileup:
	input:
		tbam=config['bamsdir'] + "{group}/recal/{group}-D.bam",
		ref=config['refGenomeDir'] + "genome.fasta",
		interval = config['refGenomeDir'] + "Regions_GATK.interval_list"
	output:
		tmpileup="results/{group}.tumor.mpileup"
	log:
		"logs/{group}/Samtools/tumor_samtools_mpileup.log"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools mpileup {input.tbam} -f {input.ref} -o {output.tmpileup} &> {log}"

#########################VarScan#########################
rule VarScan:
	input:
		nmpileup="results/{group}.normal.mpileup",
		tmpileup="results/{group}.tumor.mpileup",
		ref=config['refGenomeDir'] + "genome.fasta"	
	output:
		snp="results/{group}/VarScan/all.mpileup.output.snp.vcf",
		indel="results/{group}/VarScan/all.mpileup.output.indel.vcf"
	conda:
		"envs/VarScan.yaml"
	log:
		"logs/{group}/VarScan/varscan.log"
	shell:
		"varscan somatic {input.nmpileup} {input.tmpileup} --output-snp {output.snp} --output-indel {output.indel}  &> {log}"

rule somaticfilter_VarScan:
	input:
		snp="results/{group}/VarScan/all.mpileup.output.snp.vcf",
		indel="results/{group}/VarScan/all.mpileup.output.indel.vcf"
	output:
		snp="results/{group}/VarScan/all.mpileup.filtered.snp.vcf",
		indel="results/{group}/VarScan/all.mpileup.filtered.indel.vcf"
	log:
		"logs/{group}/VarScan/somaticfilter.log"
	conda:
		"envs/VarScan.yaml"	
	shell:
		"varscan somaticFilter {input.snp} --output-file {output.snp} | varscan somaticFilter {input.indel} --output-file {output.indel} &>log "
	
rule mpileupvcf_snp:
	input:
		snp="results/{group}/VarScan/all.mpileup.filtered.snp.vcf",
		ref=config['refGenomeDir'] + "genome.fasta"	
	output:
		snp="results/{group}/VarScan/all.vcf.filtered.snp.vcf.gz",
		
	shell:
		"bgzip {input.snp} > {output.snp}"	
		
rule mpileupvcf_indels:
	input:
		indel="results/{group}/VarScan/all.mpileup.filtered.indel.vcf",
		ref=config['refGenomeDir'] + "genome.fasta"	
	output:
		indel="results/{group}/VarScan/all.vcf.filtered.indel.vcf.gz"
	shell:
		"bgzip {input.indel} > {output.indel}"	

##https://github.com/scchess/Varscan2VCF
rule varscan2vcf_snp:
		input:
			"results/{group}/VarScan/all.mpileup.filtered.snp.vcf",
		output:
			"results/{group}/VarScan/all.mpileup.filtered.snp.fixed.vcf",	
		shell:
			"python /media/rafael/WORK/WES_pipeline/Somatic/scripts/VarScan2VCF.py {input} > {output}"
rule varscan2vcf_indel:
		input:
			"results/{group}/VarScan/all.mpileup.filtered.indel.vcf",
		output:
			"results/{group}/VarScan/all.mpileup.filtered.indel.fixed.vcf",	
		shell:
			"python /media/rafael/WORK/WES_pipeline/Somatic/scripts/VarScan2VCF.py {input} > {output}"

rule bgzip_snp:
	input:
		snp="results/{group}/VarScan/all.mpileup.filtered.snp.fixed.vcf",
		ref=config['refGenomeDir'] + "genome.fasta"	
	output:
		snp="results/{group}/VarScan/all.mpileup.filtered.snp.fixed.vcf.gz",
		
	shell:
		"bgzip -f {input.snp} > {output.snp}; tabix -p vcf {output.snp}"	
		
rule bgzip_indels:
	input:
		indel="results/{group}/VarScan/all.mpileup.filtered.indel.fixed.vcf",
		ref=config['refGenomeDir'] + "genome.fasta"	
	output:
		indel= "results/{group}/VarScan/all.mpileup.filtered.indel.fixed.vcf.gz",
	shell:
		"bgzip -f {input.indel} > {output.indel}; tabix -p vcf {output.indel}"
	
rule merge_VS_snp_indel:
	input:
		snp="results/{group}/VarScan/all.mpileup.filtered.snp.fixed.vcf.gz",
		indel="results/{group}/VarScan/all.mpileup.filtered.indel.fixed.vcf.gz"
	output:
		vcf="results/{group}/VarScan/filtered/all.filtered.vcf.gz"
	log:
		"logs/{group}/VarScan/merge_VS_snp_indel.log"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools merge {input.snp} {input.indel} --force-samples > {output}"	
		
#################################################################RTG-ANALYSIS

rule RTG_baseline:
    input:
    	ref=  config['refGenomeDir'] + "GRCh38_hs38d1.sdf/",
        baseline= ""
    output:
        directory= "results/RTG/SDF"
    conda:
        "envs/rtg.yaml"
    log:
        "logs/662266/RTG/RTG.log"
    shell:
    	"""
        rtg genomesim -o SDF --max-length INT --min-length INT -n INT
        
        """    

rule RTG_vcfeval:
    input:
        vcfmutect="results/662266/Mutect2/filtered/all.filtered.vcf.gz",
        vcfvarscan="results/662266/VarScan/filtered/all.filtered.vcf.gz",
        vcfsstrelka="results/662266/strelka/filtered/all.filtered.vcf.gz",
        ref=  config['refGenomeDir'] + "GRCh38_hs38d1.sdf/",
        baseline= ""
    output:
        directory= "results/662266/RTG"
    conda:
        "envs/rtg.yaml"
    log:
        "logs/662266/RTG/RTG.log"
    shell:
    	"""
       # rtg popsim/samplesim(baseline)
        rtg vcfeval -b {input.baseline} -c {input.vcfmutect} -c {input.vcfvarscan} -c {vcfstrelka} -o {output.directory} -t SDF
        rtg rocplot
        """
	

################################################################## NORMALIZATION AND MERGE ####################################################################################

rule Normalization:
    input:
         vcf= "results/{group}/{caller}/filtered/all.filtered.vcf.gz",
         ref= config['refGenomeDir'] + "genome.fasta",
    output:
         vcf= "results/{group}/{caller}/filtered/all.norm.filtered.vcf",
         #bcf= "results/{group}/{caller}/filtered/all.norm.filtered.bcf",
    log:
        "logs/{group}/{caller}/Normalization/Norm.log"
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools norm {input.vcf} -f {input.ref} -o {output.vcf} &> {log}"#; bcftools index {output.vcf}"
        #"bcftools norm {input.vcf} | bcftools norm -Ov --check-ref w -f {input.ref} | bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' > {output.vcf}  &> {log}"
#ruleorder: index_norm > Merge_Norms > phase_WhatsHap

rule filter_PASS:
	input:
		"results/{group}/{caller}/filtered/all.norm.filtered.vcf",
	output:
		"results/{group}/{caller}/filtered/all.norm.filtered_2.vcf"
	log:
		"logs/{group}/{caller}/filter_PASS.log"
	conda:
		"envs/bcftools.yaml"     
	shell:
		"bcftools view -f 'PASS,.' {input} > {output}"

rule norm_gzip:
    input:
    	"results/{group}/{caller}/filtered/all.norm.filtered_2.vcf",
    output:
    	"results/{group}/{caller}/filtered/all.norm.filtered_2.vcf.gz"
    log:
    	"logs/{group}/{caller}/norm_gzip.log"
    conda:
       "envs/bcftools.yaml"     
    shell:
        "bgzip -c {input} > {output}"

rule index_norm:
    input:
       "results/{group}/{caller}/filtered/all.norm.filtered_2.vcf.gz"
    output:
       "results/{group}/{caller}/filtered/all.norm.filtered_2.vcf.gz.csi"
    log:
    	"logs/{group}/{caller}/index_norm.log"
    conda:
       "envs/bcftools.yaml"     
    shell:
       "bcftools index {input} -o {output} &> {log}"

rule Merge_Norms:
    input:
        vcfs= lambda wildcards: familias_Combine[wildcards.group], 
        idx = lambda wildcards: familias_Combine_idx[wildcards.group],
        #vcfs= "results/{group}/{caller}/filtered/all.norm.filtered.vcf.gz", 
        #idx = "results/{group}/{caller}/filtered/all.norm.filtered.vcf.gz.csi", 
        #ref= config['refGenomeDir'] + "genome.fasta",
        #intervals = config['refGenomeDir'] +  "Regions.bed",
    output:
        "results/{group}/Combination/0000.vcf.gz",
        "results/{group}/Combination/0000.vcf.gz.tbi",
        "results/{group}/Combination/0001.vcf.gz",
        "results/{group}/Combination/0001.vcf.gz.tbi",
        "results/{group}/Combination/0002.vcf.gz",
        "results/{group}/Combination/0002.vcf.gz.tbi",
        "results/{group}/Combination/README.txt",
        "results/{group}/Combination/sites.txt",
                
    log:
        "logs/{group}/Combination/all.log"
    conda:
        "envs/bcftools.yaml"
    params:
        #lambda wildcards, input: [f" -V {v}" for v in input["vcfs"]]
        directory("results/{group}/Combination/")
    shell:
        "bcftools isec -p {params} -Oz --nfiles=3 {input.vcfs} &> {log}"
       


rule Merge_Norms_Strelka_VarScan:
    input:
        vcfs= lambda wildcards: familias_Combine_2[wildcards.group], 
        idx = lambda wildcards: familias_Combine_idx_2[wildcards.group],
        #vcfs= "results/{group}/{caller}/filtered/all.norm.filtered.vcf.gz", 
        #idx = "results/{group}/{caller}/filtered/all.norm.filtered.vcf.gz.csi", 
        #ref= config['refGenomeDir'] + "genome.fasta",
        #intervals = config['refGenomeDir'] +  "Regions.bed",
    output:
        "results/{group}/Combination_2/0000.vcf.gz",
        "results/{group}/Combination_2/0000.vcf.gz.tbi",
        "results/{group}/Combination_2/0001.vcf.gz",
        "results/{group}/Combination_2/0001.vcf.gz.tbi",
        "results/{group}/Combination_2/README.txt",
        "results/{group}/Combination_2/sites.txt",
                
    log:
        "logs/{group}/Combination/all_2.log"
    conda:
        "envs/bcftools.yaml"
    params:
        #lambda wildcards, input: [f" -V {v}" for v in input["vcfs"]]
        directory("results/{group}/Combination_2/")
    shell:
        "bcftools isec -p {params} -Oz --nfiles=2 {input.vcfs} &> {log}"
       

#######################################################################################################################################################################################################
   ######################################################################## Corregir los Multiple Nucleotide Variants #############################################################################
#######################################################################################################################################################################################################
 
rule bcftools_csq:
    input:
        vcf= "results/{group}/Combination/0000.vcf.gz",
        
        ref= config['refGenomeDir'] + "genome.fasta",
        gff= config['refGenomeDir'] + "Homo_sapiens.GRCh38.105.chr.gff3.gz"
    output:
        "results/{group}/Merge/csq_phased_all.vcf.gz",
    log:
        "logs/{group}/bcf_csq/CSQ.log"
    conda:
        "envs/bcftools.yaml"
    threads:
        32
    shell:
        "bcftools csq -f {input.ref} -g {input.gff} {input.vcf} -Oz -o {output} --phase R --ncsq 1 &> {log}"

#############################################################################ANNOTATION################################################################################################
  

rule VEP_custom:
    input:
        vcf = "results/{group}/Merge/csq_phased_all.vcf.gz",
        ref = config['refGenomeDir'] + "genome.fasta",
    output:
        vcf = "results/{group}/Merge/all_VEP_custom.vcf.gz",
    log:
        "logs/{group}/vep/annotate_custom.log",
    conda:
        "envs/VEP.yaml"
    params:
        clinvar = "--custom /media/rafael/DATA/genomes/references/VEP/ClinVar/clinvar_20220403.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN",
        gnomad = "--custom /media/rafael/DATA/genomes/references/VEP/gnomAD/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz,gnomAD,vcf,exact,0,non_cancer_AC,non_cancer_AN,non_cancer_AF,non_cancer_nhomalt",
        IARC = "--custom /media/rafael/DATA/genomes/references/VEP/IARC/IARC_other.bed.gz,IARC_indel,bed,overlap,0",
        COSMIC = "--custom /media/rafael/DATA/genomes/references/VEP/COSMIC/COSMIC_sorted.bed.gz,COSMICHotspot,bed,overlap,0",
        LOVD_MSH2 = "--custom /media/rafael/DATA/genomes/references/VEP/LOVD/LOVD_MSH2_other.bed.gz,LOVD_MSH2_indel,bed,overlap,0",
        LOVD_MSH6 = "--custom /media/rafael/DATA/genomes/references/VEP/LOVD/LOVD_MSH6_other.bed.gz,LOVD_MSH6_indel,bed,overlap,0",
        LOVD_APC = "--custom /media/rafael/DATA/genomes/references/VEP/LOVD/LOVD_APC_other.bed.gz,LOVD_APC_indel,bed,overlap,0",
        LOVD_NF1 = "--custom /media/rafael/DATA/genomes/references/VEP/LOVD/LOVD_NF1_other.bed.gz,LOVD_NF1_indel,bed,overlap,0",
        CADD = "--custom /media/rafael/DATA/genomes/references/VEP/CADD/whole_genome_SNVs_CADD.vcf.gz,CADD,vcf,exact,0,RAW,PHRED",
        Human_106 ="--custom /home/rafael/.vep/homo_sapiens/106_GRCh38,RefSeq,SIFT,dbSNP,COSMIC,ClinVar,gnomAD",       
        dbNFSP4 = "--custom /media/rafael/DATA/genomes/references/VEP/dbNSFP/dbNSFP4.0a.gz,GERP++_RS,phastCons100way_vertebrate,phyloP100way_vertebrate,phastCons30way_mammalian,phastCons100way_vertebrate",
        dbscSNV = "--custom /media/rafael/DATA/genomes/references/VEP/dbscSNV/dbscSNV1.1_GRCh38.txt.gz,GRCh38",
        DisGeNET = "--custom /media/rafael/DATA/references/VEP/DisGeNET/all_variant_disease_pmid_associations_2020.tsv.gz,disease,rsid",
        dbSNP = "--custom media/rafael/DATA/genomes/references/dbSNPs/GCF_000001405.39.gz,dbSNP,vcf,,,COMMON"
    shell:
        "/home/rafael/Aplicaciones/VEP/ensembl-vep/./vep --fork 32 --species homo_sapiens --everything --keep_csq --assembly GRCh38 --force_overwrite --offline --vcf --compress_output bgzip --cache {params.clinvar} {params.gnomad} {params.IARC} {params.LOVD_APC} {params.LOVD_MSH2} {params.LOVD_MSH6} {params.LOVD_NF1} {params.CADD} {params.COSMIC} --fasta {input.ref} --input_file {input.vcf} --output_file {output.vcf} &> {log}"


rule VEP_custom_tsv:
    input:
        vcf = "results/{group}/Merge/csq_phased_all.vcf.gz",
        ref = config['refGenomeDir'] + "genome.fasta",
    output:
        tsv = "results/{group}/Merge/all_VEP_custom.tsv.gz",
    log:
        "logs/{group}/vep/annotate_custom.log",
    conda:
        "envs/VEP.yaml"
    threads:32
    params:
        clinvar = "--custom /media/rafael/DATA/genomes/references/VEP/ClinVar/clinvar_20220403.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN",
        gnomad = "--custom /media/rafael/DATA/genomes/references/VEP/gnomAD/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz,gnomAD,vcf,exact,0,non_cancer_AC,non_cancer_AN,non_cancer_AF,non_cancer_nhomalt",
        IARC = "--custom /media/rafael/DATA/genomes/references/VEP/IARC/IARC_other.bed.gz,IARC_indel,bed,overlap,0",
        COSMIC = "--custom /media/rafael/DATA/genomes/references/VEP/COSMIC/COSMIC_sorted.bed.gz,COSMICHotspot,bed,overlap,0",
        LOVD_MSH2 = "--custom /media/rafael/DATA/genomes/references/VEP/LOVD/LOVD_MSH2_other.bed.gz,LOVD_MSH2_indel,bed,overlap,0",
        LOVD_MSH6 = "--custom /media/rafael/DATA/genomes/references/VEP/LOVD/LOVD_MSH6_other.bed.gz,LOVD_MSH6_indel,bed,overlap,0",
        LOVD_APC = "--custom /media/rafael/DATA/genomes/references/VEP/LOVD/LOVD_APC_other.bed.gz,LOVD_APC_indel,bed,overlap,0",
        LOVD_NF1 = "--custom /media/rafael/DATA/genomes/references/VEP/LOVD/LOVD_NF1_other.bed.gz,LOVD_NF1_indel,bed,overlap,0",
        CADD = "--custom /media/rafael/DATA/genomes/references/VEP/CADD/whole_genome_SNVs_CADD.vcf.gz,CADD,vcf,exact,0,RAW,PHRED",
        Human_106 ="--custom /home/rafael/.vep/homo_sapiens/106_GRCh38,RefSeq,SIFT,dbSNP,COSMIC,ClinVar,gnomAD",       
        dbNFSP4 = "--custom /media/rafael/DATA/genomes/references/VEP/dbNSFP/dbNSFP4.0a.gz,GERP++_RS,phastCons100way_vertebrate,phyloP100way_vertebrate,phastCons30way_mammalian,phastCons100way_vertebrate",
        dbscSNV = "--custom /media/rafael/DATA/genomes/references/VEP/dbscSNV/dbscSNV1.1_GRCh38.txt.gz,GRCh38",
        DisGeNET = "--custom /media/rafael/DATA/references/VEP/DisGeNET/all_variant_disease_pmid_associations_2020.tsv.gz,disease,rsid",
        dbSNP = "--custom media/rafael/DATA/genomes/references/dbSNPs/GCF_000001405.39.gz,dbSNP,vcf,,,COMMON"
    shell:
        "/home/rafael/Aplicaciones/VEP/ensembl-vep/./vep --fork 32 --species homo_sapiens --everything --assembly GRCh38 --pick --no_stats --keep_csq --force_overwrite --offline --tab --compress_output bgzip --cache {params.clinvar} {params.gnomad} {params.IARC} {params.LOVD_APC} {params.LOVD_MSH2} {params.LOVD_MSH6} {params.LOVD_NF1} {params.CADD} {params.COSMIC} --fasta {input.ref} --input_file {input.vcf} --output_file {output.tsv} &> {log}"


## Alternative Annotation tools
#rule snpeff:
#    input:
#        calls= "results/{group}/Merge/all_VEP_custom.vcf.gz",
         
#    output:
#        calls=  "results/{group}/Merge/all_VEP_custom_SNEPFF.vcf.gz",  # annotated calls (vcf, bcf, or vcf.gz)
#        stats= "results/{group}/Merge/all_VEP_custom_SNEPFF.html",  # summary statistics (in HTML), optional
#        csvstats= "results/{group}/Merge/all_VEP_custom_SNEPFF.csv" # summary statistics in CSV, optional
#    log:
#        "logs/{group}/snpeff/snpeff.log"
    
#    resources:
#        mem_mb=4096
#    params:
#        db="GRCh38.99"
#    threads: 
#        32
#    shell:
#        "java -jar /home/rafael/SnpEff/snpEff.jar -c /home/rafael/SnpEff/snpEff.config -v {params.db} {input.calls} -csvStats {output.csvstats} -stats {output.stats} > {output.calls} 2> {log}"


#rule ANNOVAR: 
#     input: 
#         vcf = "results/{group}/{caller}/annotated/all_VEP_custom_SNEPFF.vcf.gz",
#     params: 
#          humanDB = "/home/rafael/Aplicaciones/annovar/humandb_38/",
#          version = "hg38",
#          ANNOVAR = "/home/rafael/Aplicaciones/annovar/", 
#          output = "results/{group}/annotated/ANNOVAR",
#     log: 
#         "logs/{group}/{caller}/ANNOVAR/Annovar.log"
#     output:
#        "results/{group}/{caller}/annotated/all_VEP_custom_SNEPFF_ANNOVAR.hg38_multianno.vcf", 
#     shell:
#        """
#          {params.ANNOVAR}/table_annovar.pl {input.vcf} {params.humanDB} -buildver {params.version} -out {params.output} \
#          -remove -protocol refGene -operation g  -nastring . -vcfinput --withzyg &> {log}
#        """    



rule bcftools_query:
    input:
        vcf= "results/{group}/Merge/csq_phased_all.vcf.gz",
        
    output:
        "results/{group}/Merge/query_all.tsv",
    log:
        "logs/{group}/bcf_query/query.log"
    
    threads:
        32
    shell:
        "/home/rafael/Aplicaciones/bcftools/bcftools-1.15.1/bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%TBCSQ\n]' {input.vcf} > {output} &> {log}"
        
    
rule bcftools_query_format:
    input:
        #vcf= "results/{group}/Merge/all_DNM2.vcf.gz",# con trios completos si
        vcf= "results/{group}/Merge/csq_phased_all.vcf.gz", 
    output:
        "results/{group}/Merge/query_FORMAT_all.tsv",
    log:
        "logs/{group}/bcf_query/query_format.log"
    
    threads:
        32
    shell:
       # "/home/rafael/Aplicaciones/bcftools/bcftools-1.15.1/bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%AD\t%DP\t%GQ\t%PGT\t%PID\t%PL\t%PS\t%DNM\t%VA\t%VAF\n]' {input.vcf} > {output}"
        "/home/rafael/Aplicaciones/bcftools/bcftools-1.15.1/bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%AD\t%DP\t%GQ\t%PGT\t%PID\t%PL\t%PS\n]' {input.vcf} > {output}"
    

rule vcf_to_tsv:
    input:
        #"results/{group}/Merge/all_DNM2.vcf.gz", # con trios completos si
        "results/{group}/Merge/csq_phased_all.vcf.gz",
    output:
        report(
            "results/{group}/Merge/calls.tsv.gz",
            caption="../report/calls.rst",
            category="Calls",
        ),
    log:
        "logs/{group}/vcf-to-tsv.log",
    conda:
        "envs/rbt.yaml"
    shell:
        "(bcftools view --apply-filters PASS --output-type u {input} | "
        "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
        "gzip > {output}) 2> {log}"


rule InterVar:
    input:
        "results/{group}/Merge/all_VEP_custom.vcf.gz"
    output:
        "results/{group}/InterVar/Test"
    shell:
        "python /home/rafael/Aplicaciones/InterVar/Intervar.py -b hg38 -i {input} --input_type=VCF_m --table_annovar=/home/rafael/Aplicaciones/InterVar/table_annovar.pl  --convert2annovar=/home/rafael/Aplicaciones/InterVar/convert2annovar.pl --annotate_variation=/home/rafael/Aplicaciones/InterVar/annotate_variation.pl --database_locat=/home/rafael/Aplicaciones/InterVar/humandb  --database_intervar=/home/rafael/Aplicaciones/InterVar/intervardb -o {output} &> {log}"

rule Final_filter:
    input:
    	query = "results/{group}/Merge/query_all.tsv",
    	vep = "results/{group}/Merge/all_VEP_custom.tsv.gz",
    	format = "results/{group}/Merge/query_FORMAT_all.tsv",
    	acmg = config['list_ACMG'],
    	allcancer = config['list_AllCancer'],
    	blodcancer = config['list_ALL_MDS_AMLCancer'],
    	childhood = config['list_Childhood_CPS_UOG'],
       
    output:
        directory("results/{group}/Filter_tables/")
    
    threads:
        32
    shell:
        "scripts/Final_filter.R --query {input.query} --vepcustom {input.vep} --format {input.format} --acmg {input.acmg} --allcancer {input.allcancer} --blodcancer {input.blodcancer} --childhood {input.childhood} --output {output} &> {log}"


#########################Maftools#############
##########DESCOMPRIMIR

rule descompress:
	input:
		vcf = "results/{group}/Merge/all_VEP_custom.vcf.gz",
		tsv = "results/{group}/Merge/all_VEP_custom.tsv.gz",
	output:
		vcf = "results/{group}/Merge/all_VEP_custom.vcf",
		tsv = "results/{group}/Merge/all_VEP_custom.tsv",
		#dir = directory("results/{group}/Merge/descompress"),
		#actualdir = directory ("results/{group}/Merge")
	log:
		"logs/{group}/vcf2maf/descompress.log"
	shell:
		"""		 
		 gzip -dk {input.vcf} {input.tsv} &> {log}
		 
		"""
		
		# gzip -dk {input.vcf} {input.tsv} > {output.vcf} {output.tsv} &> {log}
		
########VCF2MAF
##################ONE SAMPLE#################

rule vcf2maf:
		input:
			vcf2maf= "/home/rafael/Aplicaciones/mskcc-vcf2maf-754d68a/vcf2maf.pl",
			vcf = "results/{group}/Merge/all_VEP_custom.vcf",
			ref = config['refGenomeDir'] + "genome.fasta",
		output:
			maf="results/{group}/Merge/maf/all_VEP_custom.maf"
		log:
			"logs/{group}/vcf2maf/vcf2maf.log"
		params:
			normal = "{group}-H",
			tumor = "{group}-D"
				
		shell:
			"""
			perl {input.vcf2maf} --input-vcf {input.vcf} --output-maf {output.maf} --ref-fasta {input.ref} --tumor-id {params.tumor} --normal-id {params.normal} --inhibit-vep --ncbi-build GRCh38
			"""
			
###########MAFTOOLS

rule maftools_onesample:
	input:
		directory= "/media/rafael/WORK/WES_pipeline/Somatic/results/{group}/"  #take care you must give the absolute path to the Rscript
	
	output:
		"results/{group}/ResultsMaftools/plotmafSummary"
	
	log:
		"logs/{group}/Resultsmaftoolsonesample"
		
	shell:
		"""
		Rscript scripts/Maftoolsonesample.R --directory {input.directory} &> {log}
		"""
rule maftools_cohort:
	input:
		directory= "/media/rafael/WORK/WES_pipeline/Somatic/results/"  #take care you must give the absolute path to the Rscript
	output:
		"results/ResultsMaftoolsCohort/plotmafSummary"
	log:
		"logs/Resultsmaftoolscohort"
	shell:
		"""
		Rscript scripts/Maftoolscohort.R --directory {input.directory} &> {log}
		"""

##################################################################
#############Gistic2 - CNVs#################################

rule index_fasta:
	input:
		ref = config['refGenomeDir'] + "genome.fasta"
	output:
		indexref = "/results/ref.fasta.fai"
	log:
		 "logs/indexfasta.log"
	shell:
		"""
		samtools faidx {input.ref} 
		"""

rule VarScan_sh:
	input:
		ref = config['refGenomeDir'] + "genome.fasta",
		normal = config['bamsdir'] + "{group}/recal/{group}-H.bam",
		tumor = config['bamsdir'] + "{group}/recal/{group}-D.bam",
		varscanpath = "/home/rafael/Aplicaciones/VarScan.v2.3.9.jar",
	#dir = "/media/rafael/WORK/WES_pipeline/Somatic/results/{group}"
	
	output:
	#dir = directory("/media/rafael/WORK/WES_pipeline/Somatic/results/{group}"),
		copynumber = "/media/rafael/WORK/WES_pipeline/Somatic/results/{group}/output.copynumber"
	
	shell:
		"""
		cd /media/rafael/WORK/WES_pipeline/Somatic/results/{wildcards.group}
		samtools mpileup -E -q 1 -f {input.ref} {input.normal} {input.tumor} | java -jar {input.varscanpath} copynumber varScan --mpileup 1 #this command doesn't allow the redirection to a log file
		echo "VarScan copynumber" {wildcards.group} "is done"
		cd
		"""


#Manual command to check if samtools and java are okey
#samtools mpileup -E -f /media/rafael/DATA/Genome_GRCh38_100/resources/genome.fasta /media/rafael/Elements/BAMs_Jorge/662266/recal/662266-D.bam /media/rafael/Elements/BAMs_Jorge/662266/recal/662266-H.bam | java -jar /home/rafael/Aplicaciones/VarScan.v2.3.9.jar copynumber varScan --mpileup 1	
		
rule VarScan_copycaller:
	input:
		path = "/home/rafael/Aplicaciones/VarScan.v2.3.9.jar",
		copynumber = "/media/rafael/WORK/WES_pipeline/Somatic/results/{group}/output.copynumber"
	output:
		copynumbercalledgc = "results/{group}/Gistic2/varScan.copynumber.caller" #GC adjustment information
	
	shell:
		"""
		cd /media/rafael/WORK/WES_pipeline/Somatic/results/{wildcards.group}/Gistic2
		java -jar {input.path} copyCaller {input.copynumber} --output-file varScan.copynumber.caller --output-homdel-file varScan.copynumber.called.homdel
		echo "VarScan copynumber" {wildcards.group} "is done"
		cd


        """
#  java -jar {input.path} copyCaller {input.copynumber} --output-file varScan.copynumber.called > {output.copynumbercalledgc}

rule DNAcopy:
	input:
		input = "/media/rafael/WORK/WES_pipeline/Somatic/results/{group}/output.copynumber"
	output:
		seq = "/media/rafael/WORK/WES_pipeline/Somatic/results/{group}/cbs_P1_pvalue"
	log:
	    "logs/{group}/Gistic2/DNAcopy.log"
	shell:
		"Rscript scripts/DNAcopymodificado.R --input {input.input} --output {output.seq}"

from os import path

import pandas as pd

gistic_root = ["/home/rafael/Aplicaciones/GISTIC2/"]
def _mcr_env(mcr_root=None, mcr_ver='v83'):
    if not mcr_root:
        return ""

    # LD path.
    ld_paths = ["{mcr_root}/{mcr_ver}/runtime/glnxa64",
                "{mcr_root}/{mcr_ver}/bin/glnxa64",
                "{mcr_root}/{mcr_ver}/sys/os/glnxa64",
                "${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"]

    ld_path = "LD_LIBRARY_PATH=" + ":".join(ld_paths)
    ld_path = ld_path.format(mcr_root=mcr_root, mcr_ver=mcr_ver)

    # XAPPLRESDIR path.
    xappl_path = ("XAPPLRESDIR={mcr_root}/{mcr_ver}"
                  "/MATLAB_Component_Runtime/{mcr_ver}/X11/app-defaults")
    xappl_path = xappl_path.format(mcr_root=mcr_root, mcr_ver=mcr_ver)

    return ld_path + " " + xappl_path


# Determine gistic and mcr root paths. If both are not given,
# we assume that gistic is available in PATH and that the
# LD_LIBRARY_PATH path has already been configured for the MCR.
#gistic_root = "/home/rafael/Aplicaciones/GISTIC2/".get("gistic_root", "")
#mcr_root = "/home/rafael/Aplicaciones/GISTIC2/".get("mcr_root", None)


mcr_root = path.join("/home/rafael/Aplicaciones/GISTIC2/", "MATLAB_Compiler_Runtime")





rule GISTIC2:
	input:
		gisticpath = "/home/rafael/Aplicaciones/GISTIC2/gistic2",
		ref =  "/home/rafael/Aplicaciones/GISTIC2/refgenefiles/hg18.mat",
		seg = "/media/rafael/WORK/WES_pipeline/Somatic/results/cbs.csv"
	output:
		amp_qplot="results/Gistic2/amp_qplot.pdf"
		#alllesionsfile= "results/{group}/Gistic2/all_lesions.conf.txt",
		#amplificationgenesfile= "results/{group}/Gistic2/amp_genes.conf.txt",
		#deletiongenesfile= "results/{group}/Gistic2/del_genes.conf.txt",
		#gisticfile= "results/{group}/Gistic2/scores.gistic",
		#segmentedcopynumber= "results/{group}/Gistic2/raw_copy_number.pdf"
	params:
		gistic_root="/home/rafael/Aplicaciones/GISTIC2/",
		mcr_env=_mcr_env(mcr_root),
		outputdir = directory("results/Gistic2"),
	log:
		"logs/Gistic2/gistic2.log"
	shell:
		"{input.gisticpath} -b {params.outputdir} -seg {input.seg} -refgene {input.ref} -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme &> {log}"
			
			
###################MSI-SENSOR PRO########################
###Version 1.2.0 is mandatory######
rule msisensorpro:
	input:
		ref =  config['refGenomeDir'] + "genome.fasta",
		tumor= config['bamsdir'] + "{group}/recal/{group}-D.bam",
		normal= config['bamsdir'] + "{group}/recal/{group}-H.bam",
		bed= config['refGenomeDir'] + "Regions.bed"
	output:
		microsatfile = "results/{group}/Msisensor/microsat.list",
		microsat = "results/{group}/Msisensor/report",
		microsatbed = "results/{group}/Msisensor/bedreport"
	log:
		"logs/{group}/msisensor/msisensor.log"
		
	shell:
		"""
		msisensor-pro scan -d {input.ref} -o {output.microsatfile}
		msisensor-pro msi -d {output.microsatfile} -n {input.normal} -t {input.tumor} -e {input.bed} -o {output.microsatbed}
		msisensor-pro msi -d {output.microsatfile} -n {input.normal} -t {input.tumor} -o {output.microsat}
		"""
	
		


