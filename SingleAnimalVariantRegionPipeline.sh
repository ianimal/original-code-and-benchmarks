#!/bin/bash

#---------------------------------------------
# Script for variant calling 
#
# Original author: marlies dolezal (may-june 2013)
# Alterations: christine baes (june, july, august 2013)
# Modifications: eric fritz-waters (september, october 2013)
#---------------------------------------------

#---------------------------------------------
# DIRECTORIES
#---------------------------------------------
FILE=$1   
region=$2
bull=`echo "$FILE" | sed 's/.bam//g'`
temp=${bull}.temp/
log=${bull}.logs/
results=${bull}.results/

mkdir $temp
mkdir $log
mkdir $results

#---------------------------------------------
# PARAMETER FILES
#---------------------------------------------
reference=/home/eric.fritz/christine/UMD31.sorted.fa
dbSNP_orig=/home/eric.fritz/christine/Bos_taurus.vcf
dictionary=/home/eric.fritz/christine/UMD31.sorted.dict
dbSNP_sort=/home/eric.fritz/christine/Bos_taurus.sort.vcf
dbSNP_work=/home/eric.fritz/christine/Bos_taurus.dbSNP.vcf

#---------------------------------------------
# PROGRAMME VERSIONS
#---------------------------------------------
PICARD=/home/eric.fritz/picard/
GATK=/home/eric.fritz/GATK2.7/
SAMTOOLS=/home/eric.fritz/tophat2/
BCFTOOLS=/home/eric.fritz/tophat2/bcftools/
vcfutils=/home/eric.fritz/tophat2/bcftools/vcfutils.pl 
platypus=/home/eric.fritz/Platypus_0.2.3/
JAVADIR=/home/eric.fritz/java1.7/bin/


#---------------------------------------------
# FUNCTIONS
#---------------------------------------------

function GATKpipeline () {
	runTimeLog=${log}${region}.runtime.txt
	RIGHT_NOW=$(date +"%s")
	echo "Beginning: $FUNCNAME" $RIGHT_NOW >> ${runTimeLog}
	
	function extractRegion () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME " $RIGHT_NOW >> ${runTimeLog}
		${JAVADIR}java -Xmx3g -jar ${GATK}GenomeAnalysisTK.jar \
		-T PrintReads \
		-R ${reference} \
		-I ${FILE} \
		-o ${temp}${region}.bam \
		-L ${region} \
		-allowPotentiallyMisencodedQuals
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME " ${RIGHT_NOW} >> ${runTimeLog}
	}
	extractRegion

	function RealignerTargetCreator () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" $RIGHT_NOW >> ${runTimeLog}
		${JAVADIR}java -Xmx3g -jar ${GATK}GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-R ${reference} \
		-I ${temp}${region}.bam \
		-known ${dbSNP_work} \
		-o ${temp}${region}.RTCintervals.list \
		-L ${region} \
		-l INFO \
		-log ${log}${region}.RealignerTargetCreator.log \
		-allowPotentiallyMisencodedQuals
		#attention -o file must be named .list or .bed for IndelRealigner
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
	}
	RealignerTargetCreator

	function IndelRealigner () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx3g -jar ${GATK}GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ${reference} \
		-I ${temp}${region}.bam \
		-targetIntervals ${temp}${region}.RTCintervals.list \
		-known ${dbSNP_work} \
		-L ${region} \
		--consensusDeterminationModel USE_SW \
		-o ${temp}${region}.realigned.bam \
		-l INFO \
		-log ${log}${region}.IndelRealigner.log \
		-allowPotentiallyMisencodedQuals
		#@ -known Non-indel variants in these files will be ignored. so we can use
		#batch download of all variants from NCBI or ensembl
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
	}
	IndelRealigner

	function BaseRecalibrator () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME Part 1" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx3g -jar ${GATK}GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-R ${reference} \
		-L ${region} \
		-I ${temp}${region}.realigned.bam \
     	-knownSites ${dbSNP_work} \
		-o ${temp}${region}.realigned.grp \
		-l INFO \
		-allowPotentiallyMisencodedQuals #\
		#--plot_pdf_file ${Reg}${bull}.pre_recal.pdf
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME Part 1" ${RIGHT_NOW} >> ${runTimeLog}
		
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME Part 2" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx3g -jar ${GATK}GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-R ${reference} \
		-L ${region} \
		-I ${temp}${region}.realigned.bam \
		-BQSR ${temp}${region}.realigned.grp \
     	-knownSites ${dbSNP_work} \
		-o ${temp}${region}.realigned.recal.grp \
		-l INFO \
		-allowPotentiallyMisencodedQuals #\
		#--plot_pdf_file ${Reg}${bull}.post_recal.pdf
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME Part 2" ${RIGHT_NOW} >> ${runTimeLog}
			
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME Part 3" ${RIGHT_NOW} >> ${runTimeLog}		
		${JAVADIR}java -Xmx3g -jar ${GATK}GenomeAnalysisTK.jar \
		-T PrintReads \
		-R ${reference} \
		-L ${region} \
		-I ${temp}${region}.realigned.bam \
		-BQSR ${temp}${region}.realigned.grp \
		-o  ${temp}${region}.realigned.recal.bam \
		-l INFO \
		-log ${log}${region}.realigned.recal.bam.log \
		-allowPotentiallyMisencodedQuals
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME Part 3" ${RIGHT_NOW} >> ${runTimeLog}
	}
	BaseRecalibrator

	function ReduceReads () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx3g -jar ${GATK}GenomeAnalysisTK.jar \
		-T ReduceReads \
		-R ${reference} \
		-L ${region} \
		-I ${temp}${region}.realigned.recal.bam \
		-o ${temp}${region}.realigned.recal.reduced.bam \
		-log ${log}${region}.realigned.recal.reduced.bam.log \
		-allowPotentiallyMisencodedQuals
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
	}
	#ReduceReads

	function UnifiedGenotyper () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx3g -jar ${GATK}GenomeAnalysisTK.jar \
		-T UnifiedGenotyper \
		-R ${reference} \
		-L ${region} \
		-I ${temp}${region}.realigned.recal.bam \
		-o ${results}${region}.UG.raw.vcf \
		--dbsnp $dbSNP_work \
		-out_mode EMIT_VARIANTS_ONLY \
		-stand_call_conf 30.0 \
		-stand_emit_conf 30.0 \
		--genotype_likelihoods_model BOTH \
		-l INFO  \
		-log ${log}${region}.UnifiedGenotyper.log \
		-allowPotentiallyMisencodedQuals
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
	}
	UnifiedGenotyper

	function HaplotypeCaller () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx3g -jar ${GATK}GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R ${reference} \
		-L ${region} \
		-I ${temp}${region}.realigned.recal.bam \
		--dbsnp $dbSNP_work \
		-stand_call_conf 30.0 \
		-stand_emit_conf 30.0 \
		-o ${results}${region}.HC.raw.vcf \
		-l INFO \
		-log ${log}${region}.HaplotypeCaller.log \
		-allowPotentiallyMisencodedQuals
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
	}
	#HaplotypeCaller

	function mpileup () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 1st" ${RIGHT_NOW} >> ${runTimeLog}
		${SAMTOOLS}samtools mpileup -f ${reference} -r ${region} -A -C50 -E -D -S -u ${temp}${region}.bam | \
		${BCFTOOLS}bcftools view -cvgb - > ${temp}${region}.raw.bcf 2> ${temp}${region}.raw.bcf.log 
		${BCFTOOLS}bcftools view ${temp}${region}.raw.bcf | ${vcfutils} varFilter -D800 - > ${results}${region}.flt.vcf
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 1st" ${RIGHT_NOW} >> ${runTimeLog} 
		
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 2nd" ${RIGHT_NOW} >> ${runTimeLog}
		${SAMTOOLS}samtools mpileup -f ${reference} -r ${region} -A -C50 -E -D -S -u ${temp}${region}.realigned.bam | \
		${BCFTOOLS}bcftools view -cvgb - > ${temp}${region}.raw2.bcf 2> ${temp}${region}.raw2.bcf.log 
		${BCFTOOLS}bcftools view ${temp}${region}.raw2.bcf | ${vcfutils} varFilter -D800 - > ${results}${region}.realigned.flt.vcf 
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 2nd" ${RIGHT_NOW} >> ${runTimeLog}
		
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
		${SAMTOOLS}samtools mpileup -f ${reference} -r ${region} -A -C50 -E -D -S -u ${temp}${region}.realigned.recal.bam | \
		${BCFTOOLS}bcftools view -cvgb - > ${temp}${region}.raw3.bcf 2> ${temp}${region}.raw3.bcf.log 
		${BCFTOOLS}bcftools view ${temp}${region}.raw3.bcf | ${vcfutils} varFilter -D800 - > ${results}${region}.realigned.recal.flt.vcf 
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
	}
	mpileup
		
	function platypus () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 1st" ${RIGHT_NOW} >> ${runTimeLog}
		python ${platypus}Platypus.py callVariants --refFile=${reference} --bamFiles=${temp}${region}.bam --output=${results}${region}.platypus.vcf
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 1st" ${RIGHT_NOW} >> ${runTimeLog}
		
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 2nd" ${RIGHT_NOW} >> ${runTimeLog}
		python ${platypus}Platypus.py callVariants --refFile=${reference} --bamFiles=${temp}${region}.realigned.bam --output=${results}${region}.realigned.platypus.vcf 
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 2nd" ${RIGHT_NOW} >> ${runTimeLog}
		
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
		python ${platypus}Platypus.py callVariants --refFile=${reference} --bamFiles=${temp}${region}.realigned.recal.bam --output=${results}${region}.realigned.recal.platypus.vcf 
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
	}
	platypus
	
	RIGHT_NOW=$(date +"%s")
	echo "ENDING: $FUNCNAME " ${RIGHT_NOW} >> ${runTimeLog}
}
GATKpipeline

