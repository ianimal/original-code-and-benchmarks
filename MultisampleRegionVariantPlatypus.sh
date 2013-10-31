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
region=$1
nameRegion=`echo "$region" | sed -e 's/:/./g'`
workDir=/work/reecygroup/christine/data/
runTimeLog=${workDir}${nameRegion}.platypus.log.txt

#---------------------------------------------
# PARAMETER FILES
#---------------------------------------------
reference=/work/reecygroup/christine/UMD31.fasta
dbSNP_orig=/work/reecygroup/christine/Bos_taurus.vcf
dictionary=/work/reecygroup/christine/UMD31.dict
dbSNP_sort=/work/reecygroup/christine/Bos_taurus.sort.vcf
dbSNP_work=/work/reecygroup/christine/Bos_taurus.dbSNP.vcf

#---------------------------------------------
# PROGRAM VERSIONS
#---------------------------------------------
PICARD=/work/reecygroup/picard/
GATK=/work/reecygroup/variant.calling.suite/GATK2.7/
SAMTOOLS=/work/reecygroup/variant.calling.suite/tophat2/
BCFTOOLS=/work/reecygroup/variant.calling.suite/tophat2/bcftools/
vcfutils=/work/reecygroup/variant.calling.suite/tophat2/bcftools/vcfutils.pl 
JAVADIR=/work/reecygroup/variant.calling.suite/java1.7/bin/
platypus=/work/reecygroup/variant.calling.suite/Platypus_0.2.3/

#---------------------------------------------
# FUNCTIONS
#---------------------------------------------

function GATKpipeline () {
	RIGHT_NOW=$(date +"%s")
	echo "Beginning: $FUNCNAME" $RIGHT_NOW >> ${runTimeLog}
		
	function platypus () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 1st" ${RIGHT_NOW} >> ${runTimeLog}
		python ${platypus}Platypus.py callVariants --refFile=${reference} \
		--bamFiles=\
		${workDir}66522.realigned.recal.bam,\
   		${workDir}66524.realigned.recal.bam,\
   		${workDir}68268.realigned.recal.bam,\
   		${workDir}68269.realigned.recal.bam,\
   		${workDir}68270.realigned.recal.bam,\
   		${workDir}68271.realigned.recal.bam,\
   		${workDir}68272.realigned.recal.bam,\
   		${workDir}68273.realigned.recal.bam,\
   		${workDir}68274.realigned.recal.bam,\
   		${workDir}68275.realigned.recal.bam,\
   		${workDir}68276.realigned.recal.bam,\
   		${workDir}68277.realigned.recal.bam,\
   		${workDir}68278.realigned.recal.bam,\
   		${workDir}68279.realigned.recal.bam,\
   		${workDir}68280.realigned.recal.bam,\
   		${workDir}68281.realigned.recal.bam,\
   		${workDir}68282.realigned.recal.bam,\
   		${workDir}68283.realigned.recal.bam,\
   		${workDir}68284.realigned.recal.bam,\
   		${workDir}68285.realigned.recal.bam,\
   		${workDir}68286.realigned.recal.bam,\
   		${workDir}68287.realigned.recal.bam,\
   		${workDir}68288.realigned.recal.bam,\
   		${workDir}68289.realigned.recal.bam,\
   		${workDir}68290.realigned.recal.bam,\
   		${workDir}68291.realigned.recal.bam,\
   		${workDir}68292.realigned.recal.bam,\
   		${workDir}68293.realigned.recal.bam,\
   		${workDir}68294.realigned.recal.bam,\
   		${workDir}68295.realigned.recal.bam,\
   		${workDir}68296.realigned.recal.bam,\
   		${workDir}68297.realigned.recal.bam,\
   		${workDir}68298.realigned.recal.bam,\
   		${workDir}68299.realigned.recal.bam,\
   		${workDir}68300.realigned.recal.bam,\
   		${workDir}68301.realigned.recal.bam,\
   		${workDir}68303.realigned.recal.bam,\
   		${workDir}68305.realigned.recal.bam,\
   		${workDir}68306.realigned.recal.bam,\
   		${workDir}68307.realigned.recal.bam,\
   		${workDir}68309.realigned.recal.bam,\
   		${workDir}68310.realigned.recal.bam,\
   		${workDir}68311.realigned.recal.bam,\
   		${workDir}68312.realigned.recal.bam,\
   		${workDir}68313.realigned.recal.bam,\
   		${workDir}68314.realigned.recal.bam,\
   		${workDir}68315.realigned.recal.bam,\
   		${workDir}68316.realigned.recal.bam,\
   		${workDir}68317.realigned.recal.bam,\
   		${workDir}68318.realigned.recal.bam,\
   		${workDir}68319.realigned.recal.bam,\
   		${workDir}68320.realigned.recal.bam,\
   		${workDir}68321.realigned.recal.bam,\
   		${workDir}68322.realigned.recal.bam,\
   		${workDir}68323.realigned.recal.bam,\
   		${workDir}68324.realigned.recal.bam,\
   		${workDir}68325.realigned.recal.bam,\
   		${workDir}68326.realigned.recal.bam,\
   		${workDir}68327.realigned.recal.bam,\
   		${workDir}68328.realigned.recal.bam,\
   		${workDir}68329.realigned.recal.bam,\
   		${workDir}68330.realigned.recal.bam,\
   		${workDir}68331.realigned.recal.bam,\
   		${workDir}68332.realigned.recal.bam,\
   		${workDir}9726.realigned.recal.bam,\
   		${workDir}9741.realigned.recal.bam\
   		--regions=${region} \
		--output=${workDir}${nameRegion}.platypus.vcf
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
	}
	platypus
	bzip2 ${workDir}${nameRegion}.platypus.vcf
	
	RIGHT_NOW=$(date +"%s")
	echo "ENDING: $FUNCNAME " ${RIGHT_NOW} >> ${runTimeLog}
}
GATKpipeline

