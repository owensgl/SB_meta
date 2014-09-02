#!/bin/bash

#USAGE: ./SB_gatk_pipe.sh <plate identifier> 
#Sept 2013 GLO
#Pipeline for using GATK and BWA to call snps
plate="$1"
Gpath='/home/owens/bin'
Hpath='/home/owens/SB'
rawdata='/home/owens/raw_data'
bwa='/home/owens/bin/bwa-0.7.9a'
demultiplex='/home/owens/bin/GBS_fastq_Demultiplexer_v8.GO.pl'
barcodes="/home/owens/SB/Barcodes.$plate.txt"
ref='/home/owens/ref/Gasterosteus_aculeatus.BROADS1.73.dna.toplevel.fa'
stampyref='/home/owens/ref/Gasterosteus_aculeatus.BROADS1.73.dna.toplevel'
cleandata="/home/owens/SB/$plate.clean_data_paired"
trimmeddata="/home/owens/SB/$plate.trimmed_data_paired"
unpaired="/home/owens/SB/$plate.trimmed_data_unpaired"
sam="/home/owens/SB/$plate.sam"
bam="/home/owens/SB/$plate.bam"
log="/home/owens/SB/$plate.log"
gvcf="/home/owens/SB/gvcf"
tabix='/home/owens/bin/tabix-0.2.6'
vcftools='/home/owens/bin/vcftools_0.1.12a/bin'
picardtools='/home/owens/bin/picard-tools-1.114'
stampy='/home/owens/bin/stampy-1.0.23/stampy.py'
trim="/home/owens/bin/Trimmomatic-0.32"
project='sb1'
#Prerequisites:
#Index the reference for GATK and BWA
#Change the variables to fit your file structure and ensure the folders exist.

echo "Starting GBS processing pipeline on $plate."
#Make directories if they don't exist
if [ ! -d "$sam" ]; then
	mkdir $sam
fi
if [ ! -d "$bam" ]; then
	mkdir $bam
fi
if [ ! -d "$log" ]; then
	mkdir $log
fi
if [ ! -d "$gvcf" ]; then
        mkdir $gvcf
fi
if [ ! -d "$cleandata" ]; then
        mkdir $cleandata
fi
if [ ! -d "$trimmeddata" ]; then
        mkdir $trimmeddata
fi
if [ ! -d "$unpaired" ]; then
        mkdir $unpaired
fi

#Demultiplex
#perl $demultiplex $barcodes $rawdata/"$plate"_R1.fastq $rawdata/"$plate"_R2.fastq $cleandata/$project


#ls $cleandata | grep -v "nobar" | sed s/_R1//g | sed s/_R2//g | sed s/.fastq// | uniq  > $Hpath/Samplelist.$plate.txt
#Trim the data using Trimmomatic. Removes bad reads and illumina adapter contamination.
#while read prefix
#do
#java -jar $trim/trimmomatic-0.32.jar PE -phred33 $cleandata/"$prefix"_R1.fastq $cleandata/"$prefix"_R2.fastq $trimmeddata/"$prefix"_R1.fastq $unpaired/"$prefix"_unR1.fastq $trimmeddata/"$prefix"_R2.fastq $unpaired/"$prefix"_unR2.fastq ILLUMINACLIP:$trim/adapters/TruSeq3-PE.fa:2:30:10:8:T SLIDINGWINDOW:4:15 MINLEN:36
#done < Samplelist.$plate.txt
ls  $trimmeddata | grep -v "nobar" | sed s/_R1//g | sed s/_R2//g | sed s/.fastq// | uniq  > $Hpath/Samplelist.$plate.txt

###Align using BWA. Turn from sam to bam. Sort by coordinate and add read group data.
while read prefix
do
	$bwa/bwa aln -t 8 $ref $trimmeddata/"$prefix"_R1.fastq 1> $trimmeddata/"$prefix"_R1.sai
	$bwa/bwa aln -t 8 $ref $trimmeddata/"$prefix"_R2.fastq 1> $trimmeddata/"$prefix"_R2.sai
	$bwa/bwa sampe $ref $trimmeddata/"$prefix"_R1.sai $trimmeddata/"$prefix"_R2.sai $trimmeddata/"$prefix"_R1.fastq $trimmeddata/"$prefix"_R2.fastq 1> $sam/$prefix.sam 2> $log/$prefix.bwasampe.log
	samtools view -Sb $sam/$prefix.sam > $bam/$prefix.bam	
	$stampy -g $stampyref -h $stampyref -t8 --bamkeepgoodreads -M $bam/$prefix.bam -o $bam/$prefix.stampy.bam 2> $log/$prefix.stampy.log
	java -jar $picardtools/CleanSam.jar INPUT=$bam/$prefix.stampy.bam OUTPUT=$bam/$prefix.clean.bam 2> $log/$prefix.cleansam.log
	java -jar $picardtools/SortSam.jar INPUT=$bam/$prefix.clean.bam OUTPUT=$bam/$prefix.sort.bam SORT_ORDER=coordinate 2> $log/$prefix.sortsam.log
	java -jar $picardtools/AddOrReplaceReadGroups.jar I=$bam/$prefix.sort.bam O= $bam/$prefix.sortrg.bam SORT_ORDER=coordinate RGID=$prefix RGLB=$project RGPL=ILLUMINA RGPU=$project RGSM=$prefix CREATE_INDEX=True 2> $log/$prefix.addRG.log
	rm $trimmeddata/"$prefix"_R1.sai
	rm $trimmeddata/"$prefix"_R2.sai
	rm $sam/$prefix.sam
	rm $bam/$prefix.bam
	rm $bam/$prefix.stampy.bam
	rm $bam/$prefix.clean.bam
	rm $bam/$prefix.sort.bam
done < Samplelist.$plate.txt


#Make bam.list for GATK
ls -d $bam/*.* | grep sortrg.bam  > $Hpath/bamlist.$plate.list

#identify local indels

java -Xmx2g -jar $Gpath/GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R $ref \
   -I $Hpath/bamlist.$plate.list \
   -nt 8 \
   -log $log/$plate.RealignerTargetCreator.log \
   -o $Hpath/$plate.realign.intervals

#Realign around local indels
while read prefix
do
java -Xmx4g -jar $Gpath/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $ref \
	-I $bam/$prefix.sortrg.bam \
	-targetIntervals $Hpath/$plate.realign.intervals \
	-o $bam/$prefix.realign.bam \
	-log $log/$plate.$prefix.IndelRealigner.log 

#Call GATK HaplotypeCaller
java -Xmx18g -jar $Gpath/GenomeAnalysisTK.jar \
	-nct 8 \
	-l INFO \
	-R $ref \
	-log $log/$plate.$prefix.HaplotypeCaller.log \
	-T HaplotypeCaller \
	-I  $bam/$prefix.realign.bam \
	--emitRefConfidence GVCF \
	--max_alternate_alleles 2 \ 
	-variant_index_type LINEAR \
	-variant_index_parameter 128000 \
	-o $gvcf/$prefix.GATK.gvcf.vcf
done < $Hpath/Samplelist.$plate.txt

#Make input list for GATK GenotypeGVCFs
tmp=""
while read prefix
do
        tmp="$tmp --variant $gvcf/$prefix.GATK.gvcf.vcf"
done < $Hpath/Samplelist.$plate.txt

#Genotype all gvcf together into one vcf file
java -Xmx18g -jar $Gpath/GenomeAnalysisTK.jar \
	-nt 8 \
	-l INFO \
	-R $ref \
	-log $log/$plate.GenotypeGVCFs.log \
	-T GenotypeGVCFs \
	$tmp \
	-o $Hpath/SB.GATK.total.vcf \
	-inv \
	--max_alternate_alleles 4


exit

