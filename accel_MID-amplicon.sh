#!/bin/bash
# MID analysis using fgbio
# S Sandhu 160911

set -e
set -x

for f in *_R1_001.fastq.gz
do
fq1=$f
fq2=${f%_R1*}_R2_001.fastq.gz
OUT=${fq1%%_L001*}

anno=/tools/snpEff
fgbio=/home/sandhu/tools/fgbio/target/scala-2.11/fgbio-0.1.4-SNAPSHOT.jar
picard=/seq/tools/picard/build/libs/picard.jar
trim=/tools/trimmomatic36

dbnsfp=/tools/snpEff/data/GRCh37.75/dbNSFP2.9.txt.gz
cosmic=/tools/snpEff/data/GRCh37.75/CosmicCodingMuts.vcf
cosmicNonCod=/tools/snpEff/data/GRCh37.75/CosmicNonCodingVariants.vcf
clinvar=/tools/snpEff/data/GRCh37.75/clinvar/clinvar_20170516.vcf

REF="/seq/refgenomes/Homo_sapiens_assembly19broad.fasta"
BEDDIR="/home/smrtanalysis/bedfiles"
GATK="java -Xmx64g -jar /tools/GenomeAnalysisTK-2014.4-2-g9ad6aa8/GenomeAnalysisTK.jar"
MASTER="${BEDDIR}/18-2132_EGFR_MID_Masterfile.txt"

# make a (nonmerged) bed from  Masterfile
awk '{print $1,$2,$3,$4}' OFS="\t" $MASTER > nonmerged_targets.bed

# make a merged bedfile from nonmerged
bedtools merge -nms -i nonmerged_targets.bed | sed 's/;.*//' > merged_targets.tmp
awk '{print $1,$2,$3,"+",$4}' OFS="\t" merged_targets.tmp > merged_targets.bed

### overlap-removed non-merged target BED ###
# get all overlapping regions for target bedfile

bedtools intersect \
    -a nonmerged_targets.bed -b nonmerged_targets.bed |
    awk 'a[$1FS$2FS$3]++' OFS="\t" > overlapped_regions.bed

bedtools subtract \
    -a nonmerged_targets.bed -b overlapped_regions.bed \
 > nonmerged_noolaps_targets.tmp
awk '{print $1,$2,$3,"+",$4}' OFS="\t" nonmerged_noolaps_targets.tmp > nonmerged_noolaps_targets.bed

#mv overlapped_regions.bed bed
BED=merged_targets.bed
NOOLAPBED=nonmerged_noolaps_targets.bed


########## PRE-MID analysis & POST-MID with fgbio ##############
#generate index sequence by keep 1st 5 bases of each read
java -jar ${trim}/trimmomatic-0.36.jar SE -threads 20 $fq2 $OUT.R2.MID.fq.gz \
 CROP:10

# MID fastq file
mid=$OUT.R2.MID.fq.gz

######## crop 1st 10bp MID
java -jar $trim/trimmomatic-0.36.jar SE -threads 20 $fq2 $OUT.trimd.R2.fq.gz \
    HEADCROP:10

######## Adapter Trim, & crop 1st 5bp MID
#java -jar $trim/trimmomatic-0.36.jar PE -threads 20 $fq1 $fq2 $OUT.trimd.R1.fq.gz \
 #   $OUT_unpairedR1.fq.gz $OUT.trimd.R2.fq.gz $OUT_unpairedR2.fq.gz \
  #  ILLUMINACLIP:${trim}/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:10


#echo "***********bwa align 2s to hg19 and collect picard metrics********"
#align using bwa with verbosity 1 (to print errors only), mark sec alignments and add RG
bwa mem $REF $fq1 $OUT.trimd.R2.fq.gz -M -t 16 -v 1 > $OUT.sam

# Sort, Add RG for nopclip bam
java -jar $picard AddOrReplaceReadGroups I=$OUT.sam O=$OUT.preMID.bam \
 SO=coordinate VALIDATION_STRINGENCY=STRICT RGID=Accel RGLB=MID RGSM=$OUT RGPL=Illumina \
 RGPU=Miseq CREATE_INDEX=TRUE

# query sort for primerclip
java -jar $picard SortSam I=$OUT.sam O=$OUT.qsort.sam SO=queryname
    
#PrimerCLIP before MID
#/home/sandhu/.local/bin/primerclip $MASTER $OUT.qsort.sam preM${OUT}.preMID.pclip.sam
#/home/smrtanalysis/.local/bin/primerclip $MASTER $OUT.qsort.sam preM${OUT}.preMID.pclip.sam

# Sort
#java -jar $picard AddOrReplaceReadGroups I=preM${OUT}.preMID.pclip.sam O=$OUT.preMID.pclip.bam \
# SO=coordinate VALIDATION_STRINGENCY=STRICT RGID=Accel RGLB=MID RGSM=$OUT RGPL=Illumina \
# RGPU=Miseq CREATE_INDEX=TRUE

## variant call preMID primer trimmed bam
lofreq call --call-indels -f $REF $OUT.preMID.bam -l $BED -o $OUT.preMID.noPclip.lf.vcf

#bedtools to get per amplicon coverage
coverageBed -abam $OUT.bamAnnotdWumi.bam -b $BED > $OUT.bamAnnotdWumi.cov
coverageBed -abam $OUT.preMID.bam -b $BED -d > $OUT.preMID.covd
coverageBed -abam $OUT.preMID.bam -b $BED > $OUT.preMID.cov

### MID processing ########
# annotate bam with MIDs from the I2 file
java -Xmx64g -jar $fgbio AnnotateBamWithUmis -i $OUT.preMID.bam -f $mid -o $OUT.bamAnnotdWumi.bam

#RevertSam to sanitise
java -jar $picard RevertSam I=$OUT.bamAnnotdWumi.bam O=$OUT.sanitised.bam \
    SANITIZE=true REMOVE_DUPLICATE_INFORMATION=false REMOVE_ALIGNMENT_INFORMATION=false

#SetMateInfo to bring MQ tags
java -jar $fgbio SetMateInformation -i $OUT.sanitised.bam -o $OUT.sanitised.setmateinfo.bam

#sort by queryname
java -jar $picard SortSam I=$OUT.sanitised.setmateinfo.bam \
    O=$OUT.sanitised.setmateinfo.qsort.bam SO=queryname

#Group Reads by UMI
java -jar $fgbio GroupReadsByUmi -s adjacency --edits 1 \
    -i $OUT.sanitised.setmateinfo.qsort.bam -o $OUT.groupdbyMID.bam

#Make consensus N molecules to make consensus Out is coordinate sorted
java -jar $fgbio CallMolecularConsensusReads -M 3 \
    -i $OUT.groupdbyMID.bam -o $OUT.M3.consensus.bam -r $OUT.M3.notused4consensus.bam

java -jar $fgbio CallMolecularConsensusReads \
    -i $OUT.groupdbyMID.bam -o $OUT.consensus.bam -r $OUT.notused4consensus.bam

#bamtofastq
java -Xmx128g -jar $picard SamToFastq I=$OUT.M3.consensus.bam \
    F=$OUT.consensus.R1.fq F2=$OUT.consensus.R2.fq FU=$OUT.unpaired.fq

#Realign
bwa mem $REF $OUT.consensus.R1.fq $OUT.consensus.R2.fq -t 24 -M > ${OUT}_consensus.sam
#bwa mem $REF $OUT.consensus.R1.fq $OUT.consensus.R2.fq -t 24 -M -I 100,50,50,500 > ${OUT}_consensus.sam
#bwa mem $REF $OUT.pclip.consensus.R1.fq $OUT.pclip.consensus.R2.fq -t 24 -M > ${OUT}_pclip.consensus.sam

## qsort for primer-clip
java -jar $picard SortSam I=${OUT}_consensus.sam \
    O=$OUT.consensus.qsort.sam SO=queryname

#PrimerCLIP
/home/smrtanalysis/.local/bin/primerclip $MASTER $OUT.consensus.qsort.sam ${OUT}.fgbio.pclip.sam
#/home/smrtanalysis/.local/bin/primerclip $MASTER $OUT.pclip.consensus.qsort.sam pclip${OUT}.pclip.fgbio.pclip.sam

#sort and Add read groups
java -Xmx64g -jar $picard AddOrReplaceReadGroups I=${OUT}.fgbio.pclip.sam \
    O=$OUT.fgbio.bam RGID=MID-Amp RGLB=$OUT RGSM=NA12878 RGPL=Illumina \
    RGPU=MiSeq SO=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=STRICT

# make BED file for picard
cd=$PWD
cp $BED $cd/BED.bed
samtools view -H $OUT.fgbio.bam > header.txt
cat header.txt BED.bed > BED.picard.bed
pBED=BED.picard.bed

# Collect TargetedPCRMetrics from the TWO BAMs #######
java -Xmx65g -jar $picard CollectTargetedPcrMetrics I=$OUT.preMID.bam \
    O=$OUT.preMID.targetPCRmetrics.txt TI=$pBED AI=$pBED R=$REF \
    PER_TARGET_COVERAGE=$OUT.preMID.perTargetCov.txt

java -Xmx64g -jar $picard CollectTargetedPcrMetrics I=$OUT.fgbio.bam \
    O=$OUT.fgbio.targetPCRmetrics.txt TI=$pBED AI=$pBED R=$REF \
    PER_TARGET_COVERAGE=$OUT.fgbio.perTargetCov.txt

#BEDtools
coverageBed -abam $OUT.fgbio.bam -b $BED > $OUT.fgbio.cov
coverageBed -abam $OUT.fgbio.bam -b $BED -d > $OUT.fgbio.covd

coverageBed -abam $OUT.fgbio.bam -b $NOOLAPBED \
    > $OUT.noolap.fgbio.cov
coverageBed -abam $OUT.fgbio.bam -b $NOOLAPBED \
    -d > $OUT.noolap.fgbio.covd

# pre MID cov covd files with noolap BED
coverageBed -abam $OUT.preMID.bam -b $NOOLAPBED \
    > $OUT.noolap.preMID.cov
coverageBed -abam $OUT.preMID.bam -b $NOOLAPBED \
    -d > $OUT.noolap.preMID.covd

############# BQSR fgbio bam      ##########
#Create target interval for Indelrealigner
$GATK  -T RealignerTargetCreator -R $REF \
    -known $anno/data/GRCh37.75/1000G_phase1.indels.b37.vcf \
    -known $anno/data/GRCh37.75/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -I $OUT.fgbio.bam -L $BED -o ${OUT}.forIndelRealigner.intervals

#GATK Indel Realignment
$GATK  -T IndelRealigner -R $REF -I $OUT.fgbio.bam \
    -known $anno/data/GRCh37.75/1000G_phase1.indels.b37.vcf \
    -known $anno/data/GRCh37.75/Mills_and_1000G_gold_standard.indels.b37.vcf \
    --targetIntervals ${OUT}.forIndelRealigner.intervals \
    -L $BED -o ${OUT}.realigned.bam

$GATK  -T BaseRecalibrator -R $REF -I $OUT.fgbio.bam \
    --knownSites $anno/data/GRCh37.75/dbsnp141.vcf -nct 12 \
    --knownSites $anno/data/GRCh37.75/1000Genome_extra.vcf \
    --knownSites $anno/data/GRCh37.75/1000G_phase1.snps.high_confidence.b37.vcf \
    --knownSites $anno/data/GRCh37.75/hapmap_3.3.b37.vcf \
    -L $BED -o ${OUT}.recal_data.table

#generate Recalibrated bam
$GATK  -T PrintReads -R $REF -I ${OUT}.realigned.bam \
    -BQSR ${OUT}.recal_data.table -o ${OUT}.fgbio.bqsrCal.bam

#########################################
#### BQSR pre_MID bam####################
#Create target interval for Indelrealigner
$GATK  -T RealignerTargetCreator -R $REF \
    -known $anno/data/GRCh37.75/1000G_phase1.indels.b37.vcf \
    -known $anno/data/GRCh37.75/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -I $OUT.preMID.bam -L $BED -o ${OUT}.forIndelRealigner.intervals

#GATK Indel Realignment
$GATK  -T IndelRealigner -R $REF -I $OUT.preMID.bam \
    -known $anno/data/GRCh37.75/1000G_phase1.indels.b37.vcf \
    -known $anno/data/GRCh37.75/Mills_and_1000G_gold_standard.indels.b37.vcf \
    --targetIntervals ${OUT}.forIndelRealigner.intervals \
    -L $BED -o ${OUT}.preMID.realigned.bam

$GATK  -T BaseRecalibrator -R $REF -I $OUT.preMID.bam \
    --knownSites $anno/data/GRCh37.75/dbsnp141.vcf -nct 12 \
    --knownSites $anno/data/GRCh37.75/1000Genome_extra.vcf \
    --knownSites $anno/data/GRCh37.75/1000G_phase1.snps.high_confidence.b37.vcf \
    --knownSites $anno/data/GRCh37.75/hapmap_3.3.b37.vcf \
    -L $BED -o ${OUT}.recal_data.table

#generate Recalibrated bam
$GATK  -T PrintReads -R $REF -I ${OUT}.preMID.realigned.bam \
    -BQSR ${OUT}.recal_data.table -o ${OUT}.preMID.bqsrCal.bam
################################################################
###############################################################

#low freq variant calls

lofreq call --call-indels -f $REF ${OUT}.fgbio.bqsrCal.bam -l $BED \
    -o $OUT.fgbio.bqsr.lf.vcf
lofreq call --call-indels -f $REF ${OUT}.preMID.bqsrCal.bam -l $BED \
    -o $OUT.preMID.bqsr.lf.vcf

lofreq call --call-indels -f $REF ${OUT}.fgbio.bam -l $BED \
    -o $OUT.fgbio.NObqsr.lf.vcf

#lofreq vcf merge overlapping indels
/seq/tools/lofreq_star-2.1.2/lofreq2_indel_ovlp.py $OUT.fgbio.bqsr.lf.vcf \
    > $OUT.fgbio.bqsr.lf.mergeIndels.vcf
/seq/tools/lofreq_star-2.1.2/lofreq2_indel_ovlp.py $OUT.preMID.bqsr.lf.vcf \
    > $OUT.preMID.bqsr.lf.mergeIndels.vcf

#germline var calling GATK HC
$GATK -T HaplotypeCaller -I $OUT.fgbio.bqsrCal.bam -R $REF -L $BED -o $OUT.fgbio.bqsr.gatkHC.vcf
$GATK -T HaplotypeCaller -I $OUT.preMID.bqsrCal.bam -R $REF -L $BED -o $OUT.preMID.bqsr.gatkHC.vcf

#filter for qual
#java -Xmx64g -jar $anno/SnpSift.jar filter "(QUAL>=30)" $OUT.fgbio.lf.mergeIndels.vcf \
 #   > $OUT.fgbio.filter.vcf

#java -Xmx64g -jar $anno/SnpSift.jar annotate $cosmic $OUT.preMID.filter.vcf > $OUT.preMID.filter.cosmic.vcf
java -Xmx64g -jar $anno/SnpSift.jar annotate $cosmic $OUT.fgbio.bqsr.lf.mergeIndels.vcf \
    > $OUT.fgbio.lf.cosmic.vcf
java -Xmx64g -jar $anno/SnpSift.jar annotate $cosmic $OUT.fgbio.bqsr.gatkHC.vcf \
    > $OUT.fgbio.gatkHC.cosmic.vcf

#java -Xmx64g -jar $anno/SnpSift.jar annotate $clinvar $OUT.preMID.filter.cosmic.vcf > $OUT.preMID.filter.cosmic.clinvar.vcf
java -Xmx64g -jar $anno/SnpSift.jar annotate $clinvar $OUT.fgbio.lf.cosmic.vcf \
    > $OUT.fgbio.lf.cosmic.clinvar.vcf
java -Xmx64g -jar $anno/SnpSift.jar annotate $clinvar $OUT.fgbio.gatkHC.cosmic.vcf \
    > $OUT.fgbio.gatkHC.cosmic.clinvar.vcf

# dbNSFP annotation for polyphen and pathogenic score
#java -Xmx64g -jar $anno/SnpSift.jar dbNSFP -db $dbnsfp $OUT.preMID.filter.cosmic.clinvar.vcf > $OUT.preMID.filter.cosmic.clinvar.dbnsfp.vcf
java -Xmx64g -jar $anno/SnpSift.jar dbNSFP -db $dbnsfp $OUT.fgbio.lf.cosmic.clinvar.vcf \
    > $OUT.fgbio.lf.cosmic.clinvar.dbnsfp.vcf
java -Xmx64g -jar $anno/SnpSift.jar dbNSFP -db $dbnsfp $OUT.fgbio.gatkHC.cosmic.clinvar.vcf \
    > $OUT.fgbio.gatkHC.cosmic.clinvar.dbnsfp.vcf

#final vars from clinvar, cosmic &dbnsfp

java -Xmx64g -jar $anno/SnpSift.jar extractFields $OUT.fgbio.lf.cosmic.clinvar.dbnsfp.vcf \
 CHROM POS ID REF ALT QUAL DP AF SB DP4 "CLNDN" "CLNSIG" "CLNVI" "CLNVC" "GENEINFO" "ORIGIN" "RS" "MC" "dbNSFP_LRT_pred" \
    "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_MutationTaster_pred" "AF_ESP" "AF_TGP" > $OUT.fgbio.lf.finalvars.txt

java -Xmx64g -jar $anno/SnpSift.jar extractFields $OUT.fgbio.gatkHC.cosmic.clinvar.dbnsfp.vcf \
 CHROM POS ID REF ALT QUAL DP AF SB DP4 "CLNDN" "CLNSIG" "CLNVI" "CLNVC" "GENEINFO" "ORIGIN" "RS" "MC" "dbNSFP_LRT_pred" \
    "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_MutationTaster_pred" "AF_ESP" "AF_TGP" > $OUT.fgbio.gatkHC.finalvars.txt

done

######  clean up ######
rm *.MID.fq.gz
#rm unpaired*gz
rm *var.vcf
rm *cosmic.vcf
#rm *sus.bam
#rm *MID.bam
rm *.sanitise*.bam
#rm *notused4consensus.bam
#rm *.realigned.ba[mi]
rm *table

#############################################################
##########   summarize coverage results and reporting ########
for f in *covd
do
awk '{sum+=$7}END{m=(sum/NR); b=m*0.2; print m, b}' $f \
    > $f_covd.tmp
awk 'BEGIN{n=0}NR==FNR{m=$1;b=$2;next}{if($7>=b)n++}END{print m,(n/FNR*100.0)}' \
 OFS="\t" $f_covd.tmp $f \
    > ${f%.covd}_covMetrics.txt
done
#awk 'BEGIN{n=0}NR==FNR{m=$1;b=$2;next}{if($6>=b)n++}END{print m,b,(n/FNR*100.0)}'
# Summarize PCR metrics
# new PICARD rs2
for f in *targetPCRmetrics.txt
do
    awk -v n=${f%%.target*} 'NR==1{print n,$5,$11,$14*100,$22*100.0,($17+$18)/$7*100.0}' \
        OFS="\t" <(head -n8 $f | tail -n1) \
        > ${f%%.txt}_summary.txt
    f2=${f%%.target*}_covMetrics.txt
    paste ${f%%.txt}_summary.txt $f2 > ${f2%%_cov*}_combined_cov_metrics.txt

done

echo "SampleID    Total_Reads    #UQ_Reads_Aligned    %UQ_Reads_Aligned    %Bases_OnTarget_Aligned     %Bases_OnTarget_Total    Mean_Coverage    %Coverage_Uniformity" \
    > final_metrics_report.txt
cat *_combined_cov_metrics.txt >> final_metrics_report.txt

rm *_summary.txt
rm *_covMetrics.txt
rm *cov_metrics.txt
rm *fp.vcf
rm *mic.vcf *var.vcf
