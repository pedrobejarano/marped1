## Author: Mario Ruiz Velazquez and Pedro Bejarano Diaz
## Contact: marioruizvelazquez@gmail.com and pedro.bejarano@jerez.es
##Date: November 2019

#$ -S /bin/bash
#$ -N chip_sample_processing
#$ -V
#$ -cwd
#$ -j yes
#$ -o chip_sample_processing

#! /bin/bash

## Reading parameters

CHIP=$1
WD=$2
NUMCHIP=$3
NUMSAM=$4
PROMOTER=$5
OUTPUT=$6
RSCRIPT_1=$7
RSCRIPT_2=$8
SAMPLEDIR=$9

##Access chip folder

cd $WD/samples/chip/chip$CHIP

##Quality control and mapping

if [ -e chip${CHIP}_2.fastq ]
   then
     fastqc chip${CHIP}_1.fastq
     fastqc chip${CHIP}_2.fastq

     bowtie2 -x $WD/genome/index -1 chip${CHIP}_1.fastq -2 chip${CHIP}_2.fastq -S chip${CHIP}.sam
  else
     fastqc chip${CHIP}.fastq

     bowtie2 -x $WD/genome/index -U chip${CHIP}.fastq -S chip${CHIP}.sam
fi

## Transcript assembly

cd $WD/samples/chip/chip${CHIP}

samtools view -S -b chip${CHIP}.sam > chip${CHIP}.bam
rm chip${CHIP}.sam
samtools sort chip${CHIP}.bam -o chip_sorted_${CHIP}.bam
rm chip$CHIP.bam
samtools index chip_sorted_${CHIP}.bam


## Sincronisation poit trough blackboard

echo "chip${SAM_ID} DONE" >> $WD/logs/blackboard.txt

DONE_CHIP=$(wc -l $WD/logs/blackboard.txt | awk '{ print$1  }')

echo "Chip samples have been processed"

if [ $DONE_CHIP -eq $NUMSAM ]
then
   qsub -N callpeak -o $WD/logs/callpeak $WD/../calling_peaks.sh $WD $NUMCHIP $PROMOTER $OUTPUT $RSCRIPT_1 $RSCRIPT_2 $SAMPLEDIR
fi

