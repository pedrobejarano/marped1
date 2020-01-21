
## Author: Mario Ruiz Velazquez and Pedro Bejarano Diaz
## Contact: marioruizvelazquez@gmail.com and pedro.bejarano@jerez.es
##Date: November 2019

#$ -S /bin/bash
#$ -N input_sample_processing
#$ -V
#$ -cwd
#$ -j yes
#$ -o input_sample_processing

#! /bin/bash

##Reading parameters

INPUT=$1
WD=$2
NUMINPUT=$3
NUMSAM=$4
PROMOTER=$5
OUTPUT=$6
RSCRIPT_1=$7
RSCRIPT_2=$8
SAMPLEDIR=$9

## Access input folder

cd $WD/samples/input/input${INPUT}

##Quality control  and mapping

if [ -e input${INPUT}_2.fastq ]
   then
     fastqc input${INPUT}_1.fastq
     fastqc input${INPUT}_2.fastq

     bowtie2 -x $WD/genome/index -1 input${INPUT}_1.fastq -2 input${INPUT}_2.fastq -S input${INPUT}.sam

   else
     fastqc input${INPUT}.fastq

     bowtie2 -x $WD/genome/index -U input${INPUT}.fastq -S input${INPUT}.sam
fi


## Transcript assembly

cd $WD/samples/input/input$INPUT

samtools view -S -b input$INPUT.sam > input$INPUT.bam
rm input$INPUT.sam
samtools sort input$INPUT.bam -o input_sorted_${INPUT}.bam
rm input$INPUT.bam
samtools index input_sorted_${INPUT}.bam

## Sinch point through blackboard

echo "input$INPUT DONE" >> $WD/logs/blackboard.txt

DONE_INPUT=$(wc -l $WD/logs/blackboard.txt | awk '{ print$1 }')

echo "Input samples have been processed"


if [ $DONE_INPUT -eq $NUMSAM ]
then
   qsub -N callpeak -o $WD/logs/callpeak $WD/../calling_peaks.sh $WD $NUMINPUT $PROMOTER $OUTPUT $RSCRIPT_1 $RSCRIPT_2 $SAMPLEDIR
   echo "submitting calling_peaks file"
fi

