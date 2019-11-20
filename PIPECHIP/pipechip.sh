## This pipeline analysis Chip-seq data

## Author: Mario Ruiz Velazquez and Pedro Bejarano Diaz
##Contact: marioruizvelazquez@gmail.com
## Date: November 2019

#! /bin/bash

if [ $# -eq 0 ]
   then
    echo "This pipeline analysisi Chip Seq data."
    echo "Usage: piperna <param_file>"
    echo ""
    echo "param_file: File with the parameters spefification. Please, check params.sh for an example"
    echo ""
    echo "Enjoy!"

    exit 0
fi


## Parameters loading

PARAMS=$1

WD=$(grep working_directory: $PARAMS | awk '{ print$2 }')
NUMSAM=$(grep number_of_samples: $PARAMS | awk '{ print$2 }')
GENOME=$(grep genome: $PARAMS | awk '{ print$2 }')
ANNOTATION=$(grep annotation: $PARAMS | awk '{ print$2 }')
NUMCHIP=$(grep chip_num: $PARAMS | awk '{ print$2 }')
NUMINPUT=$(grep input_num: $PARAMS | awk '{ print$2 }')

SAMPLES_CHIP=( )
I=0

while [ $I -lt $NUMCHIP ]
do
   SAMPLES_CHIP[$I]=$(grep sra_chip$(($I+1)): $PARAMS | awk '{ print$2 }')
((I++))
done

SAMPLES_INPUT=( )
I=0

while [ $I -lt $NUMINPUT ]
do
   SAMPLES_INPUT[$I]=$(grep sra_input$(($I+1)): $PARAMS | awk '{ print$2 }')
   ((I++))
done


##Printing variable values

echo WORKING_DIRECTORY=$WD
echo NUMBER_OF_SAMPLES=$NUMSAM
echo GENOME=$GENOME
echo ANNOTATION=$ANNOTATION
echo NUMBER_CHIP_SAMPLES=$NUMCHIP
echo NUMBER_INPUT_SAMPLES=$NUMINPUT

I=0

while [ $I -lt $NUMCHIP ]
do
   echo chip_$((I + 1)) = ${SAMPLES_CHIP[$I]}
   ((I++))
done

I=0

while [ $I -lt $NUMINPUT ]
do
   echo input_$((I + 1)) = ${SAMPLES_INPUT[$I]}
   ((I++))
done


##Generate the working directory

mkdir $WD
cd $WD
mkdir genome annotation samples results logs

cd samples
mkdir chip input
cd chip

I=1

while [ $I -le $NUMCHIP ]
do
   mkdir chip$I
   ((I++))
done

cd ../input

I=1

while [ $I -le $NUMINPUT ]
do
     mkdir input$I
     ((I++))
done

## Download genome of reference

cd $HOME

cd $WD/genome
wget -O genome.fa.gz $GENOME
gunzip genome.fa.gz


## Download annotation

cd $HOME

cd $WD/annotation
wget -O annotation.gtf.gz $ANNOTATION
gunzip annotation.gtf.gz

## Building reference index

cd $WD/genome
bowtie2-build genome.fa index

## Download samples

cd $HOME

cd $WD/samples/chip

I=0

while [ $I -lt $NUMCHIP ]
do
   cd chip$((I+1))
   fastq-dump --split-files ${SAMPLES_CHIP[$I]}
   cd $WD/samples/chip
   sleep 1m ##Wait 1 minute
   ((I++))
done

cd $WD/samples/input
I=0

while [ $I -lt $NUMINPUT ]
do
   cd input$((I+1))
   fastq-dump --split-files ${SAMPLES_INPUT[$I]}
   cd $WD/samples/input
   ((I++))
done

## Chip processing

I=1
while [ $I -le $NUMCHIP ]
do
   qsub -N chip$I -o $WD/logs/chip$I /home/marped/Desktop/tareas/marped1/PIPECHIP/chip_sample_processing.sh $I $WD ${NUMCHIP} ${SAMPLES_CHIP[(($I-1))]}
   ((I++))

done

## Input processing

I=1
while [ $I -le $NUMINPUT ]
do
   qsub -N input$I -o $WD/logs/input$I /home/Desktop/tareas/marped1/PIPECHIP/input_sample_processing.sh $I $WD ${NUMINPUT} ${SAMPLES_INPUT[(($I-1))]}
   ((I++))
done

