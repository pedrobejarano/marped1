## This pipeline analysis Chip-seq data

## Author: Mario Ruiz Velazquez and Pedro Bejarano Diaz
##Contact: marioruizvelazquez@gmail.com or pedro.bejarano@jerez.es
## Date: November 2019

#! /bin/bash

if [ $# -eq 0 ]
   then
    echo "This pipeline analyses Chip Seq data."
    echo "Usage: pipechip <param_file>"
    echo ""
    echo "param_file: File with the parameters specification. Please, check params.sh in the test folder for an example"
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
PROMOTER=$(grep promoter: $PARAMS | awk '{ print$2  }')
RSCRIPT_1=$(grep rscript_1: $PARAMS | awk '{ print$2  }')
RSCRIPT_2=$(grep rscript_2: $PARAMS | awk '{ print$2  }')
OUTPUT=$(grep output: $PARAMS | awk '{ print$2  }')
SAMPLEDIR=$(grep sample_dir: $PARAMS | awk '{ print$2  }')

SAMPLES_CHIP=( )
I=0

while [ $I -lt $NUMCHIP ]
do
   SAMPLES_CHIP[$I]=$(grep chip_$(($I+1)): $PARAMS | awk '{ print$2 }')
   ((I++))
done

SAMPLES_INPUT=( )
I=0

while [ $I -lt $NUMINPUT ]
do
   SAMPLES_INPUT[$I]=$(grep input_$(($I+1)): $PARAMS | awk '{ print$2 }')
   ((I++))
done


##Printing variable values

echo "printing loaded variables"


echo WORKING_DIRECTORY=$WD
echo NUMBER_OF_SAMPLES=$NUMSAM
echo GENOME=$GENOME
echo ANNOTATION=$ANNOTATION
echo NUMBER_CHIP_SAMPLES=$NUMCHIP
echo NUMBER_INPUT_SAMPLES=$NUMINPUT
echo PROMOTER=$PROMOTER
echo OUTPUT=$OUTPUT
echo RSCRIPT_1=$RSCRIPT_1
echo RSCRIPT_2=$RSCRIPT_2
echo SAMPLEDIR=$SAMPLEDIR

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


##Generating the working directory

echo "Generating the needed folders in the" $WD

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

## Downloading reference genome

while true; do
     read -p "Is your reference genome already downloaded?" yn
     case $yn in

        [Yy]* ) echo "Genome must be saved in the folder $WD/test and unzipped"

        ## Copying the genome.fa file

        cd $WD/genome
        cp $GENOME genome.fa; break;;

        [Nn]* ) echo "The reference genome will be downloaded from ensemble"

        cd $WD/genome
        wget -O genome.fa.gz $GENOME
        gunzip genome.fa.gz; break;;
     esac
done

## Downloading annotation

while true; do
     read -p "Is your genome annotation already downloaded?" yn
     case $yn in

        [Yy]* ) echo "Annotation must be saved in the folder $WD/test and unzipped"

        ## Copying the genome.fa file

        cd $WD/annotation
        cp $GENOME annotation.gtf; break;;

        [Nn]* ) echo "The reference genome will be downloaded from ensemble"

        cd $WD/annotation
        wget -O annotation.gtf.gz $GENOME
        gunzip annotation.gtf.gz; break;;
     esac
done

echo "The genome and the annotation are already downloaded and saved in the genome and annotation folder"

## Building reference index

echo "Building the reference index"

cd $WD/genome
bowtie2-build genome.fa index

echo "Index built"

## Download/copy samples into the $WD/samples and name changing

while true; do
     read -p "Are your samples already downloaded?" yn
     case $yn in

        [Yy]* ) echo "Samples must be at $WD/../test"

        ## Copying the CHIP samples

        cd $WD/samples/chip
        I=0

        while [ $I -lt $NUMCHIP ]
        do
        cd $WD/samples/chip
        cp ${SAMPLES_CHIP[$I]} chip$(($I+1))/chip$(($I+1)).fastq
        cd chip$(($I+1))
        ((I++))
        done

        ## Copying the INPUT samples

        cd $WD/samples/input
        I=0

        while [ $I -lt $NUMINPUT ]
        do
        cd $WD/samples/input
        cp ${SAMPLES_INPUT[$I]} input$(($I+1))/input$(($I+1)).fastq
        cd input$(($I+1))
        ((I++))
        done; break;;

        [Nn]* ) echo "Please, make sure you added the right SRR in the params file. Your samples will be downloaded from the ncbi"

        ##Downloading CHIP samples
        cd $WD/samples/chip

        I=0
        J=1
        K=1
        while [ $I -lt $NUMCHIP ]
        do
        cd chip$((I+1))
        fastq-dump --split-files ${SAMPLES_CHIP[$I]}
        if [ -e ${SAMPLES_CHIP[$I]}_2.fastq ]
        then
           mv ${SAMPLES_CHIP[$I]}_2.fastq chip${J}_2.fastq
           mv ${SAMPLES_CHIP[$I]}_1.fastq chip${J}_1.fastq

        else
           mv ${SAMPLES_CHIP[$I]}_1.fastq chip$K.fastq
        fi
        ((I++))
        ((J++))
        ((K++))

        cd $WD/samples/chip
        sleep 15s ##Wait 15 seconds
        done

 ##Downloading INPUT samples

        cd $WD/samples/input

        I=0
        J=1
        K=1
        while [ $I -lt $NUMINPUT ]
        do
        cd input$((I+1))
        fastq-dump --split-files ${SAMPLES_INPUT[$I]}
        if [ -e ${SAMPLES_INPUT[$I]}_2.fastq ]
        then
           mv ${SAMPLES_INPUT[$I]}_2.fastq input${J}_2.fastq
           mv ${SAMPLES_INPUT[$I]}_1.fastq input${J}_1.fastq
        else
           mv ${SAMPLES_INPUT[$I]}_1.fastq input$K.fastq
        fi
        ((I++))
        ((J++))
        ((K++))
        cd $WD/samples/input
        sleep 30s ##Wait 30s
        done; break;;
     esac
done

echo "Your chip and input samples are already downloaded"

## Chip processing

I=1

while [ $I -le $NUMCHIP ]
do
   qsub -N chip$I -o $WD/logs/chip$I $WD/../chip_sample_processing.sh $I $WD $NUMCHIP $NUMSAM $PROMOTER $OUTPUT $RSCRIPT_1 $RSCRIPT_2 $SAMPLEDIR
   ((I++))

done

## Input processing

I=1

while [ $I -le $NUMINPUT ]
   do
   qsub -N input$I -o $WD/logs/input$I $WD/../input_sample_processing.sh $I $WD $NUMINPUT $NUMSAM $PROMOTER $OUTPUT $RSCRIPT_1 $RSCRIPT_2 $SAMPLEDIR
   ((I++))
done

