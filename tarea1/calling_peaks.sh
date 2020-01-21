## Author: Pedro J. Bejarano Diaz and Mario Ruiz Velazquez
## Date: December 2019
## Contact: pedro.bejarano@jerez.es and marioruizvelazquez@gmail.com

#$ -S /bin/bash
#$ -N calling_peaks
#$ -V
#$ -cwd
#$ -j yes
#$ -o calling_peaks

#! /bin/bash

## Loading parameters

WD=$1
NC=$2
PROMOTER=$3
OUTPUT=$4
RSCRIPT_1=$5
RSCRIPT_2=$6
SAMPLEDIR=$7

#Callpeak function

cd $WD/results

## Macs2 for creating the peakAnnotation file
echo "Finding the peaks"

I=1

while [ $I -le $NC ]
do
   macs2 callpeak -t $WD/samples/chip/chip$I/chip_sorted_${I}.bam -c $WD/samples/input/input$I/input_sorted_${I}.bam -n 'peaks_'$I --outdir . -f BAM
   ((I++))
done

## Rscript for the processing of peaks

if [ $NC -eq 2 ]
then
        cd $WD/results
        cp $SAMPLEDIR/$RSCRIPT_2 $WD/results
        Rscript $RSCRIPT_2 peaks_1_peaks.narrowPeak peaks_2_peaks.narrowPeak $PROMOTER $WD/results
elif [ $NC -eq 1 ]
then
        cd $WD/results
        cp $SAMPLEDIR/$RSCRIPT_1 $WD/results
        Rscript $RSCRIPT_1 peaks_1_peaks.narrowPeak $PROMOTER $WD/results
fi

echo "The plots and the target genes will be saved in the folder results"

## HOMER for finding motifs

cd $WD/results

I=1

while [ $I -le $NC ]
do
   findMotifsGenome.pl peaks_${I}_summits.bed tair10 ./HOMER_$I -size 65
   ((I++))
done

echo "DONE"


