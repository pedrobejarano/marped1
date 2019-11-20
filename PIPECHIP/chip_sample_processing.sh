## Author: Mario Ruiz Velazquez and Pedro Bejarano Diaz
## Contact: marioruizvelazquez@gmail.com
##Date: November 2019

#$ -S /bin/bash
#$ -N chip_sample_processing
#$ -V
#$ -cwd
#$ -j yes
#$ -o chip_sample_processing


## Reading parameters

CHIP_ID=$1
WD=$2
NUMCHIP=$3
SRA=$4

##Access chip folder

cd $WD/samples/chip${CHIP_ID}

##Quality control and mapping

if [ -e ${SRA}_2.fastq ]
   then
     fastqc ${SRA}_1.fastq
     fastqc ${SRA}_2.fastq

     bowtie2 -x ../../genome/index -1 ${SRA}_1.fastq -2 ${SRA}_2.fastq -S chip${CHIP_ID}.sam

   else
     fastqc ${SRA}_1.fastq

     bowtie2 -x ../../genome/index -U ${SRA}_1.fastq -S chip${CHIP_ID}.sam
fi
## Transcript assembly

samtools view -S -b chip${CHIP_ID}.sam > chip${CHIP_ID}.bam
rm chip${CHIP_ID}.sam
samtools sort chip${CHIP_ID}.bam -o chip_sorted_${CHIP_ID}.bam
rm chip${CHIP_ID}.bam
samtools index chip_sorted_${CHIP_ID}.bam


## Synchronisation point through blackboard

echo "chip${SAM_ID} DONE" >> $WD/logs/blackboard

DONE_SAMPLES=$(wc -l $WD/logs/blackboard)

I=0
if [ $DONE_SAMPLES -eq $NUMCHIP ]
then
   qsub -N input_sample_processing -o $WD/logs/input_sample_processing input_sample_processing.sh PARAMETROS
   ((I++))
fi

