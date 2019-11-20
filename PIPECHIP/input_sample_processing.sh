## Author: Mario Ruiz Velazquez and Pedro Bejarano Diaz
## Contact: marioruizvelazquez@gmail.com
##Date: November 2019

#$ -S /bin/bash
#$ -N input_sample_processing
#$ -V
#$ -cwd
#$ -j yes
#$ -o input_sample_processing

##Reading input parameters

INPUT=$1
WD=$2
NUMCHIP=$3
SRA=$4

## Access input folder

cd $WD/samples/input${INPUT}

##Quality control and mapping

if [ -e {$SRA}_2.fastq ]
   then
     fastqc {$SRA}_1.fastq
     fastqc {$SRA}_2.fastq

     bowtie2 -x ../../genome/index -1 {$SRA}_1.fastq -2 {$SRA}_2.fastq -S input${INPUT}.sam

   else
     fastqc {$SRA}_1.fastq

     bowtie2 -x ../../genome/index -U {$SRA}_1.fastq -S input${INPUT}.sam
fi
## Transcript assembly

samtools view -S -b chip${INPUT}.sam > chip${INPUT}.bam
rm chip${INPUT}.sam
samtools sort chip${INPUT}.bam -o chip_sorted_${INPUT}.bam
rm chip${INPUT}.bam
samtools index chip_sorted_${INPUT}.bam

