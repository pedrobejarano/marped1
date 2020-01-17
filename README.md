This scripts are done by Pedro J. Bejarano-Diaz and Mario Ruiz-Velazquez for that users that require to perform a ChIP-seq analysis (identify the binding site of transcriptions factors in a concrete genome) with any kind of organism and number of samples. It may help you to reach your conclusion in less time due to their automatization.

Take time to read this README. Here you will find the information to premorí the analysis successfully.


*NECESSARY TOOLS*


- R (https://cran.r-project.org/)
- Bowtie 2
- Fastqc
- Samtools
- MACS2
- HOMER and genome data (tair10)
- ChIPseeker (R package)
- Annotation package of the organism of study (TxDb.Athaliana.BioMart.plantsmart28 by default)
- clusterProfiler (R package)
- .db data of the organism of study (org.At.tair.db by default)
- path view (R package)


*SAMPLE REQUIREMENTS*


- If your experiment has the same number of chip samples than input samples, you can use this program by default. For different conditions, you have to modify some parameters of call_peaks.sh.

- This program is configured for Arabidopsis thaliana samples by default, but you can modify some parameters to use other organisms. In that case, you will have to change the annotation package, the .db data of peak_processiong_onerep.R and peak_processing_tworep.R, and the genome data for Arabidopsis thaliana (tair10) in findMotifsGenome.pl function for HOMER.


- If your experiment has 1 or 2 replicas, you can use this program by default. If your experiment have more than 2 replicas, modify peak_processing_onerep.R and peak_processing_tworep.R following the next instructions:

  1º Add both args and input file you have, adding new arguments. For example, if you have another .narrowPeak, add:

  input.file.name3 <- args[[4]].

  2º In the section called "reading peaks", read as many files as input.file.name. For example, for input_file_name3, name another peaks3.

  3º Add the new peaks in your function intersect.

 In this case, you will also need to change the script "execute_scriptR.sh" adding the new args that you have or deleting the args that you do not use. You will need to change pipechip.sh and you will need to add new parameters in params.txt.


*HOW WORKS*


This program has 6 bash scripts and 2 R scripts.

The instruction to execute is: “bash pipe chip.sh params.txt”.

	This instruction will execute pipe chip.sh, that create the working directory and copy/download the genome file, the annotation file, and the samples files. Then, this script will launch the next scripts to a High Performance Computing (HPC) with “qsub” instruction. The script pipe chip.sh will launch one script for each sample, due to samples processing is a parallel task.

Once executed, chip_processing.sh and input_processing.sh will be executed for the sample processing. This two scripts, with have the same instructions to process the chip samples and the input samples, respectively, will do a quality study using fastqc and it will map the short sequences with reference genome using bowtie2. This program CAN PROCESS paired data and unpaired data. Next, the script will generate the .bam file using samtools and it will use a synchronization point through blackboards to launch one call_peaks.sh (next script) for each replica in the right moment.

Next script is call_peaks.sh. One call_peaks will be launched for each replicas (for each chip and input samples). This script use MACS2 to create a peakAnnotation file though sorted .bam files, that are generated with chip_processing.sh and input_processing.sh. Also, this script will find motifs in the genome using HOMER tool. This process can take a while, be patient. Then, if the number of replicas of your experiment is 1, call_peaks.sh will launch execute_Rscript_onerep.sh. If the number of replicas of your experiment is 2, call_peaks.sh will launch execute_Rscript_tworep.sh. In other case (number of replicas is higher than 2), you must edit the scripts like we tell at *sample requirements*.
