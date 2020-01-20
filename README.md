This scripts are done by Pedro J. Bejarano-Diaz and Mario Ruiz-Velazquez for that users that require to perform a ChIP-seq analysis (identify the binding site of transcriptions factors in a concrete genome) with any kind of organism and number of samples. It may help you to reach your conclusion in less time due to their automatization.

Take time to read this README. Here you will find the information to premorí the analysis successfully.


### NECESSARY TOOLS


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


### SAMPLE REQUIREMENTS


- If your experiment has the same number of chip samples than input samples, you can use this program by default. For different conditions, you have to modify some parameters of call_peaks.sh.

- This program is configured for Arabidopsis thaliana samples by default, but you can modify some parameters to use other organisms. In that case, you will have to change the annotation package, the .db data of peak_processiong_onerep.R and peak_processing_tworep.R, and the genome data for Arabidopsis thaliana (tair10) in findMotifsGenome.pl function for HOMER.


- If your experiment has 1 or 2 replicas, you can use this program by default. If your experiment have more than 2 replicas, modify peak_processing_onerep.R and peak_processing_tworep.R following the next instructions:

  1º Add both args and input file you have, adding new arguments. For example, if you have another .narrowPeak, add:

  input.file.name3 <- args[[4]].

  2º In the section called "reading peaks", read as many files as input.file.name. For example, for input_file_name3, name another peaks3.

  3º Add the new peaks in your function intersect.

 In this case, you will also need to change the script "execute_scriptR.sh" adding the new args that you have or deleting the args that you do not use. You will need to change pipechip.sh and you will need to add new parameters in params.txt.


### USAGE AND PARAMS FILE


First of all, you must to create a parameters file in order to create those with they may find necessary.

   -‘Working_directory’: This parameter create the directory where the analysis will be done. We suggest to create a folder called PIPECHIP where you can paste the scripts and a folder called TEST where params file can be added.

   -‘Number_of_samples’: This parameter refers to the total number of samples (both chips and inputs).

   -‘Genome’: There might be 2 possible options depending on having already downloaded or not the reference genome:

      -If you have already download your genome, you must paste it into the TEST folder inside PIPECHIP. Then, the param genome will be the path to the reference genome inside the TEST folder.

      -If you haven’t download your genome yet, you must paste in the genome param, the link to the reference genome. 

   -‘annotation’: There might be 2 possible options depending on having already download or not the annotation:

      -If you have already download the annotation, you must paste it into the TEST folder inside PIPECHIP. Then, the param genome will be the path to the annotation inside the TEST folder.

      -If you haven’t download the annotation yet, you must pase in the genome param the link to the reference genome. 

   -‘chip_num’: This parameter refers to the number of chip samples.

   -‘input_num’: This parameter refers to the number of input samples.

Regarding the samples params, there might be 2 possible options:

      -If you have already download the samples, you should paste them as well into the TEST folder and add in this parameter the path to this folder.

‘chip_1’: paste the path to the chip 1 sample.
‘chip_2’: paste the path to the chip 2 sample.
‘input_1’: paste the path to the input 1 sample.
‘input_2’: paste the path to the input 2 sample.

      -If you haven’t download the samples yet, you must paste here the SRR accession number from NCBI.
‘chip_1’: add the right SRR.
‘chip_2’: add the right SRR.
‘input_1’: add the right SRR.
‘input_2’: add the right SRR.

In the case that you have more than 2 chip and 2 input samples, you must add extra parameter following the previous order and make sure you write them as chip_x and input_x, x referring to the sample number.

   -‘sample_dir’: paste the path to the refers to the TEST folder where you have pasted the params file along with the genome, annotation and samples, given the case you have them already downloaded. 

   -‘promoter’: This parameter refers to the length of the promoter for the processing of picks obtained using macs2 function.

   -‘output’: This parameter refers to the directory where you want to save the peak processing file.



### HOW WORKS



This program have 4 scripts. In this section we describe each script so as you know how the analysis work.

   -‘pipechip.sh’. This is the main script, and will launch the others. It create the working directory with the samples folders and subfolders. After that it will prompt a question whether you have or not already downloaded the genome. If yes, it will copy it from the TEST folder to the correct one inside the working directory. If not, it will download it from the link given in the params file. Then, the script will prompt another question on whether you have or not already download the annotation. If yes, it will copy it from the TEST folder into the correct one inside the working directory. If not, it will download it from the link given in the params file.Once having downloaded both the genome and annotation it will create the index. After that the script will prompt a third message asking if the samples are downloaded or not. If yes, it will copy them from the TEST folder into the correct one, taking into account chip and inputs samples (__you shall make sure you have added the input and chip samples correctly in the params file__). If not, it will download them using the SRR given in the params file. Once this is done, the script will launch the chip and input processing scripts.
   

   -‘chip_sample_processing.sh’. This script processes the samples. The .sam, .bam, sorted.bam and sorted.bam.bai files of the chip samples will be created.

   -‘input_sample_processing.sh’. See ‘chip_sample_processing.sh’ but this processes the input samples.

   -‘calling_peaks.sh’. This scripts in launched once the others are done.



### ANALYSIS



This pipeline performs a ChIP-Seq analysis along with a functional analysis of the peaks obtained in order to reveal which molecular functions are related to these peaks and a motif finding to identify the family to which the transcription factor (TF) belongs.

Once the peaks are obtained, they will be saved as peaks_n_peak.narrowPeak, among other files. These narrowPeak files will provide you the vital information for the analysis. In order to extract all the information from it, the script `peak_analysis.R` will be launched to the queue. Please, make sure depending on the number of samples, one script or another will be launched (if you have more than 3 samples, thus 3 different .narrowPeak files, check the comments on the Rscript).

To clarify to which family does the TF belong, a HOMER analysis using the findMotifsGenome.pl will be perform. It will perform as many analysis as summit.bed files there are. Each analysis will be saved into a folder named HOMER_x (x as the number of analysis depending on the number of peaks) inside the RESULTS folder.



## EXAMPLE GIVEN



For further information about this tool you may like to perform the analysis using the example samples. You may find them in the NCBI accession number __GSE115358__. The genome and the annotation used were obtained from ensembl plants.
