## Script de R para obtener los genes dianas de un factor de transcripción. Para ello
## se usa el fichero narrowPeak.

## Autor: Mario Ruiz Velázquez y Pedro Jesús Bejarano Díaz
## Fecha: Diciembre 2019

## Loading arguments

args <- commandArgs(trailingOnly = TRUE)

input.file.name <- args[[1]]
promoter.length <- as.numeric(args[[2]])
directory <- args[[3]]

## Setting working directory

setwd(dir = directory)

library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28

library(clusterProfiler)
library(org.At.tair.db)

library("pathview")


## Peak analysis

## Reading peak file

peaks <- readPeakFile(peakfile = input.file.name,header=FALSE)
head(peaks)

## Defining the region that is considered as a promoter.
## Normaly de region contaions a 1000 pb upstream and downstream the TSS

promoter <- getPromoters(TxDb=txdb,
                         upstream=-promoter.length,
                         downstream=promoter.length)


## Checking the number of genes from de A.Thaliana genome. It should have 33602 genes

genes <- as.data.frame(genes(txdb))
genes_names <- genes$gene_id
length(genes_names)


## Annotating peaks

peakanno <- annotatePeak(peak = peaks,
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb)

## Binding sites in specific regions of the genome

plotAnnoPie(peakanno)
vennpie(peakanno)

plotAnnoBar(peakanno)

## Distribution of genomic loci relative to TSS
plotDistToTSS(peakanno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")


## Converting annotation to data frame and writing a table with target genes

annotation_dataframe <- as.data.frame(peakanno)
target_genes <- annotation_dataframe$geneId[annotation_dataframe$annotation == "Promoter"]

write(x = target_genes, file="target_genes.txt")

## GO TERMS ENRICHMENT 

##Reading the target genes

gene.set <- read.table(file = "target_genes.txt", header = F, as.is = T)[[1]]
length(gene.set)

## GO Enrichment analysis of a gene set

ego <- enrichGO(gene = gene.set, OrgDb = org.At.tair.db, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.01, universe = genes_names, keyType = "TAIR")
ego.res <- as.data.frame(ego)
head(ego.res)

dotplot(ego)
barplot(ego, showCategory=10 )
emapplot(ego, showCategory = 10)

## KEGG ENRICHMENT 


kk <- enrichKEGG(gene = gene.set, organism = "ath", universe = genes_names)
kk
kk.res <- as.data.frame(kk)
head(kk.res)

kk.ID<-kk.res$ID

my.universe <- rep(0,length(genes_names))
names(my.universe) <- genes_names
my.universe[gene.set] <- 1
my.universe


pathway.res <- pathview(gene.data = my.universe, pathway.id = kk.ID, species = "ath", gene.idtype = "TAIR")

