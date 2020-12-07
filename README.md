# Metagenomic-project
Songjun's master degree project about functional metagenomic 

## Description
This project is a workflow to process the raw metagenomic dataset for two annotation servers: MG-RAST and GhostKOALA, and integrate the final abundance number into KEGG pathwy.

The workflow include two parts:

The first part is a snakemake pipeline, which performs quality assessment and filtering of raw fastq files to generate filtered fasta files and amino acid files that are ready for annotation servers.

Example: 
<img src="./example_workflow.svg">

Second part consists of a python script and a R shiny app, which generates annotated results from GhostKOALA in KEGG website and visualizes them in KEGG pathways by using the shiny app based on PathView package.

usage: 

```
./group.py -m metadata.txt -o output.csv
```

where metadata.txt contains the experimental setup (see example below). The result will be printed to  output.csv 

A metadata file is also necessary for the python script. It is used for grouping results from annotation server. A example file is attached in metadata.txt.

** Example metadata file:**

top	bottom	root
11.txt	12.txt	13.txt
14.txt	15.txt	16.txt
17.txt	18.txt	19.txt
20.txt	21.txt	22.txt
23.txt	24.txt	25.txt
26.txt	27.txt	28.txt
29.txt	30.txt	31.txt
32.txt	33.txt	34.txt
35.txt	36.txt	37.txt

### Dependencies
Required software and version:
* Python/3.7.3, R/4.0.1, snakemake/5.26.1, Fastqc/0.11.8, Mulitqc/1.8, PEAR/0.9.1, seqtk/1.0, BLAST/2.5.0+, SPAdes/3.9.1, bwa/0.7.17




