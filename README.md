Metagenomic-project: Songjun's master degree project about functional metagenomic 

raw data are fast.gz files,
each sample has two files: R1_001 and R2_001, which represents forward and reverse reads.

Required software: snakemake/5.26.1, Fastqc/0.11.8, Mulitqc/1.8, PEAR/0.9.1, seqtk/1.0, BLAST/2.5.0+, SPAdes/3.9.1, bwa/0.7.17

The workflow include two steps:
First part is a snakemake pipeline, which translate raw fastq files into filtered fasta files and faa files that are ready for annotation servers.
Specific steps can be visulized in example_workflow.svg.

Second part is a python script and a R shiny app, which generates annotated results from GhostKOALA in KEGG website and visualizes them in kegg pathway by using pathview package.
A metadata file is also necessary for the python script. It is used for grouping results from annotation server. A exapmle file is attached in metadata.txt.
usage: ./group.py metadata.txt output.csv
