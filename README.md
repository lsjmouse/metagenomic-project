# metagenomic-project
# Songjun's master degree project about functional metagenomic 

# raw data are fast.gz files
# each sample has two files: R1_001 and R2_001, which represents forward and reverse reads.

# The workflow include two steps
# First part is a snakemake pipeline, which translate raw fastq files into filtered fasta files and faa files that are ready for annotation servers
# Specific steps including Pear, tBlastx, seqtk, metaSpades and several custom python scripts.


# Second part is a python script and a R shiny app, which generates annotated results from GhostKOALA in KEGG website and visualizes them in kegg pathway by using pathview package
# A metadata file is also necessary for the python script. It is used for grouping results from annotation server. A exapmle file is attached in this github.
# usage: ./group.py metadata.txt output.csv
