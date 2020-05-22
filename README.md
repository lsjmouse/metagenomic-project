# metagenomic-project
# Songjun's master degree project about functional metagenomic 

# raw data are fast.gz files
# each sample has two files: R1_001 and R2_001, which represents forward and reverse reads.

# step1 Fastqc and Multiqc
# run Fastqc for each files and use multiqc to generate them into one html file
# version: FastQC/0.11.8  MultiQC/1.8
for i in *fastq.gz; do fastqc $i -t 6 -o ../../Analysis/FastQC/ ;done
cd /proj/snic2019-30-58/Songjun/Analysis/FastQC
multiqc .

# step2 Pear
# use pear to match overlapped sequences
# version: pear/0.9.10
for f in *R1_001.fastq.gz; do "file is $f"; p=${f%_*_*}; pear -f ${p}_R1_001.fastq.gz -r ${p}_R2_001.fastq.gz -o ../../Analysis/Pear/${p}.pear.fastq.gz; done;
# the output for each sample has 4 files, which are assembled.fastq, discard.fastq, unassembled.forward.fastq and unassembled.reverse.fastq

# step3 organize output files from Pear
cd Pear
mkdir assembled_fastq
mkdir unassembled_fastq
mkdir discarded_fastq
mv *.assembled.fastq assembled_fastq/
mv *.unassembled.forward.fastq unassembled_fastq/
mv *.unassembled.reverse.fastq unassembled_fastq/
mv *.discarded.fastq discarded_fastq

# step4 count assembled percentage and generate a report
cd assembled_fastq
for i in *assembled.fastq; do  p=${i%.assembled.fastq};
a=$(< "$i" wc -l); b=$(< "../unassembled_fastq/${p}.unassembled.forward.fastq" wc -l); 
c=$(< "../discarded_fastq/${p}.discarded.fastq" wc -l); d=$(($a + $b + $c)); e=$(bc<<<"scale=2; $a/$d"); f=$(bc<<<"scale=2; $b/$d") ;
echo $p assembled percent = $e ',' unassembled percent = $f >> count_report;  done  


# step5 use seqtk to transfer fastq into fasta
# make the assembled sequence prepared for tblastx
# version: seqtk/1.2-r101
cd assembled_fastq
for i in *fastq; do seqtk seq -A $i > ${i}.fasta


# step6 bridging unassembled sequence together
# use my customized python script to bridging unassembled forward and reverse sequences, with 50N between them.
# since these bridging sequences are used for blastx, so the output file is fasta file.
cd unassembled_fastq
for f in *forward.fastq; do p=${f%.forward.fastq}; ./pear_add_n.py ${p}.forward.fastq ${p}.reverse.fastq ${p}.combine.fasta ;done
mv *.combine.fasta briding_sequence/

# step7 makeblastdb
# the file "CH4_database.fasta" is the relevant function sequence downloaded from KEGG website.
# use this fasta file as the blast database for tblastx
# version blast/2.9.0+
makeblastdb -in CH4_database.fasta -input_type fasta -dbtype nuca -out Data/database

# step8 tblastx
# run tblastx for both assembled sequences that generated from Pear and bridging sequences that combined by customized python.
# the evalue here is setted as 1e-5 for this project
# version blast/2.9.0+
cd assembled_fastq
for i in *.fasta; do tblastx -query $i -out Analysis/tblastx/assembled_sequence_results/${i}.tblastx -db Data/database/CH4_database -outfmt 6 -evalue 1e-5
cd bridging_sequence
for i in *.fasta; do tblastx -query $i -out Analysis/tblastx/unassembled_sequence_results/${i}.tblastx -db Data/database/CH4_database -outfmt 6 -evalue 1e-5

# step9 get hit sequences for tblastx results
# pick up the sequences id that have hit in tblastx, and use these ids to find their original fastq sequence
# the sequences get from this step are the filtered sequences by tblastx, which can be used for downstream analysis
cd unassembled_sequence_results
for i in *.tblastx; do p=${i%.pear.fastq.gz.unassembled.combine.fasta.tblastx}; cat $i|cut -f 1|uniq >>../filter_list/${p}.filtered.list;done
cd assembled_sequence_results
for i in *.tblastx; do p=${i%.pear.fastq.gz.assembled.fastq.fasta.tblastx}; cat $i|cut -f 1| uniq >> ../filter_list/${p}.filtered.list; done
cd filter_list
for i in *.list; do p=${i%.filtered.list}; seqtk subseq ../../Data/Sequences_fastq.gz/${p}_R1_001.fastq.gz $i >../filtered_fastq/${p}_R1_001.filtered.fastq;
seqtk subseq ../../Data/Sequences_fastq.gz/${p}_R2_001.fastq.gz $i >../filtered_fastq/${p}_R2_001.filtered.fastq; done 


# step10 meta assembly
# use metaspades to make assembly for these filtered fastq
# version: SPAdes v3.9.1 [metaSPAdes mode]
cd filtered_fastq
for i in *R1_001.filtered.fastq; do p=${i%_R1_001.filtered.fastq}; metaspades.py -1 ${p}_R1_001.filtered.fastq -2 ${p}_R2_001.filtered.fastq -o ../filtered_assembly/$p; done
