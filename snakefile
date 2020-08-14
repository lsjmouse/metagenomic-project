sample_names, = glob_wildcards("Data/Sequences_fastq.gz/{sample}_R1_001.fastq.gz")

'''
rule all:
    input: expand("{sample}.fasta", id=sample_names)
'''

rule FastQC_for_rawdata:
    input:
        r1 = "Data/Sequences_fastq.gz/{sample}_R1_001.fastq.gz"
        r2 = "Data/Sequences_fastq.gz/{sample}_R2_001.fastq.gz"
    shell:
        """
        fastqc {input} -t 6 -o Analysis/1.FastQC/
        """

rule multiqc:
    shell:
        """
        multiqc Analysis/1.FastQC/
        """

rule Pear:
    input:
        r1 = "Data/Sequences_fastq.gz/{sample}_R1_001.fastq.gz"
        r2 = "Data/Sequences_fastq.gz/{sample}_R2_001.fastq.gz"
    output:
        outfile = "Analysis/2.Pear/{sample}.pear.fastq.gz"
    shell:
        """
        pear -f {input.r1} -r {input.r2} -o {output.outfile}
        """

rule mkdir_for_pear_output:
    shell:
        """
        mkdir Analysis/2.Pear/assembled_fastq Analysis/2.Pear/unassembled_fastq Analysis/2.Pear/discarded_fastq
        """

rule mv_pear_output:
    input:
        assembled_fastq = "Analysis/2.Pear/{sample}.pear.fastq.gz.assembled.fastq"
        unassembled_r1 = "Analysis/2.Pear/{sample}.pear.fastq.gz.unassembled.forward.fastq"
        unassembled_r2 = "Analysis/2.Pear/{sample}.pear.fastq.gz.unassembled.reverse.fastq"
        discarded_fastq = "Analysis/2.Pear/{sample}.pear.fastq.gz.discarded.fastq"
    shell:
        """
        mv {input.assembled_fastq} Analysis/2.Pear/assembled_fastq;
        mv {input.unassembled_r1} Analysis/2.Pear/unassembled_fastq;
        mv {input.unassembled_r2} Analysis/2.Pear/unassembled_fastq;
        mv {input.discarded_fastq} Analysis/2.Pear/discarded_fastq
        """

rule count_pear_percent:
    input:
        assembled_fastq = "Analysis/2.Pear/assembled_fastq/{sample}.pear.fastq.gz.assembled.fastq"
        unassembled_r1 = "Analysis/2.Pear/unassembled_fastq/{sample}.pear.fastq.gz.unassembled.forward.fastq"
        discarded_fastq = "Analysis/2.Pear/discarded_fastq/{sample}.pear.fastq.gz.discarded.fastq"
    output:
        count_report = "Analysis/2.Pear/count_report"
    shell:
        """
        a=$(< "{input.assembled_fastq}" wc -l); b=$(< "{input.unassembled_r1}" wc -l); c=$(< "{input.discarded_fastq}" wc -l); d=$(($a + $b + $c)); e=$(bc<<<"scale=2; $a/$d"); f=$(bc<<<"scale=2; $b/$d");
        echo ${sample} assembled percent = $e ',' unassembled percent = $f >> {output.count_report}
        """

rule mkdir_for_fasta:
    shell:
        """
        mkdir Analysis/2.Pear/assembled_fasta
        """

rule translate_assembled_fastq_to_fasta:
    input:
        assembled_fastq = "Analysis/2.Pear/assembled_fastq/{sample}.pear.fastq.gz.assembled.fastq"
    output:
        assembled_fasta = "Analysis/2.Pear/assembled_fasta/{sample}.pear.fastq.gz.assembled.fasta"
    shell:
        """
        seqtk seq -A {input.assembled_fastq} > {output.assembled_fasta}
        """

# the CH4_fasta is dowanloaded from KEGG website and created by customized R script

rule make_blast_database:
    input:
        CH4_fasta = "Data/database/CH4_database.fasta"

    shell:
        """
        makeblastdb -in {input.CH4_fasta} -input_type fasta -dbtype nucl -out Data/database
        """

rule mkdir_for_aseembled_tblastx_result:
    shell:
        """
        mkdir Analysis/3.tblastx/assembled_sequence_results
        """


# for now the evalue set is 1e-05, it can be changed due to different demands

rule run_tblastx_for_assembled_fasta:
    input:
        assembled_fasta = "Analysis/2.Pear/assembled_fasta/{sample}.pear.fastq.gz.assembled.fasta"
    output:
        assembled_tblastx_result = "Analysis/3.tblastx/assembled_sequence_results/{sample}.pear.fastq.gz.assembled.fasta.tblastx"
    shell:
        """
        tblastx -query {input.assembled_fasta} -out {output.assembled_tblastx_result} -db Data/database/CH4_database -outfmt 6 -evalue 1e-5
        """

rule mkdir_for_brdiging_sequence:
    shell:
        """
        mkdir Analysis/2.Pear/bridging_sequence
        """

rule add_50_N_to_unassembled_sequence:
    input:
        unassembled_r1 = "Analysis/2.Pear/unassembled_fastq/{sample}.pear.fastq.gz.unassembled.forward.fastq"
        unassembled_r2 = "Analysis/2.Pear/unassembled_fastq/{sample}.pear.fastq.gz.unassembled.reverse.fastq"
    output:
        combiend_sequence = "Analysis/2.Pear/bridging_sequence/{sample}.pear.fastq.gz.unassembled.combine.fasta"
    shell:
        """
        ./Script/pear_add_n.py {input.unassembled_r1} {input.unassembled_r2} {output.combiend_sequence}
        """

rule mkdir_for_unaseembled_tblastx_result:
    shell:
        """
        mkdir Analysis/3.tblastx/unassembled_sequence_results
        """

rule run_tblastx_for_unassembled_fasta:
    input:
        combiend_sequence = "Analysis/2.Pear/bridging_sequence/{sample}.pear.fastq.gz.unassembled.combine.fasta"
    output:
        unassembled_tblastx_result = "Analysis/3.tblastx/unassembled_sequence_results/{sample}.pear.fastq.gz.unassembled.combine.fasta.tblastx"
    shell:
        """
        tblastx -query {input.combiend_sequence} -out {output.unassembled_tblastx_result} -db Data/database/CH4_database -outfmt 6 -evalue 1e-5
        """

rule mkdir_for_filtered_sequence:
    shell:
        """
        mkdir Analysis/4.assembly/filtered_id_list Analysis/4.assembly/filtered_fastq
        """

rule get_hitted_sequence_id_list:
    input:
        assembled_tblastx_result = "Analysis/3.tblastx/assembled_sequence_results/{sample}.pear.fastq.gz.assembled.fasta.tblastx"
        unassembled_tblastx_result = "Analysis/3.tblastx/unassembled_sequence_results/{sample}.pear.fastq.gz.unassembled.combine.fasta.tblastx"
    output:
        id_list = "Analysis/4.assembly/filtered_id_list/{sample}.filtered.list"
    shell:
        """
        cat {input.assembled_tblastx_result}|cut -f 1|uniq >>{output.id_list};
        cat {input.unassembled_tblastx_result}|cut -f 1|uniq >>{output.id_list}
        """

rule get_hitted_sequence_fastq:
    input:
        id_list = "Analysis/4.assembly/filtered_id_list/{sample}.filtered.list"
        r1 = "Data/Sequences_fastq.gz/{sample}_R1_001.fastq.gz"
        r2 = "Data/Sequences_fastq.gz/{sample}_R2_001.fastq.gz"
    output:
        hitted_forward_fastq = "Analysis/4.assembly/filtered_fastq/{sample}_R1_001.filtered.fastq"
        hitted_reverse_fastq = "Analysis/4.assembly/filtered_fastq/{sample}_R2_001.filtered.fastq"
    shell:
        """
        seqtk subseq {input.r1} {input.id_list} >{output.hitted_forward_fastq};
        seqtk subseq {input.r2} {input.id_list} >{output.hitted_reverse_fastq}
        """

rule metaspade:
    input:
        hitted_forward_fastq = "Analysis/4.assembly/filtered_fastq/{sample}_R1_001.filtered.fastq"
        hitted_reverse_fastq = "Analysis/4.assembly/filtered_fastq/{sample}_R2_001.filtered.fastq"
    output:
        "Analysis/4.assembly/filtered_assembly/{sample}"
    shell:
        """
        metaspade.py -1 {input.hitted_forward_fastq} -2 {input.hitted_reverse_fastq} -o {output}
        """

rule mkdir_for_mapping:
    shell:
        """
        mkdir Analysis/5.mapping
        """

#use software:bwa for mapping, need to get it prepared in the directory

rule mapping:
    input:
        scaffolds = "Analysis/4.assembly/filtered_assembly/{sample}/scaffolds.fasta"
        hitted_forward_fastq = "Analysis/4.assembly/filtered_fastq/{sample}_R1_001.filtered.fastq"
        hitted_reverse_fastq = "Analysis/4.assembly/filtered_fastq/{sample}_R2_001.filtered.fastq"
    output:
        "Analysis/5.mapping/{sample}.sam"
    shell:
        """
        ./bwa index {input.scaffolds};
        ./bwa mem {input.scaffolds} {input.hitted_forward_fastq} {input.hitted_reverse_fastq} >{output}
        """

rule mkdir_for_annotation:
    shell:
        """
        mkdir Analysis/6.annotation/pear_result Analysis/6.annotation/filtered_fasta
        """

rule get_hitted_sequence_fasta:
    input:
        hitted_forward_fastq = "Analysis/4.assembly/filtered_fastq/{sample}_R1_001.filtered.fastq"
        hitted_reverse_fastq = "Analysis/4.assembly/filtered_fastq/{sample}_R2_001.filtered.fastq"
    output:
        "Analysis/6.annotation/filtered_fasta/{sample}.filtered.fasta"
    shell:
        """
        pear -f {input.hitted_forward_fastq} -r {input.hitted_reverse_fastq} -o Analysis/6.annotation/pear_result;
        seqtk seq -A Analysis/6.annotation/pear_result/{sample}.assembled.fastq >> {output};
        ./Script/pear_add_n.py Analysis/6.annotation/pear_result/{sample}.unassembled.forward.fastq Analysis/6.annotation/pear_result/{sample}.unassembled.reverse.fastq Analysis/6.annotation/pear_result/{sample}.unassembled.combined.fasta;
        cat Analysis/6.annotation/pear_result/{sample}.unassembled.combined.fasta >> {output}
        """
#now we get filtered fasta files for annotation
#for MG-RAST, these filtered fasta files can be submitted directly, just need to prepare a metadata file manually
#for GhostKOALA, these filtered fasta files need to be translated into amino acid sequence(faa file)

rule get_amino_acid_sequence:
    input:
        "Analysis/6.annotation/filtered_fasta/{sample}.filtered.fasta"
    output:
        "Analysis/6.annotation/filtered_fasta/{sample}.filtered.fasta.faa"
    shell:
        """
        ./Script/translate.py {input} {output}
        """
