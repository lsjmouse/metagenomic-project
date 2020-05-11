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
        mv {input.assembled_fastq} Analysis/2.Pear/assembled_fastq
        mv {input.unassembled_r1} Analysis/2.Pear/unassembled_fastq
        mv {input.unassembled_r2} Analysis/2.Pear/unassembled_fastq
        mv {input.discarded_fastq} Analysis/2.Pear/discarded_fastq
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

# also how to get this CH4_fasta, maybe ask joel about the R script
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
        unassembled_tblastx_result = "Analysis/3.tblastx/assembled_sequence_results/{sample}.pear.fastq.gz.unassembled.combine.fasta.tblastx"
    shell:
        """
        tblastx -query {input.combiend_sequence} -out {output.unassembled_tblastx_result} -db Data/database/CH4_database -outfmt 6 -evalue 1e-5
        """
