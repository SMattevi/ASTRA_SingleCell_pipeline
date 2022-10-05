##########################
#### alignment ###########
##########################

#Alignment performed with hisat using the "--no-spliced-alignment" and the paird end option for ATAC
rule alignment_atac:
    input:
        R1=expand("results/atac/barcodecorrection/{sample}_R1.fastq.gz", sample=config["lanes"],sep=","),
        R2=expand("results/atac/barcodecorrection/{sample}_R2.fastq.gz", sample=config["lanes"],sep=",")
    threads: 
        config["threads_num"]
    params: 
        index_gen = config["hisat_index"],
    output:
        temp("results/atac/alignment/atac.sam")
    log: "logs/atac_hisat.log"
    conda:
        "envs/hisat.yml"
    shell:
        """ R1=$(echo {input.R1})
            R1new=$(echo $R1 | sed 's/ /,/g ')
            R2=$(echo {input.R2})
            R2new=$(echo $R2 | sed 's/ /,/g ')
            hisat2 -x {params.index_gen} \
            -1 $R1new \
            -2 $R2new \
            -S {output} \
            --no-spliced-alignment \
            -p {threads} \
            --summary-file {log}"""

#Alignment performed with hisat for RNAseq
rule alignment_gex:
    input:
        "results/gex/umitools_extr/merged_R2_extracted.fastq.gz",
    output:
        temp("results/gex/alignment/gex.sam")
    threads: 
        config["threads_num"]
    params: 
        config["hisat_index"]
    log: "logs/gex_hisat.log"
    conda:
        "envs/hisat.yml"
    shell:
        """ hisat2 -x {params} \
            -U {input} \
            -S {output} \
            -p {threads} \
            --summary-file {log}"""

#Alignment of pseudo-bulk 
rule alignment_bulk_gex:
    input:
        "results/gex/fastq/merged_R2.fastq.gz",
    output:
        temp("results/gex_bulk/alignment/gex_bulk.sam")
    threads: 
        config["threads_num"]
    params: 
        config["hisat_index"]
    log: "logs/bulkgex_hisat.log"
    conda:
        "envs/hisat.yml"
    shell:
        """ hisat2 -x {params} \
            -U {input} \
            -S {output} \
            -p {threads} \
            --summary-file {log}"""
