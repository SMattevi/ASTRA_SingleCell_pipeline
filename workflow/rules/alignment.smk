##########################
#### alignment ###########
##########################

#Alignment performed with hisat using the "--no-spliced-alignment" and the paird end option for ATAC
rule alignment_atac:
    input:
        R1=expand("results_{sample_id}/atac/barcodecorrection/{sample}_R1.fastq.gz", sample=config["lanes"],sep=",",sample_id=config["sample_name"]),
        R2=expand("results_{sample_id}/atac/barcodecorrection/{sample}_R2.fastq.gz", sample=config["lanes"],sep=",",sample_id=config["sample_name"])
    threads: 
        config["threads_num"]
    params: 
        index_gen = config["hisat_index"],
    output:
        temp("results_{sample_id}/atac/alignment/atac.sam")
    log: "results_{sample_id}/logs/atac_hisat.log"
    conda:
        "../envs/hisat.yml"
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
        "results_{sample_id}/gex/umitools_extr/merged_R2_extracted.fastq.gz",
    output:
        temp("results_{sample_id}/gex/alignment/gex.sam")
    threads: 
        config["threads_num"]
    params: 
        config["hisat_index"]
    log: "results_{sample_id}/logs/gex_hisat.log"
    conda:
        "../envs/hisat.yml"
    shell:
        """ hisat2 -x {params} \
            -5 30 \
            -U {input} \
            -S {output} \
            -p {threads} \
            --summary-file {log}"""

#Alignment of pseudo-bulk 
rule alignment_bulk_gex:
    input:
        "results_{sample_id}/gex/fastq/merged_R2.fastq.gz",
    output:
        temp("results_{sample_id}/gex_bulk/alignment/gex_bulk.sam")
    threads: 
        config["threads_num"]
    params: 
        config["hisat_index"]
    log: "results_{sample_id}/logs/bulkgex_hisat.log"
    conda:
        "../envs/hisat.yml"
    shell:
        """ hisat2 -x {params} \
            -5 30 \
            -U {input} \
            -S {output} \
            -p {threads} \
            --summary-file {log}"""

rule alignment_exome:
    input:
        R1=expand("{path}/{sample}_1.fastq.gz", path=config["path_exome"],sample=config["fastqs_exome"],sep=","),
        R2=expand("{path}/{sample}_2.fastq.gz", path=config["path_exome"],sample=config["fastqs_exome"],sep=",")
    threads: 
        config["threads_num"]
    params: 
        index_gen = config["hisat_index"]
    output:
        temp("results_{sample_id}/exome/alignment/exome.sam")
    log: "results_{sample_id}/logs/exome.log"
    conda:
        "../envs/hisat.yml"
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
            --summary-file {log} """
