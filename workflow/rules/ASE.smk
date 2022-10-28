##########################
#### Read group add ######
##########################

#functions 

def get_file_names(wildcards):
    ck_output = checkpoints.split_bam.get(**wildcards).output[0]
    global SMP
    SMP, = glob_wildcards(os.path.join(ck_output, "{sample}.bam"))
    return expand("results/{tech}/data_by_clusters/{SAMPLE}.bam", SAMPLE=SMP,tech=wildcards.tec)

#modify bam header-> add read group needed for ASEReadCounter
rule GATK_AddorRep:
    input:
        "results/{tec}/alignment/{tec}.positionsort.bam"
    output:
        "results/{tec}/alignment/gatkgroup.sorted_in.bam"
    conda: "../envs/gatk.yml"
    shell:
        """ gatk AddOrReplaceReadGroups -I {input} -O {output} -RGLB DNA -RGPL ILLUMINA -RGPU {wildcards.tec} -RGSM {wildcards.tec}_sample -VALIDATION_STRINGENCY SILENT """

rule corr_barcode:
    input:
        bam="results/{tec}/alignment/gatkgroup.sorted_in.bam"
    output:
        "results/{tec}/alignment/gatkgroup.sorted.bam"
    conda: "../envs/samtools.yml"
    shell:
        """ if [ {wildcards.tec} == gex ]
        then 
            samtools view -h {input.bam} | awk -F"_" '{{if ($1 ~ /^@/) print $0; else print $2":"$0;}}' > barcode_corr.sam
            samtools view barcode_corr.sam -o {input.bam} -@{threads}
            rm barcode_corr.sam
        fi 
        mv {input} {output}
        samtools index {output}"""

checkpoint split_bam:
    input: 
        cluster_file="results/{tec}/features/cluster_{tec}.tsv",
        bam="results/{tec}/alignment/gatkgroup.sorted.bam"
    conda:
        "../envs/sinto.yml"
    output:
        directory("results/{tec}/data_by_clusters/")
    threads: config["threads_num"]
    conda: "../envs/sinto.yml"
    shell:
        """ sinto filterbarcodes -b {input.bam} -p {threads} -c {input.cluster_file} --outdir {output} --barcode_regex "[^:]*" """

##########################
#### ASE GATK ############
##########################

#Allelic specific expression: gatk ASEReadCounter over the phased vcf for each cluster bam
rule ASEReadCount_cluster_sc:
    input:
        vcf= "results/merged_vcf/chr{chrom}QC.vcf.gz",
        bam = get_file_names
    output:
        directory("results/{tec}/ASE{chrom}")
    conda: "../envs/ASE.yml"
    params:
        fa=config["genome_fa"]
    shell:
        """ mkdir -p "results/{wildcards.tec}/ASE{wildcards.chrom}"
        for C in {input.bam}; do 
        samtools index $C
        gatk ASEReadCounter -R {params.fa}  -I $C -V {input.vcf}  -O {output}/$(basename "$C" | cut -d. -f1).table -L {wildcards.chrom}; done"""

rule index_bulk_bam:
    input:
        "results/gex_bulk/alignment/gatkgroup.sorted.bam"
    output:
        "results/gex_bulk/alignment/gatkgroup.sorted.bam.bai"
    conda: "../envs/samtools.yml"
    shell:
        """ samtools index {input}"""

#Allelic specific expression: gatk ASEReadCounter over the phased vcf for each cluster bam
rule ASEReadCount_cluster_bulk:
    input:
        vcf= "results/merged_vcf/chr{chrom}QC.vcf.gz",
        bam = "results/gex_bulk/alignment/gatkgroup.sorted.bam",
        bai = "results/gex_bulk/alignment/gatkgroup.sorted.bam.bai",
        fa=config["genome_fa"]
    output:
        directory("results/gex_bulk/ASE{chrom}_bulk")
    conda: "../envs/gatk.yml"
    params:
        o = "results/gex_bulk/ASE{chrom}_bulk"
    shell:
        """ mkdir -p {params.o}
        gatk ASEReadCounter -R {input.fa}  -I {input.bam} -V {input.vcf}  -O {output}/gatkgroup.sorted.table -L {wildcards.chrom}"""
