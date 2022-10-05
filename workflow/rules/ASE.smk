##########################
#### Read group add ######
##########################

#functions 

def get_file_names(wildcards):
    ck_output = checkpoints.split_bam.get(**wildcards).output[0]
    global SMP
    SMP, = glob_wildcards(os.path.join(ck_output, "cluster_{sample}.bam"))
    return expand("results/{tech}/data_by_clusters/cluster_{SAMPLE}.bam", SAMPLE=SMP,tech=wildcards.tec)

#modify bam header-> add read group needed for ASEReadCounter
rule GATK_AddorRep:
    input:
        "results/{tec}/alignment/{tec}.positionsort.bam"
    output:
        "results/{tec}/alignment/gatkgroup.sorted.bam"
    conda: "../envs/gatk.yml"
    shell:
        """ gatk AddOrReplaceReadGroups -I {input} -O {output} -RGLB DNA -RGPL ILLUMINA -RGPU {wildcards.tec} -RGSM {wildcards.tec}_sample -VALIDATION_STRINGENCY SILENT """

#split bam by clusters            
checkpoint split_bam:
    input: 
        cluster_file="results/{tec}/features/cluster_{tec}.tsv",
        bam="results/{tec}/alignment/gatkgroup.sorted.bam"
    conda:
        "../envs/samtools.yml"
    output:
        directory("results/{tec}/data_by_clusters/")
    shell:
        """ mkdir -p {output}
        if [ {wildcards.tec} == gex ]
        then 
            samtools view -h {input.bam} | awk -F"_" '{{if ($1 ~ /^@/) print $0; else print $2":"$0;}}' > barcode_corr.sam
            samtools view barcode_corr.sam -o {input.bam} -@{threads}
            rm barcode_corr.sam
        fi
        perl  workflow/scripts/split_bam2clusters.pl --cluster_file {input.cluster_file} --bam_file {input.bam} --output_dir {output} 
        for i in {output}/*.bam; do samtools index $i -@{threads}; done """

##########################
#### ASE GATK ############
##########################

#Allelic specific expression: gatk ASEReadCounter over the phased vcf for each cluster bam
rule ASEReadCount_cluster_sc:
    input:
        vcf= "results/merged_vcf/chr{chrom}QC.vcf.gz",
        bam = get_file_names,
        fa=config["genome_fa"]
    output:
        directory("results/{tec}/ASE{chrom}")
    conda: "../envs/gatk.yml"
    params:
        o = "results/{tec}/ASE{chrom}"
    shell:
        """ mkdir -p {params.o}
        for C in {input.bam}; do gatk ASEReadCounter -R {input.fa}  -I $C -V {input.vcf}  -O {output}/$(basename "$C" | cut -d. -f1).table -L {wildcards.chrom}; done"""

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
