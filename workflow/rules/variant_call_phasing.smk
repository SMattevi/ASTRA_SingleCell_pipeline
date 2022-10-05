##########################
#### variant calling ####
##########################
            
#variant calling performed with strelka using the "GermlineWorkflow" and the "--exome" option
rule variant_calling:
    input:
        bam="results/{tec}/alignment/{tec}.positionsort.bam",
        fa=config["genome_fa"]
    output: 
        "results/{tec}/variant_calling/strelka/results/variants/variants.vcf.gz"
    params:
        config["strelka_path"]
    shell:
        """ {params}/bin/configureStrelkaGermlineWorkflow.py \
        --bam {input.bam} \
        --referenceFasta {input.fa} \
        --runDir results/{wildcards.tec}/variant_calling/strelka \
        --exome

        results/{wildcards.tec}/variant_calling/strelka/runWorkflow.py  -m local -j 15 --quiet """

#QC: extract only het variants, with AF and DP as specified in config file
rule QC_VCF:
    input:
        "results/{tec}/variant_calling/strelka/results/variants/variants.vcf.gz"
    output:
        final="results/{tec}/variant_calling/strelka/results/variants/variantsQC.vcf.gz",
	    initial="results/{tec}/variant_calling/strelka/results/variants/variantsPASS.vcf.gz",
	    vcftbi="results/{tec}/variant_calling/strelka/results/variants/variantsQC.vcf.gz.tbi"
    conda:
        "envs/samtools.yml"
    params:
        AF=config["AF"],
        DP=config["DP"]
    shell:
        """ bcftools view {input} -f PASS -O z -o {output.initial}
	    bcftools view {output.initial} -i 'GT=="het" && AF>{params.AF} && DP>{params.DP}' -m2 -M2 -O z -o {output.final}
	    tabix {output.final} """

rule prephasing_WHATSHAP:
    input:
        vcf=expand("results/{tec}/variant_calling/strelka/results/variants/variantsQC.vcf.gz",tec=config["tech"]),
        bam=expand("results/{tec}/alignment/{tec}.positionsort.bam",tec=config["tech"])
    output:
        v="results/{tec}/prephasing/pre_phased.vcf.gz",
        t="results/{tec}/prephasing/pre_phased.vcf.gz.tbi"
    params: 
        config["genome_fa"]
    conda:
        "envs/whatshap.yml"
    shell:
        """ mkdir -p results/{wildcards.tec}/prephasing
        whatshap phase -o results/{wildcards.tec}/prephasing/pre_phased.vcf --reference={params} {input.vcf} {input.bam} --ignore-read-groups
        bgzip results/{wildcards.tec}/prephasing/pre_phased.vcf
        tabix {output.v} """

rule merge_tec_vcf:
    input:
        pp=expand("results/{tec}/variant_calling/strelka/results/variants/variantsQC.vcf.gz",tec=config["tech"]),
    output:
        v="results/merged_vcf/variantsQC.vcf.gz",
        t="results/merged_vcf/variantsQC.vcf.gz.tbi",
    conda:
        "envs/samtools.yml"
    shell:
        """ bcftools concat {input.pp} -a | bcftools norm -d all -O z -o {output.v} 
        tabix {output.v} """

rule merge_tec_vcf_prephased:
    input:
        np=expand("results/{tec}/prephasing/pre_phased.vcf.gz",tec=config["tech"])
    output:
        vp="results/phased/pre_phased.vcf.gz",
        tp="results/phased/pre_phased.vcf.gz.tbi",
    conda:
        "envs/samtools.yml"
    shell:
        """ bcftools concat {input.np} -a | bcftools norm -d all -O z -o {output.vp} 
        tabix {output.vp} """

rule divide_chr_VCF:
    input:
        "results/merged_vcf/variantsQC.vcf.gz"
    output:
        cv="results/merged_vcf/chr{chrom}QC.vcf.gz",
        ct="results/merged_vcf/chr{chrom}QC.vcf.gz.tbi",
    conda:
        "envs/samtools.yml"
    shell:
        """ bcftools view {input} --regions {wildcards.chrom} -O z -o {output.cv}
        tabix {output.cv} """

##########################
######## phasing ########
##########################
        
#phase the quality controlled called SNPs with shapeit4
rule phasing_SHAPEIT4:
    input:
        vcftbi="results/phased/pre_phased.vcf.gz.tbi",
        vcf="results/phased/pre_phased.vcf.gz",
        ref_vcf=expand("{path}/{prefix}{chrom}.{extension}", path=config["path_ALLvcf"],prefix=config["prefix_ALLvcf"],extension=config["extension_ALLvcf"],  allow_missing=True)
    output: 
        "results/phased/chr{chrom}_phased.vcf"
    conda:
        "envs/shapeit.yml"
    threads: config["threads_num"]
    params: 
        sex=config["sex"]
    shell:
        """  if [ {wildcards.chrom} == "X" ]
        then
            if [ {params.sex} == "male" ]
            then
                shapeit4 --input {input.vcf} --region X:10000-2781479,X:155701382-156030895 --output {output} --thread {threads} --reference {input.ref_vcf}
            else
                shapeit4 --input {input.vcf} --region X --output {output} --thread {threads} --reference {input.ref_vcf}
            fi
        else 
            shapeit4 --input {input.vcf} --region {wildcards.chrom} --output {output} --thread {threads} --reference {input.ref_vcf}
        fi """

rule bgzip_and_indexing:
    input:
        "results/phased/chr{chrom}_phased.vcf"
    output: 
        vcf="results/phased/chr{chrom}_phased.vcf.gz",
        vcftbi="results/phased/chr{chrom}_phased.vcf.gz.tbi"
    conda:
        "envs/samtools.yml"
    shell:
        """ bgzip {input} 
            tabix {output.vcf} """
 
