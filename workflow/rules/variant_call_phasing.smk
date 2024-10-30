##########################
#### variant calling ####
##########################
            
#variant calling performed with strelka using the "GermlineWorkflow" and the "--exome" option
rule variant_calling:
    input:
        bam="results_{sample_id}/{tec}/alignment/{tec}.positionsort.bam",
        fa=config["genome_fa"]
    output: 
        "results_{sample_id}/{tec}/variant_calling/strelka/results/variants/variants.vcf.gz"
    params:
        config["strelka_path"]
    conda:
        "../envs/py2.yml"
    shell:
        """ {params}/bin/configureStrelkaGermlineWorkflow.py \
        --bam {input.bam} \
        --referenceFasta {input.fa} \
        --runDir results_{wildcards.sample_id}/{wildcards.tec}/variant_calling/strelka \
        --exome

        results_{wildcards.sample_id}/{wildcards.tec}/variant_calling/strelka/runWorkflow.py  -m local -j 15 --quiet """

#QC: extract only het variants, with AF and DP as specified in config file
rule QC_VCF:
    input:
        "results_{sample_id}/{tec}/variant_calling/strelka/results/variants/variants.vcf.gz"
    output:
        final="results_{sample_id}/{tec}/variant_calling/strelka/results/variants/variantsQC.vcf.gz",
        initial="results_{sample_id}/{tec}/variant_calling/strelka/results/variants/variantsPASS.vcf.gz",
        inter=temp("results_{sample_id}/{tec}/variant_calling/strelka/results/variants/variantsrename.vcf.gz"),
        vcftbi="results_{sample_id}/{tec}/variant_calling/strelka/results/variants/variantsQC.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    params:
        AF=config["AF"],
        DP=config["DP"],
        samp=config["sample_name"],
        samplefile="results_{sample_id}/sample.txt"
    shell:
        """ echo "SAMPLE1 {params.samp}"> {params.samplefile}
        bcftools view {input} -f PASS -O z -o {output.initial}
	    bcftools view {output.initial} -i 'GT=="het" && AF>{params.AF} && DP>{params.DP}' -m2 -M2 -O z -o {output.inter}
	    bcftools reheader -s {params.samplefile} {output.inter} -o {output.final}
        tabix {output.final} """

rule prephasing_WHATSHAP:
    input:
        vcf=myinput,
        bam=expand("results_{sample_id}/{tec}/alignment/{tec}.positionsort.bam",tec=config["tech"],sample_id=config["sample_name"])
    output:
        v="results_{sample_id}/{tec}/prephasing/pre_phased.vcf.gz",
        t="results_{sample_id}/{tec}/prephasing/pre_phased.vcf.gz.tbi"
    params: 
        config["genome_fa"]
    conda:
        "../envs/whatshap.yml"
    shell:
        """ mkdir -p results_{wildcards.sample_id}/{wildcards.tec}/prephasing
        whatshap phase -o results_{wildcards.sample_id}/{wildcards.tec}/prephasing/pre_phased.vcf --reference={params} {input.vcf} {input.bam} --ignore-read-groups
        bgzip results_{wildcards.sample_id}/{wildcards.tec}/prephasing/pre_phased.vcf
        tabix {output.v} """

rule merge_tec_vcf:
    input:
        pp=myinput
    output:
        v="results_{sample_id}/merged_vcf/variantsQC.vcf.gz",
        t="results_{sample_id}/merged_vcf/variantsQC.vcf.gz.tbi",
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input.pp} -a | bcftools norm -d all -O z -o {output.v} 
        tabix {output.v} """

rule merge_tec_vcf_prephased:
    input:
        np=expand("results_{sample_id}/{tec}/prephasing/pre_phased.vcf.gz",tec=config["tech"],sample_id=config["sample_name"])
    output:
        vp="results_{sample_id}/phased/pre_phased.vcf.gz",
        tp="results_{sample_id}/phased/pre_phased.vcf.gz.tbi",
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input.np} -a | bcftools norm -d all -O z -o {output.vp} 
        tabix {output.vp} """

# rule divide_chr_VCF:
#     input:
#         "results_{sample_id}/merged_vcf/variantsQC.vcf.gz"
#     output:
#         cv="results_{sample_id}/merged_vcf/{chrom}QC.vcf.gz",
#         ct="results_{sample_id}/merged_vcf/{chrom}QC.vcf.gz.tbi",
#     conda:
#         "../envs/samtools.yml"
#     shell:
#         """ bcftools view {input} --regions {wildcards.chrom} -O z -o {output.cv}
#         tabix {output.cv} """

##########################
######## exome ########
##########################

rule HaplotypeCaller_e:
    input:
        "results_{sample_id}/exome/recalibration/exome.recal.bam"
    output:
        initial=temp("results_{sample_id}/exome/haplotypeCaller/exome.g.vcf.gz"),
        final="results_{sample_id}/exome/haplotypeCaller/exome.vcf.gz"
    conda: "../envs/gatk.yml"
    params: 
        fa=config["genome_fa"]
    shell: 
        """ gatk --java-options "-Xmx4g" HaplotypeCaller -R {params.fa} -I {input} -O {output.initial} -ERC GVCF
        gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.fa} -V {output.initial} -O {output.final} """

rule HardFiltering:
    input:
        "results_{sample_id}/exome/haplotypeCaller/exome.vcf.gz"
    output:
        final="results_{sample_id}/exome/filtration/filtered.vcf.gz",
        ind_in="results_{sample_id}/exome/filtration/indels.vcf.gz",
        snps_in="results_{sample_id}/exome/filtration/snps.vcf.gz",
        ind_fin="results_{sample_id}/exome/filtration/indels_filtered.vcf.gz",
        snps_fin="results_{sample_id}/exome/filtration/snps_filtered.vcf.gz"
    conda: "../envs/gatk.yml"
    params: 
        dbsnp=config["known_sites"],
        dict_gen=config["genome_dict"],
        sampleid="{sample_id}"
    shell:
        """ mkdir -p tmp
        gatk SelectVariants -V {input} --tmp-dir tmp -select-type SNP -O {output.snps_in}
        gatk SelectVariants -V {input} --tmp-dir tmp -select-type INDEL -O {output.ind_in}

        gatk VariantFiltration -V {output.snps_in} --tmp-dir tmp -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O {output.snps_fin}

        gatk VariantFiltration  -V {output.ind_in} --tmp-dir tmp -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O {output.ind_fin}

        gatk MergeVcfs I={output.snps_fin} I={output.ind_fin} O={output.final}

        gatk CollectVariantCallingMetrics -I {output.final} --DBSNP {params.dbsnp} -SD {params.dict_gen} -O results_{params.sampleid}/exome/filtration/metrics """
    

rule het_selection_exome:
    input: 
        "results_{sample_id}/exome/filtration/snps_filtered.vcf.gz"
    output:
        fin="results_{sample_id}/exome/filtration/snps_het.vcf.gz",
        inter=temp("results_{sample_id}/exome/filtration/snps_het_1.vcf.gz")
    params:
        samp=config["sample_name"],
        samplefile="results_{sample_id}/exome/filtration/sample.txt"
    conda: "../envs/samtools.yml"
    shell:
        """ echo "exome_{params.samp} {params.samp}"> {params.samplefile}
        bcftools view {input} -i 'GT=="het"' -m2 -M2 -O z -o {output.inter}
        bcftools reheader -s {params.samplefile} {output.inter} -o {output.fin}
        tabix {output.fin} """  

##########################
######## phasing ########
##########################
        
#phase the quality controlled called SNPs with shapeit4
rule phasing_SHAPEIT4:
    input:
        vcftbi="results_{sample_id}/phased/pre_phased.vcf.gz.tbi",
        vcf="results_{sample_id}/phased/pre_phased.vcf.gz",
        ref_vcf=expand("{path}/{prefix}{chrom}.{extension}", path=config["path_ALLvcf"],prefix=config["prefix_ALLvcf"],extension=config["extension_ALLvcf"],  allow_missing=True)
    output: 
        "results_{sample_id}/phased/{chrom}_phased.vcf"
    conda:
        "../envs/shapeit.yml"
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

rule phasing_haptreex:
    input:
        vcf="results_{sample_id}/merged_vcf/variantsQC.vcf.gz",
        bam="results_{sample_id}/gex_bulk/alignment/gex_bulk.positionsort.bam"
    output:
        "results_{sample_id}/phased/haptreex.tsv"
    params:
        gtffile=config["genome_gtf"],
        haptreex=config["haptreex_exe"],
        tec=config["tech"],
        htslib=config["htslib_path"],
        sampleid="{sample_id}"
    shell:
        """ bcftools view {input.vcf} -Ov -o temp.vcf
        export LD_LIBRARY_PATH={params.htslib}
        if [[ "{params.tec}" == *exome* ]]
        then
            {params.haptreex} -v temp.vcf -r {input.bam} -g {params.gtffile} -o {output} -d results_{params.sampleid}/exome/recalibration/exome.recal.bam
        elif [[ "{params.tec}" == *atac* ]]
        then
            {params.haptreex} -v temp.vcf -r {input.bam} -g {params.gtffile} -o {output} -d results_{params.sampleid}/atac/mapping_result/atac.positionsort.MAPQ30.bam
        else
            {params.haptreex} -v temp.vcf -r {input.bam} -g {params.gtffile} -o {output}
        fi
        rm temp.vcf """

rule bgzip_and_indexing:
    input:
        expand("results_{sample_id}/phased/{chrom}_phased.vcf",chrom=config["chromosomes_to_phase"],sample_id=config["sample_name"])
    output: 
        vcf="results_{sample_id}/phased/shapeit_whatshap.vcf.gz",
        vcftbi="results_{sample_id}/phased/shapeit_whatshap.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input} -Oz -o {output.vcf}
            tabix {output.vcf} """

rule manual_phasing:
    input:
        haptreex="results_{sample_id}/phased/haptreex.tsv",
        shapeit="results_{sample_id}/phased/shapeit_whatshap.vcf.gz",
        not_phased="results_{sample_id}/merged_vcf/variantsQC.vcf.gz",
        ase="results_{sample_id}/gex_bulk/ASE{chrom}_bulk"
    output:
        temp("results_{sample_id}/phased/manual_phasing{chrom}.tsv")
    params:
        sample=config["sample_name"]
    shell:
        """ Rscript --vanilla workflow/scripts/iterative.R -c {wildcards.chrom} -v {input.not_phased} \
        -s {input.shapeit} \
        -x {input.haptreex} \
        -a {input.ase}/gatkgroup.sorted.table \
        -n {wildcards.sample_id} \
        -o {output}
        """

rule tsv_to_vcf:
    input:
        "results_{sample_id}/phased/manual_phasing{chrom}.tsv"
    output: 
        temp("results_{sample_id}/phased/manual_refinment{chrom}.vcf")
    conda:
        "../envs/samtools.yml"
    params:
        sample=config["sample_name"]
    shell:
        """ day=$(date "+%d/%m/%4Y %T")
        echo "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=$day
##source=manual_phasing
##contig=<ID={wildcards.chrom}>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">
##INFO=<ID=CM,Number=A,Type=Float,Description="Interpolated cM position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{params.sample}" >  {output}

        less {input} | awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t.\t.\tAC=1;AF=0.5\tGT\t"$6"|"$7}}'>> {output}
        """

rule bgzip_and_indexing_man:
    input:
        expand("results_{sample_id}/phased/manual_refinment{chrom}.vcf",chrom=config["chromosomes_to_phase"],sample_id=config["sample_name"])
    output: 
        vcf="results_{sample_id}/phased/manual_refinment.vcf.gz",
        vcftbi="results_{sample_id}/phased/manual_refinment.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input} -Oz -o {output.vcf}
            tabix {output.vcf} """
