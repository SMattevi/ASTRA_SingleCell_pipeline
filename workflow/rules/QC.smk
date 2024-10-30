##########################
#### QC post alignment ###
##########################

#Sort alignment result and mark duplicate 
rule mark_dup_BAM:
    input:
        "results_{sample_id}/{tec}/alignment/{tec}.sam"
    output:
        final_markdup="results_{sample_id}/{tec}/alignment/{tec}.positionsort.bam",
        final_markdup_index="results_{sample_id}/{tec}/alignment/{tec}.positionsort.bam.bai",
        tmpsort=temp("results_{sample_id}/{tec}/{tec}.sorted.bam"),
        out0=temp("results_{sample_id}/{tec}/{tec}.positionsort0.bam"),
        fixmate=temp("results_{sample_id}/{tec}/{tec}.fixmate.bam")
    conda: "../envs/samtools.yml"
    params: config["memory"]
    threads: config["threads_num"]
    shell: 
        """ mkdir -p tmp
        samtools sort -m {params} -T results_{wildcards.sample_id}/tmp/ -@ {threads} -n -o {output.tmpsort} {input}
        
        samtools fixmate -@ {threads} -m {output.tmpsort} {output.fixmate}
        
        samtools sort -m {params} -@ {threads} -T results_{wildcards.sample_id}/tmp/ -o {output.out0} {output.fixmate}
       
        samtools markdup -@ {threads} {output.out0} {output.final_markdup}
        samtools index -@ {threads} {output.final_markdup} 
        
        rmdir results_{wildcards.sample_id}/tmp """

#QC of the aligned bam file-> params -q 30 -> used to extract cells
rule QC_BAM:
    input:
        "results_{sample_id}/{tec}/alignment/{tec}.positionsort.bam"
    threads: config["threads_num"]
    conda: "../envs/samtools.yml"
    output:
        "results_{sample_id}/{tec}/mapping_result/{tec}.positionsort.MAPQ30.bam"
    shell:
        """ samtools view -f 0x2 -b -h -q 30 -@ {threads} {input} -o {output}
        samtools index -@ {threads} {output} 
        bash workflow/scripts/createsummary.sh {input} {output} {threads} {wildcards.tec} {wildcards.sample_id} """

rule BaseRecalibrator:
    input:
        "results_{sample_id}/exome/alignment/exome.positionsort.bam"
    output: 
        bam=temp("results_{sample_id}/exome/recalibration/exome.positionsort.gatkgroup.bam"),
        table="results_{sample_id}/exome/recalibration/recal_data.table"
    conda: "../envs/gatk.yml"
    params: 
        sites=config["known_sites"], 
        fa=config["genome_fa"],
        samp="{sample_id}"
    threads: config["threads_num"]
    shell:
        """ gatk AddOrReplaceReadGroups -I {input} -O {output.bam} -RGLB DNA -RGPL ILLUMINA -RGPU exome -RGSM exome_{params.samp} -VALIDATION_STRINGENCY SILENT
        gatk BaseRecalibrator -I {output.bam} --known-sites {params.sites} -R {params.fa} -O {output.table} """

rule index:
    input:
        "results_{sample_id}/exome/recalibration/exome.positionsort.gatkgroup.bam"
    output:
        temp("results_{sample_id}/exome/recalibration/exome.positionsort.gatkgroup.bam.bai")
    conda: "../envs/samtools.yml"
    threads: config["threads_num"]
    shell: "samtools index -@ {threads} {input}"

rule BQSR:
    input:
        bam="results_{sample_id}/exome/recalibration/exome.positionsort.gatkgroup.bam",
        bai="results_{sample_id}/exome/recalibration/exome.positionsort.gatkgroup.bam.bai",
        table="results_{sample_id}/exome/recalibration/recal_data.table"
    output:
        "results_{sample_id}/exome/recalibration/exome.recal.bam"
    conda: "../envs/gatk.yml"
    shell: 
        "gatk ApplyBQSR --bqsr {input.table} -I {input.bam} -O {output}"
