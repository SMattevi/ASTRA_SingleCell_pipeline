##########################
#### QC post alignment ###
##########################

#Sort alignment result and mark duplicate 
rule mark_dup_BAM:
    input:
        "results/{tec}/alignment/{tec}.sam"
    output:
        final_markdup="results/{tec}/alignment/{tec}.positionsort.bam",
        final_markdup_index="results/{tec}/alignment/{tec}.positionsort.bam.bai",
        tmpsort=temp("results/{tec}/{tec}.sorted.bam"),
        out0=temp("results/{tec}/{tec}.positionsort0.bam"),
        fixmate=temp("results/{tec}/{tec}.fixmate.bam")
    conda: "envs/samtools.yml"
    params: config["memory"]
    threads: config["threads_num"]
    shell: 
        """ mkdir -p tmp
        samtools sort -m {params} -T tmp/ -@ {threads} -n -o {output.tmpsort} {input}
        
        samtools fixmate -@ {threads} -m {output.tmpsort} {output.fixmate}
        
        samtools sort -m {params} -@ {threads} -T tmp/ -o {output.out0} {output.fixmate}
       
        samtools markdup -@ {threads} {output.out0} {output.final_markdup}
        samtools index -@ {threads} {output.final_markdup} 
        
        rmdir tmp """

#QC of the aligned bam file-> params -q 30 -> used to extract cells
rule QC_BAM:
    input:
        "results/{tec}/alignment/{tec}.positionsort.bam"
    threads: config["threads_num"]
    conda: "envs/samtools.yml"
    output:
        "results/{tec}/mapping_result/{tec}.positionsort.MAPQ30.bam"
    shell:
        """ samtools view -f 0x2 -b -h -q 30 -@ {threads} {input} -o {output}
        samtools index -@ {threads} {output} 
        bash workflow/scripts/createsummary.sh {input} {output} {threads} {wildcards.tec} """
