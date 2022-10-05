##########################
#### ATACseq analysis ####
##########################

#creation of the fragment file for ATACseq with sinto pipeline 
rule fragment_file:
    input: 
        "results/atac/mapping_result/atac.positionsort.MAPQ30.bam"
    threads: config["threads_num"]
    output:
        temp("results/atac/summary/fragments.notsorted.tsv")
    conda: "envs/sinto.yml"
    shell:
        """ sinto fragments -b {input} -p {threads} -f {output} --barcode_regex "[^:]*" --use_chrom "(?i)" """

#sort fragments
rule sort:
    input: 
        "results/atac/summary/fragments.notsorted.tsv"
    output:
        final="results/atac/summary/atac.fragments.tsv.gz"
    conda: "envs/samtools.yml"
    threads: 
        config["threads_num"]
    params:
        "results/atac/summary/atac.fragments.tsv"
    shell:
        """ sort -k1,1 -k2,2n {input} > {params}
        bgzip -@ {threads} {params}
        tabix -p bed {output.final} """

#Call scATAC peaks with MACS2(can be changed from scATACconfig file)
rule peak_calling:
    input: 
        "results/atac/mapping_result/atac.positionsort.MAPQ30.bam"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output:
        "results/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed"
    shell:
        """ scATAC-pro -s call_peak \
            -i {input} \
            -c {params} \
            -o results/atac """

rule get_mtx:
    input: 
        frag="results/atac/summary/atac.fragments.tsv.gz",
        peaks="results/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output: 
        "results/atac/raw_matrix/MACS2/matrix.mtx"
    shell:
        """ scATAC-pro -s get_mtx \
            -i {input.frag},{input.peaks} \
            -c {params} \
            -o results/atac """

rule qc_per_barcode:
    input: 
        frag="results/atac/summary/atac.fragments.tsv.gz",
        peaks="results/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output:
        "results/atac/summary/atac.MACS2.qc_per_barcode.txt"
    shell:
        """ scATAC-pro -s qc_per_barcode \
            -i {input.frag},{input.peaks} \
            -c {params} \
            -o results/atac """

#Call cell-> FILTER on MACS2 results-> change method in scATACconfig.txt
rule call_cell:
    input: 
        raw="results/atac/raw_matrix/MACS2/matrix.mtx",
        qc="results/atac/summary/atac.MACS2.qc_per_barcode.txt"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output:
        "results/atac/filtered_matrix/MACS2/FILTER/barcodes.txt",
        "results/atac/filtered_matrix/MACS2/FILTER/matrix.rds"
    shell:
        """ scATAC-pro -s call_cell \
            -i {input.raw} \
            -c {params} \
            -o results/atac """

rule rmDoublets:
    input: 
        "results/atac/filtered_matrix/MACS2/FILTER/matrix.rds"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output:
        "results/atac/filtered_matrix/MACS2/FILTER/barcodes_doubletsRemoved.txt"
    shell:
        """ taskset 32 scATAC-pro -s rmDoublets \
            -i {input},0.03 \
            -c {params} \
            -o results/atac """

##########################
#### RNAseq analysis #####
##########################

#6. Assign reads to genes
rule featureCounting:
    input:
        bam="results/gex/alignment/gex.positionsort.bam",
        genome= config["genome_gtf"]
    output:
        outreal="results/gex/alignment/gex.positionsort.bam.featureCounts.bam",
        outsum="results/gex/alignment/gene_assigned"
    conda:
        "envs/umitools.yml"
    threads: config["threads_num"]
    shell:
        """featureCounts -a {input.genome} \
              -o {output.outsum} \
              -R BAM {input.bam} \
              -T {threads}  """

#7. Sorting and indexing
rule sortingANDindexing:
    input:
        "results/gex/alignment/gex.positionsort.bam.featureCounts.bam"
    output:
        "results/gex/alignment/assigned_sorted.bam"
    conda:
        "envs/samtools.yml"
    threads: config["threads_num"]
    shell:
        """samtools sort {input} -o {output} -@{threads}
        samtools index {output} -@{threads}"""

#8. Count UMIs per gene per cell
rule umi_counts:
    input:
        "results/gex/alignment/assigned_sorted.bam"
    output:
        "results/gex/features/counts.tsv.gz"
    conda:
        "envs/umitools.yml"
    shell:
        """ umi_tools count --per-gene \
                --gene-tag=XT --assigned-status-tag=XS \
                --per-cell -I {input} \
                -S {output} 
        """

rule clustering_scGEX:
    input:
        "results/gex/features/counts.tsv.gz"
    output:
        "results/gex/features/cluster_gex.tsv"
    singularity:
        "docker://timoast/signac"
    script:
        "scripts/gex_clustering.R"

rule assign_cluster_to_scATAC:
    input:
        "results/atac/filtered_matrix/MACS2/FILTER/matrix.rds",
        "results/gex/features/cluster_gex.tsv"
    output:
        "results/atac/features/fromgex_cluster_atac.tsv"
    singularity:
        "docker://timoast/signac"
    shell:
        """ mkdir -p results/atac/features
        mkdir -p results/plot
        Rscript workflow/scripts/atac_clustering.R """

rule clustering_scATAC:
    input: 
        "results/atac/filtered_matrix/MACS2/FILTER/matrix.rds"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output:
        "results/atac/downstream_analysis/MACS2/FILTER/cell_cluster_table.tsv"
    shell:
        """ taskset 32 scATAC-pro -s clustering \
            -i {input} \
            -c {params} \
            -o results/atac """
            
rule decision_scATAC_clustering:
    input: 
        config["which_clustering"]
    output:
        "results/atac/features/cluster_atac.tsv"
    shell:
        "mv {input} {output}"
