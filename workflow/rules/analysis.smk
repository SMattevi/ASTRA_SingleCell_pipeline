##########################
#### ATACseq analysis ####
##########################

#creation of the fragment file for ATACseq with sinto pipeline 
rule fragment_file:
    input: 
        "results_{sample_id}/atac/mapping_result/atac.positionsort.MAPQ30.bam"
    threads: config["threads_num"]
    output:
        temp("results_{sample_id}/atac/summary/fragments.notsorted.tsv")
    conda: "../envs/sinto.yml"
    shell:
        """ sinto fragments -b {input} -p {threads} -f {output} --barcode_regex "[^:]*" --use_chrom "(?i)" """

#sort fragments
rule sort:
    input: 
        "results_{sample_id}/atac/summary/fragments.notsorted.tsv"
    output:
        final="results_{sample_id}/atac/summary/atac.fragments.tsv.gz"
    conda: "../envs/samtools.yml"
    threads: 
        config["threads_num"]
    params:
        "results_{sample_id}/atac/summary/atac.fragments.tsv"
    shell:
        """ sort -k1,1 -k2,2n {input} > {params}
        bgzip -@ {threads} {params}
        tabix -p bed {output.final} """

#Call scATAC peaks with MACS2(can be changed from scATACconfig file)
rule peak_calling:
    input: 
        "results_{sample_id}/atac/mapping_result/atac.positionsort.MAPQ30.bam"
    output:
        "results_{sample_id}/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed"
    params:
        config["bl_file"]
    conda:
        "../envs/peakcall.yml"
    shell:
        """ 
        mkdir -p results_{wildcards.sample_id}/atac/peaks/MACS2/
        macs2 callpeak -t {input} --outdir results_{wildcards.sample_id}/atac/peaks/MACS2/ -n atac -f BAMPE -q 0.05 -g hs --nomodel --extsize 200 --shift -100
		
        awk '{{ if ($1>=1 && $1<=22 || $1=="X" || $1=="Y" || $1=="M") {{print $0}}}}' results_{wildcards.sample_id}/atac/peaks/MACS2/atac_peaks.narrowPeak\
            >results/atac/peaks/MACS2/tmp.narrowPeak
        mv results_{wildcards.sample_id}/atac/peaks/MACS2/tmp.narrowPeak results_{wildcards.sample_id}/atac/peaks/MACS2/atac_peaks.narrowPeak

        bedtools intersect -a results_{wildcards.sample_id}/atac/peaks/MACS2/atac_peaks.narrowPeak -b {params} -v \
            > results_{wildcards.sample_id}/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed
        """

rule get_mtx:
    input: 
        frag="results_{sample_id}/atac/summary/atac.fragments.tsv.gz",
        peaks="results_{sample_id}/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output: 
        "results_{sample_id}/atac/raw_matrix/MACS2/matrix.mtx"
    shell:
        """ scATAC-pro -s get_mtx \
            -i {input.frag},{input.peaks} \
            -c {params} \
            -o results_{wildcards.sample_id}/atac """

rule qc_per_barcode:
    input: 
        frag="results_{sample_id}/atac/summary/atac.fragments.tsv.gz",
        peaks="results_{sample_id}/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output:
        "results_{sample_id}/atac/summary/atac.MACS2.qc_per_barcode.txt"
    shell:
        """ scATAC-pro -s qc_per_barcode \
            -i {input.frag},{input.peaks} \
            -c {params} \
            -o results_{wildcards.sample_id}/atac """

#Call cell-> FILTER on MACS2 results-> change method in scATACconfig.txt
rule call_cell:
    input: 
        raw="results_{sample_id}/atac/raw_matrix/MACS2/matrix.mtx",
        qc="results_{sample_id}/atac/summary/atac.MACS2.qc_per_barcode.txt"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output:
        "results_{sample_id}/atac/filtered_matrix/MACS2/FILTER/barcodes.txt",
        "results_{sample_id}/atac/filtered_matrix/MACS2/FILTER/matrix.rds"
    shell:
        """ scATAC-pro -s call_cell \
            -i {input.raw} \
            -c {params} \
            -o results_{wildcards.sample_id}/atac """

rule rmDoublets:
    input: 
        "results_{sample_id}/atac/filtered_matrix/MACS2/FILTER/matrix.rds"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output:
        "results_{sample_id}/atac/filtered_matrix/MACS2/FILTER/barcodes_doubletsRemoved.txt"
    shell:
        """ taskset 32 scATAC-pro -s rmDoublets \
            -i {input},0.03 \
            -c {params} \
            -o results_{wildcards.sample_id}/atac """

##########################
#### RNAseq analysis #####
##########################

#6. Assign reads to genes
rule featureCounting:
    input:
        bam="results_{sample_id}/gex/alignment/gex.positionsort.bam"
    output:
        outreal="results_{sample_id}/gex/alignment/gex.positionsort.bam.featureCounts.bam",
        outsum="results_{sample_id}/gex/alignment/gene_assigned"
    params:
        genome= config["genome_gtf"]
    conda:
        "../envs/umitools.yml"
    threads: config["threads_num"]
    shell:
        """featureCounts -a {params.genome} \
              -o {output.outsum} \
              -R BAM {input.bam} \
              -T {threads}  """

#7. Sorting and indexing
rule sortingANDindexing:
    input:
        "results_{sample_id}/gex/alignment/gex.positionsort.bam.featureCounts.bam"
    output:
        "results_{sample_id}/gex/alignment/assigned_sorted.bam"
    conda:
        "../envs/samtools.yml"
    threads: config["threads_num"]
    shell:
        """samtools sort {input} -o {output} -@{threads}
        samtools index {output} -@{threads}"""

#8. Count UMIs per gene per cell
rule umi_counts:
    input:
        "results_{sample_id}/gex/alignment/assigned_sorted.bam"
    output:
        "results_{sample_id}/gex/features/counts.tsv.gz"
    conda:
        "../envs/umitools.yml"
    shell:
        """ umi_tools count --per-gene \
                --gene-tag=XT --assigned-status-tag=XS \
                --per-cell -I {input} \
                -S {output} 
        """

rule clustering_scGEX:
    input:
        "results_{sample_id}/gex/features/counts.tsv.gz"
    output:
        "results_{sample_id}/gex/features/cluster_gex.tsv"
    singularity:
        "docker://timoast/signac"
    shell:
        """ cd results_{wildcards.sample_id}
        Rscript ../workflow/scripts/gex_clustering.R"""

rule assign_cluster_to_scATAC:
    input:
        "results_{sample_id}/atac/filtered_matrix/MACS2/FILTER/matrix.rds",
        "results_{sample_id}/gex/features/cluster_gex.tsv"
    output:
        "results_{sample_id}/atac/features/fromgex_cluster_atac.tsv"
    singularity:
        "docker://timoast/signac"
    shell:
        """ mkdir -p results_{wildcards.sample_id}/atac/features
        mkdir -p results_{wildcards.sample_id}/plot
        Rscript workflow/scripts/atac_clustering.R """

rule clustering_scATAC:
    input: 
        "results_{sample_id}/atac/filtered_matrix/MACS2/FILTER/matrix.rds"
    params:
        config["config_file_scATACpro"]
    singularity:
        "docker://wbaopaul/scatac-pro:latest"
    output:
        "results_{sample_id}/atac/downstream_analysis/MACS2/FILTER/cell_cluster_table.tsv"
    shell:
        """ taskset 32 scATAC-pro -s clustering \
            -i {input} \
            -c {params} \
            -o results_{wildcards.sample_id}/atac """
            
rule decision_scATAC_clustering:
    input: 
        config["which_clustering"]
    output:
        "results_{sample_id}/atac/features/cluster_atac.tsv"
    shell:
        "mv {input} {output}"
