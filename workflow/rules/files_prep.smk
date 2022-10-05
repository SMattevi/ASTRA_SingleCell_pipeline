##########################
#### ATAC preparation ####
##########################

#Debarcode 10x scATAC FASTQ files by adding the cell barcode from R2 in from of the original read name for each read in R1 and R3 
#(https://raw.githubusercontent.com/aertslab/single_cell_toolkit/2da8cbf09474903d050ecdb073da1afd99347eee/debarcode_10x_scatac_fastqs.sh)
rule debarcoding_FASTQ:
    input:
        expand("{path}/{sample}_{num}_001.fastq.gz",path=config["input_path"], num=['R1', 'R2','R3'], allow_missing=True)
    conda: "envs/samtools.yml"
    output:
        temp("results/atac/barcodecorrection/{sample}_R1.fastq.gz"),
        temp("results/atac/barcodecorrection/{sample}_R2.fastq.gz")        
    shell:
        """bash workflow/scripts/debarcode_10x_scatac_fastqs.sh {input} results/atac/barcodecorrection/{wildcards.sample}"""

##########################
#### GEX preparation #####
##########################

#Merge the different lanes reads for R2 (assuming only 1 sample)
rule merge_samples_R2:
    input:
        expand("{sample}_R2_001.fastq.gz", sample = config["sample_prefix_lane"])
    output:
        temp("results/gex/fastq/merged_R2.fastq.gz")
    shell:
        "cat {input} > {output}"

#Merge the different lanes reads for R1 (assuming only 1 sample)
rule merge_samples_R1:
    input:
        expand("{sample}_R1_001.fastq.gz", sample = config["sample_prefix_lane"])
    output:
        "results/gex/fastq/merged_R1.fastq.gz"
    shell:
        "cat {input} > {output}"

#Cells are called from scATAC and these are used for the extraction in scRNA-> 10x give a "translation" from atac to rna barcodes (called barcode_metrix.tsv)
rule whitelist_creation_from_atac:
    input:
        atac_bar="results/atac/filtered_matrix/MACS2/FILTER/barcodes_doubletsRemoved.txt",
        barcode_info=config["barcoded_metrics"] #from 10x 
    output:
        "results/gex/umitools_extr/whitelist_atac.txt"
    shell:
        "Rscript workflow/scripts/whitelist_atac_creation.R {input.atac_bar} {input.barcode_info} {output}"

#Cells are called from scRNA and these are used for the extraction in scATAC through Seurat integration
rule whitelist_creation:
    input:
        "results/gex/fastq/merged_R1.fastq.gz"
    output:
        "results/gex/umitools_extr/whitelist_umitools.txt"
    params:
        pattern=config["pattern_umi"],
        cell_num=config["cells_number_expected"]
    conda:
        "envs/umitools.yml"
    shell:
        "umi_tools whitelist --stdin {input} --bc-pattern={params.pattern} --knee-method=density"
                " --log2stderr > {output} "

#Effective extraction of barcodes and UMIs (given the whitelist created in #2) and apposition to read name using umi_tools extract
rule extract_barcodes:
    input:
        Rfirst="results/gex/fastq/merged_R1.fastq.gz",
        Rsecond="results/gex/fastq/merged_R2.fastq.gz",
        whitelist=config["which_whitelist"]
    output:
        first="results/gex/umitools_extr/merged_R1_extracted.fastq.gz",
        second="results/gex/umitools_extr/merged_R2_extracted.fastq.gz"
    params:
        pattern=config["pattern_umi"]
    conda:
        "envs/umitools.yml"
    shell:
        """ umi_tools extract --bc-pattern={params.pattern} \
                 --stdin {input.Rfirst} \
                 --stdout {output.first} \
                 --read2-in {input.Rsecond} \
                 --read2-out {output.second} \
                 --whitelist {input.whitelist} 
        """
