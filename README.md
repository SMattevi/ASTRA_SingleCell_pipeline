# SingleCellAnalysis

Snakemake workflow for the analysis and integration of scRNA-seq and scATAC-seq data.

## Input files needed

As input are needed the gzipped fastq file of *one* sample (possibility to input multi-lane fastq). 

## Reference files preparation

Most of the scATAC-seq data analysis based on the [scATAC-pro](https://github.com/wbaopaul/scATAC-pro) tool. 
For this reason need to adapt also `config/scATACconfig.txt` file. Moreover needed annotation files available at scATAC-pro github page, here reported in the `scATAC-pro_annotation/` folder.

### Hisat2 indexes
Needed hisat2 indexes, for an example of creation look at `hisat_indexes.sh`

### Genome gtf and fa
See `hisat_indexes.sh` for an example of download. Files available on [Ensembl web site](https://www.ensembl.org/Homo_sapiens/Info/Index).

## How to run

## Results architecture