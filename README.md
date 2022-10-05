# SingleCellAnalysis

Snakemake workflow for the analysis and integration of scRNA-seq and scATAC-seq data.

## Input files needed

As input are needed the gzipped fastq file of *one* sample (possibility to input multi-lane fastq). 

## Reference files preparation

Most of the scATAC-seq data analysis based on the [scATAC-pro](https://github.com/wbaopaul/scATAC-pro) tool. 
For this reason need to adapt also `config/scATACconfig.txt` file. Moreover needed annotation files available at scATAC-pro github page, here reported in the `scATAC-pro_annotation/` folder.



## How to run

## Results architecture