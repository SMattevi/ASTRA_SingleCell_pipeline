#!/bin/bash

gtf=resources/Homo_sapiens.GRCh38.106.chr.gtf
fa=resources/Homo_sapiens.GRCh38.dna.primary_assembly.fa 

# RSEM info: https://github.com/deweylab/RSEM?tab=readme-ov-file#-compilation--installation
rsem-prepare-reference -gtf $gtf \
                            $fa \
                            human-GRCh38-106-tr

# transcript2gene 
cat $gtf \
  | awk 'BEGIN{OFS="\t";FS="\t"} $3=="transcript" {a=gensub(/gene_id "([^"]+)"; .+; transcript_id "(([^"]+))"; .+/,"\\2\t\\1","g",$9); print a}' \
  >transcript2gene.txt
