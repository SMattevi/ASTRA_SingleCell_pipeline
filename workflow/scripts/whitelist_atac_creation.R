args = commandArgs(trailingOnly=TRUE)

library(stringr)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

barcodeinfo<-read.csv(args[2], header = T)

barcodeinfo$gex_barcode<-str_replace(barcodeinfo$gex_barcode,"-1","")
barcodeinfo$atac_barcode<-str_replace(barcodeinfo$atac_barcode,"-1","")

atac_barcode<-read.csv(args[1],header=F)
atac_bar<-merge(barcodeinfo,atac_barcode,by.x=c("atac_barcode"),by.y=c("V1"))

write.table(atac_bar[,c("gex_barcode")],args[3], quote=FALSE, row.names = FALSE, col.names = FALSE)
