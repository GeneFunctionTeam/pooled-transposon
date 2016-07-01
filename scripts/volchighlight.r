volcHighlight <- function(df,drug,bcs=c(),end="",xl=NULL,yl=NULL) {
  pcol <- paste(drug,"_PValue",sep="")
  fccol <- paste(drug,"_logFC",sep="")
  plot(df[,fccol],log10(df[,pcol]),pch=19,main=paste(drug,end),xlab="Fold enrichment/depletion",ylab="P-value (log10)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,col=ifelse(row.names(df) %in% bcs,"red","gray"),xlim=xl,ylim=yl)
 }
