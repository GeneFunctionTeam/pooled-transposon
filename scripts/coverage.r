# Coverage plot
west_cf <- data.frame(lib=character(),genes=integer(),cf=integer(),stringsAsFactors=F)
genes <- c()
for(i in west_lib_names) {
genes <- unique(c(genes,as.character(maplist[[i]]$gene)))
west_cf <- rbind(west_cf,data.frame(lib=as.character(i),genes=length(as.character(unique(maplist[[i]]$gene))),cf=length(genes)))
}

barplot(west_cf$cf,
        ylim=c(0,14000),
        ylab="Unique Genes Hit",
        names.arg=(west_cf$lib),
        space=0.4,
        axes=T)
