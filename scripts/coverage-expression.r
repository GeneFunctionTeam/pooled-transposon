esc_expression <- read.table("~/work/RNA-seq_ESgeneColor-mm9.txt")
explevels <- c(0,1,2,5,4,3)
rnaseq <- data.frame(rgb=levels(esc_expression$V9),category=explevels)
names(esc_expression) <- c("chr","start","finish","gene","sth","str","a","b","rgb")
esc_expression <- merge(esc_expression,rnaseq)
esc_expression$hit <- esc_expression$gene %in% all_mappings$gene
barplot(
  by(esc_expression,esc_expression$category,function(x) nrow(x[x$hit == TRUE,])/nrow(x)),
  xlab="Gene Expression Category",
  ylab="Fraction of genes with insertions",
  space=0.4
)

