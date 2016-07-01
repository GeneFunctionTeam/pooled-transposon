makeord <- function(file) {
	df <- read.table(file,fill=T,quote="",sep="\t")
	# chr and pos duplicated in these hashup files for some reason - need to drop first two columns.
	df <- df[,c(-1,-2)]
	names(df) <- c("chr","str","pos","tot","reads","barcode","site","gene","desc")
	df$ratio <- df$reads/df$tot
	df$sites <- sapply(df$barcode, function(x,y) {sum(y == x)},y=df$barcode)
	df[order(df$chr,df$pos),]
}

sitesperbc <- function(df) {
plot(by(df$reads,df$barcode,function(x) log10(sum(x))),by(df$reads,df$barcode,length),pch=20,ylab="Number of sites",xlab="Total reads for barcode (log10)")
}

bcspersite <- function(df) {
plot(by(df$tot,df$pos,function(x) log10(sum(x))),by(df$tot,df$pos,length),pch=20,ylab="Number of barcodes",xlab="Total reads for integration site (log10)")
}

compareTopBarcodes <- function(x,y,thresh=0.75) {
sum(x[x$sites > range(x$sites)[2]*thresh,"barcode"] %in% y[y$sites > range(y$sites)[2]*thresh,"barcode"])
#list("x length"=length(x[x$sites > range(x$sites)[2]*thresh,"barcode"]), "y length"=length(y[y$sites > range(y$sites)[2]*thresh,"barcode"]), "overlap" = sum(as.character(x[x$sites > range(x$sites)[2]*thresh,"barcode"]) %in% as.character(y[y$sites > range(y$sites)[2]*thresh,"barcode"])))
}

# To plot number of sites vs read count by barcode:
# plot(by(crb$gtn,crb$barcode,length),by(crb$gtn,crb$barcode,sum),pch=20)

chooseMappings <- function(df,ratio=0.9,reads=10) {
	# For each barcode, discard any mapping that is < 10 normalised reads (reads/million).
	byobj <- by(df,df[,"barcode"], function(x,r=reads,rt=ratio) {
		new <- x[apply(x[,grep("norm",names(x))],1,function(y,rds) {sum(y > rds,na.rm=T) > 0},rds=r),]	
		if(nrow(new) == 0) {
			new <- x[apply(x[,grep("ratio",names(x))],1,function(y,rat) {sum(y > rat,na.rm=T) > 0},rat=rt),]	
		}
		new
	})
	# Make the "by" output into a data frame:
	do.call(rbind,byobj)	
}

