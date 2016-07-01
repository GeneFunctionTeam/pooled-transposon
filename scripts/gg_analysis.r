# Contains makeord function etc
source("analysis.r")
source("edge.r")
files <- list.files(path="../mappings",pattern="*.hashup.genes",full.names=T)
message("Reading and ordering mappings")
olist <- lapply(files,makeord)
save.image()
# Elements of the list that belong to the library:
# L1: 1-10: 285-294
L1 <- list("inv"=c(7:16,101,103,105,107,109,111,113,135,137,139,141,142),"scr"="../screens/L1-west-ggcol.tsv.rc")
## L2: 11-20: 295-304
L2spk <- list("inv"=c(17:26,115,117,119,121,123,125,127,129,131,133),"scr"="../screens/L2spk-west-ggcol.tsv.rc")
L2 <- list("inv"=c(17:26,115,117,119,121,123,125,127,129,131,133),"scr"="../screens/L2-west-ggcol.tsv.rc")
# L3: 21-30: 305-314
L3 <- list("inv"=c(27:36),"scr"="../screens/L3-west-ggcol.tsv.rc")
# L4: 31-40: 315-324
L4 <- list("inv"=c(37:46),"scr"="../screens/L4-west-ggcol.tsv.rc")
# L5: 41-50: 325-334 (+214, 217, 327 & 330 missing, two files for 333)
L5 <- list("inv"=c(1,47:55),"scr"="../screens/L5-west-ggcol.tsv.rc")
# L6: 51-60: 335-344 (336,342 missing, +274,281 [alt preps]; 2 files for 340)
L6 <- list("inv"=c(6,56:64),"scr"="../screens/L6-west-ggcol.tsv.rc")
# L7: 61-70: 345:354
L7 <- list("inv"=c(65:74),"scr"="../screens/L7-west-ggcol.tsv.rc")
# L8: 71-80: 355-364
L8 <- list("inv"=c(75:84),"scr"="../screens/L8-west-ggcol.tsv.rc")
# L9: 81-90: 365-374, 366 missing + 253
L9 <- list("inv"=c(2,85:93),"scr"="../screens/L9-west-ggcol.tsv.rc")
# L10: 91-101: 375-380; 275-280; 282-283 [missing/combined]
L10 <- list("inv"=c(3:5,93:99),"scr"="../screens/L10-west-ggcol.tsv.rc")
# Sims from run A54 - samples L1-7 and -8 missing, alt digest and L10 at end
L1sims <- list("inv"=c(100,102,104,106,108,110,112,134,136),"scr"="../screens/L1-sims-ggcol.tsv")
L2sims <- list("inv"=c(114,116,118,120,122,124,126,128,130,132),"scr"="../screens/L2-sims-ggcol.tsv")
L2spksims <- list("inv"=c(114,116,118,120,122,124,126,128,130,132),"scr"="../screens/L2spk-sims-ggcol.tsv")

# West samples from A54 added above

liblist <- list("L1" = L1,
		"L2spk" = L2spk,
		"L2" = L2,
		"L3" = L3,
		"L4" = L4,
#		"L5" = L5, # Excluded as only one DMSO arm present
		"L6" = L6,
		"L7" = L7,
		"L8" = L8,
		"L9" = L9,
		"L10" = L10,
#		"L1sims" = L1sims, # Excluded as only one DMSO arm present
		"L2sims" = L2sims,
		"L2spksims" = L2spksims)
#
process10k <- function(files,olist,scr) {
	inv <- olist[files]
	# Normalise the data:
	inv <- lapply(inv,function(x) {x$norm <- x$reads/sum(x$tot)*1000000; x})
	# Merge all 1k libraries:
	com <- inv[[1]]
	for (i in 2:length(inv)) {
		com <- merge(com,
                inv[[i]],
                by=c("chr","str","pos","barcode","gene","desc","site"),
                all=T,
                suffixes=c(paste(".",i-1),paste(".",i))
                )
	}

	# Add the screen data:
	com$gtn <- apply(com[,grep("norm",names(com))],1,sum,na.rm=T)
	scr <- read.table(scr,stringsAsFactors=F,head=T)
	scrcom <- merge(scr,com,by="barcode")
	# Choose mappings based on cutoffs
	choice <- chooseMappings(scrcom)
	choice$gtt <- apply(choice[,grep("tot\\.",names(choice))],1,sum,na.rm=T)
	choice
}

message("Creating mapping list")
maplist <- lapply(liblist,function(x,ol){ process10k(x$inv,ol,x$scr) },ol=olist)
save.image()
library(edgeR)
source("edge.r")
message("Starting edgeR analysis")
edgelist <- lapply(liblist,function(x) { scr <- read.table(x[["scr"]],stringsAsFactors=F,head=T); edgeets(scr)})
save.image()
dgelist <- lapply(liblist,function(x) { scr <- read.table(x[["scr"]],stringsAsFactors=F,head=T); dgeobjs(scr)})
edgeanlist <- lapply(liblist,function(x) { scr <- read.table(x[["scr"]],stringsAsFactors=F,head=T); edgean(scr)})
save.image()

resEt <- function(et,thresh=0) {
        et$table <- et$table[et$table$logFC > thresh,]
        et
}
sensEt <- function(et,thresh=0) {
        et$table <- et$table[et$table$logFC < thresh,]
        et
}

resist <- lapply(edgelist,function(x) lapply(x,resEt))
sens <- lapply(edgelist,function(x) lapply(x,sensEt))
save.image()

