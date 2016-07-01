# Contains makeord function etc
source("analysis.r")
source("edge.r")
L2spk <- list("scr"="../screens/L2spk-west-ggcol.tsv.rc")
L2spksims <- list("scr"="../screens/L2spk-sims-ggcol.tsv")
spike_list <- list(
		"L2spk" = L2spk,
		"L2spksims" = L2spksims
)
spike_edgeanlist <- lapply(liblist,function(x) { scr <- read.table(x[["scr"]],stringsAsFactors=F,head=T)
                                          edgean(scr,nsample_cutoff = 1)
                                          }
                    )
save.image()
