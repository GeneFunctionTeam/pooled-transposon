adistgrp <- function (x, y, cst = list(ins = 20, del = 1, sub = 20)) {
    which(adist(x, y, costs = cst) < 40)
}

collapsegrp <- function (a, columns) {
    consensus <- by(a[, c(1, columns)], a$group, function(x) x[which.max(apply(subset(x, 
        select = -barcode), 1, sum)), "barcode"])
    consensus <- data.frame(unlist(lapply(consensus, rbind)), 
        stringsAsFactors = F)
    totals <- data.frame(do.call(rbind, by(a[, columns], a$group, 
        function(x) apply(x, 2, sum))))
    df <- data.frame(consensus, totals)
    names(df)[1] = "barcode"
    return(df)
}


alphagrp3 <- function (d, cols) {
    if (missing(cols)) {
        cols = 2:ncol(d)
    }
    d$alpha <- as.factor(sub("(.{3}).+", "\\1", d$barcode))
    bytest2 <- by(d, d$alpha, function(x) {
        x$group <- grpvec(as.character(x$barcode))
        x
    })
    btc <- lapply(bytest2, collapsegrp, columns = cols)
    do.call(rbind, btc)
}

grpvec <- function (bcs) {
    groups <- rep(0, length(bcs))
    ord <- 1:length(bcs)
    bcdf <- data.frame(barcode = bcs, order = ord, stringsAsFactors = F)
    bcdf <- bcdf[order(nchar(bcdf$barcode)), ]
    bcs <- bcdf$barcode
    rows <- sapply(bcs, adistgrp, y = bcs)
    grp <- sapply(rows, function(x, g) g[x] <- x[1], g = groups)
    bcdf <- data.frame(bcdf, group = grp, stringsAsFactors = F)
    bcdf[order(bcdf$order), "group"]
}
