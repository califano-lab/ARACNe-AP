# ARGS
# input expression file
# input cnv file
# input tsv network
# output
args <- commandArgs(TRUE)
acro <- args[1]
reg <- args[2]
of <- args[3]

readRandomLines <- function(fn, prop, include) {
    f1 <- file(fn, open="rt")
    res <- NULL
    while(TRUE) {
        tmp <- readLines(f1, n=1e6)
        if (length(tmp)==0) break
        tmp <- strsplit(tmp, "\t")
        res1 <- as.numeric(sapply(tmp, function(x) x[3]))
        tmp <- t(sapply(tmp, function(x) x[1:2]))
        names(res1) <- paste(tmp[, 1], tmp[, 2], sep="-")
        rm(tmp)
        pos <- names(res1) %in% include
        if (prop>0) res1 <- c(res1[pos], res1[sample(which(!pos), min(length(which(!pos)), ceiling(length(res1)*prop)))])
        else res1 <- c(res1[pos], res1[sample(which(!pos), min(length(which(pos)), length(which(!pos))))])
        res <- c(res, res1[!duplicated(names(res1))])
    }
    close(f1)
    res
}   

expf <- paste(acro, "-expmatzero.dat", sep="")
ntwf <- paste(acro, "-", reg, "_4col.tsv", sep="")

library(viper)

tmp <- strsplit(readLines(expf), "\t")
regdset <- t(sapply(tmp[-1], function(x) as.numeric(x[-1])))
colnames(regdset) <- tmp[[1]][-1]
rownames(regdset) <- sapply(tmp[-1], function(x) x[1])
regul <- aracne2regulon(ntwf, regdset, gene=FALSE, format="3col", verbose=FALSE)
aracne <- cbind(rep(names(regul), sapply(regul, function(x) length(x$tfmode))), unlist(lapply(regul, function(x) names(x$tfmode)), use.names=FALSE), unlist(lapply(regul, function(x) x$tfmode), use.names=FALSE), unlist(lapply(regul, function(x) x$likelihood), use.names=FALSE))
index <- paste(aracne[, 1], aracne[, 2], sep="-")
miall <- readRandomLines(file.path(paste(acro, "-", reg, "_all", sep=""), "nobootstrap_network.txt"), 0, index)
pos <- names(miall) %in% index
aracned <- ecdf(miall[pos])
nulld <- ecdf(miall[!pos])
miall <- miall[match(index, names(miall))]
rl <- aracned(miall)/(aracned(miall)+1-nulld(miall))
rl[is.na(rl)] <- 1
aracne[, 4] <- rl
regul <- tapply(1:nrow(aracne), aracne[, 1], function(i, aracne) {
    tfmode <- as.numeric(aracne[i, 3])
    names(tfmode) <- aracne[i, 2]
    list(tfmode=tfmode, likelihood=as.numeric(aracne[i, 4]))
}, aracne=aracne)
class(regul) <- "regulon"
save(regul, file=of)
