# Script to correct expression by cnv
#
# Arguments
# tumor acronyms
# input file name

tmp <- commandArgs(TRUE)
acro <- tmp[1]
fn <- tmp[2]

if (!file.exists(paste(acro, "-expmat.dat", sep=""))) {
    load(sub("\\.dat", ".rda", fn))
    expmat[expmat == 0] <- NA
    tmp <- cbind(genes=rownames(expmat), round(expmat, 5))
    tmp <- rbind(colnames(tmp), tmp)
    cat(unlist(apply(tmp, 1, paste, collapse="\t"), use.names=FALSE), sep="\n", file=paste(acro, "-expmatna.dat", sep=""))
}
