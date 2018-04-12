# set the threshold to 0
tmp <- commandArgs(TRUE)
acro <- tmp[1]
reg <- tmp[2]
fn <- list.files(path=paste(acro, "-", reg, "_all", sep=""), pattern="miThreshold")
cat("0", file=file.path(paste(acro, "-", reg, "_all", sep=""), fn))

