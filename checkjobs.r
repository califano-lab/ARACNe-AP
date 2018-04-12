args <- commandArgs(TRUE)
acro <- args[1]
reg <- args[2]
rfn <- args[3]
aracne <- "/ifs/scratch/c2b2/ac_lab/malvarez/aracne/hparacne.jar"
rscript="/nfs/apps/R/3.3.1/bin/Rscript --vanilla"
java <- "/nfs/apps/java/1.7.0_25/bin/java"
fn <- list.files(pattern=paste("ar_", acro, "-", reg, sep=""))
tmp <- sapply(fn, function(x) {
    tmp <- readLines(x)
    tmp[length(tmp)]
})
names(tmp) <- sapply(strsplit(fn, "_"), function(x) sub(".txt", "", x[3]))
tmp <- as.numeric(names(tmp))[-grep("Total time elapsed", tmp)]
if (length(tmp)>0) {
    for (i in tmp) {
        qlog <- paste("./ar_", acro, "-", reg, "_", i, ".txt", sep="")
        command <- paste(java, " -Xmx7500M -jar ", aracne, " -e ./", acro, "-expmatna.dat -o ./", acro, "-", reg, "/ --tfs ", rfn, " --pvalue 0.00000001 --threads 4 --seed ", i, sep="")
        system(paste("echo \"", command, "\" | qsub -l mem=10000M,time=4:: -pe smp 4 -N ar_", acro, "-", reg, " -j yes -o ", qlog, " -cwd", sep=""))
    }
    qlog <- paste("./chk_", acro, "-", reg, ".txt", sep="")
    command <- paste(rscript, " checkjobs.r", acro, reg, rfn, sep=" ")
    system(paste("echo \"", command, "\" | qsub -l mem=2G,time=2:: -N chk_", acro, "-", reg, " -j yes -o ", qlog, " -cwd -hold_jid ar_", acro, "-", reg, sep=""))
}
