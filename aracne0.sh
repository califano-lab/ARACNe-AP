acro=$1
reg=$2
rfn=$3

aracne=/ifs/scratch/c2b2/ac_lab/hd2326/Scripts/java/hparacne.jar
rscript="/nfs/apps/R/3.3.1/bin/Rscript --vanilla"
java=/nfs/apps/java/1.7.0_25/bin/java

mkdir ./${acro}-${reg}_all
cp ./${acro}-${reg}/miThreshold*.txt ./${acro}-${reg}_all/
rename p1E-8 p1E0 ./${acro}-${reg}_all/miThreshold*.txt
$rscript mithreshold.r ${acro} ${reg}
command="$java -Xmx16000M -jar $aracne \
-e ./${acro}-expmatna.dat \
-o ./${acro}-${reg}_all/ \
--tfs ${rfn} \
--pvalue 1 \
--threads 4 \
--nodpi \
--nobootstrap"
$command

