#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem=2G,time=1::
#$ -S /bin/bash
#$ -N aracne
#
# Arguments:
# acronismID
# input file
# regulators
# regulator file

acro=${1}
ifn=${2}
reg=${3}
rfn=${4}
bst=200

aracne=/ifs/scratch/c2b2/ac_lab/hd2326/Scripts/java/hparacne.jar
rscript="/nfs/apps/R/3.3.1/bin/Rscript --vanilla"
java=/nfs/apps/java/1.7.0_25/bin/java

##### Expression matrix
qlog=./datna_${acro}-${reg}.txt
command="$rscript expdatna.r $acro $ifn"
echo $command | qsub -l mem=8000M,time=2:: -N dat_${acro}-${reg} -j yes -o $qlog -cwd
qlog=./datzero_${acro}-${reg}.txt
command="$rscript expdatzero.r $acro $ifn"
echo $command | qsub -l mem=8000M,time=2:: -N dat_${acro}-${reg} -j yes -o $qlog -cwd

##### Threshold
qlog=./thr_${acro}-${reg}.txt
command="$java -Xmx5500M -jar $aracne \
-e ./${acro}-expmatna.dat \
-o ./${acro}-${reg}/ \
--tfs ${rfn} \
--pvalue 0.00000001 --seed 1 \
--calculateThreshold"
echo $command | qsub -l mem=8000M,time=2:: -N thr_${acro}-${reg} -j yes -o $qlog -cwd -hold_jid dat_${acro}-${reg}

##### Run nobootstrap and noDPI
qlog=./ar_${acro}-${reg}_all.txt
qsub -l mem=24G,time=4:: -N ar_${acro}-${reg} -pe smp 4 -j yes -o $qlog -cwd -hold_jid thr_${acro}-${reg} aracne0.sh ${acro} ${reg} ${rfn}

##### Run bootstraps (they can be multithreaded with the option --threads 8)
for i in $(seq 1 $bst)
do
qlog=./ar_${acro}-${reg}_${i}.txt
command="$java -Xmx7500M -jar $aracne \
-e ./${acro}-expmatna.dat \
-o ./${acro}-${reg}/ \
--tfs ${rfn} \
--pvalue 0.00000001 \
--threads 4 \
--seed ${i}"
echo $command | qsub -l mem=10000M,time=4:: -N ar_${acro}-${reg} -pe smp 4 -j yes -o $qlog -cwd -hold_jid thr_${acro}-${reg} 
done

##### Check all jobs finished
qlog=./chk_${acro}-${reg}.txt
command="$rscript checkjobs.r ${acro} ${reg} ${rfn}"
echo $command | qsub -l mem=2G,time=2:: -N chk_${acro}-${reg} -j yes -o $qlog -cwd -hold_jid ar_${acro}-${reg}

##### Consolidate
qlog=./net_${acro}-${reg}.txt
command="$java -Xmx32000M -jar $aracne \
-o ./${acro}-${reg}/ \
--nobonferroni \
--consolidatepvalue 0.01 \
--consolidate"
echo $command | qsub -l mem=40000M,time=6:: -N net_${acro}-${reg} -j yes -o $qlog -cwd -hold_jid chk_${acro}-${reg}

##### Clean-up
qlog=./cl_${acro}-${reg}.txt
qsub -l mem=16G,time=8:: -N cl_${acro}-${reg} -j yes -o $qlog -cwd -hold_jid net_${acro}-${reg} aracne_cleanup.sh ${acro} ${reg}

