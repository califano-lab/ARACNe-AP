# Arguments
acro=${1}
reg=${2}

##### environment
rscript="/nfs/apps/R/3.3.1/bin/Rscript --vanilla"

##### getting the aracne network
cp ./${acro}-${reg}/network.txt ./${acro}-${reg}_4col.tsv
sed -i '1d' ./${acro}-${reg}_4col.tsv

##### generating the regulon object
$rscript regulon.r ${acro} ${reg} ${acro}-${reg}-regulon.rda

##### cleaning up
##rm -r ./${acro}-${reg}
##rm -r ./${acro}-${reg}_all
##rm dat_${acro}-${reg}.txt thr_${acro}-${reg}.txt ar_${acro}-${reg}_*.txt chk_${acro}-${reg}.txt net_${acro}-${reg}.txt cl_${acro}-${reg}.txt

