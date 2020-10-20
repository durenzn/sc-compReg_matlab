#/bin/bash

#edit prepare.sh to input genome
#edit config
#prepare Peak.bed in workdir
bash prepare.sh
cat config code_peakfree.m > run_code.m
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "run_code;exit"
