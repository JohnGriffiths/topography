#!/bin/bash
	
#export SUBJECTS_DIR=/scr/animals1/lsd_freesurfer/;
#/a/projects/mar004_lsd-lemon-preproc/freesurfer; 
#./condor_sulcal_pits.sh 1 /a/projects/mar004_lsd-lemon-preproc/freesurfer ../../sulcal_pits/subs_lsd.txt

local_dir=/scr/animals1/Dropbox/01_code/topography/sandbox/individual_distance/0_behavdist/

while read i; do
  	echo "running ${i}"

	j=\'\'${i}\'\'
	echo "executable = /a/software/matlab/8.2/bin/matlab" > ${local_dir}/condor/condor_${i}
	echo "arguments = -nodisplay -nosplash -nodesktop -r run('/scr/animals1/Dropbox/01_code/topography/sandbox/individual_distance/0_behavdist/DoPrimaryDist(${j})');exit; " >> ${local_dir}/condor/condor_${i} 
	#echo "arguments = -nodisplay -nosplash -nodesktop -r run('/scr/animals1/Dropbox/01_code/topography/sandbox/individual_distance/0_behavdist/DoSurfArea_HCP(${j})');exit; " >> ${local_dir}/condor/condor_${i} 
	echo "universe = vanilla" >> ${local_dir}/condor/condor_${i}
	echo "output = ${local_dir}/condor/condor_${i}.out" >> ${local_dir}/condor/condor_${i}
	echo "error = ${local_dir}/condor/condor_${i}.error" >> ${local_dir}/condor/condor_${i}
	echo "log = ${local_dir}/condor/condor_${i}.log" >> ${local_dir}/condor/condor_${i}
	echo "request_memory = 2000" >> ${local_dir}/condor/condor_${i}
	echo "request_cpus = 1" >> ${local_dir}/condor/condor_${i}
	echo "getenv = True" >> ${local_dir}/condor/condor_${i}
	echo "notification = Never" >> ${local_dir}/condor/condor_${i}
	echo "queue" >> ${local_dir}/condor/condor_${i}
	echo "" >> ${local_dir}/condor/condor_${i}

	condor_submit ${local_dir}/condor/condor_${i}

done<${1}

