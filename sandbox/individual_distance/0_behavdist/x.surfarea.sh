#!/bin/bash

# cd  /a/documents/connectome/_all/
# List=`ls -d *`

# dir=/a/documents/connectome/_all/
# ldir=/scr/animals1/Dropbox/01_code/topography/sandbox/individual_distance/0_behavdist/hcp.surfarea/
# cd $ldir
# cd ..


# for i in $List;
# do
# 	for h in L R;
# 	do
# 		echo ${i}
# 		wb_command -surface-vertex-areas ${dir}${i}/MNINonLinear/fsaverage_LR32k/${i}.${h}.midthickness.32k_fs_LR.surf.gii ${ldir}${i}.${h}.midthickness.32k_fs_LR.surf.surfarea.gii
# 	done
# done


dir=/scr/animals1/lsd_freesurfer/
cd $dir
List=`ls -d *`
ldir=/scr/animals1/Dropbox/01_code/topography/sandbox/individual_distance/0_behavdist/lsd.surfarea/
cd $ldir
cd ..

for i in $List;
do
	for h in lh rh;
	do
 		echo ${i}
 		wb_command -surface-vertex-areas ${dir}${i}/MNINonLinear/fsaverage_LR32k/${i}.${h}.midthickness.32k_fs_LR.surf.gii ${ldir}${i}.${h}.midthickness.32k_fs_LR.surf.surfarea.gii
 	done
 done

${dir}${i}/surf/${h}.pial ${ldir}${i}.${h}.midthickness.32k_fs_LR.surf.surfarea.gii