#!/bin/bash

# set up root environment
source /opt/root/root-6-32-04/bin/thisroot.sh

# the script to run
ScriptName=/home/jialwu/FASERnu_reco_MC_analysis/DistributeMCtruth/DistMCtruth.cxx

# Data directory and MC truth directory
DataDir=/home/jialwu/raid/FASERnu_MC_reco/MC24_PG_neut_in_fasernu_xin_100069
dir_NTUP=$DataDir/NTUP/
dir_jw_test=$DataDir/jw_test/

# Get the file number from the first argument of the condor submission script
nskip=36
let x=${1}+${nskip}
FileNum=$(printf "%05d" ${x})

# Find the correct path to linked_tracks.root
for dir_linked_tracks in "$DataDir/cp_data_files/cp_data_file_${FileNum}"/evt_*/
do
    if [ -f "${dir_linked_tracks}linked_tracks.root" ]; then
        
        # Extract evt_id_PlStart_End from dir_linked_tracks
        evt_id_PlStart_End=$(basename "$dir_linked_tracks")  
        
        # Extract EventID from evt_id_PlStart_End
        EventID=$(echo $evt_id_PlStart_End | awk -F'_' '{print $2}')

        # run the DistMCtruth.cxx
        root -l -b -q $ScriptName'("'${dir_NTUP}'", "'${dir_jw_test}'", "'${dir_linked_tracks}'", "'${FileNum}'", '${EventID}')'

    fi
done


