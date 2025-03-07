#!/bin/bash

DataDir=~/raid/FASERnu_MC_reco/MC24_PG_neut_in_fasernu_xin_100069
eventlist=../EventList/100069.txt

rm -f "$eventlist"
touch "$eventlist"

for FileNum in {00000..00199}
do
    # Find the correct path to linked_tracks.root
    for evt_dir in "$DataDir/cp_data_files/cp_data_file_${FileNum}"/evt_*/
    do
        if [ -f "${evt_dir}linked_tracks.root" ]; then
            evt_id_PlStart_End=$(basename "$evt_dir")  # Extract evt_id_PlStart_End from evt_dir
            echo "cp_data_file_${FileNum}/$evt_id_PlStart_End" >> "$eventlist"
        fi
    done
done


