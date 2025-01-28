#!/bin/bash

dir_pinu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_light_eposlhc_pi_10invab
dir_knu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_light_eposlhc_k_10invab_200026
dir_charmnu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_charm_p8_monash_central_10invab

declare -a dir_reco_MC=(${dir_pinu} ${dir_knu} ${dir_charmnu})
echo "${dir_reco_MC[0]}"
echo "${dir_reco_MC[1]}"
echo "${dir_reco_MC[2]}"

dir_EvtInfo_pinu=/home/jialwu/FASERnu_reco_MC_analysis/Classification/MC200025
dir_EvtInfo_knu=/home/jialwu/FASERnu_reco_MC_analysis/Classification/MC200026
dir_EvtInfo_charmnu=/home/jialwu/FASERnu_reco_MC_analysis/Classification/MC200035

declare -a dir_EvtInfo=(${dir_EvtInfo_pinu} ${dir_EvtInfo_knu} ${dir_EvtInfo_charmnu})

declare -a MCNum=(200025 200026 200035)


MC_count=${#dir_EvtInfo[@]}

    MCit=1
    # count the number of reco_evt_vtx_info_*.txt
    file_count=$(ls ${dir_EvtInfo[$MCit]}/reco_evt_vtx_info_*.txt 2>/dev/null | wc -l)
    echo $file_count

    # reco_evt_vtx_info_*.txt file loop
    for ((FileIt=0; FileIt<file_count; FileIt++))
    do
        fileNum=$(printf "%05d" ${FileIt})
        EvtVtxInfo=${dir_EvtInfo[$MCit]}/reco_evt_vtx_info_${fileNum}.txt

        # store the event info. into array line by line
        declare -a array=()
        LineIt=0
        while IFS= read -r line
        do
            [[ $line = \#* ]] && continue
            array[LineIt]=$line
            let "LineIt++"
        done < "$EvtVtxInfo"

        # event loop in cp_data_file_${FileIt}
        for ((EvtIt=0; EvtIt<LineIt; EvtIt++))
        do
            declare -a event=${array[$EvtIt]}
            cat=${event[5]}

            # check linked_tracks.root 
            dir_linked_tracks=${dir_reco_MC[$MCit]}/cp_data_files/cp_data_file_${FileIt}/evt_${event[0]}_pl${event[3]}_${event[4]}

            if [ -f "$dir_linked_tracks/linked_tracks.root" ]; then
                # found linked_tracks.root
                # copy into categories
                if [ "$cat" -eq 0 ]; then
                    # nue CC
                    echo "${dir_linked_tracks} is classified as nue CC"
                    dir_cp2cat=/home/jialwu/raid/FASERnu_MC_reco/nueCC/${MCNum[$MCit]}/
                    cp -rf $dir_linked_tracks $dir_cp2cat

                elif [ "$cat" -eq 1 ]; then
                    # numu CC
                    echo "${dir_linked_tracks} is classified as numu CC"
                    dir_cp2cat=/home/jialwu/raid/FASERnu_MC_reco/numuCC/${MCNum[$MCit]}/
                    cp -rf $dir_linked_tracks $dir_cp2cat

                elif [ "$cat" -eq 2 ]; then
                    # nutau CC
                    echo "${dir_linked_tracks} is classified as nutau CC"
                    dir_cp2cat=/home/jialwu/raid/FASERnu_MC_reco/nutauCC/${MCNum[$MCit]}/
                    cp -rf $dir_linked_tracks $dir_cp2cat
                
                elif [ "$cat" -eq 3 ]; then
                    # NC
                    echo "${dir_linked_tracks} is classified as NC"
                    dir_cp2cat=/home/jialwu/raid/FASERnu_MC_reco/NC/${MCNum[$MCit]}/
                    cp -rf $dir_linked_tracks $dir_cp2cat

                else
                    # nu-electron scatterings
                    echo "${dir_linked_tracks} is not classfied !"

                fi



            else
                # missing linked_tracks.root
                echo "${dir_linked_tracks}/linked_tracks.root does not exist."

            fi

        done

    
    done

