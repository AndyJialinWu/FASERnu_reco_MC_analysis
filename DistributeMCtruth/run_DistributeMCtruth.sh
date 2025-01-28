#!/bin/bash

dir_pinu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_light_eposlhc_pi_10invab
dir_knu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_light_eposlhc_k_10invab_200026
dir_charmnu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_charm_p8_monash_central_10invab

declare -a dir_reco_MC=(${dir_pinu} ${dir_knu} ${dir_charmnu})
echo "${dir_reco_MC[0]}"
echo "${dir_reco_MC[1]}"
echo "${dir_reco_MC[2]}"

dir_jw_test_pinu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_NC_10invab/200025/jw_test/
dir_jw_test_knu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_NC_10invab/200026/jw_test/
dir_jw_test_charmnu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_NC_10invab/200035/jw_test/

declare -a dir_jw_test=(${dir_jw_test_pinu} ${dir_jw_test_knu} ${dir_jw_test_charmnu})

MC_count=${#dir_reco_MC[@]}

# reco. MC 200025, 200026, 200035 loop
#for ((MCit=0; MCit<MC_count; MCit++))
#do
    MCit=0
    # count the number of reco_evt_vtx_info_*.txt
    file_count=$(ls ${dir_reco_MC[$MCit]}/reco_evt_vtx_info/reco_evt_vtx_info_*.txt 2>/dev/null | wc -l)
    echo $file_count

    # reco_evt_vtx_info_*.txt file loop
    for ((FileIt=30; FileIt<file_count; FileIt++))
    do
        fileNum=$(printf "%05d" ${FileIt})
        EvtVtxInfo=${dir_reco_MC[$MCit]}/reco_evt_vtx_info/reco_evt_vtx_info_${fileNum}.txt

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
            EventID=${event[0]}

            # check linked_tracks.root 
            dir_linked_tracks=${dir_reco_MC[$MCit]}/cp_data_files/cp_data_file_${FileIt}/evt_${event[0]}_pl${event[3]}_${event[4]}

            if [ -f "$dir_linked_tracks/linked_tracks.root" ]; then
                # found linked_tracks.root

                ScriptName=DistributeMCtruth.cxx
                root -l -q $ScriptName'("'${dir_linked_tracks}'", "'${dir_jw_test[$MCit]}'", '${MCit}', "'${fileNum}'", '${EventID}')'

            else
                # missing linked_tracks.root
                echo "$dir_linked_tracks/linked_tracks.root does not exist."

            fi
        done
    done

#done

