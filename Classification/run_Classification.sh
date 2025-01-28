#!/bin/bash

dir_pinu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_light_eposlhc_pi_10invab
dir_knu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_light_eposlhc_k_10invab_200026
dir_charmnu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_charm_p8_monash_central_10invab

declare -a dir_reco_MC=(${dir_pinu} ${dir_knu} ${dir_charmnu})

dir_jw_test_pinu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_NC_10invab/200025/jw_test/
dir_jw_test_knu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_NC_10invab/200026/jw_test/
dir_jw_test_charmnu=/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_NC_10invab/200035/jw_test/

declare -a dir_jw_test=(${dir_jw_test_pinu} ${dir_jw_test_knu} ${dir_jw_test_charmnu})

dir_out_pinu=/home/jialwu/FASERnu_reco_MC_analysis/Classification/MC200025
dir_out_knu=/home/jialwu/FASERnu_reco_MC_analysis/Classification/MC200026
dir_out_charmnu=/home/jialwu/FASERnu_reco_MC_analysis/Classification/MC200035

declare -a dir_output=(${dir_out_pinu} ${dir_out_knu} ${dir_out_charmnu})


MC_count=${#dir_reco_MC[@]}

    MCit=2

    file_count=$(ls ${dir_reco_MC[$MCit]}/reco_evt_vtx_info/reco_evt_vtx_info_*.txt 2>/dev/null | wc -l)
    echo $file_count

    # reco_evt_vtx_info_*.txt file loop
    for ((FileIt=0; FileIt<file_count; FileIt++))
    do
        fileNum=$(printf "%05d" ${FileIt})
        input_EvtVtxInfo=${dir_reco_MC[$MCit]}/reco_evt_vtx_info/reco_evt_vtx_info_${fileNum}.txt
        output_EvtVtxInfo=${dir_output[$MCit]}/reco_evt_vtx_info_${fileNum}.txt

        ScriptName=Classification.cxx
        root -l -q $ScriptName'("'${input_EvtVtxInfo}'", "'${output_EvtVtxInfo}'", "'${dir_jw_test[$MCit]}'", '${MCit}', "'${fileNum}'")'

    done

