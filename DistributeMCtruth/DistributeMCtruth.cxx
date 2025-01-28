#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

void DistributeMCtruth(TString dir_linked_tracks, TString dir_jw_test, int MCit, TString fileNum, int EventID){

    //TString dir_jw_test = "~/raid/FASERnu_MC_reco/MC24_Genie_NC_10invab/200025/jw_test/";
    
    TString name_jw_test;
    if(MCit == 0) name_jw_test = "FaserMC-MC24_Genie_light_eposlhc_pi_10invab-200025-" + fileNum + "-s0012-NTUP_jw_test.root";
    if(MCit == 1) name_jw_test = "FaserMC-MC24_Genie_light_eposlhc_k_10invab-200026-" + fileNum + "-s0012-NTUP_jw_test.root";
    if(MCit == 2) name_jw_test = "FaserMC-MC24_Genie_charm_p8_monash_central_10invab-200035-" + fileNum + "-s0012-NTUP_jw_test.root";

    TFile *f_jw_test = new TFile(dir_jw_test+name_jw_test, "READ");
    TTree *t_jw_test = (TTree *)f_jw_test->Get("NuMCTruth_kinematics");
    int NEvts = t_jw_test->GetEntries();
    std::cout<<"# events in jw_test = "<<NEvts<<std::endl;

    TFile *f_output = new TFile(dir_linked_tracks+"/jw_test.root", "RECREATE");
    TTree *t_output = t_jw_test->CloneTree(0);

    int m_event_id_MC;
    t_jw_test->SetBranchAddress("m_event_id_MC", &m_event_id_MC);

    for(int EvtIt=0; EvtIt<NEvts; EvtIt++){

        t_jw_test->GetEntry(EvtIt);

        if(m_event_id_MC == EventID){

            f_output->cd();
            t_output->Fill();

            break;

        }

    }

    f_output->cd();
    t_output->Write();
    f_output->Save();
    f_output->Close();
    delete f_output;

}
