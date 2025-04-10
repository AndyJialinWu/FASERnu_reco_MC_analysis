#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "m_NuHit_tree.h"
#include "m_NuMCTruth_tree.h"
#include "NuMCTruth_kinematics.C"

void DistMCtruth(TString dir_NTUP, TString dir_jw_test, TString dir_linked_tracks, TString FileNum, int EventID){

    TString NTUP_name = "FaserMC-MC24_PG_neut_in_fasernu_xin-100069-" + FileNum + "-s0013-NTUP.root";
    TString jw_test_name = "FaserMC-MC24_PG_neut_in_fasernu_xin-100069-" + FileNum + "-s0013-NTUP_jw_test.root";

    TFile *f_NTUP = new TFile(dir_NTUP+NTUP_name, "READ");
    TTree *t_hits = (TTree *)f_NTUP->Get("m_NuHit_tree");
    TTree *t_kinematics = (TTree *)f_NTUP->Get("m_NuMCTruth_tree");
    m_NuHit_tree *NTUP_hit = new m_NuHit_tree(t_hits);
    m_NuMCTruth_tree *NTUP_kinematics = new m_NuMCTruth_tree(t_kinematics);

    TString NTUP_output = dir_linked_tracks + "MCTruth_NTUP.root";
    TFile *f_NTUP_output = new TFile(NTUP_output, "RECREATE");
    TTree *t_hits_output = NTUP_hit->fChain->CloneTree(0);
    TTree *t_kinematics_output = NTUP_kinematics->fChain->CloneTree(0);

    // hits loop
    for(int hitIt=0; hitIt<NTUP_hit->fChain->GetEntries(); hitIt++){

        NTUP_hit->GetEntry(hitIt);

        if(NTUP_hit->m_event_id == EventID){
            f_NTUP_output->cd();
            t_hits_output->Fill();
        }

    }

    // tracks kinematics loop
    for(int trkIt=0; trkIt<NTUP_kinematics->fChain->GetEntries(); trkIt++){

        NTUP_kinematics->GetEntry(trkIt);

        if(NTUP_kinematics->m_event_id_MC == EventID){
            f_NTUP_output->cd();
            t_kinematics_output->Fill();
        }

    }

    f_NTUP_output->cd();
    t_hits_output->Write();
    t_kinematics_output->Write();
    f_NTUP_output->Save();
    f_NTUP_output->Close();
    delete f_NTUP_output;

    f_NTUP->Close();
    delete f_NTUP;

    TFile *f_jwtest = new TFile(dir_jw_test+jw_test_name, "READ");
    TTree *t_jwtest = (TTree *)f_jwtest->Get("NuMCTruth_kinematics");
    NuMCTruth_kinematics *jw_test = new NuMCTruth_kinematics(t_jwtest);

    TString jwtest_output = dir_linked_tracks + "jw_test.root";
    TFile *f_jwtest_output = new TFile(jwtest_output, "RECREATE");
    TTree *t_jwtest_output = jw_test->fChain->CloneTree(0);

    for(int EvtIt=0; EvtIt<jw_test->fChain->GetEntries(); EvtIt++){

        jw_test->GetEntry(EvtIt);

        if(jw_test->m_event_id_MC == EventID){

            f_jwtest_output->cd();
            t_jwtest_output->Fill();

            break;

        }

    }

    f_jwtest_output->cd();
    t_jwtest_output->Write();
    f_jwtest_output->Save();
    f_jwtest_output->Close();
    delete f_jwtest_output;

    f_jwtest->Close();
    delete f_jwtest;

}