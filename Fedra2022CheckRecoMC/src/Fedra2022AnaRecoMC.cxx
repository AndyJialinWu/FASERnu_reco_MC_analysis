#include <iostream>
#include <fstream>
#include <string>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "LinkedTracks.h"
#include "NuMCTruth_kinematics.h"

int main(){

    // to be placed in the event loop

    TString fname = "/home/jialwu/raid/FASERnu_MC_reco/NC/200025/evt_10001_pl1_94/linked_tracks.root";
    TFile *trfile = new TFile(fname, "READ");
    TTree *tracks = (TTree *)trfile->Get("tracks");

    TString jwname = "/home/jialwu/raid/FASERnu_MC_reco/NC/200025/evt_10001_pl1_94/jw_test.root";
    TFile *f_jwtest = new TFile(jwname, "READ");
    TTree *t_jwtest = (TTree *)f_jwtest->Get("NuMCTruth_kinematics");
    NuMCTruth_kinematics *jw_test = new NuMCTruth_kinematics(t_jwtest);
    jw_test->GetEntry(0); // !
    

    LinkedTracks *RecoMC = new LinkedTracks(tracks);
    std::cout<<"(vx, vy, vz) = "<<"("<<jw_test->m_vx<<", "<<jw_test->m_vy<<", "<<jw_test->m_vz<<")"<<" mm"<<std::endl;
    RecoMC->SetTrueVertex(jw_test->m_vx, jw_test->m_vy, jw_test->m_vz);
    RecoMC->GetRecoMCInfo();
    RecoMC->GetRecoAngleIP();
    RecoMC->SortBasedOnTrackID();

    trfile->Close();
    f_jwtest->Close();

    delete jw_test;
    delete RecoMC;

    delete trfile;
    delete f_jwtest;

    return 0;

}

