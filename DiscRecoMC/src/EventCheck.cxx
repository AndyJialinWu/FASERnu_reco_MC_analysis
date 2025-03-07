#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "Fedra2022AnaRecoMC.h"
#include "CheckRecoPlot.h"

#include "LinkedTracks.h"
#include "NuMCTruth_kinematics.h"
#include "m_NuHit_tree.h"
#include "m_NuMCTruth_tree.h"
#include "CalcDisc.h"

int main(){

    std::vector<std::string> *EvtNClist[3];
    
    for(int mcIt=0; mcIt<3; mcIt++){
        EvtNClist[mcIt] = new std::vector<std::string>();  // Allocate each vector
        ReadEventList(EvtNClist[mcIt], EventListFile[mcIt]);
        std::cout<<"# NC events from "<<EventListFile[mcIt]<<" = "<<EvtNClist[mcIt]->size()<<std::endl;
    }
    //PrintEventList(EvtNClist[0]);
    int mcIt = 0;
    std::cout<<"Which MC run nummber? (0:200025, 1:200026, 2:200035) "<<std::endl;
    std::cin>>mcIt;
    int evtID = 0;
    std::cout<<"Which event ID? "<<std::endl;
    std::cin>>evtID;
        
    int NEvts = EvtNClist[mcIt]->size();
    //double EventWeight = Run3ExpNC[mcIt] / NEvts;

    // Event Loop
    for(int EvtIt=0; EvtIt<NEvts; EvtIt++){

        int EventID = ExtractEventID(EvtNClist[mcIt]->at(EvtIt));
        //std::cout<<"Event ID = "<<EventID<<std::endl;
        if(EventID != evtID) continue;

        // to be placed in the event loop
        TString LTname = EventDataDir[mcIt] + EvtNClist[mcIt]->at(EvtIt) + "/linked_tracks.root";
        TFile *f_trfile = new TFile(LTname, "READ");
        TTree *t_tracks = (TTree *)f_trfile->Get("tracks");
        LinkedTracks *RecoMC = new LinkedTracks(t_tracks);
        if(t_tracks->GetEntries() == 0) break;

        TString jwname = EventDataDir[mcIt] + EvtNClist[mcIt]->at(EvtIt) + "/jw_test.root";
        TFile *f_jwtest = new TFile(jwname, "READ");
        TTree *t_jwtest = (TTree *)f_jwtest->Get("NuMCTruth_kinematics");
        NuMCTruth_kinematics *jw_test = new NuMCTruth_kinematics(t_jwtest);
        jw_test->GetEntry(0);      // !
        jw_test->GetRecoMCtruth(); // get reconstructed MC truth info.
        std::cout<<"Event "<<jw_test->m_event_id_MC<<" (vx, vy, vz) = "<<"("<<jw_test->m_vx<<", "<<jw_test->m_vy<<", "<<jw_test->m_vz<<")"<<" mm"<<std::endl;

        RecoMC->SetTrueVertex(jw_test->m_vx, jw_test->m_vy, jw_test->m_vz);
        RecoMC->GetRecoMCInfo(true);
        //RecoMC->GetRecoAngleIPdz();
        RecoMC->SortBasedOnTrackID(true);
        
        f_trfile->Close();
        f_jwtest->Close();
        delete f_trfile;
        delete f_jwtest;

        break;

    }// Event Loop

            
    

    return 0;

}