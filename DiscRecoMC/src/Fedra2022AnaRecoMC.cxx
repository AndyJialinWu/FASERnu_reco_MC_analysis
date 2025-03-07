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


int main(){

    std::vector<std::string> *EvtNClist[3];
    //std::cout<<__LINE__<<std::endl;
    for(int mcIt=0; mcIt<3; mcIt++){
        EvtNClist[mcIt] = new std::vector<std::string>();  // Allocate each vector
        ReadEventList(EvtNClist[mcIt], EventListFile[mcIt]);
        std::cout<<"# NC events from "<<EventListFile[mcIt]<<" = "<<EvtNClist[mcIt]->size()<<std::endl;
    }
    //PrintEventList(EvtNClist[0]);
    //std::cout<<__LINE__<<std::endl;

    // init CheckRecoPlot class
    CheckRecoPlot *CRP = new CheckRecoPlot();


    // MC 200025 200026 200035 Loop
    for(int mcIt=0; mcIt<3; mcIt++){
        
        int NEvts = EvtNClist[mcIt]->size();
        double EventWeight = Run3ExpNC[mcIt] / NEvts;

        // Event Loop
        for(int EvtIt=0; EvtIt<NEvts; EvtIt++){

            // to be placed in the event loop
            TString LTname = EventDataDir[mcIt] + EvtNClist[mcIt]->at(EvtIt) + "/linked_tracks.root";
            TFile *f_trfile = new TFile(LTname, "READ");
            TTree *t_tracks = (TTree *)f_trfile->Get("tracks");
            LinkedTracks *RecoMC = new LinkedTracks(t_tracks);
            if(t_tracks->GetEntries() == 0) continue;

            TString jwname = EventDataDir[mcIt] + EvtNClist[mcIt]->at(EvtIt) + "/jw_test.root";
            TFile *f_jwtest = new TFile(jwname, "READ");
            TTree *t_jwtest = (TTree *)f_jwtest->Get("NuMCTruth_kinematics");
            NuMCTruth_kinematics *jw_test = new NuMCTruth_kinematics(t_jwtest);
            jw_test->GetEntry(0); // !
            jw_test->GetRecoMCtruth(); // get reconstructed MC truth info.
            std::cout<<"Event "<<jw_test->m_event_id_MC<<" (vx, vy, vz) = "<<"("<<jw_test->m_vx<<", "<<jw_test->m_vy<<", "<<jw_test->m_vz<<")"<<" mm"<<std::endl;

            RecoMC->SetTrueVertex(jw_test->m_vx, jw_test->m_vy, jw_test->m_vz);
            RecoMC->GetRecoMCInfo(false);
            //RecoMC->GetRecoAngleIPdz();
            RecoMC->SortBasedOnTrackID(false);

            // match primary tracks reco. and truth
            size_t trueIt=0;

            // reco. tracks loop
            for(size_t recoIt=0; recoIt<RecoMC->trackID.size(); recoIt++){
                
                if(RecoMC->trackID.at(recoIt) < 19999){
                    // primary tracks

                    if(RecoMC->IsPrimTrack(recoIt)){
                        // cut by IsPrimTrack function

                        RecoMC->n_ch++;

                        CRP->h2_IP_dz[0]->Fill(RecoMC->dz.at(recoIt)/1000., RecoMC->IP.at(recoIt), EventWeight);

                        while(jw_test->trackID.at(trueIt) < RecoMC->trackID.at(recoIt)){
                            trueIt++;
                            if(trueIt >= jw_test->trackID.size()) break;
                        }
                        if(trueIt >= jw_test->trackID.size()) break;

                        if(jw_test->trackID.at(trueIt) == RecoMC->trackID.at(recoIt)){
                            // track ID matches

                            CRP->h1_dtan_theta->Fill(std::tan(RecoMC->theta.at(recoIt)) - std::tan(jw_test->theta.at(trueIt)), EventWeight);
                            CRP->h1_dphi->Fill(RecoMC->phi.at(recoIt) - jw_test->phi.at(trueIt), EventWeight);
                            CRP->h2_nseg_pz->Fill(jw_test->p4.at(trueIt).Pz(), RecoMC->NSeg.at(recoIt), EventWeight);
                            CRP->h2_dphi_pz->Fill(jw_test->p4.at(trueIt).Pz(), RecoMC->phi.at(recoIt) - jw_test->phi.at(trueIt), EventWeight);
                            CRP->h2_dtan_theta_pz->Fill(jw_test->p4.at(trueIt).Pz(), std::tan(RecoMC->theta.at(recoIt)) - std::tan(jw_test->theta.at(trueIt)), EventWeight);

                            double ptrue = RecoMC->ptrue.at(recoIt);
                            double pcoord = RecoMC->pmag_coord.at(recoIt);
                            double pang = RecoMC->pmag_ang.at(recoIt);
                            double pharuhi = RecoMC->pmag_haruhi.at(recoIt);

                            if(pharuhi > 0){
                                CRP->h2_preco_ptrue[0]->Fill(ptrue, pharuhi, EventWeight);
                            }
                            else{
                                CRP->h2_preco_ptrue[0]->Fill(ptrue, -5, EventWeight);
                                CRP->h1_nseg_precoFailed->Fill(RecoMC->NSeg.at(recoIt), EventWeight);
                            }

                            if(pcoord > 0){
                                CRP->h2_preco_ptrue[1]->Fill(ptrue, pcoord, EventWeight);
                            }
                            else{
                                CRP->h2_preco_ptrue[1]->Fill(ptrue, -5, EventWeight);
                                CRP->h1_nseg_precoFailed->Fill(RecoMC->NSeg.at(recoIt), EventWeight);
                            }

                            if(pang > 0){
                                CRP->h2_preco_ptrue[2]->Fill(ptrue, pang, EventWeight);
                            }
                            else{
                                CRP->h2_preco_ptrue[2]->Fill(ptrue, -5, EventWeight);
                                CRP->h1_nseg_precoFailed->Fill(RecoMC->NSeg.at(recoIt), EventWeight);
                            }
                            

                        }

                    }

                }

            }

            CRP->h1_n_ch[0]->Fill(RecoMC->n_ch, EventWeight);
            CRP->h1_n_ch[1]->Fill(jw_test->n_ch, EventWeight);

            for(size_t recoIt=0; recoIt<RecoMC->trackID.size(); recoIt++){
                if(RecoMC->trackID.at(recoIt) >= 20000) {
                    // secondary tracks
                    
                    CRP->h2_IP_dz[1]->Fill(RecoMC->dz.at(recoIt)/1000., RecoMC->IP.at(recoIt), EventWeight);
                }
            }

            f_trfile->Close();
            f_jwtest->Close();

            delete jw_test;
            delete RecoMC;
            
            delete f_trfile;
            delete f_jwtest;

        }// Event Loop

    }// MC 200025 200026 200035 Loop

    CRP->StoreHist2ROOT("Figures");
    delete CRP;
   
    return 0;

}

