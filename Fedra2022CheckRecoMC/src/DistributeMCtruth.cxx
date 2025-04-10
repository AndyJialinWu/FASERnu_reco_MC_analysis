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
    //CheckRecoPlot *CRP = new CheckRecoPlot();

    std::ofstream EvtID_MisMatch;
    EvtID_MisMatch.open("EvtID_MisMatch.txt");

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
            if(RecoMC->eventID != jw_test->m_event_id_MC){
                EvtID_MisMatch<<"Event ID does NOT match for "<<LTname<<"\n";
                continue;
            }

            std::string NTUPname = EvtIDtoMCtruthNTUPfname(jw_test->m_event_id_MC, mcIt);
            std::cout<<NTUPname<<std::endl;
            TFile *f_NTUP = new TFile(NTUPname.c_str(), "READ");
            TTree *t_hits = (TTree *)f_NTUP->Get("m_NuHit_tree");
            TTree *t_kinematics = (TTree *)f_NTUP->Get("m_NuMCTruth_tree");
            m_NuHit_tree *NTUP_hit = new m_NuHit_tree(t_hits);
            m_NuMCTruth_tree *NTUP_kinematics = new m_NuMCTruth_tree(t_kinematics);

            TString NTUP_output = EventDataDir[mcIt] + EvtNClist[mcIt]->at(EvtIt) + "/MCTruth_NTUP.root";
            TFile *f_NTUP_output = new TFile(NTUP_output, "RECREATE");
            TTree *t_hits_output = NTUP_hit->fChain->CloneTree(0);
            TTree *t_kinematics_output = NTUP_kinematics->fChain->CloneTree(0);
            
            // hits loop
            for(int hitIt=0; hitIt<NTUP_hit->fChain->GetEntries(); hitIt++){

                NTUP_hit->GetEntry(hitIt);

                if(NTUP_hit->m_event_id == jw_test->m_event_id_MC){
                    f_NTUP_output->cd();
                    t_hits_output->Fill();
                }

            }

            // tracks kinematics loop
            for(int trkIt=0; trkIt<NTUP_kinematics->fChain->GetEntries(); trkIt++){

                NTUP_kinematics->GetEntry(trkIt);

                if(NTUP_kinematics->m_event_id_MC == jw_test->m_event_id_MC){
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

/*
            RecoMC->GetRecoAngleIPdz();
            RecoMC->SortBasedOnTrackID();

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

                            if(pcoord > 0 && pang > 0){
                                CRP->h2_preco_ptrue->Fill(ptrue, std::min(pcoord, pang), EventWeight);
                            }
                            else if(pcoord*pang < 0){
                                CRP->h2_preco_ptrue->Fill(ptrue, std::max(pcoord, pang), EventWeight);
                            }
                            else{
                                CRP->h2_preco_ptrue->Fill(ptrue, -5, EventWeight); // both failed
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
*/
            f_trfile->Close();
            f_jwtest->Close();
            f_NTUP->Close();

            delete jw_test;
            delete RecoMC;
            delete NTUP_hit;
            delete NTUP_kinematics;

            delete f_trfile;
            delete f_jwtest;
            delete f_NTUP;

        }// Event Loop

    }// MC 200025 200026 200035 Loop

    //CRP->StoreHist2ROOT("Figures");

    //delete CRP;
    EvtID_MisMatch.close();
    
    return 0;

}

