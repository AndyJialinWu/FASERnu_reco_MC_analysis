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

    std::vector<std::string> *EvtNClist[4]; // 0: 200025, 1: 200026, 2: 200035, 3: 100069
    
    for(int mcIt=0; mcIt<4; mcIt++){
        EvtNClist[mcIt] = new std::vector<std::string>();   // Allocate each vector
        ReadEventList(EvtNClist[mcIt], EventListFile[mcIt]);
        std::cout<<"# NC events from "<<EventListFile[mcIt]<<" = "<<EvtNClist[mcIt]->size()<<std::endl;
    }
    //PrintEventList(EvtNClist[0]);

    TFile *f_disc = new TFile("PhysicsNTUP_ML.root", "RECREATE");
    TTree *t_disc[4][2];    // 0: 200025, 1: 200026, 2: 200035, 3: 100069; 0: TrainAndTest, 1: Validation
    for(int mcIt=0; mcIt<4; mcIt++){
        t_disc[mcIt][0] = new TTree(Form("disc%d_TrainTest", MonteCarloID[mcIt]), Form("disc%d_TrainTest", MonteCarloID[mcIt]));
        t_disc[mcIt][1] = new TTree(Form("disc%d_Validation", MonteCarloID[mcIt]), Form("disc%d_Validation", MonteCarloID[mcIt]));
        Branch(t_disc[mcIt][0]);   // TrainAndTest
        Branch(t_disc[mcIt][1]);   // Validation
    }

    // MC 200025 200026 200035 100069 Loop
    for(int mcIt=0; mcIt<4; mcIt++){
        
        int NEvts = EvtNClist[mcIt]->size();
        //double EventWeight = Run3ExpNC[mcIt] / NEvts;

        // Event Loop
        for(int EvtIt=0; EvtIt<NEvts; EvtIt++){

            //std::cout<<__LINE__<<std::endl;

            EventInit();
            Discriminators *disc = new Discriminators();

            //std::cout<<__LINE__<<std::endl;

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

            //std::cout<<__LINE__<<std::endl;

            TString NTUPname = EventDataDir[mcIt] + EvtNClist[mcIt]->at(EvtIt) + "/MCTruth_NTUP.root";
            TFile *f_NTUP = new TFile(NTUPname, "READ");
            TTree *t_kinematics = (TTree *)f_NTUP->Get("m_NuMCTruth_tree");
            m_NuMCTruth_tree *NTUP_kinematics = new m_NuMCTruth_tree(t_kinematics);

            //std::cout<<__LINE__<<std::endl;

            RecoMC->SetTrueVertex(jw_test->m_vx, jw_test->m_vy, jw_test->m_vz);
            RecoMC->GetRecoMCInfo(false);
            //RecoMC->GetRecoAngleIPdz();
            RecoMC->SortBasedOnTrackID(false);

            // match tracks reco. and truth
            size_t trueIt=0;
            NTUP_kinematics->GetEntry(trueIt);

            //std::cout<<__LINE__<<std::endl;

            // reco. tracks loop: trackID ascending order
            for(size_t recoIt=0; recoIt<RecoMC->trackID.size(); recoIt++){

                if(RecoMC->IsPrimTrack(recoIt)){
                    // cut by IsPrimTrack function

                    while(NTUP_kinematics->m_track_id < RecoMC->trackID.at(recoIt)){
                        trueIt++;
                        if(trueIt >= NTUP_kinematics->fChain->GetEntries()) break;
                        NTUP_kinematics->GetEntry(trueIt);
                    }

                    if(trueIt >= NTUP_kinematics->fChain->GetEntries()) break;

                    if(NTUP_kinematics->m_track_id == RecoMC->trackID.at(recoIt)){
                        // track ID matches

                        disc->px_reco.push_back(RecoMC->pmag_reco.at(recoIt) * std::sin(RecoMC->theta.at(recoIt)) * std::cos(RecoMC->phi.at(recoIt)));
                        disc->py_reco.push_back(RecoMC->pmag_reco.at(recoIt) * std::sin(RecoMC->theta.at(recoIt)) * std::sin(RecoMC->phi.at(recoIt)));
                        disc->pz_reco.push_back(RecoMC->pmag_reco.at(recoIt) * std::cos(RecoMC->theta.at(recoIt)));
                        disc->pmag_reco.push_back(RecoMC->pmag_reco.at(recoIt));
                        disc->phi_reco.push_back(RecoMC->phi.at(recoIt));
                        disc->theta_reco.push_back(RecoMC->theta.at(recoIt));

                        disc->pmag_ang.push_back(RecoMC->pmag_ang.at(recoIt));
                        disc->pmag_coord.push_back(RecoMC->pmag_coord.at(recoIt));
                        disc->pmag_haruhi.push_back(RecoMC->pmag_haruhi.at(recoIt));
                        
                        TVector3 p3_true(NTUP_kinematics->m_px/1000., NTUP_kinematics->m_py/1000., NTUP_kinematics->m_pz/1000.); // MeV to GeV
                        disc->px_true.push_back(NTUP_kinematics->m_px/1000.);                                                    // MeV to GeV
                        disc->py_true.push_back(NTUP_kinematics->m_py/1000.);                                                    // MeV to GeV
                        disc->pz_true.push_back(NTUP_kinematics->m_pz/1000.);                                                    // MeV to GeV
                        disc->pmag_true.push_back(p3_true.Mag());
                        disc->phi_true.push_back(p3_true.Phi());
                        disc->theta_true.push_back(p3_true.Theta());
                        
                        disc->TrackID.push_back(NTUP_kinematics->m_track_id);
                        disc->PDG.push_back(NTUP_kinematics->m_pdg_id);
                        disc->nseg.push_back(RecoMC->NSeg.at(recoIt));
                        disc->TrackLength.push_back(RecoMC->TrackLength.at(recoIt));
                        disc->dz.push_back(RecoMC->dz.at(recoIt)); // um
                        disc->IP.push_back(RecoMC->IP.at(recoIt)); // um
                        disc->PID_start.push_back(RecoMC->PID.at(recoIt).at(0));
                        disc->PID_end.push_back(RecoMC->PID.at(recoIt).back());
                        disc->itrk.push_back(RecoMC->itrk.at(recoIt));
                        disc->MaxGap.push_back(RecoMC->MaxGap.at(recoIt));
                        disc->theta_RMS.push_back(RecoMC->theta_RMS.at(recoIt));
                        disc->MaxKinkAngle.push_back(RecoMC->MaxKinkAngle.at(recoIt));
                        disc->IsPartCatChanged.push_back(RecoMC->IsPartCatChanged.at(recoIt));
                        disc->IsTrkIdChanged.push_back(RecoMC->IsTrkIdChanged.at(recoIt));

                    }

                }

            }// reco. tracks loop

            //std::cout<<__LINE__<<std::endl;

            // event level info.
            disc->mcID = MonteCarloID[mcIt];
            disc->EventID = jw_test->m_event_id_MC;
            Nu_PDG = jw_test->m_Nu_PDG;
            Nu_px = jw_test->m_Nu_px;
            Nu_py = jw_test->m_Nu_py;
            Nu_pz = jw_test->m_Nu_pz;
            Nu_e = jw_test->m_Nu_e;
            if(mcIt < 3){
                // NC
                FileNum = GenieEvtIDtoFileNum(jw_test->m_event_id_MC, mcIt);
            }
            else{
                // Bkg
                FileNum = ParticleGunEventListToFileNum(EvtNClist[mcIt]->at(EvtIt), mcIt);
            }
            

            // calculate discriminators
            disc->CalcDisc();

            //std::cout<<__LINE__<<std::endl;

            PassDiscAddressToTTreeAddress(disc);
            f_disc->cd();
            if(EvtIt < TrainTestValiRatio*NEvts){
                t_disc[mcIt][0]->Fill();
            } 
            else{
                t_disc[mcIt][1]->Fill();
            }
            
            
            
            //std::cout<<__LINE__<<std::endl;

            f_trfile->Close();
            f_jwtest->Close();
            f_NTUP->Close();

            //std::cout<<__LINE__<<std::endl;

            delete jw_test;
            delete RecoMC;
            delete NTUP_kinematics;
            delete disc;

            //std::cout<<__LINE__<<std::endl;
            
            delete f_trfile;
            delete f_jwtest;
            delete f_NTUP;

            //std::cout<<__LINE__<<std::endl;

        }// Event Loop

    }// MC 200025 200026 200035 Loop

    //std::cout<<__LINE__<<std::endl;

    f_disc->cd();
    for(int mcIt=0; mcIt<4; mcIt++){
        t_disc[mcIt][0]->Write();
        t_disc[mcIt][1]->Write();
    }
    f_disc->Save();
    f_disc->Close();
    delete f_disc;

    //std::cout<<__LINE__<<std::endl;
   
    return 0;

}

