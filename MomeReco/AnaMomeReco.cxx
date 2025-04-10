#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "AnaMomeReco.h"
#include "TPDGCode.h"
#include "discNC_TrainTest.C"
//#include "reco.C"

void AnaMomeReco(){

    TFile *f_disc = new TFile("PhysicsNTUP_ML.root", "READ");
    TTree *t_disc = (TTree *)f_disc->Get("reco");
    discNC_TrainTest *disc = new discNC_TrainTest(t_disc);
    int NrecoNC = disc->fChain->GetEntries();
    std::cout<<"# reco NC events = "<<NrecoNC<<std::endl;

    Hist_Init();

    // Event Loop
    for(int EvtIt=0; EvtIt<NrecoNC; EvtIt++){

        disc->GetEntry(EvtIt);

        if(disc->n_ch < 3) continue;

        // Track Loop
        for(int TrkIt=0; TrkIt<disc->pmag_reco->size(); TrkIt++){

            double p_haruhi = disc->pmag_haruhi->at(TrkIt);
            double p_ang = disc->pmag_ang->at(TrkIt);
            double p_coord = disc->pmag_coord->at(TrkIt);
            double p_true = disc->pmag_true->at(TrkIt);
            double p_reco = disc->pmag_reco->at(TrkIt);

            if(p_reco > 6000){
                // failure 1

                if(p_coord > 0 && p_coord < 1000){
                    p_reco = p_coord;
                }
                else if(p_ang > 0){
                    p_reco = p_ang;
                }
                else{
                    p_reco = 3;
                }

                h1_MaxGap[0]->Fill(disc->MaxGap->at(TrkIt));
                
                if(std::abs(disc->PDG->at(TrkIt)) != 11){
                    // exclude electron
                    h2_thetaRMS_ptrue[2]->Fill(p_true, disc->theta_RMS->at(TrkIt)*1e3);
                    h2_TrackLength_ptrue[2]->Fill(p_true, disc->TrackLength->at(TrkIt));
                    h2_MaxKinkAngle_ptrue[2]->Fill(p_true, disc->MaxKinkAngle->at(TrkIt)*1e3);
                }

                if(p_true > 20){
                    // event display
                    std::cout<<"failure 1 MC = "<<disc->mcID<<"; "
                    <<"Event ID = "<<disc->EventID<<"; "
                    <<"itrk = "<<disc->itrk->at(TrkIt)<<"; "
                    <<"TrackID = "<<disc->TrackID->at(TrkIt)<<"; "
                    <<"PDG = "<<disc->PDG->at(TrkIt)<<"; "
                    <<"p_true = "<<p_true<<" GeV/c; "
                    <<std::endl;
                }
                
            }

            if(p_reco < 4){
                // potential failure 2
                if(std::max(p_coord, p_ang) > p_reco){
                    p_reco = std::max(p_coord, p_ang);
                }
            }

            if(std::abs(disc->PDG->at(TrkIt)) != 11){
                // exclude electron
                h2_preco_ptrue->Fill(disc->pmag_true->at(TrkIt), p_reco);

                if(p_true > 10 && p_haruhi < 4){
                    // failure 2
                    h1_MaxGap[1]->Fill(disc->MaxGap->at(TrkIt));
                    h2_thetaRMS_ptrue[1]->Fill(p_true, disc->theta_RMS->at(TrkIt)*1e3);
                    h2_TrackLength_ptrue[1]->Fill(p_true, disc->TrackLength->at(TrkIt));
                    h2_MaxKinkAngle_ptrue[1]->Fill(p_true, disc->MaxKinkAngle->at(TrkIt)*1e3);

                    if(p_true > 20){
                        // event display
                        std::cout<<"failure 2 MC = "<<disc->mcID<<"; "
                        <<"Event ID = "<<disc->EventID<<"; "
                        <<"itrk = "<<disc->itrk->at(TrkIt)<<"; "
                        <<"TrackID = "<<disc->TrackID->at(TrkIt)<<"; "
                        <<"PDG = "<<disc->PDG->at(TrkIt)<<"; "
                        <<"p_true = "<<p_true<<" GeV/c; "
                        <<std::endl;
                    }
                }
    
                if(p_true > 10 && std::abs(p_true-p_reco) < 0.30*p_true){
                    // success
                    h2_thetaRMS_ptrue[0]->Fill(p_true, disc->theta_RMS->at(TrkIt)*1e3);
                    h2_TrackLength_ptrue[0]->Fill(p_true, disc->TrackLength->at(TrkIt));
                    h2_MaxKinkAngle_ptrue[0]->Fill(p_true, disc->MaxKinkAngle->at(TrkIt)*1e3);
                }
            }
            

            if(std::abs(disc->PDG->at(TrkIt)) == 11){
                // electron 
                h2_preco_ptrue_cat[0]->Fill(p_true, p_reco);
                h1_IP[2]->Fill(disc->IP->at(TrkIt));
                h2_thetaRMS_ptrue[3]->Fill(p_true, disc->theta_RMS->at(TrkIt)*1e3);
                h2_TrackLength_ptrue[3]->Fill(p_true, disc->TrackLength->at(TrkIt));
                h2_MaxKinkAngle_ptrue[3]->Fill(p_true, disc->MaxKinkAngle->at(TrkIt)*1e3);
                h2_ptrue_dz[2]->Fill(disc->dz->at(TrkIt), p_true);
            }
            else if(std::abs(disc->PDG->at(TrkIt)) == 13){
                h2_preco_ptrue_cat[1]->Fill(p_true, p_reco);
            }
            else if(std::abs(disc->PDG->at(TrkIt)) == 211){
                h2_preco_ptrue_cat[2]->Fill(p_true, p_reco);
            }
            else if(std::abs(disc->PDG->at(TrkIt)) == 321){
                h2_preco_ptrue_cat[3]->Fill(p_true, p_reco);
            }
            else if(std::abs(disc->PDG->at(TrkIt)) == 2212){
                h2_preco_ptrue_cat[4]->Fill(p_true, p_reco);
            }

            if(disc->TrackID->at(TrkIt) < 20000){
                h1_IP[0]->Fill(disc->IP->at(TrkIt));
                h2_IP_dz[0]->Fill(disc->dz->at(TrkIt), disc->IP->at(TrkIt));
                h2_ptrue_dz[0]->Fill(disc->dz->at(TrkIt), p_true);
            }
            else{
                h1_IP[1]->Fill(disc->IP->at(TrkIt));
                h2_IP_dz[1]->Fill(disc->dz->at(TrkIt), disc->IP->at(TrkIt));
                h2_ptrue_dz[1]->Fill(disc->dz->at(TrkIt), p_true);
            }


        } // Track Loop
    
    }// Event Loop

    StoreHist2ROOT();

}
