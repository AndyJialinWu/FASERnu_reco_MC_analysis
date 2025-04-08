#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "plot_disc.h"
#include "PhysicsNTUP.C"

void plot_disc(){


    TFile *f_disc = new TFile("PhysicsNTUP_ML.root", "READ");
    TTree *t_TrainTest[4];                       // 0: 200025, 1: 200026, 2: 200035, 3: 100069
    PhysicsNTUP *disc[4];                        // 0: 200025, 1: 200026, 2: 200035, 3: 100069
    for(int RunIt=0; RunIt<4; RunIt++){

        t_TrainTest[RunIt] = (TTree*)f_disc->Get( Form("disc%d_TrainTest", RunNumber[RunIt]) );
        SetNumRecoNC(t_TrainTest[RunIt], RunIt); // sequence matters ! Place this function before setting the branch address !

        disc[RunIt] = new PhysicsNTUP(t_TrainTest[RunIt]);

    }

    wgt_calc();
    Hist_Init();
    //std::ofstream em_shower;
    //em_shower.open("HighErgElectron.txt");

    // MC Run Loop
    for(int RunIt=0; RunIt<4; RunIt++){

        int SigBkgIdx = -1;
        if(RunIt < 3){// nu NC from 0: 200025, 1: 200026, 2: 200035
            SigBkgIdx = 0;
        }
        else{         // Bkg from 3: 100069
            SigBkgIdx = 1;
        }

        // Event Loop
        for(int EvtIt=0; EvtIt<disc[RunIt]->fChain->GetEntries(); EvtIt++){

            disc[RunIt]->GetEntry(EvtIt);   // Get the event

            h1_n_ch[SigBkgIdx]->Fill(double(disc[RunIt]->n_ch), weight[RunIt]);
            //std::cout<<"RunNumber = "<<RunNumber[RunIt]<<"; n_ch = "<<disc[RunIt]->n_ch<<"; weight = "<<weight[RunIt]<<std::endl;
            
            // Track Loop
            for(int TrkIt=0; TrkIt<disc[RunIt]->pmag_reco->size(); TrkIt++){
                
                const int PDG = disc[RunIt]->PDG->at(TrkIt);
                if( std::abs(PDG) != kElectron ){// hadron
                    h2_preco_ptrue[0]->Fill(disc[RunIt]->pmag_true->at(TrkIt), disc[RunIt]->pmag_reco->at(TrkIt), weight[RunIt]);
                }
                else{// electron
                    h2_preco_ptrue[1]->Fill(disc[RunIt]->pmag_true->at(TrkIt), disc[RunIt]->pmag_reco->at(TrkIt), weight[RunIt]);

                    //if(disc[RunIt]->pmag_true->at(TrkIt) > 50){ // p_e_true > 50 GeV/c
                    //    em_shower<<"RunNumber = "<<disc[RunIt]->mcID<<"; FileNum = "<<*(disc[RunIt]->FileNum)<<"; EventID = "<<disc[RunIt]->EventID<<"\n";
                    //}
                }
                
                Fill_h2_PartCat_n_ch(PDG, disc[RunIt]->n_ch, SigBkgIdx, weight[RunIt]);
    
            }// end of track loop

            if(disc[RunIt]->n_ch < 3) continue;
    
            h1_dphi_max_deg[SigBkgIdx][0]->Fill(disc[RunIt]->dphi_max_reco, weight[RunIt]);
            h1_dphi_max_deg[SigBkgIdx][1]->Fill(disc[RunIt]->dphi_max_true, weight[RunIt]);
            
            h1_dphi_sum_deg[SigBkgIdx][0]->Fill(disc[RunIt]->dphi_sum_reco, weight[RunIt]);
            h1_dphi_sum_deg[SigBkgIdx][1]->Fill(disc[RunIt]->dphi_sum_true, weight[RunIt]);     
    
            h1_DeltaPhiMET_deg[SigBkgIdx][0]->Fill(disc[RunIt]->DeltaPhiMET_reco, weight[RunIt]);
            h1_DeltaPhiMET_deg[SigBkgIdx][1]->Fill(disc[RunIt]->DeltaPhiMET_true, weight[RunIt]);
    
            h1_pmag_had_vis[SigBkgIdx][0]->Fill(disc[RunIt]->pmag_had_vis_reco, weight[RunIt]);
            h1_pmag_had_vis[SigBkgIdx][1]->Fill(disc[RunIt]->pmag_had_vis_true, weight[RunIt]);
    
            h1_p_had_hardest[SigBkgIdx][0]->Fill(disc[RunIt]->p3_hardest_reco, weight[RunIt]);
            h1_p_had_hardest[SigBkgIdx][1]->Fill(disc[RunIt]->p3_hardest_true, weight[RunIt]);
    
            h1_pTmiss_mag[SigBkgIdx][0]->Fill(disc[RunIt]->pTmiss_mag_reco, weight[RunIt]);
            h1_pTmiss_mag[SigBkgIdx][1]->Fill(disc[RunIt]->pTmiss_mag_true, weight[RunIt]);
            
            h1_pTabs_sum[SigBkgIdx][0]->Fill(disc[RunIt]->pTabs_sum_reco, weight[RunIt]);
            h1_pTabs_sum[SigBkgIdx][1]->Fill(disc[RunIt]->pTabs_sum_true, weight[RunIt]);
    
            h1_InvThetaCh[SigBkgIdx][0]->Fill(disc[RunIt]->InvThetaCh_reco, weight[RunIt]);
            h1_InvThetaCh[SigBkgIdx][1]->Fill(disc[RunIt]->InvThetaCh_true, weight[RunIt]);
    
            h1_tan_theta_hardest[SigBkgIdx][0]->Fill(disc[RunIt]->tan_theta_hardest_reco, weight[RunIt]);
            h1_tan_theta_hardest[SigBkgIdx][1]->Fill(disc[RunIt]->tan_theta_hardest_true, weight[RunIt]);

            h1_nTrkTanThetaLeq0point1[SigBkgIdx][0]->Fill(disc[RunIt]->nTrkTanThetaLeq0point1_reco, weight[RunIt]);
            h1_nTrkTanThetaLeq0point1[SigBkgIdx][1]->Fill(disc[RunIt]->nTrkTanThetaLeq0point1_true, weight[RunIt]);

            h2_pTmissReco_pTmissTrue[SigBkgIdx]->Fill(disc[RunIt]->pTmiss_mag_true, disc[RunIt]->pTmiss_mag_reco, weight[RunIt]);
    
        }// end of event loop

    }// end of MC Run Loop
    
    //em_shower.close();

    StoreHist2ROOT("disc_plots.root");

    TString FigAddress = "Figures/discriminators/";
    StoreHist2PDF(FigAddress);

    f_disc->Close();
    delete f_disc;

}
