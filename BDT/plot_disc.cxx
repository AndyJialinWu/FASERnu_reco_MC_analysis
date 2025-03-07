#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "discNC_TrainTest.C"
#include "plot_disc.h"

void plot_disc(){


    TFile *f_disc = new TFile("PhysicsNTUP_ML.root", "READ");
    TTree *t_disc[2]; // 0: NC, 1: Bkg
    t_disc[0] = (TTree *)f_disc->Get("discNC_TrainTest");
    t_disc[1] = (TTree *)f_disc->Get("discBkg_TrainTest");
    discNC_TrainTest *disc[2]; // 0: NC, 1: Bkg
    disc[0] = new discNC_TrainTest(t_disc[0]);
    disc[1] = new discNC_TrainTest(t_disc[1]);

    wgt_calc();
    Hist_Init();

    // signal: NC, background: Bkg Loop
    for(int i=0; i<2; i++){

        for(int EvtIt=0; EvtIt<disc[i]->fChain->GetEntries(); EvtIt++){

            disc[i]->GetEntry(EvtIt);
    
            if(disc[i]->n_ch < 3) continue;
    
            int wgt_idx;
            if(disc[i]->mcID == 200025) wgt_idx = 0;
            if(disc[i]->mcID == 200026) wgt_idx = 1;
            if(disc[i]->mcID == 200035) wgt_idx = 2;
            if (disc[i]->mcID == 100069) wgt_idx = 3;
    
            // Track Loop
            for(int TrkIt=0; TrkIt<disc[i]->pmag_reco->size(); TrkIt++){
    
                h2_preco_ptrue->Fill(disc[i]->pmag_true->at(TrkIt), disc[i]->pmag_reco->at(TrkIt));    
    
            }   
    
            h1_n_ch[i]->Fill(disc[i]->n_ch, weight[wgt_idx]);
    
            h1_dphi_max_deg[i][0]->Fill(disc[i]->dphi_max_reco, weight[wgt_idx]);
            h1_dphi_max_deg[i][1]->Fill(disc[i]->dphi_max_true, weight[wgt_idx]);
            
            h1_dphi_sum_deg[i][0]->Fill(disc[i]->dphi_sum_reco, weight[wgt_idx]);
            h1_dphi_sum_deg[i][1]->Fill(disc[i]->dphi_sum_true, weight[wgt_idx]);     
    
            h1_DeltaPhiMET_deg[i][0]->Fill(disc[i]->DeltaPhiMET_reco, weight[wgt_idx]);
            h1_DeltaPhiMET_deg[i][1]->Fill(disc[i]->DeltaPhiMET_true, weight[wgt_idx]);
    
            h1_pmag_had_vis[i][0]->Fill(disc[i]->pmag_had_vis_reco, weight[wgt_idx]);
            h1_pmag_had_vis[i][1]->Fill(disc[i]->pmag_had_vis_true, weight[wgt_idx]);
    
            h1_p_had_hardest[i][0]->Fill(disc[i]->p3_hardest_reco, weight[wgt_idx]);
            h1_p_had_hardest[i][1]->Fill(disc[i]->p3_hardest_true, weight[wgt_idx]);
    
            h1_pTmiss_mag[i][0]->Fill(disc[i]->pTmiss_mag_reco, weight[wgt_idx]);
            h1_pTmiss_mag[i][1]->Fill(disc[i]->pTmiss_mag_true, weight[wgt_idx]);
            
            h1_pTabs_sum[i][0]->Fill(disc[i]->pTabs_sum_reco, weight[wgt_idx]);
            h1_pTabs_sum[i][1]->Fill(disc[i]->pTabs_sum_true, weight[wgt_idx]);
    
            h1_InvThetaCh[i][0]->Fill(disc[i]->InvThetaCh_reco, weight[wgt_idx]);
            h1_InvThetaCh[i][1]->Fill(disc[i]->InvThetaCh_true, weight[wgt_idx]);
    
            h1_tan_theta_hardest[i][0]->Fill(disc[i]->tan_theta_hardest_reco, weight[wgt_idx]);
            h1_tan_theta_hardest[i][1]->Fill(disc[i]->tan_theta_hardest_true, weight[wgt_idx]);
    
    
        }// end of Track Loop

    }// end of signal: NC, background: Bkg Loop
    // Event Loop
    

    StoreHist2ROOT();

    TString FigAddress = "Figures/discriminators/";
    StoreHist2PDF(FigAddress);

    f_disc->Close();
    delete f_disc;

}
