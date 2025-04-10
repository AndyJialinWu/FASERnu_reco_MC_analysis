#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "reco.C"
#include "plot_disc.h"

void plot_disc(){


    TFile *f_disc = new TFile("PhysicsNTUP.root", "READ");
    TTree *t_disc = (TTree *)f_disc->Get("reco");
    reco *disc = new reco(t_disc);
    int NrecoNC = disc->fChain->GetEntries();
    std::cout<<"# reco NC events = "<<NrecoNC<<std::endl;

    wgt_calc();
    Hist_Init();

    // Event Loop
    for(int EvtIt=0; EvtIt<NrecoNC; EvtIt++){

        disc->GetEntry(EvtIt);

        if(disc->n_ch < 3) continue;

        int wgt_idx;
        if(disc->mcID == 200025) wgt_idx = 0;
        if(disc->mcID == 200026) wgt_idx = 1;
        if(disc->mcID == 200035) wgt_idx = 2;

        // Track Loop
        for(int TrkIt=0; TrkIt<disc->pmag_reco->size(); TrkIt++){

            h2_preco_ptrue->Fill(disc->pmag_true->at(TrkIt), disc->pmag_reco->at(TrkIt));
            
            int nplates = disc->PID_end->at(TrkIt) - disc->PID_start->at(TrkIt) + 1;
            double TrkEff = double(disc->nseg->at(TrkIt)) / double(nplates);

            if(disc->pmag_reco->at(TrkIt) > 6000){ 
                // failure 1
/*
                std::cout<<"P_reco = 7000 GeV/c; P_true = "<<disc->pmag_true->at(TrkIt)<<" GeV/c; "
                <<"P_coord = "<<disc->pmag_coord->at(TrkIt)<<" GeV/c; "
                <<"P_ang = "<<disc->pmag_ang->at(TrkIt)<<" GeV/c; "
                <<"# seg. = "<<disc->nseg->at(TrkIt)<<"; "
                <<"# plates = "<<nplates<<"; "
                <<"Track Eff = "<<TrkEff<<"; "
                <<"Track ID = "<<disc->TrackID->at(TrkIt)<<"; "
                //<<"TraLen = "<<disc->TrackLength->at(TrkIt)<<" mm; "
                <<"PDG = "<<disc->PDG->at(TrkIt)<<"; "
                <<std::endl;
*/              
                h2_ptrue_Nplates[0]->Fill(nplates, disc->pmag_true->at(TrkIt));
                h2_TrkEff_Nplates[0]->Fill(nplates, TrkEff);
                h2_TrkEff_ptrue[0]->Fill(disc->pmag_true->at(TrkIt), TrkEff);
            }

            if(disc->pmag_reco->at(TrkIt) < 4 && disc->pmag_true->at(TrkIt) > 10){
                // failure 2
                h2_ptrue_Nplates[1]->Fill(nplates, disc->pmag_true->at(TrkIt));
                h2_TrkEff_Nplates[1]->Fill(nplates, TrkEff);
                h2_TrkEff_ptrue[1]->Fill(disc->pmag_true->at(TrkIt), TrkEff);
            }

        }

        if(disc->p3_hardest_reco > 1000) continue;
        if(disc->pTabs_sum_reco > 40) continue;
        if(disc->pTmiss_mag_reco > 40) continue;

        h1_n_ch->Fill(disc->n_ch, weight[wgt_idx]);

        h1_dphi_max_deg[0]->Fill(disc->dphi_max_reco, weight[wgt_idx]);
        h1_dphi_max_deg[1]->Fill(disc->dphi_max_true, weight[wgt_idx]);
        
        h1_dphi_sum_deg[0]->Fill(disc->dphi_sum_reco, weight[wgt_idx]);
        h1_dphi_sum_deg[1]->Fill(disc->dphi_sum_true, weight[wgt_idx]);     

        h1_DeltaPhiMET_deg[0]->Fill(disc->DeltaPhiMET_reco, weight[wgt_idx]);
        h1_DeltaPhiMET_deg[1]->Fill(disc->DeltaPhiMET_true, weight[wgt_idx]);

        h1_pmag_had_vis[0]->Fill(disc->pmag_had_vis_reco, weight[wgt_idx]);
        h1_pmag_had_vis[1]->Fill(disc->pmag_had_vis_true, weight[wgt_idx]);

        h1_p_had_hardest[0]->Fill(disc->p3_hardest_reco, weight[wgt_idx]);
        h1_p_had_hardest[1]->Fill(disc->p3_hardest_true, weight[wgt_idx]);

        h1_pTmiss_mag[0]->Fill(disc->pTmiss_mag_reco, weight[wgt_idx]);
        h1_pTmiss_mag[1]->Fill(disc->pTmiss_mag_true, weight[wgt_idx]);
        
        h1_pTabs_sum[0]->Fill(disc->pTabs_sum_reco, weight[wgt_idx]);
        h1_pTabs_sum[1]->Fill(disc->pTabs_sum_true, weight[wgt_idx]);

        h1_InvThetaCh[0]->Fill(disc->InvThetaCh_reco, weight[wgt_idx]);
        h1_InvThetaCh[1]->Fill(disc->InvThetaCh_true, weight[wgt_idx]);

        h1_tan_theta_hardest[0]->Fill(disc->tan_theta_hardest_reco, weight[wgt_idx]);
        h1_tan_theta_hardest[1]->Fill(disc->tan_theta_hardest_true, weight[wgt_idx]);


    }

    StoreHist2ROOT();

    TString FigAddress = "Figures/discriminators/";
    StoreHist2PDF(FigAddress);

    f_disc->Close();
    delete f_disc;

}
