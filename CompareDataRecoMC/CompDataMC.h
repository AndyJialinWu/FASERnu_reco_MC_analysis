#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"

#include "PhysicsNTUP.C"

// 0: MC, 1: Data
TH1D *h1_n_ch[2];
TH1D *h1_nTrkTanThetaLeq0point1[2];
TH1D *h1_pmag_had_vis[2];
TH1D *h1_dphi_max_deg[2];
TH1D *h1_p_had_hardest[2];
TH1D *h1_InvThetaCh[2];
TH1D *h1_pTmiss_mag[2];
TH1D *h1_pTabs_sum[2];
TH1D *h1_DeltaPhiMET_deg[2];
TH1D *h1_dphi_sum_deg[2];
TH1D *h1_tan_theta_hardest[2];

TString str_cat[2] = {"MC", "Data"};

void HistInit(){

    // MC Data Loop
    for(int i=0; i<2; i++){

        h1_n_ch[i] = new TH1D("h1_n_ch_"+str_cat[i], "Charged Tracks Multiplicity", 40, 0, 40);
        h1_n_ch[i]->GetXaxis()->SetTitle("n_{ch}");
        h1_n_ch[i]->GetYaxis()->SetTitle("# Events");

        h1_dphi_max_deg[i] = new TH1D("h1_dphi_max_deg_"+str_cat[i], "Largest Neighbouring Azimuthal Gap", 50, 0, 200);
        h1_dphi_max_deg[i]->GetXaxis()->SetTitle("#Delta #phi_{max} [#circ]");
        h1_dphi_max_deg[i]->GetYaxis()->SetTitle("# Events");

        h1_pmag_had_vis[i] = new TH1D("h1_pmag_had_vis_"+str_cat[i], "Sum of Visible Hadronic Momentum Magnitude", 140, 0, 1400);
        h1_pmag_had_vis[i]->GetXaxis()->SetTitle("|#vec{p}_{had, v}| [GeV/c]");
        h1_pmag_had_vis[i]->GetYaxis()->SetTitle("# Events");

        h1_p_had_hardest[i] = new TH1D("h1_p_had_hardest_"+str_cat[i], "Hardest Track Momentum", 100, 0, 1000);
        h1_p_had_hardest[i]->GetXaxis()->SetTitle("|#vec{p}_{hardest}| [GeV/c]");
        h1_p_had_hardest[i]->GetYaxis()->SetTitle("# Events");

        h1_InvThetaCh[i] = new TH1D("h1_InvThetaCh_"+str_cat[i], "Inverse Sum of Charged Hadron Track Polar Angles", 100, 0, 2000);
        h1_InvThetaCh[i]->GetXaxis()->SetTitle("#sum |1/#theta_{ch}| [rad^{-1}]");
        h1_InvThetaCh[i]->GetYaxis()->SetTitle("# Events");

        h1_DeltaPhiMET_deg[i] = new TH1D("h1_DeltaPhiMET_deg_"+str_cat[i], "Track-MET-Angle", 50, 0, 200);
        h1_DeltaPhiMET_deg[i]->GetXaxis()->SetTitle("#Delta #phi_{MET} [#circ]");
        h1_DeltaPhiMET_deg[i]->GetYaxis()->SetTitle("# Events");

        h1_dphi_sum_deg[i] = new TH1D("h1_dphi_sum_deg_"+str_cat[i], "Sum of Neighbouring Azimuthal Gaps", 90, 0, 360);
        h1_dphi_sum_deg[i]->GetXaxis()->SetTitle("#sum_{i} #Delta #phi_{i} [#circ]");
        h1_dphi_sum_deg[i]->GetYaxis()->SetTitle("# Events");

        h1_pTmiss_mag[i] = new TH1D("h1_pTmiss_mag_"+str_cat[i], "Missing Transverse Momentum", 40, 0, 40);
        h1_pTmiss_mag[i]->GetXaxis()->SetTitle("|#vec{p}_{T, miss}| [GeV/c]");
        h1_pTmiss_mag[i]->GetYaxis()->SetTitle("# Events");

        h1_pTabs_sum[i] = new TH1D("h1_pTabs_sum_"+str_cat[i], "Sum of Magnitudes of Transverse Momenta", 40, 0, 40);
        h1_pTabs_sum[i]->GetXaxis()->SetTitle("#sum |#vec{p}_{T}| [GeV/c]");
        h1_pTabs_sum[i]->GetYaxis()->SetTitle("# Events");

        h1_tan_theta_hardest[i] = new TH1D("h1_tan_theta_hardest_"+str_cat[i], "tan#theta of Hardest Hadronic Track", 100, 0, 0.5);
        h1_tan_theta_hardest[i]->GetXaxis()->SetTitle("tan#theta_{hardest}");
        h1_tan_theta_hardest[i]->GetYaxis()->SetTitle("# Events");

        h1_nTrkTanThetaLeq0point1[i] = new TH1D("h1_nTrkTanThetaLeq0point1_"+str_cat[i], "Number of Tracks with tan#theta < 0.1", 20, 0, 20);
        h1_nTrkTanThetaLeq0point1[i]->GetXaxis()->SetTitle("# tracks with tan#theta < 0.1");
        h1_nTrkTanThetaLeq0point1[i]->GetYaxis()->SetTitle("# Events");

    } // End of MC Data Loop

}


void HistFill(PhysicsNTUP *ntup, int i){

    // Fill histograms
    h1_n_ch[i]->Fill(ntup->n_ch);
    h1_dphi_max_deg[i]->Fill(ntup->dphi_max_reco);
    h1_pmag_had_vis[i]->Fill(ntup->pmag_had_vis_reco);
    h1_p_had_hardest[i]->Fill(ntup->p3_hardest_reco);
    h1_InvThetaCh[i]->Fill(ntup->InvThetaCh_reco);
    h1_DeltaPhiMET_deg[i]->Fill(ntup->DeltaPhiMET_reco);
    h1_dphi_sum_deg[i]->Fill(ntup->dphi_sum_reco);
    h1_pTmiss_mag[i]->Fill(ntup->pTmiss_mag_reco);
    h1_pTabs_sum[i]->Fill(ntup->pTabs_sum_reco);
    h1_tan_theta_hardest[i]->Fill(ntup->tan_theta_hardest_reco);
    h1_nTrkTanThetaLeq0point1[i]->Fill(ntup->nTrkTanThetaLeq0point1_reco);

} // End of HistFill


void StoreHist2ROOT(TString RootFileName){

    TFile *f_plot_disc = new TFile("Figures/" + RootFileName, "RECREATE");
    f_plot_disc->cd();

    for(int i=0; i<2; i++){
        
        h1_n_ch[i]->Write();
        h1_dphi_max_deg[i]->Write();
        h1_pmag_had_vis[i]->Write();
        h1_p_had_hardest[i]->Write();
        h1_InvThetaCh[i]->Write();
        h1_DeltaPhiMET_deg[i]->Write();
        h1_dphi_sum_deg[i]->Write();
        h1_pTmiss_mag[i]->Write();
        h1_pTabs_sum[i]->Write();
        h1_tan_theta_hardest[i]->Write();
        h1_nTrkTanThetaLeq0point1[i]->Write();

    }

    f_plot_disc->Save();
    f_plot_disc->Close();
    delete f_plot_disc;

};

TCanvas *cvs;
TLegend *lg;
void CvsLg_Init(){
    // initialize the canvas and legend
    cvs = new TCanvas("canvas", "canvas", 800, 600);
    lg = new TLegend(.7, .7, .9, .9);
    lg->SetFillStyle(0);
};

void DrawDataMC(TCanvas *cvs, TLegend *lg, TH1D *h1_mc, TH1D *h1_data, TString fig_name, TString category){

    cvs->Clear();
    lg->Clear();

    cvs->cd();
    h1_mc->SetStats(false);
    h1_data->SetStats(false);
    h1_mc->SetLineColor(kRed);
    h1_data->SetLineColor(kBlue);

    // scale the mc to the data
    double scale = h1_data->Integral() / h1_mc->Integral();
    h1_mc->Scale(scale);
    
    if(h1_mc->GetMaximum() < h1_data->GetMaximum()){
        h1_data->DrawCopy("e1");
        h1_mc->Draw("HIST same");
    }
    else{
        h1_mc->Draw("HIST");
        h1_data->DrawCopy("e1 same");
    }
    
    lg->SetHeader(category);
    lg->AddEntry(h1_mc, "MC Neutral Hadrons");
    lg->AddEntry(h1_data, "FeedBack Data");
    lg->Draw("same");

    cvs->SetGrid();
    cvs->SaveAs(fig_name);

}

void StoreHist2PDF(TString FigAddress){

    CvsLg_Init();

    TString category = "BDT Input";
    DrawDataMC(cvs, lg, h1_n_ch[0], h1_n_ch[1], FigAddress+"h1_n_ch.pdf", category);
    DrawDataMC(cvs, lg, h1_dphi_max_deg[0], h1_dphi_max_deg[1], FigAddress+"h1_dphi_max_deg.pdf", category);
    DrawDataMC(cvs, lg, h1_dphi_sum_deg[0], h1_dphi_sum_deg[1], FigAddress+"h1_dphi_sum_deg.pdf", category);
    DrawDataMC(cvs, lg, h1_DeltaPhiMET_deg[0], h1_DeltaPhiMET_deg[1], FigAddress+"h1_DeltaPhiMET_deg.pdf", category);
    DrawDataMC(cvs, lg, h1_pmag_had_vis[0], h1_pmag_had_vis[1], FigAddress+"h1_pmag_had_vis.pdf", category);
    DrawDataMC(cvs, lg, h1_p_had_hardest[0], h1_p_had_hardest[1], FigAddress+"h1_p_had_hardest.pdf", category);
    DrawDataMC(cvs, lg, h1_InvThetaCh[0], h1_InvThetaCh[1], FigAddress+"h1_InvThetaCh.pdf", category);
    DrawDataMC(cvs, lg, h1_pTmiss_mag[0], h1_pTmiss_mag[1], FigAddress+"h1_pTmiss_mag.pdf", category);
    DrawDataMC(cvs, lg, h1_pTabs_sum[0], h1_pTabs_sum[1], FigAddress+"h1_pTabs_sum.pdf", category);
    DrawDataMC(cvs, lg, h1_tan_theta_hardest[0], h1_tan_theta_hardest[1], FigAddress+"h1_tan_theta_hardest.pdf", category);
    DrawDataMC(cvs, lg, h1_nTrkTanThetaLeq0point1[0], h1_nTrkTanThetaLeq0point1[1], FigAddress+"h1_nTrkTanThetaLeq0point1.pdf", category);

}

