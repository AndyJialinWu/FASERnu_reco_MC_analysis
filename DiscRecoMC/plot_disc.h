#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Rtypes.h"

const double IntLumi_run3 = 250.00;                    // fb^{-1}
const double eff_FV = 638.00 / 730.00;
const double wgt_NC_200025 = 50173.00 * eff_FV / 10e3; // unit: event/fb^{-1}
const double wgt_NC_200026 = 66355.00 * eff_FV / 10e3;
const double wgt_NC_200035 = 13555.00 * eff_FV / 10e3;
const double Run3ExpNC[3] = {IntLumi_run3 * wgt_NC_200025, IntLumi_run3 * wgt_NC_200026, IntLumi_run3 * wgt_NC_200035};
const int NumRecoNC[3] = {8073, 4221, 1179};

double weight[3];
void wgt_calc(){
    for(int mcIt=0; mcIt<3; mcIt++){
        weight[mcIt] = Run3ExpNC[mcIt]/NumRecoNC[mcIt];
    }
};

const float t_plate = 1.1e3;                  // tungsten plate thickness in um
const float t_film = 0.34e3;                  // film thickness in um
const int kMaxTrks = 79;

TH2D *h2_preco_ptrue;
TH2D *h2_TrkEff_Nplates[2]; // 0: failure 1; 1: failure 2;
TH2D *h2_ptrue_Nplates[2];
TH2D *h2_TrkEff_ptrue[2];
const TString str_failure[2] = {"failure1", "failure2"};

TH1D *h1_n_ch;
// 0:reco, 1:truth
TH1D *h1_dphi_max_deg[2]; 
TH1D *h1_dphi_sum_deg[2];
TH1D *h1_DeltaPhiMET_deg[2];

TH1D *h1_pmag_had_vis[2];
TH1D *h1_p_had_hardest[2];
TH1D *h1_pTmiss_mag[2];
TH1D *h1_pTabs_sum[2];

TH1D *h1_InvThetaCh[2];
TH1D *h1_tan_theta_hardest[2];

const TString str_cat[2] = {"reco", "true"};
const int LineColor[2] = {kBlue, kRed};

void Hist_Init(){

    h2_preco_ptrue = new TH2D("h2_preco_ptrue", "P_{reco} versus P_{true}", 1000, 0, 1000, 8000, 0, 8000);
    h2_preco_ptrue->GetXaxis()->SetTitle("P_{true} [GeV/c]");
    h2_preco_ptrue->GetYaxis()->SetTitle("P_{reco} [GeV/c]");
        
    h1_n_ch = new TH1D("h1_n_ch", "Charged Tracks Multiplicity", 40, 0, 40);
    h1_n_ch->GetXaxis()->SetTitle("n_{ch}");
    h1_n_ch->GetYaxis()->SetTitle("# Events");

    for(int i=0; i<2; i++){

        h2_ptrue_Nplates[i] = new TH2D("h2_ptrue_Nplates_"+str_failure[i], "P_{true} versus # plates of a Track "+str_failure[i], 100, 0, 100, 200, 0, 200);
        h2_ptrue_Nplates[i]->GetXaxis()->SetTitle("# plates");
        h2_ptrue_Nplates[i]->GetYaxis()->SetTitle("P_{true} [GeV/c]");

        h2_TrkEff_Nplates[i] = new TH2D("h2_TrkEff_Nplates_"+str_failure[i], "Track Efficiency versus # plates of a Track "+str_failure[i], 100, 0, 100, 10, 0, 1);
        h2_TrkEff_Nplates[i]->GetXaxis()->SetTitle("# plates");
        h2_TrkEff_Nplates[i]->GetYaxis()->SetTitle("Track Efficiency");

        h2_TrkEff_ptrue[i] = new TH2D("h2_TrkEff_ptrue_"+str_failure[i], "Track Efficiency versus P_{true} "+str_failure[i], 200, 0, 200, 10, 0, 1);
        h2_TrkEff_ptrue[i]->GetXaxis()->SetTitle("P_{true} [GeV/c]");
        h2_TrkEff_ptrue[i]->GetYaxis()->SetTitle("Track Efficiency");

        h1_dphi_max_deg[i] = new TH1D("h1_dphi_max_deg_"+str_cat[i], "Largest Neighbouring Azimuthal Gap", 50, 0, 200);
        h1_dphi_max_deg[i]->GetXaxis()->SetTitle("#Delta #phi_{max} [#circ]");
        h1_dphi_max_deg[i]->GetYaxis()->SetTitle("# Events");
        h1_dphi_max_deg[i]->SetLineColor(LineColor[i]);

        h1_pmag_had_vis[i] = new TH1D("h1_pmag_had_vis_"+str_cat[i], "Sum of Visible Hadronic Momentum Magnitude", 140, 0, 1400);
        h1_pmag_had_vis[i]->GetXaxis()->SetTitle("|#vec{p}_{had, v}| [GeV/c]");
        h1_pmag_had_vis[i]->GetYaxis()->SetTitle("# Events");
        h1_pmag_had_vis[i]->SetLineColor(LineColor[i]);

        h1_p_had_hardest[i] = new TH1D("h1_p_had_hardest_"+str_cat[i], "Hardest Track Momentum", 100, 0, 1000);
        h1_p_had_hardest[i]->GetXaxis()->SetTitle("|#vec{p}_{hardest}| [GeV/c]");
        h1_p_had_hardest[i]->GetYaxis()->SetTitle("# Events");
        h1_p_had_hardest[i]->SetLineColor(LineColor[i]);

        h1_InvThetaCh[i] = new TH1D("h1_InvThetaCh_"+str_cat[i], "Inverse Sum of Charged Hadron Track Polar Angles", 100, 0, 2000);
        h1_InvThetaCh[i]->GetXaxis()->SetTitle("#sum |1/#theta_{ch}| [rad^{-1}]");
        h1_InvThetaCh[i]->GetYaxis()->SetTitle("# Events");
        h1_InvThetaCh[i]->SetLineColor(LineColor[i]);

        h1_DeltaPhiMET_deg[i] = new TH1D("h1_DeltaPhiMET_deg_"+str_cat[i], "Track-MET-Angle", 50, 0, 200);
        h1_DeltaPhiMET_deg[i]->GetXaxis()->SetTitle("#Delta #phi_{MET} [#circ]");
        h1_DeltaPhiMET_deg[i]->GetYaxis()->SetTitle("# Events");
        h1_DeltaPhiMET_deg[i]->SetLineColor(LineColor[i]);

        h1_dphi_sum_deg[i] = new TH1D("h1_dphi_sum_deg_"+str_cat[i], "Sum of Neighbouring Azimuthal Gaps", 90, 0, 360);
        h1_dphi_sum_deg[i]->GetXaxis()->SetTitle("#sum_{i} #Delta #phi_{i} [#circ]");
        h1_dphi_sum_deg[i]->GetYaxis()->SetTitle("# Events");
        h1_dphi_sum_deg[i]->SetLineColor(LineColor[i]);

        h1_pTmiss_mag[i] = new TH1D("h1_pTmiss_mag_"+str_cat[i], "Missing Transverse Momentum", 40, 0, 40);
        h1_pTmiss_mag[i]->GetXaxis()->SetTitle("|#vec{p}_{T, miss}| [GeV/c]");
        h1_pTmiss_mag[i]->GetYaxis()->SetTitle("# Events");
        h1_pTmiss_mag[i]->SetLineColor(LineColor[i]);

        h1_pTabs_sum[i] = new TH1D("h1_pTabs_sum_"+str_cat[i], "Sum of Magnitudes of Transverse Momenta", 40, 0, 40);
        h1_pTabs_sum[i]->GetXaxis()->SetTitle("#sum |#vec{p}_{T}| [GeV/c]");
        h1_pTabs_sum[i]->GetYaxis()->SetTitle("# Events");
        h1_pTabs_sum[i]->SetLineColor(LineColor[i]);

        h1_tan_theta_hardest[i] = new TH1D("h1_tan_theta_hardest_"+str_cat[i], "tan#theta of Hardest Hadronic Track", 100, 0, 0.5);
        h1_tan_theta_hardest[i]->GetXaxis()->SetTitle("tan#theta_{hardest}");
        h1_tan_theta_hardest[i]->GetYaxis()->SetTitle("# Events");
        h1_tan_theta_hardest[i]->SetLineColor(LineColor[i]);

    }

};


void StoreHist2ROOT(){

    TFile *f_plot_disc = new TFile("disc_plots.root", "RECREATE");
    f_plot_disc->cd();

    h2_preco_ptrue->Write();
    h1_n_ch->Write();

    for(int i=0; i<2; i++){
        h1_dphi_max_deg[i]->Write();
        h1_pmag_had_vis[i]->Write();
        h1_p_had_hardest[i]->Write();
        h1_InvThetaCh[i]->Write();
        h1_DeltaPhiMET_deg[i]->Write();
        h1_dphi_sum_deg[i]->Write();
        h1_pTmiss_mag[i]->Write();
        h1_pTabs_sum[i]->Write();
        h1_tan_theta_hardest[i]->Write();

        h2_ptrue_Nplates[i]->Write();
        h2_TrkEff_Nplates[i]->Write();
        h2_TrkEff_ptrue[i]->Write();
    }

    f_plot_disc->Save();
    f_plot_disc->Close();
    delete f_plot_disc;

}

TCanvas *cvs;
TLegend *lg;
void CvsLg_Init(){
    // initialize the canvas and legend
    cvs = new TCanvas("canvas", "canvas", 800, 600);
    lg = new TLegend(.7, .7, .9, .9);
    lg->SetFillStyle(0);
};

void DrawDisc_RecoTrue(TCanvas *cvs, TLegend *lg, TH1D *h1_reco, TH1D *h1_true, TString fig_name){

    cvs->Clear();
    lg->Clear();

    cvs->cd();
    h1_reco->SetStats(false);
    h1_true->SetStats(false);
    
    if(h1_reco->GetMaximum() < h1_true->GetMaximum()){
        h1_true->DrawCopy("HIST");
        h1_reco->DrawCopy("HIST same");
    }
    else{
        h1_reco->DrawCopy("HIST");
        h1_true->DrawCopy("HIST same");
    }
    
    lg->SetHeader("#it{#nu} NC");
    lg->AddEntry(h1_reco, "reco");
    lg->AddEntry(h1_true, "true");
    lg->Draw("same");

    cvs->SetGrid();
    cvs->SaveAs(fig_name);
}

void StoreHist2PDF(TString FigAddress){

    CvsLg_Init();

    DrawDisc_RecoTrue(cvs, lg, h1_dphi_max_deg[0], h1_dphi_max_deg[1], FigAddress+"h1_dphi_max_deg.pdf");
    DrawDisc_RecoTrue(cvs, lg, h1_dphi_sum_deg[0], h1_dphi_sum_deg[1], FigAddress+"h1_dphi_sum_deg.pdf");
    DrawDisc_RecoTrue(cvs, lg, h1_DeltaPhiMET_deg[0], h1_DeltaPhiMET_deg[1], FigAddress+"h1_DeltaPhiMET_deg.pdf");

    DrawDisc_RecoTrue(cvs, lg, h1_pmag_had_vis[0], h1_pmag_had_vis[1], FigAddress+"h1_pmag_had_vis.pdf");
    DrawDisc_RecoTrue(cvs, lg, h1_p_had_hardest[0], h1_p_had_hardest[1], FigAddress+"h1_p_had_hardest.pdf");
    DrawDisc_RecoTrue(cvs, lg, h1_pTmiss_mag[0], h1_pTmiss_mag[1], FigAddress+"h1_pTmiss_mag.pdf");
    DrawDisc_RecoTrue(cvs, lg, h1_pTabs_sum[0], h1_pTabs_sum[1], FigAddress+"h1_pTabs_sum.pdf");

    DrawDisc_RecoTrue(cvs, lg, h1_InvThetaCh[0], h1_InvThetaCh[1], FigAddress+"h1_InvThetaCh.pdf");
    DrawDisc_RecoTrue(cvs, lg, h1_tan_theta_hardest[0], h1_tan_theta_hardest[1], FigAddress+"h1_tan_theta_hardest.pdf");

}
