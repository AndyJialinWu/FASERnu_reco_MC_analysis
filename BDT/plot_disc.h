#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Rtypes.h"

const double IntLumi_run3 = 250.00;                    // fb^{-1}
const double IntLumi_F222_F223 = 9.523 + 28.9082;      // fb^{-1}
const double eff_FV = 638.00 / 730.00;
const double eff_n_ch[4] = {0.6549, 0.7963, 0.8267, 0.7073};
const double wgt_NC_200025 = 50173.00 * eff_FV / 10e3; // unit: event/fb^{-1}
const double wgt_NC_200026 = 66355.00 * eff_FV / 10e3;
const double wgt_NC_200035 = 13555.00 * eff_FV / 10e3;
const double wgt_Bkg = 4000.00 * eff_FV;               // unit: event/fb^{-1}
const double Run3ExpNC[4] = {IntLumi_run3 * wgt_NC_200025, IntLumi_run3 * wgt_NC_200026, IntLumi_run3 * wgt_NC_200035, IntLumi_run3 * wgt_Bkg};
const double F222_F223ExpNC[4] = {IntLumi_F222_F223 * wgt_NC_200025, IntLumi_F222_F223 * wgt_NC_200026, IntLumi_F222_F223 * wgt_NC_200035, IntLumi_F222_F223 * wgt_Bkg};
const int NumRecoNC[4] = {4285, 2733, 799, 12128};

double weight[4];
void wgt_calc(){
    for(int mcIt=0; mcIt<4; mcIt++){
        weight[mcIt] = F222_F223ExpNC[mcIt] * eff_n_ch[mcIt] / NumRecoNC[mcIt];
    }
};

const float t_plate = 1.1e3;                  // tungsten plate thickness in um
const float t_film = 0.34e3;                  // film thickness in um

TH2D *h2_preco_ptrue;

TH1D *h1_n_ch[2];
// first index: 0:NC, 1: Bkg; second index: 0:reco, 1:truth
TH1D *h1_dphi_max_deg[2][2]; 
TH1D *h1_dphi_sum_deg[2][2];
TH1D *h1_DeltaPhiMET_deg[2][2];

TH1D *h1_pmag_had_vis[2][2];
TH1D *h1_p_had_hardest[2][2];
TH1D *h1_pTmiss_mag[2][2];
TH1D *h1_pTabs_sum[2][2];

TH1D *h1_InvThetaCh[2][2];
TH1D *h1_tan_theta_hardest[2][2];

const TString str_cat[2][2] = {{"NC_reco", "NC_true"}, {"Bkg_reco", "Bkg_true"}};
const int LineColor[2][2] = {{kBlue, kRed}, {kGreen, kMagenta}};

void Hist_Init(){

    h2_preco_ptrue = new TH2D("h2_preco_ptrue", "P_{reco} versus P_{true}", 3500, 0, 3500, 3500, 0, 3500);
    h2_preco_ptrue->GetXaxis()->SetTitle("P_{true} [GeV/c]");
    h2_preco_ptrue->GetYaxis()->SetTitle("P_{reco} [GeV/c]");

    for(int i=0; i<2; i++){

        h1_n_ch[i] = new TH1D("h1_n_ch_"+str_cat[i][0], "Charged Tracks Multiplicity", 40, 0, 40);
        h1_n_ch[i]->GetXaxis()->SetTitle("n_{ch}");
        h1_n_ch[i]->GetYaxis()->SetTitle("# Events");

        for(int j=0; j<2; j++){
            h1_dphi_max_deg[i][j] = new TH1D("h1_dphi_max_deg_"+str_cat[i][j], "Largest Neighbouring Azimuthal Gap", 50, 0, 200);
            h1_dphi_max_deg[i][j]->GetXaxis()->SetTitle("#Delta #phi_{max} [#circ]");
            h1_dphi_max_deg[i][j]->GetYaxis()->SetTitle("# Events");
            h1_dphi_max_deg[i][j]->SetLineColor(LineColor[i][j]);

            h1_pmag_had_vis[i][j] = new TH1D("h1_pmag_had_vis_"+str_cat[i][j], "Sum of Visible Hadronic Momentum Magnitude", 140, 0, 1400);
            h1_pmag_had_vis[i][j]->GetXaxis()->SetTitle("|#vec{p}_{had, v}| [GeV/c]");
            h1_pmag_had_vis[i][j]->GetYaxis()->SetTitle("# Events");
            h1_pmag_had_vis[i][j]->SetLineColor(LineColor[i][j]);

            h1_p_had_hardest[i][j] = new TH1D("h1_p_had_hardest_"+str_cat[i][j], "Hardest Track Momentum", 100, 0, 1000);
            h1_p_had_hardest[i][j]->GetXaxis()->SetTitle("|#vec{p}_{hardest}| [GeV/c]");
            h1_p_had_hardest[i][j]->GetYaxis()->SetTitle("# Events");
            h1_p_had_hardest[i][j]->SetLineColor(LineColor[i][j]);

            h1_InvThetaCh[i][j] = new TH1D("h1_InvThetaCh_"+str_cat[i][j], "Inverse Sum of Charged Hadron Track Polar Angles", 100, 0, 2000);
            h1_InvThetaCh[i][j]->GetXaxis()->SetTitle("#sum |1/#theta_{ch}| [rad^{-1}]");
            h1_InvThetaCh[i][j]->GetYaxis()->SetTitle("# Events");
            h1_InvThetaCh[i][j]->SetLineColor(LineColor[i][j]);

            h1_DeltaPhiMET_deg[i][j] = new TH1D("h1_DeltaPhiMET_deg_"+str_cat[i][j], "Track-MET-Angle", 50, 0, 200);
            h1_DeltaPhiMET_deg[i][j]->GetXaxis()->SetTitle("#Delta #phi_{MET} [#circ]");
            h1_DeltaPhiMET_deg[i][j]->GetYaxis()->SetTitle("# Events");
            h1_DeltaPhiMET_deg[i][j]->SetLineColor(LineColor[i][j]);

            h1_dphi_sum_deg[i][j] = new TH1D("h1_dphi_sum_deg_"+str_cat[i][j], "Sum of Neighbouring Azimuthal Gaps", 90, 0, 360);
            h1_dphi_sum_deg[i][j]->GetXaxis()->SetTitle("#sum_{i} #Delta #phi_{i} [#circ]");
            h1_dphi_sum_deg[i][j]->GetYaxis()->SetTitle("# Events");
            h1_dphi_sum_deg[i][j]->SetLineColor(LineColor[i][j]);

            h1_pTmiss_mag[i][j] = new TH1D("h1_pTmiss_mag_"+str_cat[i][j], "Missing Transverse Momentum", 40, 0, 40);
            h1_pTmiss_mag[i][j]->GetXaxis()->SetTitle("|#vec{p}_{T, miss}| [GeV/c]");
            h1_pTmiss_mag[i][j]->GetYaxis()->SetTitle("# Events");
            h1_pTmiss_mag[i][j]->SetLineColor(LineColor[i][j]);

            h1_pTabs_sum[i][j] = new TH1D("h1_pTabs_sum_"+str_cat[i][j], "Sum of Magnitudes of Transverse Momenta", 40, 0, 40);
            h1_pTabs_sum[i][j]->GetXaxis()->SetTitle("#sum |#vec{p}_{T}| [GeV/c]");
            h1_pTabs_sum[i][j]->GetYaxis()->SetTitle("# Events");
            h1_pTabs_sum[i][j]->SetLineColor(LineColor[i][j]);

            h1_tan_theta_hardest[i][j] = new TH1D("h1_tan_theta_hardest_"+str_cat[i][j], "tan#theta of Hardest Hadronic Track", 100, 0, 0.5);
            h1_tan_theta_hardest[i][j]->GetXaxis()->SetTitle("tan#theta_{hardest}");
            h1_tan_theta_hardest[i][j]->GetYaxis()->SetTitle("# Events");
            h1_tan_theta_hardest[i][j]->SetLineColor(LineColor[i][j]);
        }

    }

};


void StoreHist2ROOT(){

    TFile *f_plot_disc = new TFile("Figures/disc_plots.root", "RECREATE");
    f_plot_disc->cd();

    h2_preco_ptrue->Write();
    

    for(int i=0; i<2; i++){
        
        h1_n_ch[i]->Write();

        for(int j=0; j<2; j++){
            h1_dphi_max_deg[i][j]->Write();
            h1_pmag_had_vis[i][j]->Write();
            h1_p_had_hardest[i][j]->Write();
            h1_InvThetaCh[i][j]->Write();
            h1_DeltaPhiMET_deg[i][j]->Write();
            h1_dphi_sum_deg[i][j]->Write();
            h1_pTmiss_mag[i][j]->Write();
            h1_pTabs_sum[i][j]->Write();
            h1_tan_theta_hardest[i][j]->Write();
        }

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

void DrawDisc_RecoTrue(TCanvas *cvs, TLegend *lg, TH1D *h1_reco, TH1D *h1_true, TString fig_name, TString category){

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
    
    lg->SetHeader(category);
    lg->AddEntry(h1_reco, "reco");
    lg->AddEntry(h1_true, "true");
    lg->Draw("same");

    cvs->SetGrid();
    cvs->SaveAs(fig_name);
}

void StoreHist2PDF(TString FigAddress){

    CvsLg_Init();

    TString category[2] = {"NC", "Bkg"};

    for(int i=0; i<2; i++){

        DrawDisc_RecoTrue(cvs, lg, h1_dphi_max_deg[i][0], h1_dphi_max_deg[i][1], FigAddress+"h1_dphi_max_deg_" + category[i] + ".pdf", category[i]);
        DrawDisc_RecoTrue(cvs, lg, h1_dphi_sum_deg[i][0], h1_dphi_sum_deg[i][1], FigAddress+"h1_dphi_sum_deg_" + category[i] + ".pdf", category[i]);
        DrawDisc_RecoTrue(cvs, lg, h1_DeltaPhiMET_deg[i][0], h1_DeltaPhiMET_deg[i][1], FigAddress+"h1_DeltaPhiMET_deg_" + category[i] + ".pdf", category[i]);

        DrawDisc_RecoTrue(cvs, lg, h1_pmag_had_vis[i][0], h1_pmag_had_vis[i][1], FigAddress+"h1_pmag_had_vis_" + category[i] + ".pdf", category[i]);
        DrawDisc_RecoTrue(cvs, lg, h1_p_had_hardest[i][0], h1_p_had_hardest[i][1], FigAddress+"h1_p_had_hardest_" + category[i] + ".pdf", category[i]);
        DrawDisc_RecoTrue(cvs, lg, h1_pTmiss_mag[i][0], h1_pTmiss_mag[i][1], FigAddress+"h1_pTmiss_mag_" + category[i] + ".pdf", category[i]);
        DrawDisc_RecoTrue(cvs, lg, h1_pTabs_sum[i][0], h1_pTabs_sum[i][1], FigAddress+"h1_pTabs_sum_" + category[i] + ".pdf", category[i]);

        DrawDisc_RecoTrue(cvs, lg, h1_InvThetaCh[i][0], h1_InvThetaCh[i][1], FigAddress+"h1_InvThetaCh_" + category[i] + ".pdf", category[i]);
        DrawDisc_RecoTrue(cvs, lg, h1_tan_theta_hardest[i][0], h1_tan_theta_hardest[i][1], FigAddress+"h1_tan_theta_hardest_" + category[i] + ".pdf", category[i]);

    }

    

}
