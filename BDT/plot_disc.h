#include <iostream>

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TPDGCode.h"
#include "Rtypes.h"

const int RunNumber[4] = {200025, 200026, 200035, 100069};
const double IntLumi_run3 = 250.00;                    // fb^{-1}
const double IntLumi_F222_F223 = 9.523 + 28.9082;      // fb^{-1}
const double eff_FV = 638.00 / 730.00;
double eff_n_ch[4] = {0.6549, 0.7963, 0.8267, 0.7073}; // n_ch >= 3
const double wgt_NC_200025 = 50173.00 * eff_FV / 10e3; // unit: event/fb^{-1}
const double wgt_NC_200026 = 66355.00 * eff_FV / 10e3;
const double wgt_NC_200035 = 13555.00 * eff_FV / 10e3;
const double wgt_Bkg = 4000.00 * eff_FV;               // unit: event/fb^{-1}
const double Run3ExpNC[4] = {IntLumi_run3 * wgt_NC_200025, IntLumi_run3 * wgt_NC_200026, IntLumi_run3 * wgt_NC_200035, IntLumi_run3 * wgt_Bkg};
const double F222_F223ExpNC[4] = {IntLumi_F222_F223 * wgt_NC_200025, IntLumi_F222_F223 * wgt_NC_200026, IntLumi_F222_F223 * wgt_NC_200035, IntLumi_F222_F223 * wgt_Bkg};
int NumRecoNC[4] = {0, 0, 0, 0};
double weight[4];

void SetNumRecoNC(TTree *tree, int RunIt){

    NumRecoNC[RunIt] = tree->GetEntries();

    float n_ch;
    int n_ch_geq3 = 0;
    tree->SetBranchAddress("n_ch", &n_ch);
    for(Long64_t ievt=0; ievt<tree->GetEntries(); ievt++){
        tree->GetEntry(ievt);
        if(n_ch >= 3){
            n_ch_geq3++;
        }
    }// loop over events

    eff_n_ch[RunIt] = double(n_ch_geq3) / double(tree->GetEntries());
    std::cout<<"RunNumber: "<<RunNumber[RunIt]<<", NumRecoNC: "<<NumRecoNC[RunIt]<<", eff_n_ch_geq3: "<<eff_n_ch[RunIt]<<std::endl;

};

void wgt_calc(){

    for(int RunIt=0; RunIt<4; RunIt++){
        weight[RunIt] = F222_F223ExpNC[RunIt] / double(NumRecoNC[RunIt]);
    }

};

const float t_plate = 1.1e3;                  // tungsten plate thickness in um
const float t_film = 0.34e3;                  // film thickness in um

TH2D *h2_preco_ptrue[2];                      // 0: hadrons, 1: electrons
const TString str_hadron_electron[2] = {"hadron", "electron"};

TH2D *h2_pTmissReco_pTmissTrue[2];            // 0: nu NC signal, 1: Bkg
const TString str_pTmiss_cat[2] = {"NC", "Bkg"};

TH1D *h1_n_ch[2];
// first index: 0: nu NC signal, 1: Bkg; second index: 0:reco, 1:truth
TH1D *h1_dphi_max_deg[2][2]; 
TH1D *h1_dphi_sum_deg[2][2];
TH1D *h1_DeltaPhiMET_deg[2][2];

TH1D *h1_pmag_had_vis[2][2];
TH1D *h1_p_had_hardest[2][2];
TH1D *h1_pTmiss_mag[2][2];
TH1D *h1_pTabs_sum[2][2];

TH1D *h1_InvThetaCh[2][2];
TH1D *h1_tan_theta_hardest[2][2];

TH1D *h1_nTrkTanThetaLeq0point1[2][2];     // first index: 0: nu NC signal, 1: Bkg; second index: 0:reco, 1:truth

TH2D *h2_PartCat_n_ch[2];                  // 0: nu NC signal, 1: Bkg
const TString str_PartCat[5] = {"#pi^{#pm}", "K^{#pm}", "p + #bar{p}", "e^{#pm}", "others"};
const TString str_n_ch_cat[2] = {"NC", "Bkg"};

const TString str_cat[2][2] = {{"NC_reco", "NC_true"}, {"Bkg_reco", "Bkg_true"}};
const int LineColor[2][2] = {{kBlue, kRed}, {kGreen, kMagenta}};

TH1D *h1_BDT_output[2];                    // 0: nu NC signal, 1: Bkg
TH1D *h1_BDT_output_BackCum[2];            // backward cumulative distribution
const TString str_BDT_output[2] = {"NC", "Bkg"};

void Hist_Init(){

    for(int i=0; i<2; i++){

        h1_BDT_output[i] = new TH1D("h1_BDT_output_"+str_BDT_output[i], "BDT Output of "+str_BDT_output[i], 200, -1, 1);
        h1_BDT_output[i]->GetXaxis()->SetTitle("BDT Output");
        h1_BDT_output[i]->GetYaxis()->SetTitle("# Events");

        h2_preco_ptrue[i] = new TH2D("h2_preco_ptrue_"+str_hadron_electron[i], "P_{reco} versus P_{true} of "+str_hadron_electron[i], 3500, 0, 3500, 3500, 0, 3500);
        h2_preco_ptrue[i]->GetXaxis()->SetTitle("P_{true} [GeV/c]");
        h2_preco_ptrue[i]->GetYaxis()->SetTitle("P_{reco} [GeV/c]");

        h1_n_ch[i] = new TH1D("h1_n_ch_"+str_cat[i][0], "Charged Tracks Multiplicity", 40, 0, 40);
        h1_n_ch[i]->GetXaxis()->SetTitle("n_{ch}");
        h1_n_ch[i]->GetYaxis()->SetTitle("# Events");

        h2_PartCat_n_ch[i] = new TH2D("h2_PartCat_n_ch_"+str_n_ch_cat[i], "Particle Category versus Charged Tracks Multiplicity of "+str_n_ch_cat[i], 40, 0, 40, 5, 0, 5);
        h2_PartCat_n_ch[i]->GetXaxis()->SetTitle("n_{ch}");
        h2_PartCat_n_ch[i]->GetYaxis()->SetTitle("Particle Category");
        for(int PartCatIt=0; PartCatIt<5; PartCatIt++){
            h2_PartCat_n_ch[i]->GetYaxis()->SetBinLabel(PartCatIt+1, str_PartCat[PartCatIt]);
        }

        h2_pTmissReco_pTmissTrue[i] = new TH2D("h2_pTmissReco_pTmissTrue_"+str_pTmiss_cat[i], "P_{T, miss}^{reco} versus P_{T, miss}^{true} of "+str_pTmiss_cat[i], 100, 0, 100, 100, 0, 100);
        h2_pTmissReco_pTmissTrue[i]->GetXaxis()->SetTitle("P_{T, miss}^{true} [GeV/c]");
        h2_pTmissReco_pTmissTrue[i]->GetYaxis()->SetTitle("P_{T, miss}^{reco} [GeV/c]");

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

            h1_nTrkTanThetaLeq0point1[i][j] = new TH1D("h1_nTrkTanThetaLeq0point1_"+str_cat[i][j], "Number of Tracks with tan#theta < 0.1", 20, 0, 20);
            h1_nTrkTanThetaLeq0point1[i][j]->GetXaxis()->SetTitle("# tracks with tan#theta < 0.1");
            h1_nTrkTanThetaLeq0point1[i][j]->GetYaxis()->SetTitle("# Events");
            h1_nTrkTanThetaLeq0point1[i][j]->SetLineColor(LineColor[i][j]);
        }

    }

};

void Fill_h2_PartCat_n_ch(int PDG, float n_ch, int SigBkgIdx, double wgt){

    if(std::abs(PDG) == kPiPlus){
        h2_PartCat_n_ch[SigBkgIdx]->Fill(n_ch, 0.5, wgt);
    }
    else if(std::abs(PDG) == kKPlus){
        h2_PartCat_n_ch[SigBkgIdx]->Fill(n_ch, 1.5, wgt);
    }
    else if(std::abs(PDG) == kProton){
        h2_PartCat_n_ch[SigBkgIdx]->Fill(n_ch, 2.5, wgt);
    }
    else if(std::abs(PDG) == kElectron){
        h2_PartCat_n_ch[SigBkgIdx]->Fill(n_ch, 3.5, wgt);
    }
    else{
        h2_PartCat_n_ch[SigBkgIdx]->Fill(n_ch, 4.5, wgt);
    }

};


void StoreHist2ROOT(TString RootFileName){

    TFile *f_plot_disc = new TFile("Figures/" + RootFileName, "RECREATE");
    f_plot_disc->cd();

    for(int i=0; i<2; i++){
        
        h1_n_ch[i]->Write();
        h2_preco_ptrue[i]->Write();
        h2_PartCat_n_ch[i]->Write();
        h2_pTmissReco_pTmissTrue[i]->Write();

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
            h1_nTrkTanThetaLeq0point1[i][j]->Write();
        }

    }

    f_plot_disc->Save();
    f_plot_disc->Close();
    delete f_plot_disc;

};

void BDT_StoreHist2ROOT(TString RootFileName){

    TFile *f_plot_disc = new TFile("Figures/" + RootFileName, "RECREATE");
    f_plot_disc->cd();

    for(int i=0; i<2; i++){
        
        h1_BDT_output[i]->Write();
        h1_BDT_output_BackCum[i]->Write();

        h1_n_ch[i]->Write();
        h1_DeltaPhiMET_deg[i][0]->Write();
        h1_dphi_max_deg[i][0]->Write();
        h1_dphi_sum_deg[i][0]->Write();
        h1_tan_theta_hardest[i][0]->Write();
        h1_InvThetaCh[i][0]->Write();
        h1_pmag_had_vis[i][0]->Write();
        h1_p_had_hardest[i][0]->Write();
        h1_pTmiss_mag[i][0]->Write();
        h1_pTabs_sum[i][0]->Write();
        h1_nTrkTanThetaLeq0point1[i][0]->Write();

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

};

void BDT_DrawDisc_SigBkg(TCanvas *cvs, TLegend *lg, TH1D *h1_sig, TH1D *h1_bkg, TString fig_name, TString LgHeader, bool logy=false){

    cvs->Clear();
    lg->Clear();

    cvs->cd();
    h1_sig->SetStats(false);
    h1_bkg->SetStats(false);
    h1_sig->SetLineColor(kRed);
    h1_bkg->SetLineColor(kBlue);
    
    if(h1_sig->GetMaximum() < h1_bkg->GetMaximum()){
        h1_bkg->DrawCopy("HIST");
        h1_sig->DrawCopy("HIST same");
    }
    else{
        h1_sig->DrawCopy("HIST");
        h1_bkg->DrawCopy("HIST same");
    }
    
    lg->SetHeader(LgHeader);
    lg->AddEntry(h1_sig, "NC Signal");
    lg->AddEntry(h1_bkg, "Background");
    lg->Draw("same");

    if(logy){
        h1_sig->SetMinimum(0);
        h1_bkg->SetMinimum(0);
        cvs->SetLogy();
    }
    cvs->SetGrid();
    cvs->SaveAs(fig_name);

};

void StoreHist2PDF(TString FigAddress){

    CvsLg_Init();

    TString category[2] = {"NC", "Bkg"};

    // Loop over signal and background
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
        DrawDisc_RecoTrue(cvs, lg, h1_nTrkTanThetaLeq0point1[i][0], h1_nTrkTanThetaLeq0point1[i][1], FigAddress+"h1_nTrkTanThetaLeq0point1_" + category[i] + ".pdf", category[i]);

    }// end of loop over signal and background

};

void BDT_StoreHist2PDF(TString FigAddress){

    CvsLg_Init();

    TString category = "After BDT Cut";

    BDT_DrawDisc_SigBkg(cvs, lg, h1_BDT_output[0], h1_BDT_output[1], FigAddress+"h1_BDT_output" + ".pdf", "", true);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_BDT_output_BackCum[0], h1_BDT_output_BackCum[1], FigAddress+"h1_BDT_output_BackCum" + ".pdf", "", true);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_n_ch[0], h1_n_ch[1], FigAddress+"h1_n_ch" + ".pdf", category);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_dphi_max_deg[0][0], h1_dphi_max_deg[1][0], FigAddress+"h1_dphi_max_deg" + ".pdf", category);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_dphi_sum_deg[0][0], h1_dphi_sum_deg[1][0], FigAddress+"h1_dphi_sum_deg" + ".pdf", category);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_DeltaPhiMET_deg[0][0], h1_DeltaPhiMET_deg[1][0], FigAddress+"h1_DeltaPhiMET_deg" + ".pdf", category);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_pmag_had_vis[0][0], h1_pmag_had_vis[1][0], FigAddress+"h1_pmag_had_vis" + ".pdf", category);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_p_had_hardest[0][0], h1_p_had_hardest[1][0], FigAddress+"h1_p_had_hardest" + ".pdf", category);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_pTmiss_mag[0][0], h1_pTmiss_mag[1][0], FigAddress+"h1_pTmiss_mag" + ".pdf", category);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_pTabs_sum[0][0], h1_pTabs_sum[1][0], FigAddress+"h1_pTabs_sum" + ".pdf", category);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_InvThetaCh[0][0], h1_InvThetaCh[1][0], FigAddress+"h1_InvThetaCh" + ".pdf", category);
    BDT_DrawDisc_SigBkg(cvs, lg, h1_tan_theta_hardest[0][0], h1_tan_theta_hardest[1][0], FigAddress+"h1_tan_theta_hardest" + ".pdf", category);
    
};


TGraph *g_BDT_nuNC[3]; // 0: efficiency, 1: purity, 2: significance
const TString str_BDT_nuNC[3] = {"Efficiency", "Purity", "Significance"};

void Graph_Init(){

    const int Color[3] = {kBlue, kRed, kGreen};

    for(int i=0; i<3; i++){
        g_BDT_nuNC[i] = new TGraph();
        g_BDT_nuNC[i]->SetNameTitle("g_BDT_nuNC_"+str_BDT_nuNC[i], "BDT Cut NC Signal " + str_BDT_nuNC[i]);
        g_BDT_nuNC[i]->GetXaxis()->SetTitle("BDT Cut");
        g_BDT_nuNC[i]->GetYaxis()->SetTitle(str_BDT_nuNC[i]);
        g_BDT_nuNC[i]->SetMarkerStyle(20);
        g_BDT_nuNC[i]->SetMarkerSize(0.5);
        g_BDT_nuNC[i]->SetMarkerColor(Color[i]);
        g_BDT_nuNC[i]->SetLineColor(Color[i]);
    }

};

void BDT_Calc_Eff_Pur_Sign(){

    for(int bin=1; bin<= h1_BDT_output[0]->GetNbinsX(); bin++){
        
        double eff = h1_BDT_output_BackCum[0]->GetBinContent(bin) / h1_BDT_output[0]->Integral();
        double pur = h1_BDT_output_BackCum[0]->GetBinContent(bin) / (h1_BDT_output_BackCum[0]->GetBinContent(bin) + h1_BDT_output_BackCum[1]->GetBinContent(bin));
        double sign = h1_BDT_output_BackCum[0]->GetBinContent(bin) / std::sqrt(h1_BDT_output_BackCum[0]->GetBinContent(bin) + h1_BDT_output_BackCum[1]->GetBinContent(bin));

        g_BDT_nuNC[0]->SetPoint(bin-1, h1_BDT_output[0]->GetBinCenter(bin), eff);
        g_BDT_nuNC[1]->SetPoint(bin-1, h1_BDT_output[0]->GetBinCenter(bin), pur);   
        g_BDT_nuNC[2]->SetPoint(bin-1, h1_BDT_output[0]->GetBinCenter(bin), sign);

    }

};

void BDT_StoreGraph2ROOT(TString RootFileName){

    TFile *f_plot_disc = new TFile("Figures/" + RootFileName, "UPDATE");
    f_plot_disc->cd();

    for(int i=0; i<3; i++){
        g_BDT_nuNC[i]->Write();
    }

    f_plot_disc->Save();
    f_plot_disc->Close();
    delete f_plot_disc;

};

void BDT_StoreGraph2PDF(TString FigAddress){

    CvsLg_Init();
    const int Color[3] = {kBlue, kRed, kGreen};

    for(int i=0; i<3; i++){

        g_BDT_nuNC[i]->SetMarkerStyle(20);
        g_BDT_nuNC[i]->SetMarkerSize(0.5);
        g_BDT_nuNC[i]->SetMarkerColor(Color[i]);
        g_BDT_nuNC[i]->SetLineColor(Color[i]);
        g_BDT_nuNC[i]->GetXaxis()->SetTitle("BDT Cut");
        g_BDT_nuNC[i]->GetXaxis()->SetLimits(-1, 1);
        g_BDT_nuNC[i]->GetHistogram()->SetMinimum(0);
        if(i<2){
            g_BDT_nuNC[i]->GetHistogram()->SetMaximum(1.1);
            g_BDT_nuNC[i]->GetYaxis()->SetTitle("Efficiency (Purity)");
        }
        else{
            g_BDT_nuNC[i]->GetHistogram()->SetMaximum(7);
            g_BDT_nuNC[i]->GetYaxis()->SetTitle("Significance");
        }

    }
    cvs->cd();
    lg->SetHeader("NC Signal");
    g_BDT_nuNC[0]->Draw("APL");
    g_BDT_nuNC[1]->Draw("PL same");
    lg->AddEntry(g_BDT_nuNC[0], str_BDT_nuNC[0]);
    lg->AddEntry(g_BDT_nuNC[1], str_BDT_nuNC[1]);
    lg->Draw("same");
    cvs->SetGrid();
    cvs->SaveAs(FigAddress+"g_BDT_nuNC_eff_pur" + ".pdf");

    cvs->Clear();
    lg->Clear();
    cvs->cd();
    lg->SetHeader("NC Signal");
    g_BDT_nuNC[2]->Draw("APL");
    lg->AddEntry(g_BDT_nuNC[2], str_BDT_nuNC[2]);
    lg->Draw("same");
    cvs->SetGrid();
    cvs->SaveAs(FigAddress+"g_BDT_nuNC_sign" + ".pdf");

};


