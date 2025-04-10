#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TString.h"
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

TH2D *h2_MomStatus_nseg; // tracks level
TH1D *h1_SuccessRate;    // tracks level
TH2D *h2_preco_ptrue;    // tracks level successful reco
TH2D *h2_preso_ptrue;    // tracks level successful reco
TH2D *h2_nsegm_ptrue;    // tracks level successful reco
TH2D *h2_preso_nsegm;    // tracks level successful reco

TH2D *h2_SuccessRate_NuE;// event level

void Hist_Init(){

    h2_MomStatus_nseg = new TH2D("h2_MomStatus_nseg", "Momentum Measurement Status versus # seg.", 100, 0, 100, 3, 0, 3);
    h2_MomStatus_nseg->GetXaxis()->SetTitle("nseg");
    h2_MomStatus_nseg->GetYaxis()->SetTitle("Momentum Measurement Status");
    h2_MomStatus_nseg->GetYaxis()->SetBinLabel(1, "Success");
    h2_MomStatus_nseg->GetYaxis()->SetBinLabel(2, "Failure1");
    h2_MomStatus_nseg->GetYaxis()->SetBinLabel(3, "Failure2");

    h1_SuccessRate = new TH1D("h1_SuccessRate", "Success Rate versus # seg.", 100, 0, 100);
    h1_SuccessRate->GetXaxis()->SetTitle("nseg");
    h1_SuccessRate->GetYaxis()->SetTitle("Success Rate");
    h1_SuccessRate->GetYaxis()->SetRangeUser(0, 1);

    h2_preco_ptrue = new TH2D("h2_preco_ptrue", "Successful p_{reco} versus p_{true}", 100, 0, 1000, 100, 0, 1000);
    h2_preco_ptrue->GetXaxis()->SetTitle("p_{true} [GeV/c]");
    h2_preco_ptrue->GetYaxis()->SetTitle("p_{reco} [GeV/c]");

    h2_preso_ptrue = new TH2D("h2_preso_ptrue", "Successful (p_{reco}-p_{true})/p_{true} versus p_{true}", 100, 0, 1000, 1020, -2, 100);
    h2_preso_ptrue->GetXaxis()->SetTitle("p_{true} [GeV/c]");
    h2_preso_ptrue->GetYaxis()->SetTitle("(p_{reco}-p_{true})/p_{true}");

    h2_nsegm_ptrue = new TH2D("h2_nsegm_ptrue", "Successful nseg versus p_{true}", 100, 0, 1000, 100, 0, 100);
    h2_nsegm_ptrue->GetXaxis()->SetTitle("p_{true} [GeV/c]");
    h2_nsegm_ptrue->GetYaxis()->SetTitle("nseg");

    h2_preso_nsegm = new TH2D("h2_preso_nsegm", "Successful (p_{reco}-p_{true})/p_{true} versus nseg", 100, 0, 100, 1020, -2, 100);
    h2_preso_nsegm->GetXaxis()->SetTitle("nseg");
    h2_preso_nsegm->GetYaxis()->SetTitle("(p_{reco}-p_{true})/p_{true}");

    h2_SuccessRate_NuE = new TH2D("h2_SuccessRate_NuE", "Tracks Momentum Measurement Success Rate Event By Event versus E_{#nu}", 350, 0, 3500, 100, 0, 1);
    h2_SuccessRate_NuE->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    h2_SuccessRate_NuE->GetYaxis()->SetTitle("Success Rate");
    h2_SuccessRate_NuE->GetYaxis()->SetRangeUser(0, 1);

};

void GetSuccessRate(){

    for(int i=0; i<h2_MomStatus_nseg->GetNbinsX(); i++){

        int nseg = h2_MomStatus_nseg->GetXaxis()->GetBinCenter(i);
        int n_success = h2_MomStatus_nseg->GetBinContent(i, 1);
        int n_failure1 = h2_MomStatus_nseg->GetBinContent(i, 2);
        int n_failure2 = h2_MomStatus_nseg->GetBinContent(i, 3);

        if(n_success + n_failure1 + n_failure2 > 0){
            float success_rate = float(n_success) / float(n_success + n_failure1 + n_failure2);
            h1_SuccessRate->SetBinContent(i, success_rate);
        }
        else{
            h1_SuccessRate->SetBinContent(i, 0);
        }

    }

}




void StoreHist2ROOT(TString RootFileName){

    TFile *f_plot_disc = new TFile("Figures/" + RootFileName, "RECREATE");
    f_plot_disc->cd();

    h1_SuccessRate->Write();
    h2_MomStatus_nseg->Write();
    h2_SuccessRate_NuE->Write();
    h2_preco_ptrue->Write();
    h2_preso_ptrue->Write();
    h2_nsegm_ptrue->Write();
    h2_preso_nsegm->Write();

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


void StoreHist2PDF(TString FigAddress){

    CvsLg_Init();

};




