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

TH1D *h1_BDT_output;    

void Hist_Init(){

    h1_BDT_output = new TH1D("h1_BDT_output", "BDT Output", 200, -1, 1);
    h1_BDT_output->GetXaxis()->SetTitle("BDT Output");
    h1_BDT_output->GetYaxis()->SetTitle("# Events");

}

void BDT_StoreHist2ROOT(TString RootFileName){

    TFile *f_plot_disc = new TFile("Figures/RealData/" + RootFileName, "RECREATE");
    f_plot_disc->cd();

    h1_BDT_output->Write();

    f_plot_disc->Save();
    f_plot_disc->Close();
    delete f_plot_disc;
};