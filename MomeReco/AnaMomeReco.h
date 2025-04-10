#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Rtypes.h"


TH1D *h1_IP[3];                 // 0: primary, 1: secondary, 2:electron
const TString IP_label[3] = {"Primary", "Secondary", "Electron"};
const Color_t IP_color[3] = {kRed, kBlue, kGreen};

TH2D *h2_ptrue_dz[3];           // 0: primary, 1: secondary, 2: electron
const TString ptrue_dz_label[3] = {"Primary", "Secondary", "Electron"};

TH2D *h2_preco_ptrue_cat[5];    // P_reco versus P_true for 0:electron, 1:muon, 2:pion, 3:kaon, 4:proton
const TString cat_label[5] = {"Electron", "Muon", "Pion", "Kaon", "Proton"};
TH1D *h1_MaxGap[2];             // 0: faliure 1, 1: faliure 2
const TString MaxGap_label[2] = {"Failure1", "Failure2"};

TH2D *h2_IP_dz[2];              // 0: true primary, 1: secondary
const TString IP_dz_label[2] = {"Primary", "Secondary"};

TH2D *h2_preco_ptrue;           // P_reco versus P_true for P_reco algorithm development

TH2D *h2_thetaRMS_ptrue[4];     // 0: success, 1: failure 2, 2: failure 1, 3: electron
const TString thetaRMS_ptrue_label[4] = {"Success", "Failure2", "Failure1", "Electron"};

TH2D *h2_TrackLength_ptrue[4];  // 0: success, 1: failure 2, 2: failure 1, 3: electron
const TString TrackLength_ptrue_label[4] = {"Success", "Failure2", "Failure1", "Electron"};

TH2D *h2_MaxKinkAngle_ptrue[4]; // 0: success, 1: failure 2, 2: failure 1, 3: electron
const TString MaxKinkAngle_ptrue_label[4] = {"Success", "Failure2", "Failure1", "Electron"};

void Hist_Init(){

    h2_preco_ptrue = new TH2D("h2_preco_ptrue", "P_{reco} versus P_{true}", 1000, 0, 1000, 1000, 0, 1000);
    h2_preco_ptrue->GetXaxis()->SetTitle("P_{true} [GeV/c]");
    h2_preco_ptrue->GetYaxis()->SetTitle("P_{reco} [GeV/c]");

    for(int i=0; i<2; i++){
        h1_MaxGap[i] = new TH1D("h1_MaxGap_"+MaxGap_label[i], "Max Gap of "+MaxGap_label[i], 6, 0, 6);
        h1_MaxGap[i]->GetXaxis()->SetTitle("Max Gap");
        h1_MaxGap[i]->GetYaxis()->SetTitle("Entries");

        h2_IP_dz[i] = new TH2D("h2_IP_dz_"+IP_dz_label[i], "Impact Parameter versus dz of "+IP_dz_label[i]+" Tracks", 400, 0, 4000, 100, 0, 10);
        h2_IP_dz[i]->GetXaxis()->SetTitle("dz [#mu m]");       
        h2_IP_dz[i]->GetYaxis()->SetTitle("IP [#mu m]");
    }

    for(int i=0; i<3; i++){

        h1_IP[i] = new TH1D("h1_IP_"+IP_label[i], "Impact Parameter of True "+IP_label[i]+" Tracks", 100, 0, 10);
        h1_IP[i]->GetXaxis()->SetTitle("IP [#mu m]");
        h1_IP[i]->GetYaxis()->SetTitle("Entries");
        h1_IP[i]->SetLineColor(IP_color[i]);

        h2_ptrue_dz[i] = new TH2D("h2_ptrue_dz_"+ptrue_dz_label[i], "P_{true} versus dz of "+ptrue_dz_label[i]+" Tracks", 400, 0, 4000, 1000, 0, 1000);
        h2_ptrue_dz[i]->GetXaxis()->SetTitle("dz [#mu m]");
        h2_ptrue_dz[i]->GetYaxis()->SetTitle("P_{true} [GeV/c]");

    }    

    for(int i=0; i<4; i++){
        h2_thetaRMS_ptrue[i] = new TH2D("h2_thetaRMS_ptrue_"+thetaRMS_ptrue_label[i], "#theta_{RMS} versus P_{true} for "+thetaRMS_ptrue_label[i], 1000, 0, 1000, 170, 0, 17);
        h2_thetaRMS_ptrue[i]->GetXaxis()->SetTitle("P_{true} [GeV/c]");
        h2_thetaRMS_ptrue[i]->GetYaxis()->SetTitle("#theta_{RMS} [mrad]");

        h2_TrackLength_ptrue[i] = new TH2D("h2_TrackLength_ptrue_"+TrackLength_ptrue_label[i], "Track Length versus P_{true} for "+TrackLength_ptrue_label[i], 1000, 0, 1000, 150, 0, 150);
        h2_TrackLength_ptrue[i]->GetXaxis()->SetTitle("P_{true} [GeV/c]");
        h2_TrackLength_ptrue[i]->GetYaxis()->SetTitle("Track Length [mm]");

        h2_MaxKinkAngle_ptrue[i] = new TH2D("h2_MaxKinkAngle_ptrue_"+MaxKinkAngle_ptrue_label[i], "Max Kink Angle versus P_{true} for "+MaxKinkAngle_ptrue_label[i], 1000, 0, 1000, 250, 0, 25);
        h2_MaxKinkAngle_ptrue[i]->GetXaxis()->SetTitle("P_{true} [GeV/c]");
        h2_MaxKinkAngle_ptrue[i]->GetYaxis()->SetTitle("Max Kink Angle [mrad]");
    }
    

    for(int i=0; i<5; i++){
        h2_preco_ptrue_cat[i] = new TH2D("h2_preco_ptrue_"+cat_label[i], "P_{reco} versus P_{true} for "+cat_label[i], 1000, 0, 1000, 1000, 0, 1000);
        h2_preco_ptrue_cat[i]->GetXaxis()->SetTitle("P_{true} [GeV/c]");
        h2_preco_ptrue_cat[i]->GetYaxis()->SetTitle("P_{reco} [GeV/c]");
    }



}


void StoreHist2ROOT(){

    TFile *f_out = new TFile("AnaMomeReco_ML.root", "RECREATE");

    h2_preco_ptrue->Write();

    for(int i=0; i<3; i++){
        h1_IP[i]->Write();
        h2_ptrue_dz[i]->Write();
    }

    for(int i=0; i<4; i++){
        h2_thetaRMS_ptrue[i]->Write();
        h2_TrackLength_ptrue[i]->Write();
        h2_MaxKinkAngle_ptrue[i]->Write();
    }

    for(int i=0; i<2; i++){
        h1_MaxGap[i]->Write();
        h2_IP_dz[i]->Write();
        
    }

    for(int i=0; i<5; i++){
        h2_preco_ptrue_cat[i]->Write();
    }

    f_out->Close();

}

