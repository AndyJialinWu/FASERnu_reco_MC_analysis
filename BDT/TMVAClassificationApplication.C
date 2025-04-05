#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
 
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
 
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "plot_disc.h"
 
using namespace TMVA;

void TMVAClassificationApplication( TString inputFileName = "PhysicsNTUP_ML.root", double BDT_cut = 0.99 )
{
    // This loads the library
    TMVA::Tools::Instance();

    // Create the Reader object
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    float n_ch;
    float DeltaPhiMET_reco;
    float dphi_max_reco;
    float dphi_sum_reco;
    float tan_theta_hardest_reco;
    float InvThetaCh_reco;
    float pmag_had_vis_reco;
    float p3_hardest_reco;
    float pTmiss_mag_reco;
    float pTabs_sum_reco;
    float nTrkTanThetaLeq0point1_reco;

    // Add variables to the reader
    reader->AddVariable( "n_ch", &n_ch );
    reader->AddVariable( "DeltaPhiMET_reco", &DeltaPhiMET_reco );
    reader->AddVariable( "dphi_max_reco", &dphi_max_reco );
    reader->AddVariable( "dphi_sum_reco", &dphi_sum_reco ); 
    reader->AddVariable( "tan_theta_hardest_reco", &tan_theta_hardest_reco );
    reader->AddVariable( "InvThetaCh_reco", &InvThetaCh_reco );
    reader->AddVariable( "pmag_had_vis_reco", &pmag_had_vis_reco );
    reader->AddVariable( "p3_hardest_reco", &p3_hardest_reco );
    reader->AddVariable( "pTmiss_mag_reco", &pTmiss_mag_reco );
    reader->AddVariable( "pTabs_sum_reco", &pTabs_sum_reco );
    reader->AddVariable( "nTrkTanThetaLeq0point1_reco", &nTrkTanThetaLeq0point1_reco );

    // Book the MVA methods
    TString methodName = "BDTG method";
    TString weightfile = "dataset/weights/TMVAClassification_BDTG.weights.xml";
    reader->BookMVA( methodName, weightfile );

    // Preparation of input tree
    TFile *f_validation = new TFile( inputFileName, "READ" );
    TTree *t_validation[4]; // 0: 200025, 1: 200026, 2: 200035, 3: 100069
    for(int RunIt=0; RunIt<4; RunIt++){
        t_validation[RunIt] = (TTree*)f_validation->Get(Form("disc%d_Validation", RunNumber[RunIt]));
        SetNumRecoNC(t_validation[RunIt], RunIt);
    }

    // Histograms
    Hist_Init();
    wgt_calc();

    // Loop over RunNumbers
    for(int RunIt=0; RunIt<4; RunIt++){

        int SigBkgIdx = -1;
        if(RunIt < 3){
            SigBkgIdx = 0;
        }
        else{
            SigBkgIdx = 1;
        }

        t_validation[RunIt]->SetBranchAddress("n_ch", &n_ch);
        t_validation[RunIt]->SetBranchAddress("DeltaPhiMET_reco", &DeltaPhiMET_reco);
        t_validation[RunIt]->SetBranchAddress("dphi_max_reco", &dphi_max_reco);     
        t_validation[RunIt]->SetBranchAddress("dphi_sum_reco", &dphi_sum_reco);
        t_validation[RunIt]->SetBranchAddress("tan_theta_hardest_reco", &tan_theta_hardest_reco);
        t_validation[RunIt]->SetBranchAddress("InvThetaCh_reco", &InvThetaCh_reco);
        t_validation[RunIt]->SetBranchAddress("pmag_had_vis_reco", &pmag_had_vis_reco);
        t_validation[RunIt]->SetBranchAddress("p3_hardest_reco", &p3_hardest_reco);
        t_validation[RunIt]->SetBranchAddress("pTmiss_mag_reco", &pTmiss_mag_reco);
        t_validation[RunIt]->SetBranchAddress("pTabs_sum_reco", &pTabs_sum_reco);
    
        for (Long64_t ievt=0; ievt<t_validation[RunIt]->GetEntries(); ievt++){

            if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

            t_validation[RunIt]->GetEntry(ievt);

            // Retrieve the MVA values and fill into histograms
            float mva = reader->EvaluateMVA( methodName );
            //std::cout << "RunNumber: " << RunNumbers[RunIt] << ", MVA value: " << mva << std::endl;

            // Fill histograms
            if(n_ch < 3) continue;

            h1_BDT_output[SigBkgIdx]->Fill(mva, weight[RunIt]);

            if(mva > BDT_cut){

                h1_n_ch[SigBkgIdx]->Fill(n_ch, weight[RunIt]);
                h1_DeltaPhiMET_deg[SigBkgIdx][0]->Fill(DeltaPhiMET_reco, weight[RunIt]);
                h1_dphi_max_deg[SigBkgIdx][0]->Fill(dphi_max_reco, weight[RunIt]);
                h1_dphi_sum_deg[SigBkgIdx][0]->Fill(dphi_sum_reco, weight[RunIt]);
                h1_tan_theta_hardest[SigBkgIdx][0]->Fill(tan_theta_hardest_reco, weight[RunIt]);
                h1_InvThetaCh[SigBkgIdx][0]->Fill(InvThetaCh_reco, weight[RunIt]);
                h1_pmag_had_vis[SigBkgIdx][0]->Fill(pmag_had_vis_reco, weight[RunIt]);
                h1_p_had_hardest[SigBkgIdx][0]->Fill(p3_hardest_reco, weight[RunIt]);
                h1_pTmiss_mag[SigBkgIdx][0]->Fill(pTmiss_mag_reco, weight[RunIt]);
                h1_pTabs_sum[SigBkgIdx][0]->Fill(pTabs_sum_reco, weight[RunIt]);

            }

        }

    }// End of RunNumber Loop

    h1_BDT_output_BackCum[0] = (TH1D*)h1_BDT_output[0]->GetCumulative(kFALSE);
    h1_BDT_output_BackCum[0]->SetTitle("NC Signal BDT Output Backward Cumulation");

    h1_BDT_output_BackCum[1] = (TH1D*)h1_BDT_output[1]->GetCumulative(kFALSE);
    h1_BDT_output_BackCum[1]->SetTitle("Background BDT Output Backward Cumulation");

    // Store histograms to ROOT file
    BDT_StoreHist2ROOT("BDT_Figures.root");

    // Store histograms to PDF
    BDT_StoreHist2PDF("Figures/BDT/");

    // BDT Cut Efficiency, Purity, Significance
    Graph_Init();
    BDT_Calc_Eff_Pur_Sign();
    BDT_StoreGraph2ROOT("BDT_Figures.root");
    BDT_StoreGraph2PDF("Figures/BDT/");

    delete reader;
    std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;

}