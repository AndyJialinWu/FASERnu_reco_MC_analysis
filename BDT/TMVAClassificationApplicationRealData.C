#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
 
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
 
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "plot_RealData.h"
 
using namespace TMVA;

void TMVAClassificationApplicationRealData( TString inputFileName = "PhysicsNTUP_ML_SimonThor.root", double BDT_cut = 0.99 )
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
    std::string *FileNum;

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
    TTree *t_validation; // real data
    t_validation = (TTree*)f_validation->Get("disc_RealData");

    // Histograms
    Hist_Init();
    //wgt_calc();

    t_validation->SetBranchAddress("n_ch", &n_ch);
    t_validation->SetBranchAddress("DeltaPhiMET_reco", &DeltaPhiMET_reco);
    t_validation->SetBranchAddress("dphi_max_reco", &dphi_max_reco);     
    t_validation->SetBranchAddress("dphi_sum_reco", &dphi_sum_reco);
    t_validation->SetBranchAddress("tan_theta_hardest_reco", &tan_theta_hardest_reco);
    t_validation->SetBranchAddress("InvThetaCh_reco", &InvThetaCh_reco);
    t_validation->SetBranchAddress("pmag_had_vis_reco", &pmag_had_vis_reco);
    t_validation->SetBranchAddress("p3_hardest_reco", &p3_hardest_reco);
    t_validation->SetBranchAddress("pTmiss_mag_reco", &pTmiss_mag_reco);
    t_validation->SetBranchAddress("pTabs_sum_reco", &pTabs_sum_reco);
    t_validation->SetBranchAddress("nTrkTanThetaLeq0point1_reco", &nTrkTanThetaLeq0point1_reco);
    t_validation->SetBranchAddress("FileNum", &FileNum);

    std::fstream outFile;
    outFile.open("text/NC_Candidates.txt", std::ios::out);

    for (Long64_t ievt=0; ievt<t_validation->GetEntries(); ievt++){

        if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

        t_validation->GetEntry(ievt);

        // Retrieve the MVA values and fill into histograms
        // Fill histograms
        if(n_ch < 3) continue;
        float mva = reader->EvaluateMVA( methodName );
        std::cout << "FileName: " << *FileNum << ", MVA value: " << mva << std::endl;
        h1_BDT_output->Fill(mva);

        if(mva > BDT_cut){
            // Fill the NC Candidates text file
            outFile << "----------------------------------------" << "\n";
            outFile << "FileName: " << *FileNum << ", MVA value: " << mva << "\n";
            outFile << "n_ch: " << n_ch << "\n";
            outFile << "DeltaPhiMET_reco: " << DeltaPhiMET_reco << "\n";
            outFile << "dphi_max_reco: " << dphi_max_reco << "\n";
            outFile << "dphi_sum_reco: " << dphi_sum_reco << "\n";
            outFile << "tan_theta_hardest_reco: " << tan_theta_hardest_reco << "\n";
            outFile << "InvThetaCh_reco: " << InvThetaCh_reco << "\n";
            outFile << "pmag_had_vis_reco: " << pmag_had_vis_reco << "\n";
            outFile << "p3_hardest_reco: " << p3_hardest_reco << "\n";
            outFile << "pTmiss_mag_reco: " << pTmiss_mag_reco << "\n";
            outFile << "pTabs_sum_reco: " << pTabs_sum_reco << "\n";
            outFile << "nTrkTanThetaLeq0point1_reco: " << nTrkTanThetaLeq0point1_reco << "\n";
            outFile << "----------------------------------------" << "\n";
        }
        
    }


    // Store histograms
    BDT_StoreHist2ROOT("BDT_Figures_SimonThor.root");

    delete reader;
    std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;

}