#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
 
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
 
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/TMVAMultiClassGui.h"

void TMVAClassification(TString inputFileName, TString outputFileName){

    // This loads the library
    TMVA::Tools::Instance();
    std::cout << std::endl;
    std::cout << "==> Start TMVAMulticlass" << std::endl;

    // Read the input file
    TFile *f_PhysicsNTUP_ML = new TFile(inputFileName, "READ");

    // Register the training and test trees
    TTree *signalTree0 = (TTree*)f_PhysicsNTUP_ML->Get("disc200025_TrainTest");
    TTree *signalTree1 = (TTree*)f_PhysicsNTUP_ML->Get("disc200026_TrainTest");
    TTree *signalTree2 = (TTree*)f_PhysicsNTUP_ML->Get("disc200035_TrainTest");
    TTree *background0 = (TTree*)f_PhysicsNTUP_ML->Get("disc100069_TrainTest");

    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TString outfileName( outputFileName );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    // Create the factory object
    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                                "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

    // Add variables to the dataloader
    dataloader->AddVariable("n_ch", 'F');
    dataloader->AddVariable("DeltaPhiMET_reco", 'F');
    dataloader->AddVariable("dphi_max_reco", 'F');
    dataloader->AddVariable("dphi_sum_reco", 'F');
    dataloader->AddVariable("tan_theta_hardest_reco", 'F');
    dataloader->AddVariable("InvThetaCh_reco", 'F');

    dataloader->AddVariable("pmag_had_vis_reco", 'F');
    dataloader->AddVariable("p3_hardest_reco", 'F');
    dataloader->AddVariable("pTmiss_mag_reco", 'F');
    dataloader->AddVariable("pTabs_sum_reco", 'F');

    dataloader->AddVariable("nTrkTanThetaLeq0point1_reco", 'F');

    // global event weights per tree (see below for setting event-wise weights)
    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;
 
    // You can add an arbitrary number of signal or background trees
    dataloader->AddSignalTree (signalTree0, signalWeight);
    dataloader->AddSignalTree (signalTree1, signalWeight);
    dataloader->AddSignalTree (signalTree2, signalWeight);
    dataloader->AddBackgroundTree(background0, backgroundWeight);

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycuts = "n_ch >= 3"; // vertex cut
    TCut mycutb = "n_ch >= 3"; // vertex cut

    dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,"SplitMode=random:!V" );

    // Book MVA methods
    // Gradient Boosted Decision Trees
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
                        "!H:!V:NTrees=2000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=1000" );


    // Train MVAs using the set of training events
    factory->TrainAllMethods();
 
    // Evaluate all MVAs using the set of test events
    factory->TestAllMethods();
 
    // Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // Save the output  
    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;
 
    delete factory;
    delete dataloader;
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

}