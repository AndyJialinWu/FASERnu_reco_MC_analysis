#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "AnaChHadMom.h"
#include "PhysicsNTUP.C"

void AnaChHadMom(){


    TFile *f_disc = new TFile("PhysicsNTUP_ML.root", "READ");
    TTree *t_TrainTest[4];                       // 0: 200025, 1: 200026, 2: 200035, 3: 100069
    PhysicsNTUP *disc[4];                        // 0: 200025, 1: 200026, 2: 200035, 3: 100069
    for(int RunIt=0; RunIt<4; RunIt++){

        t_TrainTest[RunIt] = (TTree*)f_disc->Get( Form("disc%d_TrainTest", RunNumber[RunIt]) );
        SetNumRecoNC(t_TrainTest[RunIt], RunIt); // sequence matters ! Place this function before setting the branch address !

        disc[RunIt] = new PhysicsNTUP(t_TrainTest[RunIt]);

    }

    wgt_calc();
    Hist_Init();

    // MC Run Loop
    for(int RunIt=0; RunIt<3; RunIt++){
        // neutrino NC MC
        

        // Event Loop
        for(int EvtIt=0; EvtIt<disc[RunIt]->fChain->GetEntries(); EvtIt++){

            disc[RunIt]->GetEntry(EvtIt);   // Get the event
            int nChHadMomSuccess = 0, nChHadMomFailure1 = 0, nChHadMomFailure2 = 0;
            
            //std::cout<<"RunNumber = "<<RunNumber[RunIt]<<"; n_ch = "<<disc[RunIt]->n_ch<<"; weight = "<<weight[RunIt]<<std::endl;
            
            // Track Loop
            for(int TrkIt=0; TrkIt<disc[RunIt]->pmag_reco->size(); TrkIt++){
                
                const int PDG = disc[RunIt]->PDG->at(TrkIt);
                if( std::abs(PDG) != kElectron ){
                    // hadron

                    int MomStatus = -1;
                    if(disc[RunIt]->pmag_haruhi->at(TrkIt) > 3500){
                        // Failure1
                        MomStatus = 1;
                        nChHadMomFailure1++;
                    }
                    else if(disc[RunIt]->pmag_haruhi->at(TrkIt) < 20 && disc[RunIt]->pmag_true->at(TrkIt) > 40){
                        //Failure2
                        MomStatus = 2;
                        nChHadMomFailure2++;
                    }
                    else{
                        // Success
                        MomStatus = 0;
                        nChHadMomSuccess++;
                        h2_preco_ptrue->Fill(disc[RunIt]->pmag_true->at(TrkIt), disc[RunIt]->pmag_haruhi->at(TrkIt));
                        double p_reso = (disc[RunIt]->pmag_reco->at(TrkIt) - disc[RunIt]->pmag_true->at(TrkIt)) / disc[RunIt]->pmag_true->at(TrkIt);
                        h2_preso_ptrue->Fill(disc[RunIt]->pmag_true->at(TrkIt), p_reso);
                        h2_nsegm_ptrue->Fill(disc[RunIt]->pmag_true->at(TrkIt), disc[RunIt]->nseg->at(TrkIt));
                        h2_preso_nsegm->Fill(disc[RunIt]->nseg->at(TrkIt), p_reso);
                    }

                    h2_MomStatus_nseg->Fill(disc[RunIt]->nseg->at(TrkIt), MomStatus);
                    if(nChHadMomSuccess + nChHadMomFailure1 + nChHadMomFailure2 > 3){
                        double MomSuccessRateEventLevel = double(nChHadMomSuccess) / double(nChHadMomSuccess + nChHadMomFailure1 + nChHadMomFailure2);
                        h2_SuccessRate_NuE->Fill(disc[RunIt]->Nu_e, MomSuccessRateEventLevel);
                    }
                    

                }
    
            }// end of track loop
    
        }// end of event loop

    }// end of MC Run Loop
    

    GetSuccessRate();


    StoreHist2ROOT("ChgHadMomPlots.root");

    //TString FigAddress = "Figures/discriminators/";
    //StoreHist2PDF(FigAddress);

    f_disc->Close();
    delete f_disc;

}
