#include <iostream>
#include <cmath>

#include "CompDataMC.h"


void CompDataMC(){

    TFile *f_data = new TFile("PhysicsNTUP_ML_SimonThor.root");
    TTree *t_data = (TTree*)f_data->Get("disc_RealData");

    TFile *f_mc = new TFile("PhysicsNTUP_ML.root");
    TTree *t_mc = (TTree*)f_mc->Get("disc100069_TrainTest");

    PhysicsNTUP *data = new PhysicsNTUP(t_data);
    PhysicsNTUP *mc = new PhysicsNTUP(t_mc);

    HistInit();

    // data loop
    for(int DataIt=0; DataIt<data->fChain->GetEntries(); DataIt++){
        
        data->GetEntry(DataIt);
        if(DataIt%10==0) std::cout << "Data Entry: " << DataIt << std::endl;
        HistFill(data, 1);
        
    }

    // mc loop
    for(int MCIt=0; MCIt<mc->fChain->GetEntries(); MCIt++){
        
        mc->GetEntry(MCIt);
        if(MCIt%1000==0) std::cout << "MC Entry: " << MCIt << std::endl;

        // neutron energy cut
        if(mc->Nu_e < 80) continue;

        HistFill(mc, 0);
        
    }


    StoreHist2ROOT("Disc_Data_MC.root");
    StoreHist2PDF("Figures/BDT_input/");

}

