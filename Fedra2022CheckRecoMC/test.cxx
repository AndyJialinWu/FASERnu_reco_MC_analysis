#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

//#if ROOT_VERSION_CODE < 6 //headers should not be included in ROOT6, loaded with .pcm files
//#include "EdbSegP.h"
//#include "EdbMomentumEstimator.h"
//#include "EdbPattern.h"
//#endif

void test(){
    
    TString fname = "/home/jialwu/raid/FASERnu_MC_reco/NC/200025/evt_10001_pl1_94/linked_tracks.root";
    TFile *trfile = new TFile(fname, "READ");
    TTree *tree = (TTree *)trfile->Get("tracks");

    int NTracks = tree->GetEntries();
    std::cout<<"# volume tracks = "<<NTracks<< std::endl;

    Int_t   nseg=0;
    EdbSegP *trk=0;
    TClonesArray *seg  = new TClonesArray("EdbSegP", 79);
    tree->SetBranchAddress("nseg", &nseg);
    tree->SetBranchAddress("t.", &trk);
    tree->SetBranchAddress("s",  &seg);

    // track loop
    for(int trkIt=0; trkIt<NTracks; trkIt++){

        tree->GetEntry(trkIt);

        std::cout<<"TrackIt = "<<trkIt<<"; ";
        std::cout<<"Track ID = "<<trk->Track()<<"; ";
        std::cout<<"PDG = "<<trk->MCTrack()<<"; ";
        std::cout<<"P = "<<trk->P()<<" GeV/c; ";
        std::cout<<"# seg. = "<<nseg<<"; ";
        std::cout<<"Z0 = "<<trk->Z()<<" um;"<<std::endl;

        if(nseg < 3) continue; // # segments cut

        EdbTrackP *trkP = new EdbTrackP();

        // segment loop
        for(int segIt=0; segIt<nseg; segIt++){

            EdbSegP *s = (EdbSegP *)seg->At(segIt);
            trkP->AddSegment(s);

        }

        trkP->SetNpl();
        EdbMomentumEstimator *MomEst = new EdbMomentumEstimator();
        MomEst->SetParPMS_Mag();

        float eP = MomEst->PMScoordinate(*trkP);
        std::cout<<"eP (coord. method) = "<<eP<<" GeV/c"<<std::endl;
    }


}