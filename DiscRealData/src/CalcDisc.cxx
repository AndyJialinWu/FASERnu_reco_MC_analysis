#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "LinkedTracks.h"
#include "AnaRealData.h"
#include "CalcDisc.h"

int main(){

    std::string RealDataAddress = "../AnaFeedBackFiles/linked_tracks/";
    std::vector<std::string> LinkedTracksFileNames = GetLinkedTracksFileNames(RealDataAddress);

    TFile *f_disc = new TFile("PhysicsNTUP_ML.root", "RECREATE");
    TTree *t_disc = new TTree("disc_RealData", "disc_RealData");
    Branch(t_disc);
    
    // Loop over the events (i.e. neutral vertices, one event per file)
    for(size_t EvtIt=0; EvtIt<LinkedTracksFileNames.size(); EvtIt++){

        EventInit();
        Discriminators *disc = new Discriminators();

        // to be placed in the event loop
        TString LTname = RealDataAddress + LinkedTracksFileNames[EvtIt];
        std::cout << "Processing file: " << LTname << std::endl;
        // Open the file and get the tree
        TFile *f_trfile = new TFile(LTname, "READ");
        TTree *t_tracks = (TTree *)f_trfile->Get("tracks");
        TTree *t_vertex = (TTree *)f_trfile->Get("vertex");
        LinkedTracks *RealData = new LinkedTracks(t_tracks);
        if(t_tracks->GetEntries() == 0) continue;

        // Get primary vertex info
        float vx, vy, vz;
        t_vertex->SetBranchAddress("vx", &vx);
        t_vertex->SetBranchAddress("vy", &vy);
        t_vertex->SetBranchAddress("vz", &vz);
        t_vertex->GetEntry(0);
        RealData->SetRecoVertex(vx, vy, vz);
        RealData->GetRecoMCInfo(true);

        // Loop over the tracks
        for(size_t recoIt=0; recoIt<RealData->trackID.size(); recoIt++){

            if(RealData->IsPrimTrack(recoIt)){
                // cut by IsPrimTrack function

                disc->px_reco.push_back(RealData->pmag_reco.at(recoIt) * std::sin(RealData->theta.at(recoIt)) * std::cos(RealData->phi.at(recoIt)));
                disc->py_reco.push_back(RealData->pmag_reco.at(recoIt) * std::sin(RealData->theta.at(recoIt)) * std::sin(RealData->phi.at(recoIt)));
                disc->pz_reco.push_back(RealData->pmag_reco.at(recoIt) * std::cos(RealData->theta.at(recoIt)));
                disc->pmag_reco.push_back(RealData->pmag_reco.at(recoIt));
                disc->phi_reco.push_back(RealData->phi.at(recoIt));
                disc->theta_reco.push_back(RealData->theta.at(recoIt));

                disc->pmag_ang.push_back(RealData->pmag_ang.at(recoIt));
                disc->pmag_coord.push_back(RealData->pmag_coord.at(recoIt));
                disc->pmag_haruhi.push_back(RealData->pmag_haruhi.at(recoIt));
                
                disc->TrackID.push_back(RealData->trackID.at(recoIt));
                disc->PDG.push_back(0);
                disc->nseg.push_back(RealData->NSeg.at(recoIt));
                disc->TrackLength.push_back(RealData->TrackLength.at(recoIt));
                disc->dz.push_back(RealData->dz.at(recoIt)); // um
                disc->IP.push_back(RealData->IP.at(recoIt)); // um
                disc->PID_start.push_back(RealData->PID.at(recoIt).at(0));
                disc->PID_end.push_back(RealData->PID.at(recoIt).back());
                disc->itrk.push_back(RealData->itrk.at(recoIt));
                disc->MaxGap.push_back(RealData->MaxGap.at(recoIt));
                disc->theta_RMS.push_back(RealData->theta_RMS.at(recoIt));
                disc->MaxKinkAngle.push_back(RealData->MaxKinkAngle.at(recoIt));
                disc->IsPartCatChanged.push_back(RealData->IsPartCatChanged.at(recoIt));
                disc->IsTrkIdChanged.push_back(RealData->IsTrkIdChanged.at(recoIt));

            }
        }

        // event level info.
        disc->mcID = 0;
        disc->EventID = RealData->eventID;
        Nu_PDG = -999;
        Nu_px = -999;
        Nu_py = -999;
        Nu_pz = -999;
        Nu_e = -999;

        // calculate discriminators
        disc->CalcDisc();
        PassDiscAddressToTTreeAddress(disc);
        f_disc->cd();
        t_disc->Fill();

    }// event loop

    f_disc->cd();
    t_disc->Write();
    f_disc->Close();
    delete f_disc;

    return 0;

}