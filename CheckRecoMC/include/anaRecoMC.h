#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "tracks.h"
#include "NuMCTruth_kinematics.h"
#include "line3Dfit.h"
#include "CheckRecoPlot.h"

const double IntLumi_run3 = 250.00; // fb^
const double wgt_NC_200025 = 50173.00 / 10e3; // event/fb^{-1}
const double wgt_NC_200026 = 66355.00 / 10e3;
const double wgt_NC_200035 = 13555.00 / 10e3;

const float t_plate = 1.1e3; // tungsten plate thickness in um
const float t_film = 0.34e3; // film thickness in um

const int kMaxs = 79;

// read in a list of evt_ID_plStart_plEnd
void ReadEventList(std::vector<std::string> *EvtList, std::string FilePathName){

    std::ifstream file(FilePathName);
    if (!file.is_open()) {
        std::cerr << "Unable to open file." << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        EvtList->push_back(line);
    }

    file.close();

}

void PrintEventList(std::vector<std::string> * EvtList){

    for(int EvtIt=0; EvtIt<EvtList->size(); EvtIt++){
        std::cout<<EvtList->at(EvtIt)<<std::endl;
    }

}

void Global2HitCoordinate(float &vx_hit, float &vy_hit, float &vz_hit, float m_vx_global, float m_vy_global, float m_vz_global){

    vx_hit = (m_vx_global - 10.) * 1000.; // micro-metre
    vy_hit = (m_vy_global + 21.) * 1000.;
    vz_hit = (m_vz_global + 2986.27) * 1000.;

}

bool IsPrimTrackCandidate(float vx_hit, float vy_hit, float vz_hit, float s_eX_first, float s_eY_first, float s_eZ_first, float tanThetaCut){

    float dz = s_eZ_first - vz_hit;
    bool StartWithin3films = dz < (3. * (t_plate + t_film)) && dz > 0; // start within 3 films downstream

    float dr = std::sqrt(std::pow(s_eX_first-vx_hit,2) + std::pow(s_eY_first-vy_hit,2));
    bool TransConfine = dr < (std::abs(dz)*tanThetaCut);

    return StartWithin3films && TransConfine;
    //return StartWithin3films;

}

double GetImpactParameter(TVector3 TrkDirVect, double p0, double p2, float vx_hit, float vy_hit, float vz_hit){
    // p0 and p2 are line 3D fitting parameters

    TVector3 xv(vx_hit, vy_hit, vz_hit);
    TVector3 x0(p0, p2, 0. );
    TVector3 u = TrkDirVect.Unit();

    double IP = std::sqrt( ((xv-x0).Cross(u)).Mag2() );
    return IP;

}

void GetTrackParam(int nseg, std::vector<float> &X, std::vector<float> &Y, std::vector<float> &Z, std::vector<float> &TX, std::vector<float> &TY,
                   float *s_eX, float *s_eY, float *s_eZ, float *s_eTX, float *s_eTY, float *s_eProb){
    
    for(int SegIt=0; SegIt<nseg; SegIt++){

        if(s_eProb[SegIt] < 1e-10) break;

        X.push_back(s_eX[SegIt]);
        Y.push_back(s_eY[SegIt]);
        Z.push_back(s_eZ[SegIt]);
        TX.push_back(s_eTX[SegIt]);
        TY.push_back(s_eTY[SegIt]);

        if(SegIt<nseg-1){
            float dz = s_eZ[SegIt+1] - s_eZ[SegIt];
            if(dz < 0 || dz > 3.*(t_plate + t_film)) break; // remove back-scattered and gapped part
        }

    }

}

void GetTrackMCtruth(int nseg, std::vector<int> &TrackID, std::vector<int> &PDG, std::vector<int> &EventID, std::vector<float> &pmag,
                     int *s_eTrack, int *s_eMCTrack, int *s_eMCEvt, float *s_eP, float *s_eProb){
    
    for(int SegIt=0; SegIt<nseg; SegIt++){

        if(s_eProb[SegIt] < 1e-10) break;

        TrackID.push_back(s_eTrack[SegIt]);
        PDG.push_back(s_eMCTrack[SegIt]);
        EventID.push_back(s_eMCEvt[SegIt] - 100000);
        pmag.push_back(s_eP[SegIt]);

    }

}


void ana_linked_tracks(tracks *linked_tracks, float vx_hit, float vy_hit, float vz_hit,
                       std::vector<int> &trackID, int &n_ch, std::vector<double> &tan_theta, std::vector<double> &phi_angle,
                       std::vector<int> &nseg, std::vector<float> &pmag_MCtruth, std::vector<double> &ImpactParameter, std::vector<float> &DeltaZ){

    int NTracks = linked_tracks->fChain->GetEntries();

    // loop tracks in the linked_tracks.root
    for(int TrkIt=0; TrkIt<NTracks; TrkIt++){

        linked_tracks->GetEntry(TrkIt);
        std::cout<<"TrackIt = "<<TrkIt<<"; nseg = "<<linked_tracks->nseg<<std::endl;

        if(linked_tracks->nseg < 3) continue; // # segments >= 3

        std::vector<float> X, Y, Z, TX, TY; // segment parameters
        GetTrackParam(linked_tracks->nseg, X, Y, Z, TX, TY, linked_tracks->s_eX, linked_tracks->s_eY, linked_tracks->s_eZ, linked_tracks->s_eTX, linked_tracks->s_eTY, linked_tracks->s_eProb);

        //if(!IsPrimTrackCandidate(vx_hit, vy_hit, vz_hit, X[0], Y[0], Z[0], 0.5)) continue; // start within 3 films downstream and tanTheta cut @ 0.5

        std::vector<int> TrackID, PDG, EventID; // segment MC truth
        std::vector<float> pmag; // segment MC true momentum

        GetTrackMCtruth(linked_tracks->nseg, TrackID, PDG, EventID, pmag, linked_tracks->s_eTrack, linked_tracks->s_eMCTrack, linked_tracks->s_eMCEvt, linked_tracks->s_eP, linked_tracks->s_eProb);


        for(int it=0; it<Z.size(); it++){
            //if (val < 1e-10) break;
            std::cout<<"(X, Y, Z) = "<<"("<<X[it]<<", "<<Y[it]<<", "<<Z[it]<<") um; Track ID = "<<TrackID[it]<<"; PDG = "<<PDG[it]<<"; P = "<<pmag[it]<<" GeV/c"<<std::endl;
        }
      
        // track 3D linear fitting
        const double *para_line3Dfit = line3Dfit(X, Y, Z);
/*      
        std::cout<<"p0 = "<<para_line3Dfit[0]<<std::endl;
        std::cout<<"p1 = "<<para_line3Dfit[1]<<std::endl;
        std::cout<<"p2 = "<<para_line3Dfit[2]<<std::endl;
        std::cout<<"p3 = "<<para_line3Dfit[3]<<std::endl;
*/ 
        TVector3 TrkDirVect(para_line3Dfit[1], para_line3Dfit[3], 1.00);

        double IP = GetImpactParameter(TrkDirVect, para_line3Dfit[0], para_line3Dfit[2], vx_hit, vy_hit, vz_hit); // impact parameter
        double tanTheta = std::tan(TrkDirVect.Theta());
        double phi = TrkDirVect.Phi();

        trackID.push_back(TrackID[0]); // take the first segment track ID as the whole reco. track ID
        phi_angle.push_back(phi);
        tan_theta.push_back(tanTheta);
        nseg.push_back(linked_tracks->nseg);
        pmag_MCtruth.push_back(pmag[0]); // take the first segment MC truth momentum as the whole reco. track MC truth momentum
        n_ch++;
        ImpactParameter.push_back(IP);
        DeltaZ.push_back(Z[0] - vz_hit);

//        std::cout<<"Track ID = "<<TrackID[0]<<"; tan theta = "<<tanTheta<<"; phi = "<<phi<<" rad;"<<std::endl;   
//        std::cout<<std::endl;

    } // track loop
        
}


void sort_linked_tracks(std::vector<int> &trackID, std::vector<double> &tan_theta, std::vector<double> &phi_angle, 
                        std::vector<int> &nseg, std::vector<float> &pmag, std::vector<double> &ImpactParameter,
                        std::vector<float> &DeltaZ){

    // Create an index vector
    std::vector<size_t> indices(trackID.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }

    // Sort indices based on trackID
    std::sort(indices.begin(), indices.end(), [&trackID](size_t i1, size_t i2) {
        return trackID[i1] < trackID[i2];
    });

    // Rearrange vectors based on sorted indices
    std::vector<int> sorted_trackID(trackID.size()), sorted_nseg(nseg.size()); 
    std::vector<double> sorted_tan_theta(tan_theta.size()), sorted_phi_angle(phi_angle.size()), sorted_IP(ImpactParameter.size());
    std::vector<float> sorted_pmag(pmag.size()), sorted_dz(DeltaZ.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        sorted_trackID[i] = trackID[indices[i]];
        sorted_tan_theta[i] = tan_theta[indices[i]];
        sorted_phi_angle[i] = phi_angle[indices[i]];
        sorted_nseg[i] = nseg[indices[i]];
        sorted_pmag[i] = pmag[indices[i]];
        sorted_IP[i] = ImpactParameter[indices[i]];
        sorted_dz[i] = DeltaZ[indices[i]];
    }

    trackID = sorted_trackID;
    tan_theta = sorted_tan_theta;
    phi_angle = sorted_phi_angle;
    nseg = sorted_nseg;
    pmag = sorted_pmag;
    ImpactParameter = sorted_IP;
    DeltaZ = sorted_dz;

}

void ana_jw_test(NuMCTruth_kinematics *jw_test, 
                std::vector<int> &trackID, int &n_ch, std::vector<double> &tan_theta, std::vector<double> &phi_angle,
                std::vector<double> &pz){

    jw_test->GetEntry(0);
    int NHad = jw_test->m_hadrons_track_id->size();

    // loop tracks in the jw_test.root
    for(int HadIt=0; HadIt<NHad; HadIt++){

        if(jw_test->m_hadrons_ch->at(HadIt) == 0) continue; // skip neutral hadrons

        TLorentzVector p4had(jw_test->m_hadrons_px->at(HadIt), jw_test->m_hadrons_py->at(HadIt), jw_test->m_hadrons_pz->at(HadIt), jw_test->m_hadrons_e->at(HadIt));
        double tanTheta = std::tan(p4had.Vect().Theta());
        //if(tanTheta > 0.5 || tanTheta < 0) continue; // forward angular cut
        double phi = p4had.Vect().Phi();

        double Ekin = p4had.E() - p4had.M(); // kinetic energy
        //if(Ekin < 1) continue; // kinetic energy cut
        double pmag = p4had.Vect().Mag(); // the magnitude of the momentum
        //if(pmag < 1.0) continue;

        trackID.push_back(jw_test->m_hadrons_track_id->at(HadIt));
        tan_theta.push_back(tanTheta);
        phi_angle.push_back(phi);
        if(p4had.Pz() > 1.0 && tanTheta < 0.5) n_ch++;
        pz.push_back(p4had.Pz());
/*
        std::cout<<"Track ID = "<<jw_test->m_hadrons_track_id->at(HadIt)
        <<"; tan theta = "<<tanTheta
        <<"; phi = "<<phi<<" rad"
        <<"; PDG = "<<jw_test->m_hadrons_PDG->at(HadIt)
        <<"; P = "<<pmag<<" GeV/c"
        <<"; Pz = "<<p4had.Pz()<<" GeV/c"
        <<"; m = "<<p4had.M()<<" GeV/c^2"
        <<std::endl;
*/
    }

}



void ana_MC_event(std::string EventDir, std::vector<std::string> *EventList, double EventWeight, CheckRecoPlot *CRP){
    
    int NEvts = EventList->size();
    std::cout<<"# events = "<<NEvts<<" under "<<EventDir<<std::endl;

    // loop the event
    for(int EvtIt=0; EvtIt<1; EvtIt++){

        std::string dir_linked_tracks = EventDir + EventList->at(EvtIt) + "/linked_tracks.root";
        std::string dir_jw_test = EventDir + EventList->at(EvtIt) + "/jw_test.root";

        TFile *f_linked_tracks = new TFile(dir_linked_tracks.c_str(), "READ");
        TTree *t_linked_tracks = (TTree *)f_linked_tracks->Get("tracks");
        tracks *linked_tracks = new tracks(t_linked_tracks);
        int NTracks = linked_tracks->fChain->GetEntries();
        std::cout<<"# tracks of "<<dir_linked_tracks<<" = "<<NTracks<<std::endl;

        TFile *f_jw_test = new TFile(dir_jw_test.c_str(), "READ");
        TTree *t_jw_test = (TTree *)f_jw_test->Get("NuMCTruth_kinematics");
        NuMCTruth_kinematics *jw_test = new NuMCTruth_kinematics(t_jw_test);
        jw_test->GetEntry(0);
        //int NHad = jw_test->m_hadrons_PDG->size();
        //jw_test->Show();

        Float_t vx_hit, vy_hit, vz_hit; // vertex position in hit coordinate
        Global2HitCoordinate(vx_hit, vy_hit, vz_hit, jw_test->m_vx, jw_test->m_vy, jw_test->m_vz);

        std::cout<<"vx = "<<vx_hit<<" um"<<std::endl;
        std::cout<<"vy = "<<vy_hit<<" um"<<std::endl;
        std::cout<<"vz = "<<vz_hit<<" um"<<std::endl;

        int n_ch[2] = {0, 0}; // 0: reco, 1: true
        std::vector<int> trackID[2], nseg; // nseg: the number of segments of the reconstructed tracks
        std::vector<double> tan_theta[2]; // from 3D line fit
        std::vector<double> phi_angle[2]; // from 3D line fit
        std::vector<double> ImpactParameter; // from 3D line fit
        std::vector<double> pz; // from jw_test MC truth
        std::vector<float> pmag_MCtruth; // the magnitude of the momentum from the MC truth
        std::vector<float> DeltaZ; // the z-distance between the vertex and the first segment of a track
        
//        std::cout<<__LINE__<<std::endl;
        ana_linked_tracks(linked_tracks, vx_hit, vy_hit, vz_hit, trackID[0], n_ch[0], tan_theta[0], phi_angle[0], nseg, pmag_MCtruth, ImpactParameter, DeltaZ);
//        std::cout<<__LINE__<<std::endl;
        sort_linked_tracks(trackID[0], tan_theta[0], phi_angle[0], nseg, pmag_MCtruth, ImpactParameter, DeltaZ);
//        std::cout<<__LINE__<<std::endl;
        
        for(int recoIt=0; recoIt<trackID[0].size(); recoIt++){
            std::cout<<"Track ID = "<<trackID[0].at(recoIt)
            <<"; tan theta = "<<tan_theta[0].at(recoIt)
            <<"; phi = "<<phi_angle[0].at(recoIt)<<" rad"
            <<"; nseg = "<<nseg.at(recoIt)
            <<"; P = "<<pmag_MCtruth.at(recoIt)<<" GeV/c"
            <<"; IP = "<<ImpactParameter.at(recoIt)<<" um"
            <<"; dz = "<<DeltaZ.at(recoIt)<<" um"
            <<std::endl;
        }
        

        ana_jw_test(jw_test, trackID[1], n_ch[1], tan_theta[1], phi_angle[1], pz);
//        std::cout<<"n_ch = "<<n_ch[0]<<", "<<n_ch[1]<<std::endl;

        CRP->h1_n_ch[0]->Fill(n_ch[0]);
        CRP->h1_n_ch[1]->Fill(n_ch[1]);
        
        int trueIt=0;
        for(int recoIt=0; recoIt<trackID[0].size(); recoIt++){

            while(trackID[1].at(trueIt) < trackID[0].at(recoIt)){
                trueIt++;
                if(trueIt >= trackID[1].size()) break;
            }
            if(trueIt >= trackID[1].size()) break;

            if(trackID[1].at(trueIt) == trackID[0].at(recoIt)){
                // the track ID is matched between reco and truth

                CRP->h1_dtan_theta->Fill(tan_theta[0].at(recoIt) - tan_theta[1].at(trueIt), EventWeight);
                CRP->h1_dphi->Fill(phi_angle[0].at(recoIt) - phi_angle[1].at(trueIt), EventWeight);
                CRP->h2_nseg_pz->Fill(pz.at(trueIt), nseg.at(recoIt), EventWeight);
                CRP->h2_dphi_pz->Fill(pz.at(trueIt), phi_angle[0].at(recoIt) - phi_angle[1].at(trueIt), EventWeight);
                CRP->h2_dtan_theta_pz->Fill(pz.at(trueIt), tan_theta[0].at(recoIt) - tan_theta[1].at(trueIt), EventWeight);

            }

        }
        
        f_linked_tracks->Close();
        f_jw_test->Close();
        delete f_linked_tracks;
        delete f_jw_test;

    }

}

