#include <vector>
#include <cmath>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TPDGCode.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "EdbSegP.h"
#include "EdbMomentumEstimator.h"
#include "EdbPattern.h"
#include "EdbSegCouple.h"
#include "EdbPVRec.h"

#include "line3Dfit.h"
#include "FnuMomCoord.hpp"

const double IntLumi_run3 = 250.00;           // fb^{-1}
const double eff_FV = 638.00 / 730.00;
const double wgt_NC_200025 = 50173.00 * eff_FV / 10e3; // unit: event/fb^{-1}
const double wgt_NC_200026 = 66355.00 * eff_FV / 10e3;
const double wgt_NC_200035 = 13555.00 * eff_FV / 10e3;
const double Run3ExpNC[3] = {IntLumi_run3 * wgt_NC_200025, IntLumi_run3 * wgt_NC_200026, IntLumi_run3 * wgt_NC_200035};

const float t_plate = 1.1e3;                  // tungsten plate thickness in um
const float t_film = 0.34e3;                  // film thickness in um
//const int kMaxTrks = 79;

// event based
class LinkedTracks{
    public:
        
        TTree *tree;                    //!pointer to the analyzed TTree or TChain

        //  event level
        int n_ch;                       // primary charged tracks multiplicity
        float vx_hit_reco;              // reconstructed vertex position
        float vy_hit_reco;              // reconstructed vertex position
        float vz_hit_reco;              // reconstructed vertex position 
        float vx_hit_true;              // MC truth vertex position
        float vy_hit_true;              // MC truth vertex position
        float vz_hit_true;              // MC truth vertex position
        int eventID;                    // MC event ID

        // track level  
        std::vector<int> trackID;           // MC truth track ID
        std::vector<int> PDG;               // MC truth track pdg code
        std::vector<int> itrk;              // reconstructed track iterator
        std::vector<double> theta;          // reconstructed track polar angle
        std::vector<double> phi;            // reconstructed track azimuthal angle
        std::vector<float> pmag_reco;       // reconstructed track momentum
        std::vector<float> pmag_coord;      // reconstructed track momentum by MCS the coordinate method, FEDRA2022
        std::vector<float> pmag_ang;        // reconstructed track momentum by MCS the angular method, FEDRA2022
        std::vector<float> pmag_haruhi;     // reconstructed track momentum by MCS the coordinate method, Haruhi
        std::vector<float> ptrue;           // MC truth track momentum
        std::vector<int> NSeg;              // reconstructed track the number of segments
        std::vector<double> IP;             // reconstructed impact parameter um
        std::vector<float> dz;              // reconstructed (z_FirstSeg - z_vertex) um
        std::vector<double> TrackLength;    // reconstructed track length mm
        std::vector<int> MaxGap;            // reconstructed track the largest gap between 2 consecutive segments
        std::vector<double> theta_RMS;      // reconstructed track the standard deviation of the segment polar angle
        std::vector<double> MaxKinkAngle;   // reconstructed track the largest kink angle between 2 consecutive cells by positions
        std::vector<bool> IsPartCatChanged; // reconstructed track the segment particle category changed or not
        std::vector<bool> IsTrkIdChanged;   // reconstructed track the segment track ID changed or not
                
        // segment level
        std::vector<std::vector<float>> X;              // reconstructed segment positions
        std::vector<std::vector<float>> Y;              // reconstructed segment positions
        std::vector<std::vector<float>> Z;              // reconstructed segment positions 
        std::vector<std::vector<float>> TX;             // reconstructed segment tangent polar angles
        std::vector<std::vector<float>> TY;             // reconstructed segment tangent polar angles
        std::vector<std::vector<int>>   PID;            // reconstructed segment plate ID: from 0 to 100
        std::vector<std::vector<int>>   SegPDG;         // reconstructed segment pdg code
        std::vector<std::vector<int>>   SegTrkID;       // reconstructed segment track ID


    public:
        LinkedTracks(TTree *tracks=0);
        virtual void Init(TTree *tracks);
        
        void SetTrueVertex(float vx_global_MC, float vy_global_MC, float vz_global_MC); // MC truth (global coord.) to hit coord.
        void SetRecoVertex(float vx_hit_reco, float vy_hit_reco, float vz_hit_reco);    // Real Data hit coord.
        void PrintTrueVertex();

        // set eventID, trackID, PDG, pmag, ptrue, NSeg, X, Y, Z, TX, TY
        //void SetRecoMCInfo(int FedraEvtID, std::vector<int> FedraTrkID, std::vector<int> FedraPDG, std::vector<float> FedraPmagCoord, std::vector<float> FedraPmagAng, 
        //                   std::vector<float> FedraPtrue, std::vector<int> FedraNSeg,
        //                   std::vector<std::vector<float>> FedraX, std::vector<std::vector<float>> FedraY, std::vector<std::vector<float>> FedraZ,
        //                   std::vector<std::vector<float>> FedraTX, std::vector<std::vector<float>> FedraTY); 
        
        static double CalcThetaRMS(std::vector<float> TX, std::vector<float> TY);
        static double CalcMaxKinkAngle(std::vector<float> X, std::vector<float> Y, std::vector<float> Z);
        static bool IsValueChanged(std::vector<int> vec);
        static float GetPmagReco(float ePharuhi, float ePcoord, float ePang, float ptrue, int PDG, int nseg);
        void GetRecoMCInfo(bool ToPrintTrkInfo);
        void PrintTrkInfo(EdbSegP *trk, int TrackIt, int nseg, float epCoord, float epAng, float epHaruhi);
        
        static double CalcTrackLength(std::vector<float> x, std::vector<float> y, std::vector<float> z);
        //void GetRecoAngleIPdz();  // get theta, phi, IP, dz, TrackLength

        bool IsPrimTrack(int TrackIt);

        void SortBasedOnTrackID(bool PrintTrkInfo);// sorting ascendingly based on the track ID

        virtual ~LinkedTracks();
};

LinkedTracks::LinkedTracks(TTree *tracks) : tree(0)
{

    if(tracks == 0){
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("linked_tracks.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("linked_tracks.root");
        }
        f->GetObject("tracks",tree);
    }
    Init(tracks);

};

void LinkedTracks::Init(TTree *tracks){
    // Set branch addresses and branch pointers
    if (!tracks) return;

    tree = tracks;

    n_ch = 0;
    eventID = -9999;
    vx_hit_true = -9999.;
    vy_hit_true = -9999.;
    vz_hit_true = -9999.;
    vx_hit_reco = -9999.;
    vy_hit_reco = -9999.;
    vz_hit_reco = -9999.;

    trackID.clear();
    PDG.clear();
    itrk.clear();
    theta.clear();
    phi.clear();
    pmag_reco.clear();
    pmag_haruhi.clear();
    pmag_coord.clear();
    pmag_ang.clear();
    NSeg.clear();
    ptrue.clear();
    IP.clear();
    dz.clear();
    TrackLength.clear();
    MaxGap.clear();
    theta_RMS.clear();
    MaxKinkAngle.clear();
    IsPartCatChanged.clear();
    IsTrkIdChanged.clear();
    
    
    X.clear();
    Y.clear();
    Z.clear();
    TX.clear();
    TY.clear();
    PID.clear();
    SegPDG.clear();
    SegTrkID.clear();

};



LinkedTracks::~LinkedTracks(){
    if (!tree) return;
    delete tree->GetCurrentFile();
};


void LinkedTracks::PrintTrueVertex(){
    std::cout<<"(vx_hit, vy_hit, vz_hit) = "<<"("<<vx_hit_true/1e3<<", "<<vy_hit_true/1e3<<", "<<vz_hit_true/1e3<<")"<<" mm"<<std::endl;
};

void LinkedTracks::SetTrueVertex(float vx_global_MC, float vy_global_MC, float vz_global_MC){
    
    vx_hit_true = (vx_global_MC - 10.) * 1000.; // micro-metre
    vy_hit_true = (vy_global_MC + 21.) * 1000.;
    vz_hit_true = (vz_global_MC + 2986.27) * 1000.;

    PrintTrueVertex();

};

void LinkedTracks::SetRecoVertex(float vx_hit_reco, float vy_hit_reco, float vz_hit_reco){
    
    vx_hit_true = vx_hit_reco;
    vy_hit_true = vy_hit_reco;
    vz_hit_true = vz_hit_reco;

    std::cout<<"(vx_hit, vy_hit, vz_hit) = "<<"("<<vx_hit_reco<<", "<<vy_hit_reco<<", "<<vz_hit_reco<<")"<<" um"<<std::endl;

};


/*
void LinkedTracks::SetRecoMCInfo(int FedraEvtID, std::vector<int> FedraTrkID, std::vector<int> FedraPDG, std::vector<float> FedraPmagCoord, std::vector<float> FedraPmagAng, 
                                 std::vector<float> FedraPtrue, std::vector<int> FedraNSeg,
                                 std::vector<std::vector<float>> FedraX, std::vector<std::vector<float>> FedraY, std::vector<std::vector<float>> FedraZ,
                                 std::vector<std::vector<float>> FedraTX, std::vector<std::vector<float>> FedraTY){

    eventID = FedraEvtID;

    // track level
    trackID = FedraTrkID;
    PDG = FedraPDG;
    pmag_coord = FedraPmagCoord;
    pmag_ang = FedraPmagAng;
    ptrue= FedraPtrue;
    NSeg = FedraNSeg;

    // segment level
    X = FedraX;
    Y = FedraY;
    Z = FedraZ;
    TX = FedraTX;
    TY = FedraTY;

};
*/

void LinkedTracks::PrintTrkInfo(EdbSegP *trk, int TrackIt, int nseg, float epCoord, float epAng, float epHaruhi){

    std::cout<<"itrk = "<<TrackIt<<"; ";
    std::cout<<"Track ID = "<<trk->Track()<<"; ";
    std::cout<<"PDG = "<<trk->MCTrack()<<"; ";
    std::cout<<"P = "<<trk->P()<<" GeV/c; ";
    std::cout<<"Pcoord = "<<epCoord<<" GeV/c; ";
    std::cout<<"Pang = "<<epAng<<" GeV/c; ";
    std::cout<<"Pharuhi = "<<epHaruhi<<" GeV/c; ";
    std::cout<<"# seg. = "<<nseg<<"; ";
    std::cout<<"dz = "<< (trk->Z()-vz_hit_true)/1000. <<" mm;"<<std::endl;

};

double LinkedTracks::CalcThetaRMS(std::vector<float> TX, std::vector<float> TY){

    double theta_RMS = 0.00;
    double theta_mean = 0.00;

    for(size_t segIt=0; segIt<TX.size(); segIt++){
        theta_mean += std::sqrt( std::pow(TX.at(segIt), 2) + std::pow(TY.at(segIt), 2) );
    }

    theta_mean /= TX.size();

    for(size_t segIt=0; segIt<TX.size(); segIt++){
        theta_RMS += std::pow( std::sqrt( std::pow(TX.at(segIt), 2) + std::pow(TY.at(segIt), 2) ) - theta_mean, 2);
    }

    theta_RMS = std::sqrt(theta_RMS/TX.size());

    return theta_RMS;

};

double LinkedTracks::CalcMaxKinkAngle(std::vector<float> X, std::vector<float> Y, std::vector<float> Z){

    double MaxKinkAngle = 0.00;

    for(size_t segIt=0; segIt<X.size()-2; segIt++){

        TVector3 v1(X.at(segIt+1)-X.at(segIt), Y.at(segIt+1)-Y.at(segIt), Z.at(segIt+1)-Z.at(segIt));
        TVector3 v2(X.at(segIt+2)-X.at(segIt+1), Y.at(segIt+2)-Y.at(segIt+1), Z.at(segIt+2)-Z.at(segIt+1));

        double cosTheta = v1.Dot(v2) / (v1.Mag() * v2.Mag());
        double kinkAngle = std::acos(cosTheta);

        if(std::abs(kinkAngle) > std::abs(MaxKinkAngle)) MaxKinkAngle = kinkAngle;

    }

    return MaxKinkAngle; // in radian

};

bool LinkedTracks::IsValueChanged(std::vector<int> vec){

    bool IsChanged = false;

    for(size_t i=0; i<vec.size()-1; i++){
        if(vec.at(i) != vec.at(i+1)){
            IsChanged = true;
            break;
        }
    }

    return IsChanged;

};

float LinkedTracks::GetPmagReco(float ePharuhi, float ePcoord, float ePang, float ptrue, int PDG, int nseg){
    
    float pmag_reco = 0.00;

    if(std::abs(PDG) == kElectron){
        // gamma EM shower algorithm ?
        pmag_reco = ptrue;
    }
    else{
        // hadron
        if(ePharuhi > 3500){
            // failure 1
            if(ePcoord > 0 && ePcoord < 1000){
                pmag_reco = ePcoord;
            }
            else if(ePang > 0 && ePang < 1000){
                pmag_reco = ePang;
            }
            else{
                pmag_reco = 3;
            }
        }
        else{
           
            pmag_reco = ePharuhi;
            
        }
    }

    return pmag_reco;

};

void LinkedTracks::GetRecoMCInfo(bool ToPrintTrkInfo){

    Int_t nseg = 0;
    EdbSegP *trk = nullptr;
    TClonesArray *seg  = new TClonesArray("EdbSegP", 1000);
    tree->SetBranchAddress("nseg", &nseg);
    tree->SetBranchAddress("t.", &trk);
    tree->SetBranchAddress("s",  &seg);

    //PrintTrueVertex();

    int NTracks = tree->GetEntries();
    tree->GetEntry(0);
    eventID = trk->MCEvt() - 100000;
    std::cout<<"# volume tracks in event "<<eventID<<" = "<<NTracks<<std::endl;

    // track loop
    for(int trkIt=0; trkIt<NTracks; trkIt++){

        tree->GetEntry(trkIt);
        
        if(nseg < 3) continue; // # segments cut

        EdbTrackP *trkP = new EdbTrackP();
        std::vector<float> SegPara[5];
        std::vector<int> SegInt[3]; // 0: segment Plate ID, 1: segment PDG, 2: segment TrackID

        // segment loop
        for(int segIt=0; segIt<nseg; segIt++){

            EdbSegP *s = (EdbSegP *)seg->At(segIt);
            trkP->AddSegment(s);

            SegPara[0].push_back(s->X());
            SegPara[1].push_back(s->Y());
            SegPara[2].push_back(s->Z());
            SegPara[3].push_back(s->TX());
            SegPara[4].push_back(s->TY());

            SegInt[0].push_back(s->PID());
            SegInt[1].push_back(s->MCTrack());
            SegInt[2].push_back(s->Track());

        }

        trkP->SetNpl();

        EdbMomentumEstimator *MomEst = new EdbMomentumEstimator();
        MomEst->SetParPMS_Mag();

        //std::cout<<__LINE__<<std::endl;
        float eP_coord = MomEst->PMScoordinate(*trkP);
        //MomEst->~EdbMomentumEstimator();
        delete MomEst;
        //std::cout<<"eP (coord. method) = "<<eP_coord<<" GeV/c"<<std::endl;

        //std::cout<<__LINE__<<std::endl;
        MomEst = new EdbMomentumEstimator();
        MomEst->SetParPMS_Mag();

        float eP_ang = MomEst->PMSang(*trkP);
        //float eP_ang = -99.;
        //MomEst->~EdbMomentumEstimator();
        //std::cout<<"eP (angular method) = "<<eP_ang<<" GeV/c"<<std::endl;
        delete MomEst;

        FnuMomCoord *MomEstHaruhi = new FnuMomCoord();
        float eP_haruhi = MomEstHaruhi->CalcMomentum(trkP, 0, 0);
        //std::cout<<"eP (Haruhi) = "<<eP_haruhi<<" GeV/c"<<std::endl;
        delete MomEstHaruhi;

        // track level
        trackID.push_back(trk->Track());
        itrk.push_back(trkIt);
        PDG.push_back(trk->MCTrack());
        pmag_coord.push_back(eP_coord);
        pmag_ang.push_back(eP_ang);
        pmag_haruhi.push_back(eP_haruhi);
        ptrue.push_back(trk->P());
        NSeg.push_back(nseg);

        // segment level
        X.push_back(SegPara[0]);
        Y.push_back(SegPara[1]);
        Z.push_back(SegPara[2]);
        TX.push_back(SegPara[3]);
        TY.push_back(SegPara[4]);
        PID.push_back(SegInt[0]);
        SegPDG.push_back(SegInt[1]);
        SegTrkID.push_back(SegInt[2]);

        // track level
        theta_RMS.push_back(CalcThetaRMS(SegPara[3], SegPara[4]));
        MaxGap.push_back(trkP->CheckMaxGap());
        MaxKinkAngle.push_back(CalcMaxKinkAngle(SegPara[0], SegPara[1], SegPara[2]));
        IsPartCatChanged.push_back(IsValueChanged(SegInt[1]));
        IsTrkIdChanged.push_back(IsValueChanged(SegInt[2]));

        // track 3D linear fitting
        const double *para_line3Dfit = line3Dfit(SegPara[0], SegPara[1], SegPara[2]);
        TVector3 TrkDirVect(para_line3Dfit[1], para_line3Dfit[3], 1.00);
        double ImpactParameter = GetImpactParameter(TrkDirVect, para_line3Dfit[0], para_line3Dfit[2], vx_hit_true, vy_hit_true, vz_hit_true);
        double TraLen = CalcTrackLength(SegPara[0], SegPara[1], SegPara[2]);
        
        // track level
        IP.push_back(ImpactParameter);
        dz.push_back(SegPara[2].at(0) - vz_hit_true);
        theta.push_back(TrkDirVect.Theta());
        phi.push_back(TrkDirVect.Phi());
        TrackLength.push_back(TraLen);
        pmag_reco.push_back(GetPmagReco(eP_haruhi, eP_coord, eP_ang, trk->P(), trk->MCTrack(), nseg));

        if(ToPrintTrkInfo) PrintTrkInfo(trk, trkIt, nseg, eP_coord, eP_ang, eP_haruhi);

        delete trkP;
        

    }// track loop

    delete trk;
    delete seg;

};


double LinkedTracks::CalcTrackLength(std::vector<float> x, std::vector<float> y, std::vector<float> z){

    double TraLen = 0.00;

    for(size_t segIt=0; segIt<x.size()-1; segIt++){

        TraLen += std::sqrt( std::pow(x.at(segIt+1)-x.at(segIt), 2) + std::pow(y.at(segIt+1)-y.at(segIt), 2) + std::pow(z.at(segIt+1)-z.at(segIt), 2) );
        
    }

    return TraLen/1e3; // track length in mm

};
/*
void LinkedTracks::GetRecoAngleIPdz(){

    size_t NRecoTrks = X.size();

    for(size_t trkIt=0; trkIt<NRecoTrks; trkIt++){
        
        // track 3D linear fitting
        const double *para_line3Dfit = line3Dfit(X[trkIt], Y[trkIt], Z[trkIt]);

        TVector3 TrkDirVect(para_line3Dfit[1], para_line3Dfit[3], 1.00);
        double ImpactParameter = GetImpactParameter(TrkDirVect, para_line3Dfit[0], para_line3Dfit[2], vx_hit_true, vy_hit_true, vz_hit_true);
        
        double TraLen = CalcTrackLength(X[trkIt], Y[trkIt], Z[trkIt]);
        
        IP.push_back(ImpactParameter);
        dz.push_back(Z.at(trkIt).at(0) - vz_hit_true);
        theta.push_back(TrkDirVect.Theta());
        phi.push_back(TrkDirVect.Phi());
        TrackLength.push_back(TraLen);
        
    }

};
*/

bool LinkedTracks::IsPrimTrack(int TrackIt){
    
    float dz = Z.at(TrackIt).at(0) - vz_hit_true;
    bool StartWithin3PlatesDownStream = dz < (3. * (t_plate + t_film)) && dz > 0;

    float dr = std::sqrt(std::pow(X[TrackIt].at(0)-vx_hit_true,2) + std::pow(Y[TrackIt].at(0)-vy_hit_true,2));
    bool TransConfine = dr < (dz * 0.5);

    bool tanThetaCut = std::tan(theta.at(TrackIt)) < 0.5;
    bool IPcut = IP.at(TrackIt) < 5; // um

    return StartWithin3PlatesDownStream && TransConfine && tanThetaCut && IPcut;

};


void LinkedTracks::SortBasedOnTrackID(bool PrintTrkInfo){

    int NRecoTrks = trackID.size();

    // track level  
    std::vector<int> sorted_trackID(NRecoTrks);             // MC truth track ID
    std::vector<int> sorted_PDG(NRecoTrks);                 // MC truth track pdg code
    std::vector<int> sorted_itrk(NRecoTrks);                // reconstructed track iterator
    std::vector<double> sorted_theta(NRecoTrks);            // reconstructed track polar angle
    std::vector<double> sorted_phi(NRecoTrks);              // reconstructed track azimuthal angle
    std::vector<float> sorted_pmagReco(NRecoTrks);          // reconstructed track momentum
    std::vector<float> sorted_pmagCoord(NRecoTrks);         // reconstructed track momentum by MCS the coordinate method, FEDRA2022
    std::vector<float> sorted_pmagAng(NRecoTrks);           // reconstructed track momentum by MCS the angular method, FEDRA2022
    std::vector<float> sorted_pmagHaruhi(NRecoTrks);        // reconstructed track momentum by MCS the coordinate method, Haruhi
    std::vector<float> sorted_ptrue(NRecoTrks);             // MC truth track momentum
    std::vector<int> sorted_NSeg(NRecoTrks);                // reconstructed track the number of segments
    std::vector<double> sorted_IP(NRecoTrks);               // reconstructed impact parameter
    std::vector<float> sorted_dz(NRecoTrks);                // reconstructed (z_FirstSeg - z_vertex) um
    std::vector<double> sorted_TrackLength(NRecoTrks);      // reconstructed track length mm
    std::vector<int> sorted_MaxGap(NRecoTrks);              // reconstructed track the largest gap between 2 consecutive segments
    std::vector<double> sorted_theta_RMS(NRecoTrks);        // reconstructed track the standard deviation of the segment polar angle
    std::vector<double> sorted_MaxKinkAngle(NRecoTrks);     // reconstructed track the largest kink angle between 2 consecutive cells by positions
    std::vector<bool> sorted_IsPartCatChanged(NRecoTrks);   // reconstructed track the segment particle category changed or not
    std::vector<bool> sorted_IsTrkIdChanged(NRecoTrks);     // reconstructed track the segment track ID changed or not
    

    // segment level
    std::vector<std::vector<float>> sorted_X(NRecoTrks);        // reconstructed segment positions
    std::vector<std::vector<float>> sorted_Y(NRecoTrks);        // reconstructed segment positions
    std::vector<std::vector<float>> sorted_Z(NRecoTrks);        // reconstructed segment positions 
    std::vector<std::vector<float>> sorted_TX(NRecoTrks);       // reconstructed segment tangent polar angles
    std::vector<std::vector<float>> sorted_TY(NRecoTrks);       // reconstructed segment tangent polar angles
    std::vector<std::vector<int>>   sorted_PID(NRecoTrks);      // reconstructed segment plate ID: from 0 to 100
    std::vector<std::vector<int>>   sorted_SegPDG(NRecoTrks);   // reconstructed segment pdg code
    std::vector<std::vector<int>>   sorted_SegTrkID(NRecoTrks); // reconstructed segment track ID

    // Create an index vector
    std::vector<size_t> indices(NRecoTrks);
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }

    // Sort indices based on trackID
    std::vector<int> TrkID = trackID;
    std::sort(indices.begin(), indices.end(), [&TrkID](size_t i1, size_t i2){
        return TrkID[i1] < TrkID[i2];
    });

    // Rearrange vectors based on sorted indices
    for (size_t i=0; i<indices.size(); ++i){
        // track level  
        sorted_trackID[i] = trackID[indices[i]];
        sorted_PDG[i] = PDG[indices[i]];
        sorted_itrk[i] = itrk[indices[i]];          
        sorted_theta[i] = theta[indices[i]];
        sorted_phi[i] = phi[indices[i]];
        sorted_pmagReco[i] = pmag_reco[indices[i]];        
        sorted_pmagCoord[i] = pmag_coord[indices[i]];        
        sorted_pmagAng[i] = pmag_ang[indices[i]];
        sorted_pmagHaruhi[i] = pmag_haruhi[indices[i]];
        sorted_ptrue[i] = ptrue[indices[i]];      
        sorted_NSeg[i] = NSeg[indices[i]];       
        sorted_IP[i] = IP[indices[i]];
        sorted_dz[i] = dz[indices[i]];    
        sorted_TrackLength[i] = TrackLength[indices[i]];
        sorted_MaxGap[i] = MaxGap[indices[i]];
        sorted_theta_RMS[i] = theta_RMS[indices[i]];
        sorted_MaxKinkAngle[i] = MaxKinkAngle[indices[i]];
        sorted_IsPartCatChanged[i] = IsPartCatChanged[indices[i]];
        sorted_IsTrkIdChanged[i] = IsTrkIdChanged[indices[i]];
                
        // segment level
        sorted_X[i] = X[indices[i]]; 
        sorted_Y[i] = Y[indices[i]]; 
        sorted_Z[i] = Z[indices[i]]; 
        sorted_TX[i] = TX[indices[i]];
        sorted_TY[i] = TY[indices[i]];
        sorted_PID[i] = PID[indices[i]];
        sorted_SegPDG[i] = SegPDG[indices[i]];
        sorted_SegTrkID[i] = SegTrkID[indices[i]];

        if(PrintTrkInfo){
            std::cout<<"sorted track ID = "<<sorted_trackID[i]<<"; "
            <<"itrk = "<<sorted_itrk[i]<<"; "
            <<"PDG = "<<sorted_PDG[i]<<"; "
            <<"# seg. = "<<sorted_NSeg[i]<<"; "
            <<"pmag_coord = "<<sorted_pmagCoord[i]<<" GeV/c; "
            <<"pmag_ang = "<<sorted_pmagAng[i]<<" GeV/c; "
            <<"pmag_haruhi = "<<sorted_pmagHaruhi[i]<<" GeV/c; "
            <<"ptrue = "<<sorted_ptrue[i]<<" GeV/c; "
            <<"IP = "<<sorted_IP[i]<<" um; "
            <<"dz = "<<sorted_dz[i]/1e3<<" mm; "
            <<"TraLen = "<<sorted_TrackLength[i]<<" mm; "
            <<"MaxGap = "<<sorted_MaxGap[i]<<"; "
            <<"theta_RMS = "<<sorted_theta_RMS[i]<<" rad; "
            <<"MaxKinkAngle = "<<sorted_MaxKinkAngle[i]<<" rad; "
            <<std::endl;
        }
        
    }

    trackID = sorted_trackID;
    PDG = sorted_PDG;
    itrk = sorted_itrk;
    theta = sorted_theta;
    phi = sorted_phi;
    pmag_reco = sorted_pmagReco;
    pmag_coord = sorted_pmagCoord;
    pmag_ang = sorted_pmagAng;
    pmag_haruhi = sorted_pmagHaruhi;
    ptrue = sorted_ptrue;
    NSeg = sorted_NSeg;
    IP = sorted_IP;
    dz = sorted_dz;
    TrackLength = sorted_TrackLength;
    MaxGap = sorted_MaxGap;
    theta_RMS = sorted_theta_RMS;
    MaxKinkAngle = sorted_MaxKinkAngle;
    IsPartCatChanged = sorted_IsPartCatChanged;
    IsTrkIdChanged = sorted_IsTrkIdChanged;

    X = sorted_X;
    Y = sorted_Y;
    Z = sorted_Z;
    TX = sorted_TX;
    TY = sorted_TY;
    PID = sorted_PID;
    SegPDG = sorted_SegPDG;
    SegTrkID = sorted_SegTrkID;

};

