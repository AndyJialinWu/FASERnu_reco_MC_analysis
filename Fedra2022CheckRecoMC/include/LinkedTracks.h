#include <vector>
#include <cmath>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "EdbSegP.h"
#include "EdbMomentumEstimator.h"
#include "EdbPattern.h"
#include "EdbSegCouple.h"
#include "EdbPVRec.h"

#include "line3Dfit.h"

const double IntLumi_run3 = 250.00;           // fb^
const double wgt_NC_200025 = 50173.00 / 10e3; // event/fb^{-1}
const double wgt_NC_200026 = 66355.00 / 10e3;
const double wgt_NC_200035 = 13555.00 / 10e3;

const float t_plate = 1.1e3;                  // tungsten plate thickness in um
const float t_film = 0.34e3;                  // film thickness in um
const int kMaxTrks = 79;

// event based
class LinkedTracks{
    public:
        
        TTree *tree;                //!pointer to the analyzed TTree or TChain

        //  event level
        //int n_ch;                       // primary charged tracks multiplicity
        float vx_hit_reco;              // reconstructed vertex position
        float vy_hit_reco;              // reconstructed vertex position
        float vz_hit_reco;              // reconstructed vertex position 
        float vx_hit_true;              // MC truth vertex position
        float vy_hit_true;              // MC truth vertex position
        float vz_hit_true;              // MC truth vertex position
        int eventID;                    // MC event ID

        // track level  
        std::vector<int> trackID;       // MC truth track ID
        std::vector<int> PDG;           // MC truth track pdg code
        std::vector<double> theta;      // reconstructed track polar angle
        std::vector<double> phi;        // reconstructed track azimuthal angle
        std::vector<float> pmag;        // reconstructed track momentum by MCS
        std::vector<float> ptrue;       // MC truth track momentum
        std::vector<int> NSeg;          // reconstructed track the number of segments
        std::vector<double> IP;         // reconstructed impact parameter
                
        // segment level
        std::vector<std::vector<float>> X; // reconstructed segment positions
        std::vector<std::vector<float>> Y; // reconstructed segment positions
        std::vector<std::vector<float>> Z; // reconstructed segment positions 
        std::vector<std::vector<float>> TX;// reconstructed segment tangent polar angles
        std::vector<std::vector<float>> TY;// reconstructed segment tangent polar angles




    public:
        LinkedTracks(TTree *tracks=0);
        virtual void Init(TTree *tracks);
        
        void SetTrueVertex(float vx_global_MC, float vy_global_MC, float vz_global_MC); // MC truth (global coord.) to hit coord.
        void PrintTrueVertex();

        // set eventID, trackID, PDG, pmag, ptrue, NSeg, X, Y, Z, TX, TY
        void SetRecoMCInfo(int FedraEvtID, std::vector<int> FedraTrkID, std::vector<int> FedraPDG, std::vector<float> FedraPmag, std::vector<float> FedraPtrue, std::vector<int> FedraNSeg,
                           std::vector<std::vector<float>> FedraX, std::vector<std::vector<float>> FedraY, std::vector<std::vector<float>> FedraZ,
                           std::vector<std::vector<float>> FedraTX, std::vector<std::vector<float>> FedraTY); 

        void GetRecoMCInfo();    
        void GetRecoAngleIP();  // get theta, phi, IP

        bool IsPrimTrack(int TrackIt);

        void SortBasedOnTrackID();// sorting ascendingly based on the track ID

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

    //n_ch = -9999;
    eventID = -9999;
    vx_hit_true = -9999.;
    vy_hit_true = -9999.;
    vz_hit_true = -9999.;
    vx_hit_reco = -9999.;
    vy_hit_reco = -9999.;
    vz_hit_reco = -9999.;

    trackID.clear();
    PDG.clear();
    theta.clear();
    phi.clear();
    pmag.clear();
    NSeg.clear();
    ptrue.clear();
    IP.clear();
    
    X.clear();
    Y.clear();
    Z.clear();
    TX.clear();
    TY.clear();

};



LinkedTracks::~LinkedTracks(){
    if (!tree) return;
    delete tree->GetCurrentFile();
};


void LinkedTracks::PrintTrueVertex(){
    std::cout<<"(vx, vy, vz) = "<<"("<<vx_hit_true<<", "<<vy_hit_true<<", "<<vz_hit_true<<")"<<" um"<<std::endl;
};

void LinkedTracks::SetTrueVertex(float vx_global_MC, float vy_global_MC, float vz_global_MC){
    
    vx_hit_true = (vx_global_MC - 10.) * 1000.; // micro-metre
    vy_hit_true = (vy_global_MC + 21.) * 1000.;
    vz_hit_true = (vz_global_MC + 2986.27) * 1000.;

    std::cout<<"(vx, vy, vz) = "<<"("<<vx_hit_true<<", "<<vy_hit_true<<", "<<vz_hit_true<<")"<<" um"<<std::endl;

    PrintTrueVertex();

};



void LinkedTracks::SetRecoMCInfo(int FedraEvtID, std::vector<int> FedraTrkID, std::vector<int> FedraPDG, std::vector<float> FedraPmag, std::vector<float> FedraPtrue, std::vector<int> FedraNSeg,
                                 std::vector<std::vector<float>> FedraX, std::vector<std::vector<float>> FedraY, std::vector<std::vector<float>> FedraZ,
                                 std::vector<std::vector<float>> FedraTX, std::vector<std::vector<float>> FedraTY){

    eventID = FedraEvtID;

    // track level
    trackID = FedraTrkID;
    PDG = FedraPDG;
    pmag= FedraPmag;
    ptrue= FedraPtrue;
    NSeg = FedraNSeg;

    // segment level
    X = FedraX;
    Y = FedraY;
    Z = FedraZ;
    TX = FedraTX;
    TY = FedraTY;

};

void LinkedTracks::GetRecoMCInfo(){

    Int_t nseg = 0;
    EdbSegP *trk = nullptr;
    TClonesArray *seg  = new TClonesArray("EdbSegP", 79);
    tree->SetBranchAddress("nseg", &nseg);
    tree->SetBranchAddress("t.", &trk);
    tree->SetBranchAddress("s",  &seg);

    PrintTrueVertex();

    int NTracks = tree->GetEntries();
    tree->GetEntry(0);
    eventID = trk->MCEvt() - 100000;
    std::cout<<"# volume tracks in event "<<eventID<<" = "<<NTracks<<std::endl;

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
        std::vector<float> SegPara[5];

        // segment loop
        for(int segIt=0; segIt<nseg; segIt++){

            EdbSegP *s = (EdbSegP *)seg->At(segIt);
            trkP->AddSegment(s);

            SegPara[0].push_back(s->X());
            SegPara[1].push_back(s->Y());
            SegPara[2].push_back(s->Z());
            SegPara[3].push_back(s->TX());
            SegPara[4].push_back(s->TY());

        }

        trkP->SetNpl();
        EdbMomentumEstimator *MomEst = new EdbMomentumEstimator();
        MomEst->SetParPMS_Mag();

        float eP = MomEst->PMScoordinate(*trkP);
        std::cout<<"eP (coord. method) = "<<eP<<" GeV/c"<<std::endl;

        trackID.push_back(trk->Track());
        PDG.push_back(trk->MCTrack());
        pmag.push_back(eP);
        ptrue.push_back(trk->P());
        NSeg.push_back(nseg);
        X.push_back(SegPara[0]);
        Y.push_back(SegPara[1]);
        Z.push_back(SegPara[2]);
        TX.push_back(SegPara[3]);
        TY.push_back(SegPara[4]);

    }

}

void LinkedTracks::GetRecoAngleIP(){

    int NRecoTrks = X.size();

    for(int trkIt=0; trkIt<NRecoTrks; trkIt++){
        
        // track 3D linear fitting
        const double *para_line3Dfit = line3Dfit(X[trkIt], Y[trkIt], Z[trkIt]);

        TVector3 TrkDirVect(para_line3Dfit[1], para_line3Dfit[3], 1.00);
        double ImpactParameter = GetImpactParameter(TrkDirVect, para_line3Dfit[0], para_line3Dfit[2], vx_hit_true, vy_hit_true, vz_hit_true);
        
        IP.push_back(ImpactParameter);
        theta.push_back(TrkDirVect.Theta());
        phi.push_back(TrkDirVect.Phi());
        
    }

};


bool LinkedTracks::IsPrimTrack(int TrackIt){
    
    float dz = Z[TrackIt].at(0) - vz_hit_true;
    bool StartWithin3PlatesDownStream = dz < (3. * (t_plate + t_film)) && dz > 0;

    float dr = std::sqrt(std::pow(X[TrackIt].at(0)-vx_hit_true,2) + std::pow(Y[TrackIt].at(0)-vy_hit_true,2));
    bool TransConfine = dr < (dz * 0.5);

    bool tanThetaCut = theta.at(TrackIt) < 0.5;
    bool IPcut = IP.at(TrackIt) < 10; // IP < 10 um

    return StartWithin3PlatesDownStream && TransConfine && tanThetaCut && IPcut;

};


void LinkedTracks::SortBasedOnTrackID(){

    int NRecoTrks = trackID.size();

    // track level  
    std::vector<int> sorted_trackID(NRecoTrks);       // MC truth track ID
    std::vector<int> sorted_PDG(NRecoTrks);           // MC truth track pdg code
    std::vector<double> sorted_theta(NRecoTrks);      // reconstructed track polar angle
    std::vector<double> sorted_phi(NRecoTrks);        // reconstructed track azimuthal angle
    std::vector<float> sorted_pmag(NRecoTrks);        // reconstructed track momentum by MCS
    std::vector<float> sorted_ptrue(NRecoTrks);       // MC truth track momentum
    std::vector<int> sorted_NSeg(NRecoTrks);          // reconstructed track the number of segments
    std::vector<double> sorted_IP(NRecoTrks);         // reconstructed impact parameter
            
    // segment level
    std::vector<std::vector<float>> sorted_X(NRecoTrks); // reconstructed segment positions
    std::vector<std::vector<float>> sorted_Y(NRecoTrks); // reconstructed segment positions
    std::vector<std::vector<float>> sorted_Z(NRecoTrks); // reconstructed segment positions 
    std::vector<std::vector<float>> sorted_TX(NRecoTrks);// reconstructed segment tangent polar angles
    std::vector<std::vector<float>> sorted_TY(NRecoTrks);// reconstructed segment tangent polar angles

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
        sorted_theta[i] = theta[indices[i]];
        sorted_phi[i] = phi[indices[i]];        
        sorted_pmag[i] = pmag[indices[i]];        
        sorted_ptrue[i] = ptrue[indices[i]];      
        sorted_NSeg[i] = NSeg[indices[i]];       
        sorted_IP[i] = IP[indices[i]];    
                
        // segment level
        sorted_X[i] = X[indices[i]]; 
        sorted_Y[i] = Y[indices[i]]; 
        sorted_Z[i] = Z[indices[i]]; 
        sorted_TX[i] = TX[indices[i]];
        sorted_TY[i] = TY[indices[i]];

        std::cout<<"sorted track ID = "<<sorted_trackID[i]<<"; "
        <<"PDG = "<<sorted_PDG[i]<<"; "
        <<"# seg. = "<<sorted_NSeg[i]<<"; "
        <<"pmag = "<<sorted_pmag[i]<<" GeV/c; "
        <<"ptrue = "<<sorted_ptrue[i]<<" GeV/c; "
        <<"IP = "<<sorted_IP[i]<<" um;"
        <<std::endl;
    }

    trackID = sorted_trackID;
    PDG = sorted_PDG;
    theta = sorted_theta;
    phi = sorted_phi;
    pmag = sorted_pmag;
    ptrue = sorted_ptrue;
    NSeg = sorted_NSeg;
    IP = sorted_IP;

    X = sorted_X;
    Y = sorted_Y;
    Z = sorted_Z;
    TX = sorted_TX;
    TY = sorted_TY;

};

