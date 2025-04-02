#include <vector>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "Discriminators.h"


const double TrainTestValiRatio = 0.80;

// truth info
int EventID;
int mcID;                         // 200025, 200026, 200035, 100069
int Nu_PDG;                       // neutral mother PDG code
double Nu_px, Nu_py, Nu_pz, Nu_e; // neutral mother 3-momentum

std::string FileNum;
std::vector<int> TrackID;
std::vector<int> PDG;
std::vector<int> itrk;          // reconstructed track iterator

// kinematics: default reco
std::vector<double> px_reco;
std::vector<double> py_reco;
std::vector<double> pz_reco;
std::vector<double> pmag_reco;
std::vector<double> theta_reco;
std::vector<double> phi_reco;

std::vector<double> pmag_haruhi; // Haruhi's method
std::vector<double> pmag_ang;    // Fedra built-in
std::vector<double> pmag_coord;

std::vector<double> px_true;
std::vector<double> py_true;
std::vector<double> pz_true;
std::vector<double> pmag_true;
std::vector<double> theta_true;
std::vector<double> phi_true;

// other track-level variables
std::vector<int> nseg;
std::vector<double> TrackLength;
std::vector<float> dz;              // z-distance to the primary vertex in um
std::vector<double> IP;             // impact parameter in um
std::vector<int> PID_start;         // the starting plate ID
std::vector<int> PID_end;           // the ending plate ID
std::vector<int> MaxGap;            // the largest gap between 2 consecutive segments
std::vector<double> theta_RMS;      // the standard deviation of the segment polar angle
std::vector<double> MaxKinkAngle;   // the largest kink angle between 2 consecutive cells by positions
std::vector<bool> IsPartCatChanged; // the segment particle category changed or not
std::vector<bool> IsTrkIdChanged;   // the segment track ID changed or not


// discriminators
float n_ch;
float nTrkTanThetaLeq0point1_reco;
float nTrkTanThetaLeq0point1_true;

float pmag_had_vis_reco;
float dphi_max_reco;
float p3_hardest_reco;
float InvThetaCh_reco;
float pTmiss_mag_reco;
float pTabs_sum_reco;
float DeltaPhiMET_reco;
float dphi_sum_reco;
float tan_theta_hardest_reco;

float pmag_had_vis_true;
float dphi_max_true;
float p3_hardest_true;
float InvThetaCh_true;
float pTmiss_mag_true;
float pTabs_sum_true;
float DeltaPhiMET_true;
float dphi_sum_true;
float tan_theta_hardest_true;

void Branch(TTree *tree){
    
    tree->Branch("mcID", &mcID);
    tree->Branch("FileNum", &FileNum);
    tree->Branch("Nu_PDG", &Nu_PDG);
    tree->Branch("Nu_px", &Nu_px);
    tree->Branch("Nu_py", &Nu_py);  
    tree->Branch("Nu_pz", &Nu_pz);
    tree->Branch("Nu_e", &Nu_e);
    tree->Branch("EventID", &EventID);
    tree->Branch("TrackID", &TrackID);
    tree->Branch("PDG", &PDG);
    tree->Branch("itrk", &itrk);

    tree->Branch("px_reco", &px_reco);
    tree->Branch("py_reco", &py_reco);
    tree->Branch("pz_reco", &pz_reco);
    tree->Branch("pmag_reco", &pmag_reco);
    tree->Branch("theta_reco", &theta_reco);
    tree->Branch("phi_reco", &phi_reco);
    

    tree->Branch("px_true", &px_true);
    tree->Branch("py_true", &py_true);
    tree->Branch("pz_true", &pz_true);
    tree->Branch("pmag_true", &pmag_true);
    tree->Branch("theta_true", &theta_true);
    tree->Branch("phi_true", &phi_true);

    tree->Branch("pmag_haruhi", &pmag_haruhi);
    tree->Branch("pmag_ang", &pmag_ang);
    tree->Branch("pmag_coord", &pmag_coord);
    tree->Branch("nseg", &nseg);
    tree->Branch("TrackLength", &TrackLength);
    tree->Branch("dz", &dz);
    tree->Branch("IP", &IP);
    tree->Branch("PID_start", &PID_start);
    tree->Branch("PID_end", &PID_end);
    tree->Branch("MaxGap", &MaxGap);
    tree->Branch("theta_RMS", &theta_RMS);
    tree->Branch("MaxKinkAngle", &MaxKinkAngle);
    tree->Branch("IsPartCatChanged", &IsPartCatChanged);
    tree->Branch("IsTrkIdChanged", &IsTrkIdChanged);

    tree->Branch("n_ch", &n_ch);
    tree->Branch("nTrkTanThetaLeq0point1_reco", &nTrkTanThetaLeq0point1_reco);
    tree->Branch("nTrkTanThetaLeq0point1_true", &nTrkTanThetaLeq0point1_true);

    tree->Branch("pmag_had_vis_reco", &pmag_had_vis_reco);
    tree->Branch("dphi_max_reco", &dphi_max_reco);
    tree->Branch("p3_hardest_reco", &p3_hardest_reco);
    tree->Branch("InvThetaCh_reco", &InvThetaCh_reco);
    tree->Branch("pTmiss_mag_reco", &pTmiss_mag_reco);
    tree->Branch("pTabs_sum_reco", &pTabs_sum_reco);
    tree->Branch("DeltaPhiMET_reco", &DeltaPhiMET_reco);
    tree->Branch("dphi_sum_reco", &dphi_sum_reco);
    tree->Branch("tan_theta_hardest_reco", &tan_theta_hardest_reco);

    tree->Branch("pmag_had_vis_true", &pmag_had_vis_true);
    tree->Branch("dphi_max_true", &dphi_max_true);
    tree->Branch("p3_hardest_true", &p3_hardest_true);
    tree->Branch("InvThetaCh_true", &InvThetaCh_true);
    tree->Branch("pTmiss_mag_true", &pTmiss_mag_true);
    tree->Branch("pTabs_sum_true", &pTabs_sum_true);
    tree->Branch("DeltaPhiMET_true", &DeltaPhiMET_true);
    tree->Branch("dphi_sum_true", &dphi_sum_true);
    tree->Branch("tan_theta_hardest_true", &tan_theta_hardest_true);

}

void EventInit(){

    EventID = -999;
    mcID = -999;
    Nu_PDG = -999;
    Nu_px = -999.;
    Nu_py = -999.;
    Nu_pz = -999.;
    Nu_e = -999.;
    FileNum = "";

    n_ch = -999;

    pmag_had_vis_reco = -999.;
    dphi_max_reco = -999.;
    p3_hardest_reco = -999.;
    InvThetaCh_reco = -999.;
    pTmiss_mag_reco = -999.;
    pTabs_sum_reco = -999.;
    DeltaPhiMET_reco = -999.;
    dphi_sum_reco = -999.;
    tan_theta_hardest_reco = -999.;

    pmag_had_vis_true = -999.;
    dphi_max_true = -999.;
    p3_hardest_true = -999.;
    InvThetaCh_true = -999.;
    pTmiss_mag_true = -999.;
    pTabs_sum_true = -999.;
    DeltaPhiMET_true = -999.;
    dphi_sum_true = -999.;
    tan_theta_hardest_true = -999.;

    TrackID.clear();
    PDG.clear();
    itrk.clear();
    IsPartCatChanged.clear();
    IsTrkIdChanged.clear();
    
    px_reco.clear();
    py_reco.clear();
    pz_reco.clear();
    pmag_reco.clear();
    theta_reco.clear();
    phi_reco.clear();

    pmag_ang.clear();
    pmag_coord.clear();
    pmag_haruhi.clear();

    px_true.clear();
    py_true.clear();
    pz_true.clear();
    pmag_true.clear();
    theta_true.clear();
    phi_true.clear();

    nseg.clear();
    TrackLength.clear();
    IP.clear();
    dz.clear();
    PID_start.clear();
    PID_end.clear();
    MaxGap.clear();
    theta_RMS.clear();
    MaxKinkAngle.clear();

}

void PassDiscAddressToTTreeAddress(Discriminators *disc){

    mcID = disc->mcID;
    EventID = disc->EventID;

    TrackID = disc->TrackID;
    PDG = disc->PDG;
    itrk = disc->itrk;
    px_reco = disc->px_reco;
    py_reco = disc->py_reco;
    pz_reco = disc->pz_reco;
    pmag_reco = disc->pmag_reco;
    theta_reco = disc->theta_reco;
    phi_reco = disc->phi_reco;
    px_true = disc->px_true;
    py_true = disc->py_true;
    pz_true = disc->pz_true;
    pmag_true = disc->pmag_true;
    theta_true = disc->theta_true;
    phi_true = disc->phi_true;

    pmag_ang = disc->pmag_ang;
    pmag_coord = disc->pmag_coord;
    pmag_haruhi = disc->pmag_haruhi;

    nseg = disc->nseg;
    TrackLength = disc->TrackLength;
    IP = disc->IP;
    dz = disc->dz;
    PID_start = disc->PID_start;
    PID_end = disc->PID_end;
    MaxGap = disc->MaxGap;
    theta_RMS = disc->theta_RMS;
    MaxKinkAngle = disc->MaxKinkAngle;
    IsPartCatChanged = disc->IsPartCatChanged;
    IsTrkIdChanged = disc->IsTrkIdChanged;

    n_ch = disc->n_ch;
    nTrkTanThetaLeq0point1_reco = disc->nTrkTanThetaLeq0point1_reco;
    nTrkTanThetaLeq0point1_true = disc->nTrkTanThetaLeq0point1_true;

    pmag_had_vis_reco = disc->pmag_had_vis_reco;
    dphi_max_reco = disc->dphi_max_reco;
    p3_hardest_reco = disc->p3_hardest_reco;
    InvThetaCh_reco = disc->InvThetaCh_reco;
    pTmiss_mag_reco = disc->pTmiss_mag_reco;
    pTabs_sum_reco = disc->pTabs_sum_reco;
    DeltaPhiMET_reco = disc->DeltaPhiMET_reco;
    dphi_sum_reco = disc->dphi_sum_reco;
    tan_theta_hardest_reco = disc->tan_theta_hardest_reco;

    pmag_had_vis_true = disc->pmag_had_vis_true;
    dphi_max_true = disc->dphi_max_true;
    p3_hardest_true = disc->p3_hardest_true;
    InvThetaCh_true = disc->InvThetaCh_true;
    pTmiss_mag_true = disc->pTmiss_mag_true;
    pTabs_sum_true = disc->pTabs_sum_true;
    DeltaPhiMET_true = disc->DeltaPhiMET_true;
    dphi_sum_true = disc->dphi_sum_true;
    tan_theta_hardest_true = disc->tan_theta_hardest_true;

}


