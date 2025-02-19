#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "Discriminators.h"

// truth info
int EventID;
int mcID; // 200025, 200026, 200035
std::vector<int> TrackID;
std::vector<int> PDG;

// kinematics
std::vector<double> px_reco;
std::vector<double> py_reco;
std::vector<double> pz_reco;
std::vector<double> pmag_reco;
std::vector<double> theta_reco;
std::vector<double> phi_reco;

std::vector<double> px_true;
std::vector<double> py_true;
std::vector<double> pz_true;
std::vector<double> pmag_true;
std::vector<double> theta_true;
std::vector<double> phi_true;

// discriminators
int n_ch;

double pmag_had_vis_reco;
double dphi_max_reco;
double p3_hardest_reco;
double InvThetaCh_reco;
double pTmiss_mag_reco;
double pTabs_sum_reco;
double DeltaPhiMET_reco;
double dphi_sum_reco;

double pmag_had_vis_true;
double dphi_max_true;
double p3_hardest_true;
double InvThetaCh_true;
double pTmiss_mag_true;
double pTabs_sum_true;
double DeltaPhiMET_true;
double dphi_sum_true;

void Branch(TTree *tree){
    
    tree->Branch("mcID", &mcID);
    tree->Branch("EventID", &EventID);
    tree->Branch("TrackID", &TrackID);
    tree->Branch("PDG", &PDG);

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


    tree->Branch("n_ch", &n_ch);

    tree->Branch("pmag_had_vis_reco", &pmag_had_vis_reco);
    tree->Branch("dphi_max_reco", &dphi_max_reco);
    tree->Branch("p3_hardest_reco", &p3_hardest_reco);
    tree->Branch("InvThetaCh_reco", &InvThetaCh_reco);
    tree->Branch("pTmiss_mag_reco", &pTmiss_mag_reco);
    tree->Branch("pTabs_sum_reco", &pTabs_sum_reco);
    tree->Branch("DeltaPhiMET_reco", &DeltaPhiMET_reco);
    tree->Branch("dphi_sum_reco", &dphi_sum_reco);

    tree->Branch("pmag_had_vis_true", &pmag_had_vis_true);
    tree->Branch("dphi_max_true", &dphi_max_true);
    tree->Branch("p3_hardest_true", &p3_hardest_true);
    tree->Branch("InvThetaCh_true", &InvThetaCh_true);
    tree->Branch("pTmiss_mag_true", &pTmiss_mag_true);
    tree->Branch("pTabs_sum_true", &pTabs_sum_true);
    tree->Branch("DeltaPhiMET_true", &DeltaPhiMET_true);
    tree->Branch("dphi_sum_true", &dphi_sum_true);

}

void EventInit(){

    EventID = -999;
    mcID = -999;

    n_ch = -999;

    pmag_had_vis_reco = -999.;
    dphi_max_reco = -999.;
    p3_hardest_reco = -999.;
    InvThetaCh_reco = -999.;
    pTmiss_mag_reco = -999.;
    pTabs_sum_reco = -999.;
    DeltaPhiMET_reco = -999.;
    dphi_sum_reco = -999.;

    pmag_had_vis_true = -999.;
    dphi_max_true = -999.;
    p3_hardest_true = -999.;
    InvThetaCh_true = -999.;
    pTmiss_mag_true = -999.;
    pTabs_sum_true = -999.;
    DeltaPhiMET_true = -999.;
    dphi_sum_true = -999.;

    TrackID.clear();
    PDG.clear();
    
    px_reco.clear();
    py_reco.clear();
    pz_reco.clear();
    pmag_reco.clear();
    theta_reco.clear();
    phi_reco.clear();

    px_true.clear();
    py_true.clear();
    pz_true.clear();
    pmag_true.clear();
    theta_true.clear();
    phi_true.clear();

}

void PassDiscAddressToTTreeAddress(Discriminators *disc){

    mcID = disc->mcID;
    EventID = disc->EventID;
    TrackID = disc->TrackID;
    PDG = disc->PDG;

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

    n_ch = disc->n_ch;

    pmag_had_vis_reco = disc->pmag_had_vis_reco;
    dphi_max_reco = disc->dphi_max_reco;
    p3_hardest_reco = disc->p3_hardest_reco;
    InvThetaCh_reco = disc->InvThetaCh_reco;
    pTmiss_mag_reco = disc->pTmiss_mag_reco;
    pTabs_sum_reco = disc->pTabs_sum_reco;
    DeltaPhiMET_reco = disc->DeltaPhiMET_reco;
    dphi_sum_reco = disc->dphi_sum_reco;

    pmag_had_vis_true = disc->pmag_had_vis_true;
    dphi_max_true = disc->dphi_max_true;
    p3_hardest_true = disc->p3_hardest_true;
    InvThetaCh_true = disc->InvThetaCh_true;
    pTmiss_mag_true = disc->pTmiss_mag_true;
    pTabs_sum_true = disc->pTabs_sum_true;
    DeltaPhiMET_true = disc->DeltaPhiMET_true;
    dphi_sum_true = disc->dphi_sum_true;

}


