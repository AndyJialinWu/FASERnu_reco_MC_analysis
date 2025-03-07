#include <vector>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TPDGCode.h"
#include "Math/Vector2D.h"
#include "Math/Vector3D.h"

#ifndef Discriminators_h
#define Discriminators_h

class Discriminators{
public:

    // truth info
    int EventID;
    int mcID; // 200025, 200026, 200035
    std::vector<int> TrackID;
    std::vector<int> PDG;
    std::vector<int> itrk;              // reconstructed MC track iterator

    // kinematics: default reco
    std::vector<double> px_reco;        // GeV
    std::vector<double> py_reco;
    std::vector<double> pz_reco;
    std::vector<double> pmag_reco;
    std::vector<double> theta_reco;
    std::vector<double> phi_reco;

    std::vector<double> pmag_ang;       // Fedra built-in
    std::vector<double> pmag_coord;
    std::vector<double> pmag_haruhi;    // Haruhi's method

    std::vector<double> px_true;        // GeV
    std::vector<double> py_true;
    std::vector<double> pz_true;
    std::vector<double> pmag_true;
    std::vector<double> theta_true;
    std::vector<double> phi_true;

    // other track-level variables
    std::vector<int> nseg;              // the number of segments of a track
    std::vector<double> TrackLength;    // track length in mm
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
    int n_ch;

    double pmag_had_vis_reco;
    double dphi_max_reco;
    double p3_hardest_reco;
    double InvThetaCh_reco;
    double pTmiss_mag_reco;
    double pTabs_sum_reco;
    double DeltaPhiMET_reco;
    double dphi_sum_reco;
    double tan_theta_hardest_reco;

    double pmag_had_vis_true;
    double dphi_max_true;
    double p3_hardest_true;
    double InvThetaCh_true;
    double pTmiss_mag_true;
    double pTabs_sum_true;
    double DeltaPhiMET_true;
    double dphi_sum_true;
    double tan_theta_hardest_true;

    // member functions
    Discriminators();
    ~Discriminators();

    static int ChTrk_Multi(std::vector<double> pmag);
    static double VisibleHadronicMomentum(std::vector<double> pmag);
    static double Delta_phi_max(std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> pmag, double pmag_had_vis);
    static double HardestTrackMomentum(std::vector<double> pmag);
    static double InverseSumChargedTrackAngle(std::vector<double> theta);
    static double pT_miss_mag(std::vector<double> px, std::vector<double> py);
    static double pT_abs_sum(std::vector<double> px, std::vector<double> py);
    static double TrackMETangle(std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> pmag, double pmag_had_vis);
    static double Delta_phi_sum(std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> pmag, double pmag_had_vis);
    static double TanThetaHardest(std::vector<double> pmag, std::vector<double> theta);

    void CalcDisc(); // calculate discriminators


};
#endif

Discriminators::Discriminators(){

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
    MaxGap.clear();
    theta_RMS.clear();
    MaxKinkAngle.clear();
    IsPartCatChanged.clear();
    IsTrkIdChanged.clear();
    
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

    pmag_ang.clear();
    pmag_coord.clear();
    pmag_haruhi.clear();

    nseg.clear();
    TrackLength.clear();
    IP.clear();
    dz.clear();
    PID_start.clear();
    PID_end.clear();

};

Discriminators::~Discriminators(){
    ;
};

int Discriminators::ChTrk_Multi(std::vector<double> pmag){

    int n_ch = pmag.size();
    return n_ch;

};

double Discriminators::VisibleHadronicMomentum(std::vector<double> pmag){

    double pmag_had_vis = std::accumulate(pmag.begin(), pmag.end(), 0.00);
    return pmag_had_vis;

};

double Discriminators::Delta_phi_max(std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> pmag, double pmag_had_vis){

    std::vector<double> phi;

    // hadron loop
    for(size_t HadIt=0; HadIt<pmag.size(); HadIt++){
        
        if(pmag[HadIt] <= (0.10*pmag_had_vis)){
            // count the hadron carrying at least 10% of the hadronic system energy
            continue;
        }
        
        ROOT::Math::XYVector pT_hadron(px[HadIt], py[HadIt]);
        phi.push_back(pT_hadron.Phi());
    }

    if(phi.size() <= 1){
        return 0;
    }
    else{
        // sort the phi ascendingly
        std::sort(phi.begin(), phi.end());
        double dphi_max = 0.00;

        // phi loop for the largest azimuthal gap
        for(size_t phiIt=0; phiIt<phi.size()-1; phiIt++){
            // azimuthal angle difference between the azimuthally neighboring track
            double dphi = phi[phiIt+1] - phi[phiIt];
            if(dphi < M_PI && dphi > dphi_max){
                dphi_max = dphi;
            }
        }

        double dphi = phi.back() - phi.front();
        if(dphi > M_PI){
            dphi = 2*M_PI - dphi;
            if(dphi > dphi_max){
                dphi_max = dphi;
            }
        }

        double dphi_max_deg = dphi_max * (180.00/M_PI);
        return dphi_max_deg;
    }

    return 0;

};

double Discriminators::HardestTrackMomentum(std::vector<double> pmag){

    if(pmag.size() == 0) return 0;

    double p3_hardest = *std::max_element(pmag.begin(), pmag.end());
    return p3_hardest;

};

double Discriminators::InverseSumChargedTrackAngle(std::vector<double> theta){

    double InvThetaCh = 0.00;

    // hadron loop
    for(size_t HadIt=0; HadIt<theta.size(); HadIt++){
        
        if(theta[HadIt] != 0){
            InvThetaCh += 1.00/theta[HadIt];
        }
        else{
            InvThetaCh += 800;
        }

    }

    return InvThetaCh;

};

double Discriminators::pT_miss_mag(std::vector<double> px, std::vector<double> py){

    double px_tot = std::accumulate(px.begin(), px.end(), 0.00);
    double py_tot = std::accumulate(py.begin(), py.end(), 0.00);

    double pT_mag = std::sqrt(std::pow(px_tot,2) + std::pow(py_tot,2));
    return pT_mag;

};


double Discriminators::pT_abs_sum(std::vector<double> px, std::vector<double> py){

    double pT_tot = 0.00;

    // hadron loop
    for(size_t HadIt=0; HadIt<px.size(); HadIt++){
        pT_tot += std::sqrt(std::pow(px[HadIt],2) + std::pow(py[HadIt],2));
    }

    return pT_tot;

};

double Discriminators::TrackMETangle(std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> pmag, double pmag_had_vis){

    std::vector<double> phi;
    double px_tot = 0.00;
    double py_tot = 0.00;

    for(size_t HadIt=0; HadIt<pmag.size(); HadIt++){
        
        // charged hadrons
        px_tot += px[HadIt];
        py_tot += py[HadIt];
       
        if(pmag[HadIt] > (0.10*pmag_had_vis)){
            ROOT::Math::XYVector pT_hadron(px[HadIt], py[HadIt]);
            phi.push_back(pT_hadron.Phi());
        }
        
    }

    if(phi.size() == 0){
        return 0;
    }
    else{
        ROOT::Math::XYVector pT_miss(-px_tot, -py_tot);
        double phi_miss = pT_miss.Phi();
        double DeltaPhiMET = M_PI;
        // phi loop
        for(size_t phiIt=0; phiIt<phi.size(); phiIt++){
            double dPhiMET = std::abs(phi[phiIt] - phi_miss);
            if(dPhiMET > M_PI){
                dPhiMET = 2*M_PI - dPhiMET;
            }
            if(dPhiMET < DeltaPhiMET){
                DeltaPhiMET = dPhiMET;
            }
        }

        double DeltaPhiMET_deg = DeltaPhiMET * (180.00/M_PI);
        return DeltaPhiMET_deg;
    }

    return 0;

};

double Discriminators::Delta_phi_sum(std::vector<double> px, std::vector<double> py, std::vector<double> pz, std::vector<double> pmag, double pmag_had_vis){

    std::vector<double> phi;
        
    // hadron loop
    for(size_t HadIt=0; HadIt<pmag.size(); HadIt++){

        if(pmag[HadIt] <= (0.10*pmag_had_vis)){
            // count the hadron carrying at least 10% of the hadronic system energy
            continue;
        }

        ROOT::Math::XYVector pT_hadron(px[HadIt], py[HadIt]);
        phi.push_back(pT_hadron.Phi());

    }

    if(phi.size() <= 1){
        return 0;
    }
    else{
        // sort the phi ascendingly
        std::sort(phi.begin(), phi.end());
        double dphi_sum = 0.00;

        // phi loop
        for(size_t phiIt=0; phiIt<phi.size()-1; phiIt++){
            // azimuthal angle difference between the azimuthally neighboring track
            double dphi = phi[phiIt+1] - phi[phiIt];
            if(dphi < M_PI){
                dphi_sum += dphi;
            }
            else{
                dphi_sum += 2*M_PI - dphi;
            }
        }

        double dphi_sum_deg = dphi_sum * (180.00/M_PI);
        return dphi_sum_deg;
    }

    return 0;

};


double Discriminators::TanThetaHardest(std::vector<double> pmag, std::vector<double> theta){

    if(pmag.size() == 0) return 0;

    size_t idx_hardest = std::max_element(pmag.begin(), pmag.end()) - pmag.begin();
    double tan_theta_hardest = std::tan(theta.at(idx_hardest));

    return tan_theta_hardest;

};


void Discriminators::CalcDisc(){

    n_ch = ChTrk_Multi(pmag_reco);

    pmag_had_vis_reco = VisibleHadronicMomentum(pmag_reco);
    dphi_max_reco = Delta_phi_max(px_reco, py_reco, pz_reco, pmag_reco, pmag_had_vis_reco);
    p3_hardest_reco = HardestTrackMomentum(pmag_reco);
    InvThetaCh_reco = InverseSumChargedTrackAngle(theta_reco);
    pTmiss_mag_reco = pT_miss_mag(px_reco, py_reco);
    pTabs_sum_reco = pT_abs_sum(px_reco, py_reco);
    DeltaPhiMET_reco = TrackMETangle(px_reco, py_reco, pz_reco, pmag_reco, pmag_had_vis_reco);
    dphi_sum_reco = Delta_phi_sum(px_reco, py_reco, pz_reco, pmag_reco, pmag_had_vis_reco);
    tan_theta_hardest_reco = TanThetaHardest(pmag_reco, theta_reco);

    pmag_had_vis_true = VisibleHadronicMomentum(pmag_true);
    dphi_max_true = Delta_phi_max(px_true, py_true, pz_true, pmag_true, pmag_had_vis_true);
    p3_hardest_true = HardestTrackMomentum(pmag_true);
    InvThetaCh_true = InverseSumChargedTrackAngle(theta_true);
    pTmiss_mag_true = pT_miss_mag(px_true, py_true);
    pTabs_sum_true = pT_abs_sum(px_true, py_true);
    DeltaPhiMET_true = TrackMETangle(px_true, py_true, pz_true, pmag_true, pmag_had_vis_true);
    dphi_sum_true = Delta_phi_sum(px_true, py_true, pz_true, pmag_true, pmag_had_vis_true);
    tan_theta_hardest_true = TanThetaHardest(pmag_true, theta_true);

};

