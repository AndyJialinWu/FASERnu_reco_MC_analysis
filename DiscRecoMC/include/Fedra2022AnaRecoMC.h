#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>


const std::string EventListFile[3] = {
    "/home/jialwu/FASERnu_reco_MC_analysis/CheckRecoMC/EventList/NC_200025.txt",
    "/home/jialwu/FASERnu_reco_MC_analysis/CheckRecoMC/EventList/NC_200026.txt",
    "/home/jialwu/FASERnu_reco_MC_analysis/CheckRecoMC/EventList/NC_200035.txt"
};
const std::string EventDataDir[3] = {
    "/home/jialwu/raid/FASERnu_MC_reco/NC/200025/",
    "/home/jialwu/raid/FASERnu_MC_reco/NC/200026/",
    "/home/jialwu/raid/FASERnu_MC_reco/NC/200035/"
};

const std::string MCtruthNTUPdir[3] = {
    "/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_NC_10invab/200025/NTUP/",
    "/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_NC_10invab/200026/NTUP/",
    "/home/jialwu/raid/FASERnu_MC_reco/MC24_Genie_NC_10invab/200035/NTUP/",
};
const std::string MCtruthNameFront[3] = {
    "FaserMC-MC24_Genie_light_eposlhc_pi_10invab-200025-",
    "FaserMC-MC24_Genie_light_eposlhc_k_10invab-200026-",
    "FaserMC-MC24_Genie_charm_p8_monash_central_10invab-200035-"
};
const std::string MCtruthNameRear[3] = {
    "-s0012-NTUP.root",
    "-s0012-NTUP.root",
    "-s0012-NTUP.root"
};
const int NumEvtsPerNTUPfile[3] = {1000, 1000, 500};
const int MonteCarloID[3] = {200025, 200026, 200035};


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

};


void PrintEventList(std::vector<std::string> * EvtList){

    for(size_t EvtIt=0; EvtIt<EvtList->size(); EvtIt++){
        std::cout<<EvtList->at(EvtIt)<<std::endl;
    }

};


std::string EvtIDtoMCtruthNTUPfname(int EventID, int mcIt){

    int fnum = EventID / NumEvtsPerNTUPfile[mcIt];

    // ensure that the number is formatted with a width of 5, filling unused spaces with 0
    std::ostringstream oss;
    oss << std::setw(5) << std::setfill('0') << fnum;
    std::string str_fnum = oss.str();

    std::string fname_NTUP = MCtruthNTUPdir[mcIt] + MCtruthNameFront[mcIt] + str_fnum + MCtruthNameRear[mcIt];

    return fname_NTUP;

};

/*
double MCSAngErrFunc(double *x, double *p){

  double error = TMath::Sqrt( (214.3296*x[0]/p[1]) * TMath::Power(1.00+0.038*TMath::Log(x[0]/p[1]), 2) * (1.00/TMath::Power(p[0], 2)) + p[2]);
  return error;

};


double MCSCoordErrFunc(double *x, double *p){
  
  double error = TMath::Sqrt( TMath::Power(p[1], 2) + (2.00/3.00) *  TMath::Power(x[0]*TMath::Sqrt(1.00 + (p[2]*p[2])), 3) * (0.0136*0.0136*p[0] / p[3]));
  return error;

};
*/
