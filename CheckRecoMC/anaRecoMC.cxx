#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

//#include "tracks.C"
//#include "NuMCTruth_kinematics.C"
#include "./include/anaRecoMC.h"
//#include "line3Df

int main(){

    std::vector<std::string> *EvtNC200025 = new std::vector<std::string>();
    std::vector<std::string> *EvtNC200026 = new std::vector<std::string>();
    std::vector<std::string> *EvtNC200035 = new std::vector<std::string>();

    ReadEventList(EvtNC200025, "/home/jialwu/FASERnu_reco_MC_analysis/CheckRecoMC/EventList/NC_200025.txt");
    ReadEventList(EvtNC200026, "/home/jialwu/FASERnu_reco_MC_analysis/CheckRecoMC/EventList/NC_200026.txt");
    ReadEventList(EvtNC200035, "/home/jialwu/FASERnu_reco_MC_analysis/CheckRecoMC/EventList/NC_200035.txt");
    //PrintEventList(EvtNC200025);

    const std::size_t NEvtsNC[3] = {EvtNC200025->size(), EvtNC200026->size(), EvtNC200035->size()};
    std::cout<<"# NC events from MC200025 = "<<NEvtsNC[0]<<std::endl;
    std::cout<<"# NC events from MC200026 = "<<NEvtsNC[1]<<std::endl;
    std::cout<<"# NC events from MC200035 = "<<NEvtsNC[2]<<std::endl;

    const double Run3ExpNC[3] = {IntLumi_run3 * wgt_NC_200025, IntLumi_run3 * wgt_NC_200026, IntLumi_run3 * wgt_NC_200035};

    CheckRecoPlot *CRP = new CheckRecoPlot();

    ana_MC_event("/home/jialwu/raid/FASERnu_MC_reco/NC/200025/", EvtNC200025, Run3ExpNC[0]/NEvtsNC[0], CRP);
    //ana_MC_event("/home/jialwu/raid/FASERnu_MC_reco/NC/200026/", EvtNC200026, Run3ExpNC[1]/NEvtsNC[1], CRP);
    //ana_MC_event("/home/jialwu/raid/FASERnu_MC_reco/NC/200035/", EvtNC200035, Run3ExpNC[2]/NEvtsNC[2], CRP);

    //CRP->PlotHist("/home/jialwu/FASERnu_reco_MC_analysis/CheckRecoMC/Figures");
    //CRP->StoreHist2ROOT("/home/jialwu/FASERnu_reco_MC_analysis/CheckRecoMC/Figures");

    return 0;
    
}

