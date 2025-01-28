#include <iostream>
#include <fstream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"

#include "Classification.h"

void Classification(std::string input_EvtVtxInfo, std::string output_EvtVtxInfo, TString dir_jw_test, int MCit, TString fileNum){

    std::ifstream infile(input_EvtVtxInfo);
    std::string line;

    std::ofstream reco_evt_vtx_info;
	reco_evt_vtx_info.open(output_EvtVtxInfo);

    TString name_jw_test;
    if(MCit == 0) name_jw_test = "FaserMC-MC24_Genie_light_eposlhc_pi_10invab-200025-" + fileNum + "-s0012-NTUP_jw_test.root";
    if(MCit == 1) name_jw_test = "FaserMC-MC24_Genie_light_eposlhc_k_10invab-200026-" + fileNum + "-s0012-NTUP_jw_test.root";
    if(MCit == 2) name_jw_test = "FaserMC-MC24_Genie_charm_p8_monash_central_10invab-200035-" + fileNum + "-s0012-NTUP_jw_test.root";

    TFile *f_jw_test = new TFile(dir_jw_test+name_jw_test, "READ");
    TTree *t_jw_test = (TTree *)f_jw_test->Get("NuMCTruth_kinematics");
    int NEvts = t_jw_test->GetEntries();
    std::cout<<"# events in jw_test = "<<NEvts<<std::endl;

    int m_event_id_MC, m_NC0_CC1, m_Nu_PDG, m_N_target;
    t_jw_test->SetBranchAddress("m_event_id_MC", &m_event_id_MC);
    t_jw_test->SetBranchAddress("m_NC0_CC1", &m_NC0_CC1);
    t_jw_test->SetBranchAddress("m_Nu_PDG", &m_Nu_PDG);
    t_jw_test->SetBranchAddress("m_N_target", &m_N_target);

    int EvtIt = 0;
    t_jw_test->GetEntry(EvtIt);
    
    if (infile.is_open()) {

        // looping through each line
        while (std::getline(infile, line)) {
            
            // 0: Event ID, 1: vx, 2: vy, 3: starting plate, 4: ending plate
            std::vector<double> EvtVtxInfo;

            // Remove the parentheses
            line.erase(std::remove(line.begin(), line.end(), '('), line.end());
            line.erase(std::remove(line.begin(), line.end(), ')'), line.end());

            // Create a stringstream to extract numbers
            std::stringstream ss(line);
            double num;

            // Extract numbers and store them in the vector
            while (ss >> num) {
                EvtVtxInfo.push_back(num);
            }

            // Output the extracted numbers for this line
            std::cout << "Extracted numbers: ";
            for (double n : EvtVtxInfo) {
                std::cout << n << " ";
            }
            std::cout << std::endl;

            
            while (m_event_id_MC < EvtVtxInfo[0]) {
                EvtIt++;
                t_jw_test->GetEntry(EvtIt);
            }

            if(m_event_id_MC == EvtVtxInfo[0]){
                // 0: nue CC, 1: numu CC, 2: nutau CC, 3: NC
                int class_flag = Classfier(m_NC0_CC1, m_Nu_PDG, m_N_target);
                EvtVtxInfo.push_back(class_flag);
                reco_evt_vtx_info<<"("<<EvtVtxInfo[0]<<" "<<EvtVtxInfo[1]<<" "<<EvtVtxInfo[2]<<" "<<EvtVtxInfo[3]<<" "<<EvtVtxInfo[4]<<" "<<EvtVtxInfo[5]<<")\n";
            }

        }

        // Close the file
        infile.close();
        reco_evt_vtx_info.close();

    } 
    else {
        std::cerr << "Unable to open the file." << std::endl;
    }

}