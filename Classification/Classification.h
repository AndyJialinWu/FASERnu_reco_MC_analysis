#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TPDGCode.h"

int Classfier(int m_NC0_CC1, int m_Nu_PDG, int m_N_target){

    if(m_N_target == 1){
        if(m_NC0_CC1 == 0){
            return 3; // NC
        } 
        else if(m_NC0_CC1 == 1){
            if(std::abs(m_Nu_PDG) == kNuE) return 0; // nue CC
            if(std::abs(m_Nu_PDG) == kNuMu) return 1; // numu CC
            if(std::abs(m_Nu_PDG) == kNuTau) return 2; // nutau CC
        }
        else{
            return -1;
        }

        return -1; 
    }
    else{
        return -1;
    }

    return -1;

}