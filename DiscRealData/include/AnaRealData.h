#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <regex>
#include <filesystem>

namespace fs = std::filesystem;

std::vector<std::string> GetLinkedTracksFileNames(std::string directory){

    std::vector<std::string> fileNames;
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.path().extension() == ".root") {
            fileNames.push_back(entry.path().filename().string());
        }
    }

    // Output the collected file names
    std::cout<<"Found "<<fileNames.size()<<" files in "<<directory<<std::endl;
    for (const auto& fileName : fileNames) {
        std::cout << fileName << std::endl;
    }

    return fileNames;

}





