#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include "AnaFeedBack.h"
//#include"EdbEDAUtil.h"

namespace fs = std::filesystem;

int main(){

    //EdbEDA *gEDA = new EdbEDA("F222_zone4_p489_100197.6_18577.9_v03262025.feedback");
    //gEDA->SetBeamAngle(-0.0037, 0.0000);        //(TX, TY)
    //gEDA->Run();                                // Load all data and open the GUI

    std::string FeedBackAddress = "feedback/FASERnu_CC_Candidates/numu/";
    std::vector<std::string> FeedBackFiles;

    try {
        for (const auto& entry : fs::directory_iterator(FeedBackAddress)) {
            if (entry.is_regular_file() && entry.path().extension() == ".feedback") {
                FeedBackFiles.push_back(entry.path().filename().string());
            }
        }

        // Output the collected file names
        std::cout << "Found " << FeedBackFiles.size() << " .feedback files:\n";
        for (const auto& file : FeedBackFiles) {
            std::cout << file << std::endl;
        }
    } 
    catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }

    for(size_t FileIt=0; FileIt<FeedBackFiles.size(); FileIt++){

        std::string filename = FeedBackAddress + FeedBackFiles[FileIt];
        std::cout << "Processing file: " << filename << std::endl;
        ParseFeedBack(FeedBackAddress, FeedBackFiles[FileIt]);
        
    }


    //std::string filename = "F222_zone4_p489_100197.6_18577.9_v03262025.feedback";
    //ParseFeedBack("feedback/" + filename);
    //char filename[] = "feedback/F222_zone4_p489_100197.6_18577.9_v03262025.feedback";
    //EdbPVRec * pvr = EdbEDAUtil::ReadFeedbackPVR(filename);
    //EdbDataProc *dproc = new EdbDataProc();
	//dproc->MakeTracksTree(pvr, "linked_tracks.root");

    return 0;

}
