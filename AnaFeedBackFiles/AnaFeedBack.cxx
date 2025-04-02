#include <string>
#include "AnaFeedBack.h"
//#include"EdbEDAUtil.h"

int main(){

    //EdbEDA *gEDA = new EdbEDA("F222_zone4_p489_100197.6_18577.9_v03262025.feedback");
    //gEDA->SetBeamAngle(-0.0037, 0.0000);        //(TX, TY)
    //gEDA->Run();                                // Load all data and open the GUI

    std::string filename = "F222_zone4_p489_100197.6_18577.9_v03262025.feedback";
    ParseFeedBack("feedback/" + filename);
    
    //char filename[] = "feedback/F222_zone4_p489_100197.6_18577.9_v03262025.feedback";
    //EdbPVRec * pvr = EdbEDAUtil::ReadFeedbackPVR(filename);
    //EdbDataProc *dproc = new EdbDataProc();
	//dproc->MakeTracksTree(pvr, "linked_tracks.root");


    return 0;

}
