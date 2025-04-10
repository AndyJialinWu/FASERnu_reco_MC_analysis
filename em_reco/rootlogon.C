#include "TSystem.h"
void rootlogon(){
	gSystem->Load("libvt");
	gSystem->Load("libEmath");
	gSystem->Load("libEdb");
	gSystem->Load("libEbase");
	gSystem->Load("libEdr");
	gSystem->Load("libEIO");
	gSystem->Load("libEDA");
	gSystem->Load("libScan");
	gSystem->Load("libShower");
	//gSystem->Load("FitAlignTracks_C.dll");
}
