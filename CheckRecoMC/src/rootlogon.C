#include "TSystem.h"

void rootlogon(){
	gSystem->Load("libMatrix");
	gSystem->Load("libTree");
	gSystem->Load("libHist");
	gSystem->Load("libPhysics");
	gSystem->Load("libEve");
	gSystem->Load("libGeom");
	gSystem->Load("libEve");
	gSystem->Load("libvt");
	gSystem->Load("libEphys");
	gSystem->Load("libEmath");
	gSystem->Load("libEdb");
	gSystem->Load("libEbase");
	gSystem->Load("libEdr");
	gSystem->Load("libEIO");
	gSystem->Load("libAlignment");
	gSystem->Load("libAnalysis");
	gSystem->Load("libScan");
	gSystem->Load("libDataConversion");
	gSystem->Load("libEGA");
	gSystem->Load("libEdd");
	gSystem->Load("libEMC");
	gSystem->Load("libShower");
	gSystem->Load("libEmr");
	gSystem->Load("libEDA");
};

