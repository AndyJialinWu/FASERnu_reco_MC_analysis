//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar  6 13:40:10 2025 by ROOT version 6.32.04
// from TTree discNC_TrainTest/discNC_TrainTest
// found on file: PhysicsNTUP_ML.root
//////////////////////////////////////////////////////////

#ifndef discNC_TrainTest_h
#define discNC_TrainTest_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class discNC_TrainTest {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           mcID;
   Int_t           EventID;
   vector<int>     *TrackID;
   vector<int>     *PDG;
   vector<int>     *itrk;
   vector<double>  *px_reco;
   vector<double>  *py_reco;
   vector<double>  *pz_reco;
   vector<double>  *pmag_reco;
   vector<double>  *theta_reco;
   vector<double>  *phi_reco;
   vector<double>  *px_true;
   vector<double>  *py_true;
   vector<double>  *pz_true;
   vector<double>  *pmag_true;
   vector<double>  *theta_true;
   vector<double>  *phi_true;
   vector<double>  *pmag_haruhi;
   vector<double>  *pmag_ang;
   vector<double>  *pmag_coord;
   vector<int>     *nseg;
   vector<double>  *TrackLength;
   vector<float>   *dz;
   vector<double>  *IP;
   vector<int>     *PID_start;
   vector<int>     *PID_end;
   vector<int>     *MaxGap;
   vector<double>  *theta_RMS;
   vector<double>  *MaxKinkAngle;
   vector<bool>    *IsPartCatChanged;
   vector<bool>    *IsTrkIdChanged;
   Float_t         n_ch;
   Float_t         pmag_had_vis_reco;
   Float_t         dphi_max_reco;
   Float_t         p3_hardest_reco;
   Float_t         InvThetaCh_reco;
   Float_t         pTmiss_mag_reco;
   Float_t         pTabs_sum_reco;
   Float_t         DeltaPhiMET_reco;
   Float_t         dphi_sum_reco;
   Float_t         tan_theta_hardest_reco;
   Float_t         pmag_had_vis_true;
   Float_t         dphi_max_true;
   Float_t         p3_hardest_true;
   Float_t         InvThetaCh_true;
   Float_t         pTmiss_mag_true;
   Float_t         pTabs_sum_true;
   Float_t         DeltaPhiMET_true;
   Float_t         dphi_sum_true;
   Float_t         tan_theta_hardest_true;

   // List of branches
   TBranch        *b_mcID;   //!
   TBranch        *b_EventID;   //!
   TBranch        *b_TrackID;   //!
   TBranch        *b_PDG;   //!
   TBranch        *b_itrk;   //!
   TBranch        *b_px_reco;   //!
   TBranch        *b_py_reco;   //!
   TBranch        *b_pz_reco;   //!
   TBranch        *b_pmag_reco;   //!
   TBranch        *b_theta_reco;   //!
   TBranch        *b_phi_reco;   //!
   TBranch        *b_px_true;   //!
   TBranch        *b_py_true;   //!
   TBranch        *b_pz_true;   //!
   TBranch        *b_pmag_true;   //!
   TBranch        *b_theta_true;   //!
   TBranch        *b_phi_true;   //!
   TBranch        *b_pmag_haruhi;   //!
   TBranch        *b_pmag_ang;   //!
   TBranch        *b_pmag_coord;   //!
   TBranch        *b_nseg;   //!
   TBranch        *b_TrackLength;   //!
   TBranch        *b_dz;   //!
   TBranch        *b_IP;   //!
   TBranch        *b_PID_start;   //!
   TBranch        *b_PID_end;   //!
   TBranch        *b_MaxGap;   //!
   TBranch        *b_theta_RMS;   //!
   TBranch        *b_MaxKinkAngle;   //!
   TBranch        *b_IsPartCatChanged;   //!
   TBranch        *b_IsTrkIdChanged;   //!
   TBranch        *b_n_ch;   //!
   TBranch        *b_pmag_had_vis_reco;   //!
   TBranch        *b_dphi_max_reco;   //!
   TBranch        *b_p3_hardest_reco;   //!
   TBranch        *b_InvThetaCh_reco;   //!
   TBranch        *b_pTmiss_mag_reco;   //!
   TBranch        *b_pTabs_sum_reco;   //!
   TBranch        *b_DeltaPhiMET_reco;   //!
   TBranch        *b_dphi_sum_reco;   //!
   TBranch        *b_tan_theta_hardest_reco;   //!
   TBranch        *b_pmag_had_vis_true;   //!
   TBranch        *b_dphi_max_true;   //!
   TBranch        *b_p3_hardest_true;   //!
   TBranch        *b_InvThetaCh_true;   //!
   TBranch        *b_pTmiss_mag_true;   //!
   TBranch        *b_pTabs_sum_true;   //!
   TBranch        *b_DeltaPhiMET_true;   //!
   TBranch        *b_dphi_sum_true;   //!
   TBranch        *b_tan_theta_hardest_true;   //!

   discNC_TrainTest(TTree *tree=0);
   virtual ~discNC_TrainTest();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef discNC_TrainTest_cxx
discNC_TrainTest::discNC_TrainTest(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("PhysicsNTUP_ML.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("PhysicsNTUP_ML.root");
      }
      f->GetObject("discNC_TrainTest",tree);

   }
   Init(tree);
}

discNC_TrainTest::~discNC_TrainTest()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t discNC_TrainTest::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t discNC_TrainTest::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void discNC_TrainTest::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   TrackID = 0;
   PDG = 0;
   itrk = 0;
   px_reco = 0;
   py_reco = 0;
   pz_reco = 0;
   pmag_reco = 0;
   theta_reco = 0;
   phi_reco = 0;
   px_true = 0;
   py_true = 0;
   pz_true = 0;
   pmag_true = 0;
   theta_true = 0;
   phi_true = 0;
   pmag_haruhi = 0;
   pmag_ang = 0;
   pmag_coord = 0;
   nseg = 0;
   TrackLength = 0;
   dz = 0;
   IP = 0;
   PID_start = 0;
   PID_end = 0;
   MaxGap = 0;
   theta_RMS = 0;
   MaxKinkAngle = 0;
   IsPartCatChanged = 0;
   IsTrkIdChanged = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mcID", &mcID, &b_mcID);
   fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
   fChain->SetBranchAddress("TrackID", &TrackID, &b_TrackID);
   fChain->SetBranchAddress("PDG", &PDG, &b_PDG);
   fChain->SetBranchAddress("itrk", &itrk, &b_itrk);
   fChain->SetBranchAddress("px_reco", &px_reco, &b_px_reco);
   fChain->SetBranchAddress("py_reco", &py_reco, &b_py_reco);
   fChain->SetBranchAddress("pz_reco", &pz_reco, &b_pz_reco);
   fChain->SetBranchAddress("pmag_reco", &pmag_reco, &b_pmag_reco);
   fChain->SetBranchAddress("theta_reco", &theta_reco, &b_theta_reco);
   fChain->SetBranchAddress("phi_reco", &phi_reco, &b_phi_reco);
   fChain->SetBranchAddress("px_true", &px_true, &b_px_true);
   fChain->SetBranchAddress("py_true", &py_true, &b_py_true);
   fChain->SetBranchAddress("pz_true", &pz_true, &b_pz_true);
   fChain->SetBranchAddress("pmag_true", &pmag_true, &b_pmag_true);
   fChain->SetBranchAddress("theta_true", &theta_true, &b_theta_true);
   fChain->SetBranchAddress("phi_true", &phi_true, &b_phi_true);
   fChain->SetBranchAddress("pmag_haruhi", &pmag_haruhi, &b_pmag_haruhi);
   fChain->SetBranchAddress("pmag_ang", &pmag_ang, &b_pmag_ang);
   fChain->SetBranchAddress("pmag_coord", &pmag_coord, &b_pmag_coord);
   fChain->SetBranchAddress("nseg", &nseg, &b_nseg);
   fChain->SetBranchAddress("TrackLength", &TrackLength, &b_TrackLength);
   fChain->SetBranchAddress("dz", &dz, &b_dz);
   fChain->SetBranchAddress("IP", &IP, &b_IP);
   fChain->SetBranchAddress("PID_start", &PID_start, &b_PID_start);
   fChain->SetBranchAddress("PID_end", &PID_end, &b_PID_end);
   fChain->SetBranchAddress("MaxGap", &MaxGap, &b_MaxGap);
   fChain->SetBranchAddress("theta_RMS", &theta_RMS, &b_theta_RMS);
   fChain->SetBranchAddress("MaxKinkAngle", &MaxKinkAngle, &b_MaxKinkAngle);
   fChain->SetBranchAddress("IsPartCatChanged", &IsPartCatChanged, &b_IsPartCatChanged);
   fChain->SetBranchAddress("IsTrkIdChanged", &IsTrkIdChanged, &b_IsTrkIdChanged);
   fChain->SetBranchAddress("n_ch", &n_ch, &b_n_ch);
   fChain->SetBranchAddress("pmag_had_vis_reco", &pmag_had_vis_reco, &b_pmag_had_vis_reco);
   fChain->SetBranchAddress("dphi_max_reco", &dphi_max_reco, &b_dphi_max_reco);
   fChain->SetBranchAddress("p3_hardest_reco", &p3_hardest_reco, &b_p3_hardest_reco);
   fChain->SetBranchAddress("InvThetaCh_reco", &InvThetaCh_reco, &b_InvThetaCh_reco);
   fChain->SetBranchAddress("pTmiss_mag_reco", &pTmiss_mag_reco, &b_pTmiss_mag_reco);
   fChain->SetBranchAddress("pTabs_sum_reco", &pTabs_sum_reco, &b_pTabs_sum_reco);
   fChain->SetBranchAddress("DeltaPhiMET_reco", &DeltaPhiMET_reco, &b_DeltaPhiMET_reco);
   fChain->SetBranchAddress("dphi_sum_reco", &dphi_sum_reco, &b_dphi_sum_reco);
   fChain->SetBranchAddress("tan_theta_hardest_reco", &tan_theta_hardest_reco, &b_tan_theta_hardest_reco);
   fChain->SetBranchAddress("pmag_had_vis_true", &pmag_had_vis_true, &b_pmag_had_vis_true);
   fChain->SetBranchAddress("dphi_max_true", &dphi_max_true, &b_dphi_max_true);
   fChain->SetBranchAddress("p3_hardest_true", &p3_hardest_true, &b_p3_hardest_true);
   fChain->SetBranchAddress("InvThetaCh_true", &InvThetaCh_true, &b_InvThetaCh_true);
   fChain->SetBranchAddress("pTmiss_mag_true", &pTmiss_mag_true, &b_pTmiss_mag_true);
   fChain->SetBranchAddress("pTabs_sum_true", &pTabs_sum_true, &b_pTabs_sum_true);
   fChain->SetBranchAddress("DeltaPhiMET_true", &DeltaPhiMET_true, &b_DeltaPhiMET_true);
   fChain->SetBranchAddress("dphi_sum_true", &dphi_sum_true, &b_dphi_sum_true);
   fChain->SetBranchAddress("tan_theta_hardest_true", &tan_theta_hardest_true, &b_tan_theta_hardest_true);
   Notify();
}

bool discNC_TrainTest::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void discNC_TrainTest::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t discNC_TrainTest::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef discNC_TrainTest_cxx
