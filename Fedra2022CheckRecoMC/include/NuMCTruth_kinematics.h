//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 21 12:41:12 2025 by ROOT version 6.32.04
// from TTree NuMCTruth_kinematics/NuMCTruth_kinematics
// found on file: /home/jialwu/raid/FASERnu_MC_reco/NC/200025/evt_7_pl114_214/jw_test.root
//////////////////////////////////////////////////////////

#ifndef NuMCTruth_kinematics_h
#define NuMCTruth_kinematics_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
using namespace std;

class NuMCTruth_kinematics {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           m_runnumber;
   Int_t           m_event_id_MC;
   Int_t           m_NC0_CC1;
   Int_t           m_target;
   vector<int>     *m_target_PDG;
   Int_t           m_N_target;
   Float_t         m_vx;
   Float_t         m_vy;
   Float_t         m_vz;
   Int_t           m_NuMu0_NuE1_NuTau2;
   Int_t           m_Nu_PDG;
   Double_t        m_Nu_px;
   Double_t        m_Nu_py;
   Double_t        m_Nu_pz;
   Double_t        m_Nu_e;
   Int_t           m_N_leptons;
   Int_t           m_N_photons;
   Int_t           m_N_hadrons;
   vector<int>     *m_leptons_PDG;
   vector<int>     *m_photons_PDG;
   vector<int>     *m_hadrons_PDG;
   vector<int>     *m_leptons_track_id;
   vector<int>     *m_photons_track_id;
   vector<int>     *m_hadrons_track_id;
   vector<double>  *m_leptons_px;
   vector<double>  *m_leptons_py;
   vector<double>  *m_leptons_pz;
   vector<double>  *m_leptons_e;
   vector<double>  *m_photons_px;
   vector<double>  *m_photons_py;
   vector<double>  *m_photons_pz;
   vector<double>  *m_photons_e;
   vector<double>  *m_hadrons_px;
   vector<double>  *m_hadrons_py;
   vector<double>  *m_hadrons_pz;
   vector<double>  *m_hadrons_e;
   vector<double>  *m_hadrons_ch;

   // List of branches
   TBranch        *b_m_runnumber;   //!
   TBranch        *b_m_event_id_MC;   //!
   TBranch        *b_m_NC0_CC1;   //!
   TBranch        *b_m_target;   //!
   TBranch        *b_m_target_PDG;   //!
   TBranch        *b_m_N_target;   //!
   TBranch        *b_m_vx;   //!
   TBranch        *b_m_vy;   //!
   TBranch        *b_m_vz;   //!
   TBranch        *b_m_NuMu0_NuE1_NuTau2;   //!
   TBranch        *b_m_Nu_PDG;   //!
   TBranch        *b_m_Nu_px;   //!
   TBranch        *b_m_Nu_py;   //!
   TBranch        *b_m_Nu_pz;   //!
   TBranch        *b_m_Nu_e;   //!
   TBranch        *b_m_N_leptons;   //!
   TBranch        *b_m_N_photons;   //!
   TBranch        *b_m_N_hadrons;   //!
   TBranch        *b_m_leptons_PDG;   //!
   TBranch        *b_m_photons_PDG;   //!
   TBranch        *b_m_hadrons_PDG;   //!
   TBranch        *b_m_leptons_track_id;   //!
   TBranch        *b_m_photons_track_id;   //!
   TBranch        *b_m_hadrons_track_id;   //!
   TBranch        *b_m_leptons_px;   //!
   TBranch        *b_m_leptons_py;   //!
   TBranch        *b_m_leptons_pz;   //!
   TBranch        *b_m_leptons_e;   //!
   TBranch        *b_m_photons_px;   //!
   TBranch        *b_m_photons_py;   //!
   TBranch        *b_m_photons_pz;   //!
   TBranch        *b_m_photons_e;   //!
   TBranch        *b_m_hadrons_px;   //!
   TBranch        *b_m_hadrons_py;   //!
   TBranch        *b_m_hadrons_pz;   //!
   TBranch        *b_m_hadrons_e;   //!
   TBranch        *b_m_hadrons_ch;   //!

   NuMCTruth_kinematics(TTree *tree=0);
   virtual ~NuMCTruth_kinematics();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};


NuMCTruth_kinematics::NuMCTruth_kinematics(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/jialwu/raid/FASERnu_MC_reco/NC/200025/evt_7_pl114_214/jw_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/home/jialwu/raid/FASERnu_MC_reco/NC/200025/evt_7_pl114_214/jw_test.root");
      }
      f->GetObject("NuMCTruth_kinematics",tree);

   }
   Init(tree);
};

NuMCTruth_kinematics::~NuMCTruth_kinematics()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
};

Int_t NuMCTruth_kinematics::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
};
Long64_t NuMCTruth_kinematics::LoadTree(Long64_t entry)
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
};

void NuMCTruth_kinematics::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   m_target_PDG = 0;
   m_leptons_PDG = 0;
   m_photons_PDG = 0;
   m_hadrons_PDG = 0;
   m_leptons_track_id = 0;
   m_photons_track_id = 0;
   m_hadrons_track_id = 0;
   m_leptons_px = 0;
   m_leptons_py = 0;
   m_leptons_pz = 0;
   m_leptons_e = 0;
   m_photons_px = 0;
   m_photons_py = 0;
   m_photons_pz = 0;
   m_photons_e = 0;
   m_hadrons_px = 0;
   m_hadrons_py = 0;
   m_hadrons_pz = 0;
   m_hadrons_e = 0;
   m_hadrons_ch = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("m_runnumber", &m_runnumber, &b_m_runnumber);
   fChain->SetBranchAddress("m_event_id_MC", &m_event_id_MC, &b_m_event_id_MC);
   fChain->SetBranchAddress("m_NC0_CC1", &m_NC0_CC1, &b_m_NC0_CC1);
   fChain->SetBranchAddress("m_target", &m_target, &b_m_target);
   fChain->SetBranchAddress("m_target_PDG", &m_target_PDG, &b_m_target_PDG);
   fChain->SetBranchAddress("m_N_target", &m_N_target, &b_m_N_target);
   fChain->SetBranchAddress("m_vx", &m_vx, &b_m_vx);
   fChain->SetBranchAddress("m_vy", &m_vy, &b_m_vy);
   fChain->SetBranchAddress("m_vz", &m_vz, &b_m_vz);
   fChain->SetBranchAddress("m_NuMu0_NuE1_NuTau2", &m_NuMu0_NuE1_NuTau2, &b_m_NuMu0_NuE1_NuTau2);
   fChain->SetBranchAddress("m_Nu_PDG", &m_Nu_PDG, &b_m_Nu_PDG);
   fChain->SetBranchAddress("m_Nu_px", &m_Nu_px, &b_m_Nu_px);
   fChain->SetBranchAddress("m_Nu_py", &m_Nu_py, &b_m_Nu_py);
   fChain->SetBranchAddress("m_Nu_pz", &m_Nu_pz, &b_m_Nu_pz);
   fChain->SetBranchAddress("m_Nu_e", &m_Nu_e, &b_m_Nu_e);
   fChain->SetBranchAddress("m_N_leptons", &m_N_leptons, &b_m_N_leptons);
   fChain->SetBranchAddress("m_N_photons", &m_N_photons, &b_m_N_photons);
   fChain->SetBranchAddress("m_N_hadrons", &m_N_hadrons, &b_m_N_hadrons);
   fChain->SetBranchAddress("m_leptons_PDG", &m_leptons_PDG, &b_m_leptons_PDG);
   fChain->SetBranchAddress("m_photons_PDG", &m_photons_PDG, &b_m_photons_PDG);
   fChain->SetBranchAddress("m_hadrons_PDG", &m_hadrons_PDG, &b_m_hadrons_PDG);
   fChain->SetBranchAddress("m_leptons_track_id", &m_leptons_track_id, &b_m_leptons_track_id);
   fChain->SetBranchAddress("m_photons_track_id", &m_photons_track_id, &b_m_photons_track_id);
   fChain->SetBranchAddress("m_hadrons_track_id", &m_hadrons_track_id, &b_m_hadrons_track_id);
   fChain->SetBranchAddress("m_leptons_px", &m_leptons_px, &b_m_leptons_px);
   fChain->SetBranchAddress("m_leptons_py", &m_leptons_py, &b_m_leptons_py);
   fChain->SetBranchAddress("m_leptons_pz", &m_leptons_pz, &b_m_leptons_pz);
   fChain->SetBranchAddress("m_leptons_e", &m_leptons_e, &b_m_leptons_e);
   fChain->SetBranchAddress("m_photons_px", &m_photons_px, &b_m_photons_px);
   fChain->SetBranchAddress("m_photons_py", &m_photons_py, &b_m_photons_py);
   fChain->SetBranchAddress("m_photons_pz", &m_photons_pz, &b_m_photons_pz);
   fChain->SetBranchAddress("m_photons_e", &m_photons_e, &b_m_photons_e);
   fChain->SetBranchAddress("m_hadrons_px", &m_hadrons_px, &b_m_hadrons_px);
   fChain->SetBranchAddress("m_hadrons_py", &m_hadrons_py, &b_m_hadrons_py);
   fChain->SetBranchAddress("m_hadrons_pz", &m_hadrons_pz, &b_m_hadrons_pz);
   fChain->SetBranchAddress("m_hadrons_e", &m_hadrons_e, &b_m_hadrons_e);
   fChain->SetBranchAddress("m_hadrons_ch", &m_hadrons_ch, &b_m_hadrons_ch);
   Notify();
};

bool NuMCTruth_kinematics::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
};

void NuMCTruth_kinematics::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
};
Int_t NuMCTruth_kinematics::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
};

void NuMCTruth_kinematics::Loop()
{
//   In a ROOT session, you can do:
//      root> .L NuMCTruth_kinematics.C
//      root> NuMCTruth_kinematics t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
};

#endif // #ifndef NuMCTruth_kinematics_h
