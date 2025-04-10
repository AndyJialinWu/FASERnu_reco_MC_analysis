//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 28 09:54:57 2025 by ROOT version 6.32.04
// from TTree NuMCTruth_kinematics/NuMCTruth_kinematics
// found on file: FaserMC-MC24_PG_neut_in_fasernu_xin-100069-00000-s0013-NTUP_jw_test.root
//////////////////////////////////////////////////////////

#ifndef NuMCTruth_kinematics_h
#define NuMCTruth_kinematics_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class NuMCTruth_kinematics {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           m_runnumber;
   Int_t           m_event_id_MC;
   Float_t         m_vx;
   Float_t         m_vy;
   Float_t         m_vz;
   Int_t           m_Nu_PDG;
   Double_t        m_Nu_px;
   Double_t        m_Nu_py;
   Double_t        m_Nu_pz;
   Double_t        m_Nu_e;
   std::vector<int>     *m_leptons_PDG;
   std::vector<int>     *m_photons_PDG;
   std::vector<int>     *m_hadrons_PDG;
   std::vector<int>     *m_leptons_track_id;
   std::vector<int>     *m_photons_track_id;
   std::vector<int>     *m_hadrons_track_id;
   std::vector<double>  *m_leptons_px;
   std::vector<double>  *m_leptons_py;
   std::vector<double>  *m_leptons_pz;
   std::vector<double>  *m_leptons_e;
   std::vector<double>  *m_photons_px;
   std::vector<double>  *m_photons_py;
   std::vector<double>  *m_photons_pz;
   std::vector<double>  *m_photons_e;
   std::vector<double>  *m_hadrons_px;
   std::vector<double>  *m_hadrons_py;
   std::vector<double>  *m_hadrons_pz;
   std::vector<double>  *m_hadrons_e;
   std::vector<double>  *m_hadrons_ch;
   //std::vector<double>  *m_hadrons_m;
   //std::vector<double>  *m_hadrons_Ekin;

   // List of branches
   TBranch        *b_m_runnumber;   //!
   TBranch        *b_m_event_id_MC;   //!
   TBranch        *b_m_vx;   //!
   TBranch        *b_m_vy;   //!
   TBranch        *b_m_vz;   //!
   TBranch        *b_m_Nu_PDG;   //!
   TBranch        *b_m_Nu_px;   //!
   TBranch        *b_m_Nu_py;   //!
   TBranch        *b_m_Nu_pz;   //!
   TBranch        *b_m_Nu_e;   //!
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
   //TBranch        *b_m_hadrons_m;   //!
   //TBranch        *b_m_hadrons_Ekin;   //!

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

#endif

#ifdef NuMCTruth_kinematics_cxx
NuMCTruth_kinematics::NuMCTruth_kinematics(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("FaserMC-MC24_PG_neut_in_fasernu_xin-100069-00000-s0013-NTUP_jw_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("FaserMC-MC24_PG_neut_in_fasernu_xin-100069-00000-s0013-NTUP_jw_test.root");
      }
      f->GetObject("NuMCTruth_kinematics",tree);

   }
   Init(tree);
}

NuMCTruth_kinematics::~NuMCTruth_kinematics()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NuMCTruth_kinematics::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
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
}

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
   //m_hadrons_m = 0;
   //m_hadrons_Ekin = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("m_runnumber", &m_runnumber, &b_m_runnumber);
   fChain->SetBranchAddress("m_event_id_MC", &m_event_id_MC, &b_m_event_id_MC);
   fChain->SetBranchAddress("m_vx", &m_vx, &b_m_vx);
   fChain->SetBranchAddress("m_vy", &m_vy, &b_m_vy);
   fChain->SetBranchAddress("m_vz", &m_vz, &b_m_vz);
   fChain->SetBranchAddress("m_Nu_PDG", &m_Nu_PDG, &b_m_Nu_PDG);
   fChain->SetBranchAddress("m_Nu_px", &m_Nu_px, &b_m_Nu_px);
   fChain->SetBranchAddress("m_Nu_py", &m_Nu_py, &b_m_Nu_py);
   fChain->SetBranchAddress("m_Nu_pz", &m_Nu_pz, &b_m_Nu_pz);
   fChain->SetBranchAddress("m_Nu_e", &m_Nu_e, &b_m_Nu_e);
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
   //fChain->SetBranchAddress("m_hadrons_m", &m_hadrons_m, &b_m_hadrons_m);
   //fChain->SetBranchAddress("m_hadrons_Ekin", &m_hadrons_Ekin, &b_m_hadrons_Ekin);
   Notify();
}

bool NuMCTruth_kinematics::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void NuMCTruth_kinematics::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NuMCTruth_kinematics::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NuMCTruth_kinematics_cxx
