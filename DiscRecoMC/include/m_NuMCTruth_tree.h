//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb  7 10:23:08 2025 by ROOT version 6.32.08
// from TTree m_NuMCTruth_tree/Tree containing simulated truth data
// found on file: FaserMC-MC24_Genie_light_eposlhc_k_10invab-200026-00000-s0012-NTUP.root
//////////////////////////////////////////////////////////

#ifndef m_NuMCTruth_tree_h
#define m_NuMCTruth_tree_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

// Header file for the classes stored in the TTree if any.
#include <vector>

class m_NuMCTruth_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t                m_runnumber;
   Int_t                m_event_id_MC;
   Int_t                m_track_id;
   Int_t                m_pdg_id;
   Float_t              m_vx_prod;
   Float_t              m_vy_prod;
   Float_t              m_vz_prod;
   Float_t              m_vx_decay;
   Float_t              m_vy_decay;
   Float_t              m_vz_decay;
   Double_t             m_px;
   Double_t             m_py;
   Double_t             m_pz;
   Double_t             m_energy;
   Double_t             m_kinetic_energy;
   Double_t             m_mass;
   Int_t                m_num_in_particle;
   Int_t                m_num_out_particle;
   Int_t                m_trackid_begin_in_particle;
   Int_t                m_trackid_begin_out_particle;
   Int_t                m_trackid_end_in_particle;
   Int_t                m_trackid_end_out_particle;
   std::vector<int>     *m_pdg_in_particle;
   std::vector<int>     *m_pdg_out_particle;
   std::vector<int>     *m_trackid_in_particle;
   std::vector<int>     *m_trackid_out_particle;
   Int_t                m_status;

   // List of branches
   TBranch        *b_m_runnumber;   //!
   TBranch        *b_m_event_id_MC;   //!
   TBranch        *b_m_track_id;   //!
   TBranch        *b_m_pdg_id;   //!
   TBranch        *b_m_vx_prod;   //!
   TBranch        *b_m_vy_prod;   //!
   TBranch        *b_m_vz_prod;   //!
   TBranch        *b_m_vx_decay;   //!
   TBranch        *b_m_vy_decay;   //!
   TBranch        *b_m_vz_decay;   //!
   TBranch        *b_m_px;   //!
   TBranch        *b_m_py;   //!
   TBranch        *b_m_pz;   //!
   TBranch        *b_m_energy;   //!
   TBranch        *b_m_kinetic_energy;   //!
   TBranch        *b_m_mass;   //!
   TBranch        *b_m_num_in_particle;   //!
   TBranch        *b_m_num_out_particle;   //!
   TBranch        *b_m_trackid_begin_in_particle;   //!
   TBranch        *b_m_trackid_begin_out_particle;   //!
   TBranch        *b_m_trackid_end_in_particle;   //!
   TBranch        *b_m_trackid_end_out_particle;   //!
   TBranch        *b_m_pdg_in_particle;   //!
   TBranch        *b_m_pdg_out_particle;   //!
   TBranch        *b_m_trackid_in_particle;   //!
   TBranch        *b_m_trackid_out_particle;   //!
   TBranch        *b_m_status;   //!

   m_NuMCTruth_tree(TTree *tree=0);
   virtual ~m_NuMCTruth_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif


m_NuMCTruth_tree::m_NuMCTruth_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("FaserMC-MC24_Genie_light_eposlhc_k_10invab-200026-00000-s0012-NTUP.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("FaserMC-MC24_Genie_light_eposlhc_k_10invab-200026-00000-s0012-NTUP.root");
      }
      f->GetObject("m_NuMCTruth_tree",tree);

   }
   Init(tree);
}

m_NuMCTruth_tree::~m_NuMCTruth_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t m_NuMCTruth_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t m_NuMCTruth_tree::LoadTree(Long64_t entry)
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

void m_NuMCTruth_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   m_pdg_in_particle = 0;
   m_pdg_out_particle = 0;
   m_trackid_in_particle = 0;
   m_trackid_out_particle = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("m_runnumber", &m_runnumber, &b_m_runnumber);
   fChain->SetBranchAddress("m_event_id_MC", &m_event_id_MC, &b_m_event_id_MC);
   fChain->SetBranchAddress("m_track_id", &m_track_id, &b_m_track_id);
   fChain->SetBranchAddress("m_pdg_id", &m_pdg_id, &b_m_pdg_id);
   fChain->SetBranchAddress("m_vx_prod", &m_vx_prod, &b_m_vx_prod);
   fChain->SetBranchAddress("m_vy_prod", &m_vy_prod, &b_m_vy_prod);
   fChain->SetBranchAddress("m_vz_prod", &m_vz_prod, &b_m_vz_prod);
   fChain->SetBranchAddress("m_vx_decay", &m_vx_decay, &b_m_vx_decay);
   fChain->SetBranchAddress("m_vy_decay", &m_vy_decay, &b_m_vy_decay);
   fChain->SetBranchAddress("m_vz_decay", &m_vz_decay, &b_m_vz_decay);
   fChain->SetBranchAddress("m_px", &m_px, &b_m_px);
   fChain->SetBranchAddress("m_py", &m_py, &b_m_py);
   fChain->SetBranchAddress("m_pz", &m_pz, &b_m_pz);
   fChain->SetBranchAddress("m_energy", &m_energy, &b_m_energy);
   fChain->SetBranchAddress("m_kinetic_energy", &m_kinetic_energy, &b_m_kinetic_energy);
   fChain->SetBranchAddress("m_mass", &m_mass, &b_m_mass);
   fChain->SetBranchAddress("m_num_in_particle", &m_num_in_particle, &b_m_num_in_particle);
   fChain->SetBranchAddress("m_num_out_particle", &m_num_out_particle, &b_m_num_out_particle);
   fChain->SetBranchAddress("m_trackid_begin_in_particle", &m_trackid_begin_in_particle, &b_m_trackid_begin_in_particle);
   fChain->SetBranchAddress("m_trackid_begin_out_particle", &m_trackid_begin_out_particle, &b_m_trackid_begin_out_particle);
   fChain->SetBranchAddress("m_trackid_end_in_particle", &m_trackid_end_in_particle, &b_m_trackid_end_in_particle);
   fChain->SetBranchAddress("m_trackid_end_out_particle", &m_trackid_end_out_particle, &b_m_trackid_end_out_particle);
   fChain->SetBranchAddress("m_pdg_in_particle", &m_pdg_in_particle, &b_m_pdg_in_particle);
   fChain->SetBranchAddress("m_pdg_out_particle", &m_pdg_out_particle, &b_m_pdg_out_particle);
   fChain->SetBranchAddress("m_trackid_in_particle", &m_trackid_in_particle, &b_m_trackid_in_particle);
   fChain->SetBranchAddress("m_trackid_out_particle", &m_trackid_out_particle, &b_m_trackid_out_particle);
   fChain->SetBranchAddress("m_status", &m_status, &b_m_status);
   Notify();
}

bool m_NuMCTruth_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void m_NuMCTruth_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t m_NuMCTruth_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void m_NuMCTruth_tree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L m_NuMCTruth_tree.C
//      root> m_NuMCTruth_tree t
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
}

