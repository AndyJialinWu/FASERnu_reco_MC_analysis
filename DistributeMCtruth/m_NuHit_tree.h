//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb  7 10:22:46 2025 by ROOT version 6.32.08
// from TTree m_NuHit_tree/Tree containing emulsion hits data
// found on file: FaserMC-MC24_Genie_light_eposlhc_k_10invab-200026-00000-s0012-NTUP.root
//////////////////////////////////////////////////////////

#ifndef m_NuHit_tree_h
#define m_NuHit_tree_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

// Header file for the classes stored in the TTree if any.

class m_NuHit_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           m_runnumber;
   Int_t           m_event_id;
   Float_t         m_x_start;
   Float_t         m_y_start;
   Float_t         m_z_start;
   Float_t         m_x_end;
   Float_t         m_y_end;
   Float_t         m_z_end;
   Float_t         m_slope_x;
   Float_t         m_slope_y;
   Float_t         m_slope_z;
   Int_t           m_track_id_hit;
   Int_t           m_track_id_MC;
   Int_t           m_pdg_MC;
   Double_t        m_px_MC;
   Double_t        m_py_MC;
   Double_t        m_pz_MC;
   Double_t        m_energy_MC;
   Double_t        m_kinetic_energy_MC;
   Double_t        m_mass_MC;
   Float_t         m_vx_prod_hit;
   Float_t         m_vy_prod_hit;
   Float_t         m_vz_prod_hit;
   Float_t         m_vx_decay_hit;
   Float_t         m_vy_decay_hit;
   Float_t         m_vz_decay_hit;
   Float_t         m_x_start_global;
   Float_t         m_y_start_global;
   Float_t         m_z_start_global;
   Float_t         m_x_end_global;
   Float_t         m_y_end_global;
   Float_t         m_z_end_global;
   Float_t         m_slope_x_global;
   Float_t         m_slope_y_global;
   Float_t         m_slope_z_global;
   Float_t         m_emeloss;
   Int_t           m_plate;
   Int_t           m_base;
   Int_t           m_emmodule;
   Int_t           m_film; //0-bottom, 1-top

   // List of branches
   TBranch        *b_m_runnumber;   //!
   TBranch        *b_m_event_id;   //!
   TBranch        *b_m_x_start;   //!
   TBranch        *b_m_y_start;   //!
   TBranch        *b_m_z_start;   //!
   TBranch        *b_m_x_end;   //!
   TBranch        *b_m_y_end;   //!
   TBranch        *b_m_z_end;   //!
   TBranch        *b_m_slope_x;   //!
   TBranch        *b_m_slope_y;   //!
   TBranch        *b_m_slope_z;   //!
   TBranch        *b_m_track_id_hit;   //!
   TBranch        *b_m_track_id_MC;   //!
   TBranch        *b_m_pdg_MC;   //!
   TBranch        *b_m_px_MC;   //!
   TBranch        *b_m_py_MC;   //!
   TBranch        *b_m_pz_MC;   //!
   TBranch        *b_m_energy_MC;   //!
   TBranch        *b_m_kinetic_energy_MC;   //!
   TBranch        *b_m_mass_MC;   //!
   TBranch        *b_m_vx_prod_hit;   //!
   TBranch        *b_m_vy_prod_hit;   //!
   TBranch        *b_m_vz_prod_hit;   //!
   TBranch        *b_m_vx_decay_hit;   //!
   TBranch        *b_m_vy_decay_hit;   //!
   TBranch        *b_m_vz_decay_hit;   //!
   TBranch        *b_m_x_start_global;   //!
   TBranch        *b_m_y_start_global;   //!
   TBranch        *b_m_z_start_global;   //!
   TBranch        *b_m_x_end_global;   //!
   TBranch        *b_m_y_end_global;   //!
   TBranch        *b_m_z_end_global;   //!
   TBranch        *b_m_slope_x_global;   //!
   TBranch        *b_m_slope_y_global;   //!
   TBranch        *b_m_slope_z_global;   //!
   TBranch        *b_m_emeloss;   //!
   TBranch        *b_m_plate;   //!
   TBranch        *b_m_base;   //!
   TBranch        *b_m_emmodule;   //!
   TBranch        *b_m_film;   //!

   m_NuHit_tree(TTree *tree=0);
   virtual ~m_NuHit_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif


m_NuHit_tree::m_NuHit_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("FaserMC-MC24_Genie_light_eposlhc_k_10invab-200026-00000-s0012-NTUP.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("FaserMC-MC24_Genie_light_eposlhc_k_10invab-200026-00000-s0012-NTUP.root");
      }
      f->GetObject("m_NuHit_tree",tree);

   }
   Init(tree);
}

m_NuHit_tree::~m_NuHit_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t m_NuHit_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t m_NuHit_tree::LoadTree(Long64_t entry)
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

void m_NuHit_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("m_runnumber", &m_runnumber, &b_m_runnumber);
   fChain->SetBranchAddress("m_event_id", &m_event_id, &b_m_event_id);
   fChain->SetBranchAddress("m_x_start", &m_x_start, &b_m_x_start);
   fChain->SetBranchAddress("m_y_start", &m_y_start, &b_m_y_start);
   fChain->SetBranchAddress("m_z_start", &m_z_start, &b_m_z_start);
   fChain->SetBranchAddress("m_x_end", &m_x_end, &b_m_x_end);
   fChain->SetBranchAddress("m_y_end", &m_y_end, &b_m_y_end);
   fChain->SetBranchAddress("m_z_end", &m_z_end, &b_m_z_end);
   fChain->SetBranchAddress("m_slope_x", &m_slope_x, &b_m_slope_x);
   fChain->SetBranchAddress("m_slope_y", &m_slope_y, &b_m_slope_y);
   fChain->SetBranchAddress("m_slope_z", &m_slope_z, &b_m_slope_z);
   fChain->SetBranchAddress("m_track_id_hit", &m_track_id_hit, &b_m_track_id_hit);
   fChain->SetBranchAddress("m_track_id_MC", &m_track_id_MC, &b_m_track_id_MC);
   fChain->SetBranchAddress("m_pdg_MC", &m_pdg_MC, &b_m_pdg_MC);
   fChain->SetBranchAddress("m_px_MC", &m_px_MC, &b_m_px_MC);
   fChain->SetBranchAddress("m_py_MC", &m_py_MC, &b_m_py_MC);
   fChain->SetBranchAddress("m_pz_MC", &m_pz_MC, &b_m_pz_MC);
   fChain->SetBranchAddress("m_energy_MC", &m_energy_MC, &b_m_energy_MC);
   fChain->SetBranchAddress("m_kinetic_energy_MC", &m_kinetic_energy_MC, &b_m_kinetic_energy_MC);
   fChain->SetBranchAddress("m_mass_MC", &m_mass_MC, &b_m_mass_MC);
   fChain->SetBranchAddress("m_vx_prod_hit", &m_vx_prod_hit, &b_m_vx_prod_hit);
   fChain->SetBranchAddress("m_vy_prod_hit", &m_vy_prod_hit, &b_m_vy_prod_hit);
   fChain->SetBranchAddress("m_vz_prod_hit", &m_vz_prod_hit, &b_m_vz_prod_hit);
   fChain->SetBranchAddress("m_vx_decay_hit", &m_vx_decay_hit, &b_m_vx_decay_hit);
   fChain->SetBranchAddress("m_vy_decay_hit", &m_vy_decay_hit, &b_m_vy_decay_hit);
   fChain->SetBranchAddress("m_vz_decay_hit", &m_vz_decay_hit, &b_m_vz_decay_hit);
   fChain->SetBranchAddress("m_x_start_global", &m_x_start_global, &b_m_x_start_global);
   fChain->SetBranchAddress("m_y_start_global", &m_y_start_global, &b_m_y_start_global);
   fChain->SetBranchAddress("m_z_start_global", &m_z_start_global, &b_m_z_start_global);
   fChain->SetBranchAddress("m_x_end_global", &m_x_end_global, &b_m_x_end_global);
   fChain->SetBranchAddress("m_y_end_global", &m_y_end_global, &b_m_y_end_global);
   fChain->SetBranchAddress("m_z_end_global", &m_z_end_global, &b_m_z_end_global);
   fChain->SetBranchAddress("m_slope_x_global", &m_slope_x_global, &b_m_slope_x_global);
   fChain->SetBranchAddress("m_slope_y_global", &m_slope_y_global, &b_m_slope_y_global);
   fChain->SetBranchAddress("m_slope_z_global", &m_slope_z_global, &b_m_slope_z_global);
   fChain->SetBranchAddress("m_emeloss", &m_emeloss, &b_m_emeloss);
   fChain->SetBranchAddress("m_plate", &m_plate, &b_m_plate);
   fChain->SetBranchAddress("m_base", &m_base, &b_m_base);
   fChain->SetBranchAddress("m_emmodule", &m_emmodule, &b_m_emmodule);
   fChain->SetBranchAddress("m_film", &m_film, &b_m_film); //0-bottom, 1-top
   Notify();
}

bool m_NuHit_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void m_NuHit_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t m_NuHit_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void m_NuHit_tree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L m_NuHit_tree.C
//      root> m_NuHit_tree t
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
