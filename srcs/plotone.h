//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  7 11:42:29 2020 by ROOT version 6.22/02
// from TTree events/waveform tree
// found on file: datafiles/caen1.root
// Macro from Erin/Michelle 02/20/23
//////////////////////////////////////////////////////////

#ifndef plotone_h
#define plotone_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class plotone {
public :
   TTree          *fChain;   //!pointer to the plotonelyzed TTree or TChain
   TTree          *fChain2;
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           fRun;
   Int_t           fEvent;
   vector<unsigned long> *fTicksVec;
   vector<vector<unsigned short> > *fWvfmsVec;

   Float_t art_ev;
   Float_t caen_ev;
   Float_t caenv_ev_tts;

   // List of branches
   TBranch        *b_fRun;   //!
   TBranch        *b_fEvent;   //!
   TBranch        *b_fTicksVec;   //!
   TBranch        *b_fWvfmsVec;   //!

   TBranch *b_art_ev;
   TBranch *b_caen_ev;
   TBranch *b_caenv_ev_tts;

   plotone(TTree *tree=0, TTree *tree2=0);
   virtual ~plotone();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, TTree *tree2);
   virtual void     Loop(int ichoice);
   virtual void     Average();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef plotone_cxx
plotone::plotone(TTree *tree, TTree *tree2) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("eventana_simdata_run4385.root");//caendump_run4484.root");//caenv1730dump_hist_run1291.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("eventana_simdata_run4385.root");//caendump_run4484.root");//caenv1730dump_hist_run1291.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("eventana_simdata_run4385.root:/dumpall");//caendump_run4484.root:/caenv1730dump");
      dir->GetObject("events",tree);
      dir->GetObject("nt_header",tree2);

   }
   Init(tree, tree2);
}

plotone::~plotone()
{
   if (!fChain || !fChain2) return;
   delete fChain->GetCurrentFile();
   delete fChain2->GetCurrentFile();
}

Int_t plotone::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain || !fChain2) return 0;
   return fChain->GetEntry(entry), fChain2->GetEntry(entry);
}
Long64_t plotone::LoadTree(Long64_t entry)
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

void plotone::Init(TTree *tree, TTree *tree2)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   fTicksVec = 0;
   fWvfmsVec = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fRun", &fRun, &b_fRun);
   fChain->SetBranchAddress("fEvent", &fEvent, &b_fEvent);
   fChain->SetBranchAddress("fTicksVec", &fTicksVec, &b_fTicksVec);
   fChain->SetBranchAddress("fWvfmsVec", &fWvfmsVec, &b_fWvfmsVec);
   Notify();

   // Set object pointer
   art_ev = 0;
   caen_ev = 0;
   caenv_ev_tts = 0;
   // Set branch addresses and branch pointers
   if (!tree2) return;
   fChain2 = tree2;
   fCurrent = -1;
   fChain2->SetMakeClass(1);

   fChain2->SetBranchAddress("art_ev", &art_ev, &b_art_ev);
   fChain2->SetBranchAddress("caen_ev", &caen_ev, &b_caen_ev);
   fChain2->SetBranchAddress("caenv_ev_tts", &caenv_ev_tts, &b_caenv_ev_tts);
   Notify();
}

Bool_t plotone::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void plotone::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t plotone::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef plotone_cxx
