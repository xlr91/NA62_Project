//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 18 11:42:43 2019 by ROOT version 6.06/08
// from TTree nominal/tree
// found on file: /usr/local/share/physics/labs/y2/projects/MiriamWatson/zee.root
//////////////////////////////////////////////////////////

#ifndef Example_h
#define Example_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>



class Example : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<ULong64_t> eventNumber = {fReader, "eventNumber"};
   TTreeReaderValue<UInt_t> runNumber = {fReader, "runNumber"};
   TTreeReaderArray<float> el_pt = {fReader, "el_pt"};
   TTreeReaderArray<float> el_eta = {fReader, "el_eta"};
   TTreeReaderArray<float> el_cl_eta = {fReader, "el_cl_eta"};
   TTreeReaderArray<float> el_phi = {fReader, "el_phi"};
   TTreeReaderArray<float> el_e = {fReader, "el_e"};
   TTreeReaderArray<float> el_charge = {fReader, "el_charge"};
   TTreeReaderArray<float> mu_pt = {fReader, "mu_pt"};
   TTreeReaderArray<float> mu_eta = {fReader, "mu_eta"};
   TTreeReaderArray<float> mu_phi = {fReader, "mu_phi"};
   TTreeReaderArray<float> mu_e = {fReader, "mu_e"};
   TTreeReaderArray<float> jet_pt = {fReader, "jet_pt"};
   TTreeReaderArray<float> jet_eta = {fReader, "jet_eta"};
   TTreeReaderArray<float> jet_phi = {fReader, "jet_phi"};
   TTreeReaderArray<float> jet_e = {fReader, "jet_e"};
   TTreeReaderArray<char> jet_isbtagged_77 = {fReader, "jet_isbtagged_77"};
   TTreeReaderValue<Float_t> met_met = {fReader, "met_met"};
   TTreeReaderValue<Float_t> met_phi = {fReader, "met_phi"};


   Example(TTree * /*tree*/ =0) { }
   virtual ~Example() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Example,2);

};

#endif

#ifdef Example_cxx
void Example::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t Example::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef Example_cxx
