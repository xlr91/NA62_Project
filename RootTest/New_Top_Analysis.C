#define New_Top_Analysis_cxx
// The class definition in New_Top_Analysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Top_Analysis.C")
// root> T->Process("Top_Analysis.C","some options")
// root> T->Process("Top_Analysis.C+")
//


#include "New_Top_Analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <math.h>

TH1F* TopMass;

int success_Counter = 0;

float pi = 3.1415927;

void New_Top_Analysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   TopMass = new TH1F("TopMass", "Calculated Mass of Top Pair", 200, 0., 900000.);
}

void New_Top_Analysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();  
}

Bool_t New_Top_Analysis::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);

   int nE1 = el_e.GetSize();
   int nM1 = mu_e.GetSize();
   int nJ1 = jet_e.GetSize();
   int Bcounter = 0;
   int Ecounter = 0;
   int Mcounter = 0;

   float totalEnergy = 0.;

   bool failure = false;

   // First, we confirm that there are two b-jets in our event
   for (int i = 0; i < nJ1; ++i) {
	   bool isB = (jet_isbtagged_77[i] == 1);
	   if (isB) {
		   ++Bcounter;
		   totalEnergy += jet_e[i];
	   }
   }
   if (Bcounter != 2) {
	   failure = true;
   }

   // If so, we attempt to find a lepton with an appropriate energy.
   if (failure == false) {
	   int eLocation;

	   for (int i = 0; i < nE1; ++i) {
		   // We will require our electron to have an energy of 40 GeV. Obviously the energy is
		   //distributed randomly between it and the neutrino so it will sometimes get less than half, however
		   // there is a lot of free energy in the decay, so "half" will be more than 40 Gev.
		   if (el_e[i] > 40000.) {
			   eLocation = i;
			   ++Ecounter;
		   }

	   }
	   if (Ecounter == 1) {
		   totalEnergy += (el_e[eLocation] + (*met_met));
	   }
	   // If we did not get exactly one electron, we will look for a muon.
	   else {
		   int mLocation;

		   for (int i = 0; i < nM1; ++i) {
			   if (mu_e[i] > 40000.) {
				   mLocation = i;
				   ++Mcounter;
			   }
		   }
		   if (Mcounter == 1) {
			   totalEnergy += (mu_e[mLocation] + (*met_met));
		   }
		   else {
			   failure = true;
		   }
	   }
   }
   //Now we will look for a pair of jets. We will select the pair of jets with an invariant mass closest
   //to that of the W boson, which is 80.4 GeV
   if (failure == false) {
	   float bestError = 99999999999999999;
	   float sumE = 0;
	   float error;


	   for (int i = 0; i < (nJ1 - 1); i++) {
		   for (int j = (i + 1); j < nJ1; j++) {
			   float E1 = jet_e[i];
			   float E2 = jet_e[j];
			   float PT1 = jet_pt[i];
			   float PT2 = jet_pt[j];
			   float phi1 = jet_phi[i];
			   float phi2 = jet_phi[j];
			   float eta1 = jet_eta[i];
			   float eta2 = jet_eta[j];

			   // We are adding the momentum as a a vector, so we must split it 
			   //into its components, and add them individually.
			   float sumX = PT1 * cos(phi1) + PT2 * cos(phi2);
			   float sumY = PT1 * sin(phi1) + PT2 * sin(phi2);
			   float sumZ = PT1 * sinh(eta1) + PT2 * sinh(eta2);

			   //We don't bother square rooting as it would just get squared in the next line anyway.
			   float magSumP = pow(sumX, 2) + pow(sumY, 2) + pow(sumZ, 2);
			   
			   float M = sqrt(pow((E1 + E2), 2) - magSumP);
			   // We use the squared error to ensure it is positive
			   error = pow((M - 80385), 2);

				   if (error < bestError) {
					   sumE = E1 + E2;
			   }
		   }
	   }
	   //We want our invariant mass to be within 10GeV of the W mass
	   if (error < 100000000 && sumE != 0) {
		   totalEnergy += sumE;
		}
	   else {
		   failure = true;
	   }
   }

   // If none of our other sub-functions have produced a failure, we will add the total energy to the histogram.
   if (failure == false) {
	   TopMass->Fill(totalEnergy);
	   ++success_Counter;
   }

   return kTRUE;
}

void New_Top_Analysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void New_Top_Analysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
	std::cout << "Number of candidates = " << success_Counter << std::endl;
}
