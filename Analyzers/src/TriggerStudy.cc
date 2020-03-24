#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "TriggerStudy.hh"
#include "DownstreamTrack.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "L0PrimitiveHandler.hh"
#include "TLine.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include <fstream>
using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;


Double_t TriggerStudy::RICHlims(Double_t p, Double_t lim){
	Double_t result = Frichval * sqrt(
		2 - (2/nrichval) * sqrt(
			1 + (massofpion2/(p*p))
		)
	)
	+ lim;
	return result;
}

TriggerStudy::TriggerStudy(Core::BaseAnalysis *ba) : Analyzer(ba, "TriggerStudy"){

	RequestL0Data();
	RequestL1Data();
	fTriggerMaskPNN  = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-nQX-UTMC-nMUV-nLKr30");

}

void TriggerStudy::InitOutput(){

}

void TriggerStudy::InitHist(){


	BookHisto("hFullTrigStudy", new TH1I("TriggerStudy", "Kaon_Decay_Characteristics_afterL0", 10, 0, 10));
	BookHisto("hL0PNNChar", new TH1I("L0PNNChar", "L0PNN_Decay_Characteristics", 10, 0, 10));
	BookHisto("hLKREnergy", new TH1D("LKREnergyTest", "Energy_Distro_of_LKR", 30, 0, 70000));
	BookHisto("hPMom", new TH1D("PMomTest", "Momentum_Distro_from_Somewhere", 30, 0, 70000));
	BookHisto("hEOPCalc", new TH1D("EOPTestCalc", "E/p_thing", 30, 0, 1.5));
	
	BookHisto("hTotE", new TH1D("LKREnergyTot", "Energy_Distro_of_LKR", 30, 0, 70000));
	BookHisto("hTotEoP", new TH1D("LKRTotEoP", "E/p_thing", 30, 0, 1.5));

	BookHisto("hLKREoP_base",new TH1D("LKrEoP_cuts", "Histogram of LKr Energy over Spectrometer Momentum (Data)", 500, 0, 1.2));
	BookHisto("hLKREoP_pion",new TH1D("LKrEoP_cuts", "Histogram of LKrEnergy over Spectrometer Momentum (Pion)", 500, 0, 1.2));
	BookHisto("hLKREoP_electron",new TH1D("LKrEoP", "Histogram of LKrEnergy over Spectrometer Momentum (Electron)", 500, 0, 1.2));
	BookHisto("hLKREoP_muon",new TH1D("LKrEoP_cuts", "Histogram of LKrEnergy over Spectrometer Momentum (Muon)", 500, 0, 1.2));

	BookHisto("hTOTEoP_base",new TH1D("TOTEoP", "Histogram of Total Energy over Spectrometer Momentum (Data)", 500, 0, 1.2));
	BookHisto("hTOTEoP_pion",new TH1D("TOTEoP_cuts", "Histogram of TotEnergy over Spectrometer Momentum (Pion)", 500, 0, 1.2));
	BookHisto("hTOTEoP_electron",new TH1D("TOTEoPEoP_cuts", "Histogram of TotEnergy over Spectrometer Momentum (Electron)", 500, 0, 1.2));
	BookHisto("hTOTEoP_muon",new TH1D("TOTEoPEoP_cuts", "Histogram of TotEnergy over Spectrometer Momentum (Muon)", 500, 0, 1.2));

	BookHisto("hEoP_base", new TH2D("EoP_test", "2D Histogram of Total and LKr EoP (Data)", 500, 0, 1.2, 500, 0, 1.2));
	BookHisto("hEoP_electron", new TH2D("EoP_test", "2D Histogram of LKr and Total EoP", 500, 0, 1.2, 500, 0, 1.2));
	BookHisto("hEoP_muon", new TH2D("EoP_test", "2D Histogram of LKr and Total EoP", 500, 0, 1.2, 500, 0, 1.2));

	BookHisto("hRICHring", new TH2D("RichRing", "Radius of Ring vs Particle Momentum (Data)", 500, 14000, 36000, 500, 0, 240));
	BookHisto("hRICHring_exc", new TH2D("RichRing_cuts", "Radius of ring function of particle momentum (Excluded)", 500, 14000, 36000, 500, 0, 240));

	BookHisto("hRICHMissingMass_base", new TH1D("Mass_RICH", "Reconstruction of Mass from RICH (Data)", 500, 0, 0.04));
	BookHisto("hRICHMissingMass_pion", new TH1D("Mass_RICH", "Reconstruction of Mass from RICH", 500, 0, 0.04));
	BookHisto("hRICHMissingMass_electron", new TH1D("Mass_RICH", "Reconstruction of Mass from RICH", 500, 0, 0.04));
	BookHisto("hRICHMissingMass_muon", new TH1D("Mass_RICH_cuts", "Reconstruction of Mass from RICH after selection cuts", 500, 0, 0.04));
	
	///define all useful book counters
	BookCounter("TotalEvent");
	BookCounter("PhysicsEvent");
	BookCounter("L0PNN");

	BookCounter("Kmu2Selection");
	BookCounter("K2piCounter");
	BookCounter("K3piCounter");
	BookCounter("Ke3Selection");
	BookCounter("Kmu3SelectionNoSpectrometer");
	BookCounter("Pi0Selection");
	BookCounter("Main5");
	BookCounter("0TrackSize");
	
	BookCounter("ThreeTrack");
	BookCounter("TwoTrackwLepton");
	BookCounter("notAutopass");

	BookCounter("LKrEoP_counts");
	BookCounter("LKr0Counter");
	BookCounter("LKrnon0Counter");
	BookCounter("LKrMuonExcluded");
	BookCounter("LKrElectronExcluded");
	BookCounter("LKrPionsRemaining");

	BookCounter("TotalEoP_counts");
	BookCounter("TotEoP_counts");
	BookCounter("Tot0Counter");
	BookCounter("Totnon0Counter");
	BookCounter("TotMuonExcluded");
	BookCounter("TotElectronExcluded");
	BookCounter("TotPionsRemaining");

	BookCounter("SeqElectron");
	BookCounter("SeqMuon");
	BookCounter("SeqPion");



	BookCounter("TotalRich_counts");
	BookCounter("Rich_Included");
	BookCounter("Rich_Excluded");

	BookCounter("RecoMass_pion");
	BookCounter("RecoMass_muon");
	BookCounter("RecoMass_electron");


	/// Tables
	NewEventFraction("TriggerAnalysis");
	NewEventFraction("L0PNNAnalysis");
	NewEventFraction("LKrEoPCuts");
	NewEventFraction("TotEoPCuts");
	NewEventFraction("SeqEoPCuts");
	NewEventFraction("RichCuts");
	NewEventFraction("RecoMassCuts");
	

	AddCounterToEventFraction("TriggerAnalysis", "TotalEvent");
 	AddCounterToEventFraction("TriggerAnalysis", "PhysicsEvent");
	AddCounterToEventFraction("TriggerAnalysis", "notAutopass");
 	AddCounterToEventFraction("TriggerAnalysis", "L0PNN");
	
	AddCounterToEventFraction("TriggerAnalysis", "Kmu2Selection");
	AddCounterToEventFraction("TriggerAnalysis", "K2piCounter");
	AddCounterToEventFraction("TriggerAnalysis", "K3piCounter");
	AddCounterToEventFraction("TriggerAnalysis", "Ke3Selection");
	AddCounterToEventFraction("TriggerAnalysis", "Kmu3SelectionNoSpectrometer");
	AddCounterToEventFraction("TriggerAnalysis", "Pi0Selection");
	AddCounterToEventFraction("TriggerAnalysis", "ThreeTrack");
	AddCounterToEventFraction("TriggerAnalysis", "TwoTrackwLepton");
	DefineSampleSizeCounter("TriggerAnalysis", "TotalEvent");



 	AddCounterToEventFraction("L0PNNAnalysis", "L0PNN");
	AddCounterToEventFraction("L0PNNAnalysis", "Main5");
	AddCounterToEventFraction("L0PNNAnalysis", "Kmu2Selection");
	AddCounterToEventFraction("L0PNNAnalysis", "K2piCounter");
	AddCounterToEventFraction("L0PNNAnalysis", "K3piCounter");
	AddCounterToEventFraction("L0PNNAnalysis", "Ke3Selection");
	AddCounterToEventFraction("L0PNNAnalysis", "Kmu3SelectionNoSpectrometer");
	AddCounterToEventFraction("L0PNNAnalysis", "0TrackSize");
	DefineSampleSizeCounter("L0PNNAnalysis", "L0PNN");


	///AddCounterToEventFraction("EoPCuts", "TotalEoP_counts");
	AddCounterToEventFraction("LKrEoPCuts", "LKr0Counter");
	AddCounterToEventFraction("LKrEoPCuts", "LKrnon0Counter");
	AddCounterToEventFraction("LKrEoPCuts", "LKrMuonExcluded");
	AddCounterToEventFraction("LKrEoPCuts", "LKrElectronExcluded");
	AddCounterToEventFraction("LKrEoPCuts", "LKrPionsRemaining");
	DefineSampleSizeCounter("LKrEoPCuts", "LKrnon0Counter");

	AddCounterToEventFraction("TotEoPCuts", "Tot0Counter");
	AddCounterToEventFraction("TotEoPCuts", "Totnon0Counter");
	AddCounterToEventFraction("TotEoPCuts", "TotMuonExcluded");
	AddCounterToEventFraction("TotEoPCuts", "TotElectronExcluded");
	AddCounterToEventFraction("TotEoPCuts", "TotPionsRemaining");
	DefineSampleSizeCounter("TotEoPCuts", "Totnon0Counter");

	AddCounterToEventFraction("SeqEoPCuts", "Tot0Counter");
	AddCounterToEventFraction("SeqEoPCuts", "Totnon0Counter");
	AddCounterToEventFraction("SeqEoPCuts", "SeqElectron");
	AddCounterToEventFraction("SeqEoPCuts", "SeqMuon");
	AddCounterToEventFraction("SeqEoPCuts", "SeqPion");
	DefineSampleSizeCounter("SeqEoPCuts", "Totnon0Counter");

	AddCounterToEventFraction("RichCuts", "TotalRich_counts");
	AddCounterToEventFraction("RichCuts", "Rich_Included");
	AddCounterToEventFraction("RichCuts", "Rich_Excluded");
	DefineSampleSizeCounter("RichCuts", "TotalRich_counts");
	
	AddCounterToEventFraction("RecoMassCuts", "TotalRich_counts");
	AddCounterToEventFraction("RecoMassCuts", "RecoMass_electron");
	AddCounterToEventFraction("RecoMassCuts", "RecoMass_muon");
	AddCounterToEventFraction("RecoMassCuts", "RecoMass_pion");
	DefineSampleSizeCounter("RecoMassCuts", "TotalRich_counts");
	
}

void TriggerStudy::DefineMCSimple(){
}

void TriggerStudy::StartOfRunUser(){
}

void TriggerStudy::StartOfBurstUser(){
}

void TriggerStudy::ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType){
	/// \MemberDescr
	/// \param iEvent : Special event number
	/// \param triggerType : Special trigger type (-1 if not known). For this
	/// variable to be filled, RequestL0SpecialTrigger must be called in the constructor.
	///
	/// Process method for special triggers. Called on each special trigger event after each start of burst.\n
	/// \EndMemberDescr
}

void TriggerStudy::Process(int iEvent){


	IncrementCounter("TotalEvent");
	FillHisto("hFullTrigStudy", 0); ///TotalEvent

	///Retrieve trigger information
 	L0TPData* L0Packet = GetL0Data();
 	L1TPData* L1Packet = GetL1Data();
 
 	Int_t  L0DataType     = L0Packet->GetDataType();
	Bool_t PhysicsData = L0DataType & 0x1; 	
 
	EventHeader* EvtHdr = GetEventHeader(); ///this is important for l1autopass
 	Int_t RunNumber = EvtHdr->GetRunID(); 
 	Bool_t L0TriggerOnPNN    = TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0Packet, fTriggerMaskPNN);
	Bool_t ThreeTrack = *(Bool_t*)GetOutput("FilterThreeTracks.EventSelected");

	Bool_t autopass = TriggerConditions::GetInstance()->L1TriggerAutopass(GetEventHeader()); ///L1Trigger autopass
	if (autopass == false) {
		IncrementCounter("notAutopass");
		return;
	}

	if(ThreeTrack) {
		IncrementCounter("ThreeTrack");
	}
	if(PhysicsData) {
		IncrementCounter("PhysicsEvent");
		FillHisto("hFullTrigStudy", 1); ///Physics Events
	}


	std::vector<DownstreamTrack> Tracks = *GetOutput<std::vector<DownstreamTrack>>("DownstreamTrackBuilder.Output");
	

	///Particle Selections 
		///Kmu2Selection
		NA62Analysis::UserMethods::OutputState state_2mu; // can choose name of variable
		Bool_t Kmu2Selected = *(Bool_t*)GetOutput("Kmu2Selection.EventSelected", state_2mu);
		if(Kmu2Selected) {
			IncrementCounter("Kmu2Selection");
			IncrementCounter("Main5");
			FillHisto("hL0PNNChar", 1);
			FillHisto("hL0PNNChar", 2);
			FillHisto("hFullTrigStudy", 3); ///K3PiCounter
		}

		///K2pi Selection
		NA62Analysis::UserMethods::OutputState state_2pi; // can choose name of variable
		Bool_t K2PiSelected = *(Bool_t*)GetOutput("K2piSelection.EventSelected", state_2pi);
		if(K2PiSelected) {
			IncrementCounter("K2piCounter");
			IncrementCounter("Main5");
			FillHisto("hL0PNNChar", 1);
			FillHisto("hL0PNNChar", 3);
			FillHisto("hFullTrigStudy", 4); ///K3PiCounter
		}

		/// K3PiSelectoin
		NA62Analysis::UserMethods::OutputState state_3pi; // cah choose name of variable
		Bool_t K3PiSelected = *(Bool_t*)GetOutput("K3piSelection.EventSelected", state_3pi);
		
		if(K3PiSelected) {
			IncrementCounter("K3piCounter");
			IncrementCounter("Main5");
			FillHisto("hL0PNNChar", 1);
			FillHisto("hL0PNNChar", 4);
			FillHisto("hFullTrigStudy", 5); ///K3PiCounter
			}

		///Ke3Selection
		NA62Analysis::UserMethods::OutputState state_3e; // can choose name of variable
		Bool_t Ke3Selected = *(Bool_t*)GetOutput("Ke3Selection.EventSelected", state_3e);
		if(Ke3Selected) {
			IncrementCounter("Ke3Selection");
			IncrementCounter("Main5");
			FillHisto("hL0PNNChar", 1);
			FillHisto("hL0PNNChar", 5);
			FillHisto("hFullTrigStudy", 6); ///K3PiCounter
		}
		
		///Kmu3SelectionNoSpectrometer
		NA62Analysis::UserMethods::OutputState state_3mu; // can choose name of variable
		Bool_t Kmu3Selected = *(Bool_t*)GetOutput("Kmu3SelectionNoSpectrometer.EventSelected", state_3mu);
		if(Kmu3Selected) {
			IncrementCounter("Kmu3SelectionNoSpectrometer");
			IncrementCounter("Main5");
			FillHisto("hL0PNNChar", 1);
			FillHisto("hL0PNNChar", 6);
			FillHisto("hFullTrigStudy", 7); ///K3PiCounter
		}

		///Pi0Selection
		NA62Analysis::UserMethods::OutputState state_0pi; // can choose name of variable
		Bool_t K0PiSelected = *(Bool_t*)GetOutput("Pi0Selection.EventSelected", state_0pi);
		if(K0PiSelected) {
			IncrementCounter("Pi0Selection");
			FillHisto("hFullTrigStudy", 8); ///K3PiCounter
		}
		



    if(L0TriggerOnPNN) {
		IncrementCounter("L0PNN"); 
		FillHisto("hFullTrigStudy", 2); ///L0PNN
		FillHisto("hL0PNNChar", 0);

		if (Tracks.size() != 1) {
			IncrementCounter("0TrackSize");
			return;
		}
	

		Double_t Ptrack = Tracks[0].GetMomentum();
		Double_t LKREnergy = Tracks[0].GetLKrEnergy(); 
		Double_t LKREoP = Tracks[0].GetLKrEoP(); ///useful
		Double_t TotEnergy = Tracks[0].GetLKrTotalEnergy();
		///Double_t TotEoP = Tracks[0].GetLKrTotalEoP();

		Double_t RichRing = Tracks[0].GetRICHRingRadius();
		Double_t RichMass = Tracks[0].GetRICHSingleRingTrkCentredMass();
		Double_t RichMass2 = RichMass*RichMass/1000000;

		///Double_t LowEoPLim = 0.05;
		Double_t LowEoPLim = 0.2;
		Double_t HighEoPLim = 0.9;
		Double_t LowMassLim = 0.005;
		Double_t HighMassLim = 0.015;

		if (Ptrack > 35000 || Ptrack < 15000) return;

		///Test Histograms
		FillHisto("hLKREnergy", LKREnergy);
		FillHisto("hPMom", Ptrack);
		FillHisto("hEOPCalc", LKREnergy/Ptrack);
		FillHisto("hTotE", TotEnergy);
		///FillHisto("hTotEoP",TotEoP);


		

		///EoP Cuts
		IncrementCounter("LKrEoP_counts");
		if(LKREoP == 0.0) IncrementCounter("LKr0Counter");
		else{
			IncrementCounter("LKrnon0Counter");
			FillHisto("hLKREoP_base", LKREoP);		
			if (LKREoP < LowEoPLim) {
				IncrementCounter("LKrMuonExcluded");
				FillHisto("hLKREoP_muon", LKREoP);
			}
			else if (LKREoP > HighEoPLim) {
				IncrementCounter("LKrElectronExcluded");
				FillHisto("hLKREoP_electron", LKREoP);
			}
			else{
				IncrementCounter("LKrPionsRemaining");
				FillHisto("hLKREoP_pion", LKREoP);
			}
		}



		///EoP Cuts with MUV
		Double_t LKRE = Tracks[0].GetLKrEnergy();
		Double_t MUV1E = Tracks[0].GetMUV1ClusterEnergy();
		Double_t MUV2E = Tracks[0].GetMUV2ClusterEnergy();
		Double_t TotE = LKRE + MUV1E + MUV2E;
		Double_t TotEoP = TotE / (Tracks[0].GetMomentum());
		
		IncrementCounter("TotalEoP_counts");
		if(TotEoP == 0.0) IncrementCounter("Tot0Counter");
		else{
			IncrementCounter("Totnon0Counter");
			FillHisto("hTOTEoP_base", TotEoP);		
			
			
			if (TotEoP < LowEoPLim) {
				IncrementCounter("TotMuonExcluded");
				FillHisto("hTOTEoP_muon", TotEoP);
			}
			else if (TotEoP > HighEoPLim) {
				IncrementCounter("TotElectronExcluded");
				FillHisto("hTOTEoP_electron", TotEoP);
			}
			else{
				IncrementCounter("TotPionsRemaining");
				FillHisto("hTOTEoP_pion", TotEoP);
			}

			///Sequential Eop cuts 
			///use TotalEoP_Counts and Tot0Counter
			if (LKREoP > HighEoPLim) {
				IncrementCounter("SeqElectron");
				FillHisto("hEoP_electron", LKREoP, TotEoP);
			}

			else if (TotEoP < LowEoPLim && LKREoP != 0.0) {
				IncrementCounter("SeqMuon");
				FillHisto("hEoP_muon", LKREoP, TotEoP);
			}
			else {
				IncrementCounter("SeqPion");
				FillHisto("hEoP_base", LKREoP, TotEoP);
			}
		}

		///Rich Cuts
		IncrementCounter("TotalRich_counts");
		if(RICHlims(Ptrack, lowlim) < RichRing && RichRing < RICHlims(Ptrack, highlim)){
			IncrementCounter("Rich_Included");
			FillHisto("hRICHring", Ptrack, RichRing);
		}
		else {
			IncrementCounter("Rich_Excluded");
			FillHisto("hRICHring_exc", Ptrack, RichRing);
		}
	
		///RecoMass Cuts
		FillHisto("hRICHMissingMass_base", RichMass2);
		if(RichMass2 <= LowMassLim){
			IncrementCounter("RecoMass_electron");
			FillHisto("hRICHMissingMass_electron", RichMass2);
		}
		else if(LowMassLim < RichMass2 && RichMass2 < HighMassLim ){
			IncrementCounter("RecoMass_muon");
			FillHisto("hRICHMissingMass_muon", RichMass2);
			
		}
		else {
			IncrementCounter("RecoMass_pion");
			FillHisto("hRICHMissingMass_pion", RichMass2);
		}

	}

	///labelling histograms
	Int_t ix;
	const Int_t nx1 = 9;
	const Int_t nx2 = 7;
	const char *labels1[nx1] = {"TotalEvent","PhysicsEvent", "L0PNN", "Kmu2Selection","K2piCounter",
      "K3piCounter","Ke3Selection","Kmu3Selection","Pi0Selection"};
	const char *labels2[nx2] = {"L0PNN", "Main5", "Kmu2Selection","K2piCounter",
      "K3piCounter","Ke3Selection","Kmu3Selection"};

	TH1 *MyHisto1 = fHisto.GetHisto("hFullTrigStudy");
	TH1 *MyHisto2 = fHisto.GetHisto("hL0PNNChar"); 

	for (ix=1;ix<=nx1;ix++) MyHisto1->GetXaxis()->SetBinLabel(ix,labels1[ix-1]); 
	for (ix=1;ix<=nx2;ix++) MyHisto2->GetXaxis()->SetBinLabel(ix,labels2[ix-1]); 

}

void TriggerStudy::PostProcess(){
}

void TriggerStudy::EndOfBurstUser(){
}

void TriggerStudy::EndOfRunUser(){

}

void TriggerStudy::EndOfJobUser(){

	///PDF_Files/TriggerStudy
	///Could have used iterator but just need this to work

	TCanvas *c = new TCanvas;
	
	TLegend *leglkreop = new TLegend(0.15, 0.75, 0.35, 0.85);
	TLegend *legtoteop = new TLegend(0.15, 0.75, 0.35, 0.85);
	TLegend *legrich = new TLegend(0.15, 0.75, 0.40, 0.90);
	TLegend *legmass = new TLegend(0.15, 0.75, 0.35, 0.85);
	
	fHisto.GetTH1("hFullTrigStudy")->Draw();
	///fHisto.GetTH1("hFullTrigStudy")->SetXTitle("Characteristics");
	fHisto.GetTH1("hFullTrigStudy")->SetYTitle("Number of Hits");
	fHisto.GetTH1("hFullTrigStudy")->SetStats(false);
	c->SaveAs("PDF_Files/TriggerStudy/hFullTrigStudy.pdf");

	fHisto.GetTH1("hL0PNNChar")->Draw();
	///fHisto.GetTH1("hL0PNNChar")->SetXTitle("Characteristics");
	fHisto.GetTH1("hL0PNNChar")->SetYTitle("Number of Hits");
	fHisto.GetTH1("hL0PNNChar")->SetStats(false);
	c->SaveAs("PDF_Files/TriggerStudy/hL0PNNChar.pdf");
	
	///fHisto.GetTH1("hLKREnergy")->Draw();
	///c->SaveAs("PDF_Files/TriggerStudy/hLKREnergy.pdf");

	///fHisto.GetTH1("hPMom")->Draw();
	///c->SaveAs("PDF_Files/TriggerStudy/hPMom.pdf");

	///fHisto.GetTH1("hEOPCalc")->Draw();
	///c->SaveAs("PDF_Files/TriggerStudy/hEOPCalc.pdf");

		///Combining EoP Graph
	fHisto.GetTH1("hLKREoP_base")->SetXTitle("#frac{E_{LKr}}{p_{spec}}");
	fHisto.GetTH1("hLKREoP_base")->SetYTitle("Number of Hits");
	fHisto.GetTH1("hLKREoP_base")->SetStats(false);

	fHisto.GetTH1("hLKREoP_muon")->SetLineColor(kRed);
	fHisto.GetTH1("hLKREoP_pion")->SetLineColor(kBlue);
	fHisto.GetTH1("hLKREoP_electron")->SetLineColor(kGreen);

	fHisto.GetTH1("hLKREoP_base")->Draw();
	fHisto.GetTH1("hLKREoP_muon")->Draw("same");	
	fHisto.GetTH1("hLKREoP_electron")->Draw("same");
	fHisto.GetTH1("hLKREoP_pion")->Draw("same");
	
	leglkreop -> AddEntry(fHisto.GetTH1("hLKREoP_muon"), "Muon Cuts", "l");
	leglkreop -> AddEntry(fHisto.GetTH1("hLKREoP_electron"), "Positron Cuts", "l");
	leglkreop -> AddEntry(fHisto.GetTH1("hLKREoP_pion"), "Pion Remains", "l");
	leglkreop -> Draw();
	c->SaveAs("PDF_Files/TriggerStudy/LKrEoP_comb.pdf");

	fHisto.GetTH1("hTOTEoP_base")->SetXTitle("#frac{E_{Tot}}{p_{spec}}");
	fHisto.GetTH1("hTOTEoP_base")->SetYTitle("Number of Hits");
	fHisto.GetTH1("hTOTEoP_base")->SetStats(false);

	fHisto.GetTH1("hTOTEoP_muon")->SetLineColor(kRed);
	fHisto.GetTH1("hTOTEoP_pion")->SetLineColor(kBlue);
	fHisto.GetTH1("hTOTEoP_electron")->SetLineColor(kGreen);

	fHisto.GetTH1("hTOTEoP_base")->Draw();
	fHisto.GetTH1("hTOTEoP_muon")->Draw("same");
	fHisto.GetTH1("hTOTEoP_electron")->Draw("same");
	fHisto.GetTH1("hTOTEoP_pion")->Draw("same");
	
	legtoteop -> AddEntry(fHisto.GetTH1("hTOTEoP_muon"), "Muon Cuts", "l");
	legtoteop -> AddEntry(fHisto.GetTH1("hTOTEoP_electron"), "Positron Cuts", "l");
	legtoteop -> AddEntry(fHisto.GetTH1("hTOTEoP_pion"), "Pion Remains", "l");
	legtoteop -> Draw();
	c->SaveAs("PDF_Files/TriggerStudy/TotEoP_comb.pdf");


	///fHisto.GetTH1("hTotE")->Draw();
	///c->SaveAs("PDF_Files/TriggerStudy/hTotE.pdf");

	///fHisto.GetTH1("hTotEoP")->Draw();
	///c->SaveAs("PDF_Files/TriggerStudy/hTotEoP.pdf");

	///Combining RICH Graph
	fHisto.GetTH2("hRICHring_exc")->SetMarkerColor(kRed);
	fHisto.GetTH2("hRICHring")->SetMarkerColor(kBlue);
	///fHisto.GetTH2("hRICHring_exc")->SetMarkerStyle(6);
	///fHisto.GetTH2("hRICHring")->SetMarkerStyle(6);

	fHisto.GetTH2("hRICHring")->SetXTitle("Momentum (MeV)");
	fHisto.GetTH2("hRICHring")->SetYTitle("Radius (mm)");
	fHisto.GetTH2("hRICHring")->SetStats(false);
	fHisto.GetTH2("hRICHring")->Draw();
	fHisto.GetTH2("hRICHring_exc")->Draw("same");
	legrich -> AddEntry(fHisto.GetTH2("hRICHring_exc"), "Excluded Cuts", "p");
	legrich -> AddEntry(fHisto.GetTH2("hRICHring"), "Included Cuts", "p");
	legrich -> Draw();
	c->SaveAs("PDF_Files/TriggerStudy/RICH_comb.pdf");



	///Combining Mass Graph
	fHisto.GetTH1("hRICHMissingMass_base")->SetXTitle("Mass^{2} (GeV^{2})");
	fHisto.GetTH1("hRICHMissingMass_base")->SetYTitle("Number of Hits");

	fHisto.GetTH1("hRICHMissingMass_electron")->SetLineColor(kGreen);
	fHisto.GetTH1("hRICHMissingMass_muon")->SetLineColor(kRed);
	fHisto.GetTH1("hRICHMissingMass_pion")->SetLineColor(kBlue);
	

	fHisto.GetTH1("hRICHMissingMass_base")->Draw();
	fHisto.GetTH1("hRICHMissingMass_pion")->Draw("same");
	fHisto.GetTH1("hRICHMissingMass_muon")->Draw("same");
	fHisto.GetTH1("hRICHMissingMass_electron")->Draw("same");
	legmass -> AddEntry(fHisto.GetTH1("hRICHMissingMass_electron"), "Positron Cuts", "l");
	legmass -> AddEntry(fHisto.GetTH1("hRICHMissingMass_muon"), "Muon Cuts", "l");
	legmass -> AddEntry(fHisto.GetTH1("hRICHMissingMass_pion"), "Pion Cuts", "l");
	legmass -> Draw();
	c->SaveAs("PDF_Files/TriggerStudy/Mass_comb.pdf");

	
	fHisto.GetTH2("hEoP_base")->SetYTitle("#frac{E_{Tot}}{p_{spec}}");
	fHisto.GetTH2("hEoP_base")->SetXTitle("#frac{E_{LKr}}{p_{spec}}");
	fHisto.GetTH2("hEoP_base")->SetMarkerColor(kBlue);
	fHisto.GetTH2("hEoP_electron")->SetMarkerColor(kGreen);
	fHisto.GetTH2("hEoP_muon")->SetMarkerColor(kRed);
	fHisto.GetTH2("hEoP_base")->SetStats(false);
	
	fHisto.GetTH2("hEoP_base")->Draw();
	fHisto.GetTH2("hEoP_electron")->Draw("same");
	fHisto.GetTH2("hEoP_muon")->Draw("same");
	legtoteop -> Draw();
	c->SaveAs("PDF_Files/TriggerStudy/EoP_comb.pdf");

	delete c;
	SaveAllPlots();
}

void TriggerStudy::DrawPlot(){


}

TriggerStudy::~TriggerStudy(){

}
