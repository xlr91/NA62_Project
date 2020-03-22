#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include <TCanvas.h>
#include "MCAnalyzer.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "TLegend.h"
#include "TLatex.h"
using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;


Double_t MCAnalyzer::RICHlims(Double_t p, Double_t lim){
	Double_t result = Frichval * sqrt(
		2 - (2/nrichval) * sqrt(
			1 + (massofpion2/(p*p))
		)
	)
	+ lim;
	return result;
}

MCAnalyzer::MCAnalyzer(Core::BaseAnalysis *ba) : Analyzer(ba, "MCAnalyzer"){

	RequestL0Data();
	RequestL0SpecialTrigger();
	fTriggerMaskPNN  = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-nQX-UTMC-nMUV-nLKr30");
	fPrimitiveHandler = L0PrimitiveHandler::GetInstance();
	
	fPrimitiveHandler->DeclareL0Emulators(fParent, kL0RICH, kL0NewCHOD, kL0MUV3, kL0Calo);


}

void MCAnalyzer::InitOutput(){

}

void MCAnalyzer::InitHist(){



	BookHisto("hLKREoP_base",new TH1D("LKrEoP", "Histogram of LKr Energy over Spectrometer Momentum (MC)", 500, 0, 1.2));
	BookHisto("hLKREoP_pion",new TH1D("LKrEoP_cuts", "Histogram of LKrEnergy over Spectrometer Momentum (Pion)", 500, 0, 1.2));
	BookHisto("hLKREoP_electron",new TH1D("LKrEoP_cuts", "Histogram of LKrEnergy over Spectrometer Momentum (Electron)", 500, 0, 1.2));
	BookHisto("hLKREoP_muon",new TH1D("LKrEoP_cuts", "Histogram of LKrEnergy over Spectrometer Momentum (Muon)", 500, 0, 1.2));

	BookHisto("hTOTEoP_base",new TH1D("TOTEoP", "Histogram of Total Energy over Spectrometer Momentum (MC)", 500, 0, 1.2));
	BookHisto("hTOTEoP_pion",new TH1D("TOTEoP_cuts", "Histogram of TotEnergy over Spectrometer Momentum (Pion)", 500, 0, 1.2));
	BookHisto("hTOTEoP_electron",new TH1D("TOTEoPEoP_cuts", "Histogram of TotEnergy over Spectrometer Momentum (Electron)", 500, 0, 1.2));
	BookHisto("hTOTEoP_muon",new TH1D("TOTEoPEoP_cuts", "Histogram of TotEnergy over Spectrometer Momentum (Muon)", 500, 0, 1.2));

	BookHisto("hEoP_base", new TH2D("EoP_test", "2D Histogram of Total and LKr EoP (MC)", 500, 0, 1.2, 500, 0, 1.2));
	BookHisto("hEoP_electron", new TH2D("EoP_test", "2D Histogram of LKr and Total EoP", 500, 0, 1.2, 500, 0, 1.2));
	BookHisto("hEoP_muon", new TH2D("EoP_test", "2D Histogram of LKr and Total EoP", 500, 0, 1.2, 500, 0, 1.2));

	BookHisto("hRICHring", new TH2D("RichRing", "Radius of Ring vs Particle Momentum (MC)", 300, 14000, 36000, 300, 0, 240));
	BookHisto("hRICHring_exc", new TH2D("RichRing_cuts", "Radius of ring function of particle momentum (Excluded)", 300, 14000, 36000, 300, 0, 240));

	BookHisto("hRICHMissingMass_base", new TH1D("Mass_RICH", "Reconstruction of Mass from RICH (MC)", 500, 0, 0.04));
	BookHisto("hRICHMissingMass_pion", new TH1D("Mass_RICH", "Reconstruction of Mass from RICH", 500, 0, 0.04));
	BookHisto("hRICHMissingMass_electron", new TH1D("Mass_RICH", "Reconstruction of Mass from RICH", 500, 0, 0.04));
	BookHisto("hRICHMissingMass_muon", new TH1D("Mass_RICH_cuts", "Reconstruction of Mass from RICH after selection cuts", 500, 0, 0.04));


	BookCounter("TotalEvents");
	BookCounter("PhysicsEvents");
	BookCounter("PassedL0Trigger");
	BookCounter("K3PiSelected"); ///change for others
	BookCounter("QX_ok");
	BookCounter("MUV_ok");
	BookCounter("UTMC_ok");
	BookCounter("RICH_ok");
	BookCounter("LKr30_ok");
	BookCounter("L0PNN_ok");
	BookCounter("0TrackSize");


	BookCounter("Kmu2Selection");
	BookCounter("K2piCounter");
	BookCounter("K3piCounter");
	BookCounter("Ke3Selection");
	BookCounter("Kmu3SelectionNoSpectrometer");
	BookCounter("Pi0Selection");


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

	TString name1 = "L0_Performance";
	TString name2 = "Data_Selections";

	NewEventFraction(name1);
	NewEventFraction(name2);
	AddCounterToEventFraction(name1, "TotalEvents");
	AddCounterToEventFraction(name1, "PhysicsEvents");
	AddCounterToEventFraction(name1, "PassedL0Trigger");
	AddCounterToEventFraction(name1, "QX_ok");
	AddCounterToEventFraction(name1, "MUV_ok");
	AddCounterToEventFraction(name1, "UTMC_ok");
	AddCounterToEventFraction(name1, "RICH_ok");
	AddCounterToEventFraction(name1, "LKr30_ok");
	AddCounterToEventFraction(name1, "L0PNN_ok");
	AddCounterToEventFraction(name1, "0TrackSize");
	DefineSampleSizeCounter(name1, "TotalEvents");

	AddCounterToEventFraction(name2, "L0PNN_ok");
	AddCounterToEventFraction(name2, "Kmu2Selection");
	AddCounterToEventFraction(name2, "K2piCounter");
	AddCounterToEventFraction(name2, "K3piCounter");
	AddCounterToEventFraction(name2, "Ke3Selection");
	AddCounterToEventFraction(name2, "Kmu3SelectionNoSpectrometer");
	DefineSampleSizeCounter(name2, "L0PNN_ok");


	NewEventFraction("LKrEoPCuts");
	NewEventFraction("TotEoPCuts");
	NewEventFraction("SeqEoPCuts");
	NewEventFraction("RichCuts");
	NewEventFraction("RecoMassCuts");
	

	///AddCounterToEventFraction("LKrEoPCuts", "TotalEoP_counts");
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


void MCAnalyzer::DefineMCSimple(){

}

void MCAnalyzer::StartOfRunUser(){


  myfilep.open ("MCARingp.txt");
  myfiler.open ("MCARingr.txt");
  myfilepr.open ("MCARingpr.txt");


}

void MCAnalyzer::StartOfBurstUser(){

}

void MCAnalyzer::ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType){

}

void MCAnalyzer::Process(int iEvent){

	
	if(fMCSimple.fStatus == MCSimple::kMissing){printIncompleteMCWarning(iEvent);return;}
	if(fMCSimple.fStatus == MCSimple::kEmpty){printNoMCWarning();return;}
	IncrementCounter("TotalEvents");
	
	///Initialize variables
	L0TPData* L0Packet = GetL0Data();
	Int_t  L0DataType     = L0Packet->GetDataType();
	EventHeader* EvtHdr = GetEventHeader();
	Int_t RunNumber = EvtHdr->GetRunID(); 
	Bool_t PhysicsData = L0DataType & 0x1; 
	Bool_t L0TriggerOnPNN    = TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0Packet, fTriggerMaskPNN); ///doesntwork :()


	///Performs L0 primitive check for MC
	fPrimitiveHandler -> SetData(GetL0Data(), GetRunID());
	Int_t RichTime = fPrimitiveHandler->GetTriggerTime(kL0RICH);
	Bool_t QX_ok_emu = fPrimitiveHandler -> CheckEmulatedPrimitives("QX", RichTime);
	Bool_t MUV_ok_emu = fPrimitiveHandler -> CheckEmulatedPrimitives("MUV", RichTime);
	Bool_t UTMC_ok_emu = fPrimitiveHandler -> CheckEmulatedPrimitives("UTMC", RichTime);
	Bool_t RICH_ok_emu = fPrimitiveHandler -> CheckEmulatedPrimitives("RICH", RichTime);
	///Bool_t LKr30_ok_emu = fPrimitiveHandler -> CheckEmulatedPrimitives("LKr", RichTime);
	Bool_t LKr30_ok_emu = fPrimitiveHandler -> CheckEmulatedPrimitives("E20", RichTime);
	Bool_t PNN_ok_emu = RICH_ok_emu && !QX_ok_emu && UTMC_ok_emu && !MUV_ok_emu && !LKr30_ok_emu;

	if(!QX_ok_emu) IncrementCounter("QX_ok");
	if(!MUV_ok_emu) IncrementCounter("MUV_ok");
	if(UTMC_ok_emu) IncrementCounter("UTMC_ok");
	if(RICH_ok_emu) IncrementCounter("RICH_ok");
	if(!LKr30_ok_emu) IncrementCounter("LKr30_ok");
	if(PhysicsData) IncrementCounter("PhysicsEvents");

	if(L0TriggerOnPNN) IncrementCounter("PassedL0Trigger");

	///for testing purposes
	///PNN_ok_emu = true;

	///Runs
	if(PNN_ok_emu) {
		IncrementCounter("L0PNN_ok");
		///Particle Selections 
			///Kmu2Selection
			NA62Analysis::UserMethods::OutputState state_2mu; // can choose name of variable
			Bool_t Kmu2Selected = *(Bool_t*)GetOutput("Kmu2Selection.EventSelected", state_2mu);
			if(Kmu2Selected) {
				IncrementCounter("Kmu2Selection");
			}

			///K2pi Selection
			NA62Analysis::UserMethods::OutputState state_2pi; // can choose name of variable
			Bool_t K2PiSelected = *(Bool_t*)GetOutput("K2piSelection.EventSelected", state_2pi);
			if(K2PiSelected) {
				IncrementCounter("K2piCounter");
			}

			/// K3PiSelectoin
			NA62Analysis::UserMethods::OutputState state_3pi; // cah choose name of variable
			Bool_t K3PiSelected = *(Bool_t*)GetOutput("K3piSelection.EventSelected", state_3pi);
			
			if(K3PiSelected) {
				IncrementCounter("K3piCounter");
				}

			///Ke3Selection
			NA62Analysis::UserMethods::OutputState state_3e; // can choose name of variable
			Bool_t Ke3Selected = *(Bool_t*)GetOutput("Ke3Selection.EventSelected", state_3e);
			if(Ke3Selected) {
				IncrementCounter("Ke3Selection");
			}
			
			///Kmu3SelectionNoSpectrometer
			NA62Analysis::UserMethods::OutputState state_3mu; // can choose name of variable
			Bool_t Kmu3Selected = *(Bool_t*)GetOutput("Kmu3SelectionNoSpectrometer.EventSelected", state_3mu);
			if(Kmu3Selected) {
				IncrementCounter("Kmu3SelectionNoSpectrometer");
			}

		std::vector<DownstreamTrack> Tracks = *GetOutput<std::vector<DownstreamTrack>>("DownstreamTrackBuilder.Output");
		if (Tracks.size() != 1) {
			IncrementCounter("0TrackSize");
			return;
			
		}
		

		Double_t Ptrack = Tracks[0].GetMomentum();
		Double_t LKREoP = Tracks[0].GetLKrEoP(); 
		Double_t RichRing = Tracks[0].GetRICHRingRadius();
		Double_t RichMass = Tracks[0].GetRICHSingleRingTrkCentredMass();
		Double_t RichMass2 = RichMass*RichMass/1000000;

		Double_t LowEoPLim = 0.05;
		Double_t HighEoPLim = 0.9;

		///Double_t LowestMassLim = 0.005;
		Double_t LowMassLim = 0.005;
		Double_t HighMassLim = 0.015;

		if (Ptrack > 35000 || Ptrack < 15000) return;
		///if (LKREoP > 0.1) return;
		///if(LKREoP == 0.0) return;
		
		///cout << Ptrack << RichRing << endl;
		myfilep << Ptrack << endl;
		myfiler << RichRing << endl;
		myfilepr << Ptrack << " " << RichRing << endl;

		
		
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

			else if (TotEoP < LowEoPLim) {
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
}

void MCAnalyzer::PostProcess(){


}

void MCAnalyzer::EndOfBurstUser(){

}

void MCAnalyzer::EndOfRunUser(){


}

void MCAnalyzer::EndOfJobUser(){

	TCanvas *c = new TCanvas;
	TLegend *leglkreop = new TLegend(0.15, 0.75, 0.35, 0.85);
	TLegend *legtoteop = new TLegend(0.15, 0.75, 0.35, 0.85);
	TLegend *legrich = new TLegend(0.15, 0.75, 0.40, 0.90);
	TLegend *legmass = new TLegend(0.15, 0.75, 0.35, 0.85);


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
	c->SaveAs("PDF_Files/MC_simulations/MCALKrEoP_comb.pdf");

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
	c->SaveAs("PDF_Files/MC_simulations/MCATotEoP_comb.pdf");


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
	c->SaveAs("PDF_Files/MC_simulations/MCARICH_comb.pdf");

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
	c->SaveAs("PDF_Files/MC_simulations/MCAMass_comb.pdf");

	
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
	c->SaveAs("PDF_Files/MC_simulations/MCAEoP_comb.pdf");


	///fHisto.GetTH2("hEoP_base")->SetOption("LEGO");
	
	
	SaveAllPlots();
	myfilep.close();
	myfiler.close();
	delete c;
}

void MCAnalyzer::DrawPlot(){

}

MCAnalyzer::~MCAnalyzer(){

}
