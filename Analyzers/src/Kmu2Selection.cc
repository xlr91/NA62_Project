// ---------------------------------------------------------------
//
// History:
//
// Created by Evgueni Goudzovski (eg@hep.ph.bham.ac.uk) 2016-03-25
//
// ---------------------------------------------------------------

/// \class Kmu2Selection
/// \Brief
/// Kmu2 decay selection
/// \EndBrief
/// \Detailed
/// A simple Kmu2 decay selection, currently works only at low intensity due to the absence of the
/// timing conditions in the LKr and CHANTI veto conditions.
/// Outputs the basic moniroting plots: missing mass, momentum, vertex position, etc.
/// Only "physics" L0 data type is considered.
/// The L0 trigger mask to be used can be passed from the command line (it is 0xFF by default).
/// For example, to select events with L0 trigger bit 5 up, one can call
/// \code
/// ./MyApplication ... -p "Kmu2Selection:TriggerMask=0x20"
/// \endcode
/// or equivalently, using decimal notation,
/// \code
/// ./MyApplication ... -p "Kmu2Selection:TriggerMask=32"
/// \endcode
/// \author Evgueni Goudzovski (eg@hep.ph.bham.ac.uk)
/// \EndDetailed

#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "Kmu2Selection.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "DownstreamTrack.hh"
#include "BeamParameters.hh"
#include "LAVMatching.hh"
#include "SAVMatching.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

Kmu2Selection::Kmu2Selection(Core::BaseAnalysis *ba) : Analyzer(ba, "Kmu2Selection") {
  RequestTree("Cedar",  new TRecoCedarEvent,  "Reco");
  RequestTree("CHANTI", new TRecoCHANTIEvent, "Reco");
  RequestTree("CHOD",   new TRecoCHODEvent,   "Reco");
  RequestTree("RICH",   new TRecoRICHEvent,   "Reco");
  RequestTree("LKr",    new TRecoLKrEvent,    "Reco");
  RequestTree("LAV",    new TRecoLAVEvent,    "Reco");
  RequestTree("IRC",    new TRecoIRCEvent,    "Reco");
  RequestTree("SAC",    new TRecoSACEvent,    "Reco");
  RequestL0Data();

  fReadingData = kTRUE;

  fHPhysicsEventsPerBurst = 0;

  AddParam("TriggerMask", &fTriggerMask, 0xFF);
  AddParam("MaxNBursts",  &fMaxNBursts,  5000); // max number of bins in histograms
}

Kmu2Selection::~Kmu2Selection() {}

void Kmu2Selection::InitHist() {
  fReadingData = GetIsTree();

  if (fReadingData) {
    BookHisto("hNTracks",
	      new TH1F("NTracks", "Number of tracks", 11, -0.5, 10.5));
    BookHisto("hEoP",
	      new TH1F("EoP", "Track E/p; E/p", 150, 0.0, 1.5));
    BookHisto("hZvtx",
	      new TH1F("Zvtx", "Z of track-beam axis vertex;Vertex z [m]",
		       200, 50, 250));
    BookHisto("hCDA",
	      new TH1F("CDA", "CDA of the track-beam axis vertex;CDA [mm]",
		       200, 0, 200));
    BookHisto("hZvtxCDA",
	      new TH2F("ZvtxCDA", "CDA vs Z of the track-beam axis vertex;Vertex z [m];CDA [mm]",
		       75, 50, 200, 100, 0, 200));
    BookHisto("hdTTrackCedar",
	      new TH1F("dTTrackCedar", "dTTrackCedar;#Deltat [ns]", 100, -25, 25));
    BookHisto("hMMiss2Mu",
	      new TH1F("MMiss2Mu",
		       "Squared missing mass in muon hypothesis;M_{miss}^{2}(#mu) [GeV^{2}/c^{4}]",
		       300, -0.15, 0.15));
    BookHisto("hPMMiss2Mu",
	      new TH2F("PMMiss2Mu",
		       "Squared missing mass in muon hypothesis vs momentum; Track momentum [GeV/c];M_{miss}^{2}(#mu) [GeV^{2}/c^{4}]",
		       160, 0, 80, 100, -0.15, 0.15));
    BookHisto("hPTheta",
	      new TH2F("PTheta",
		       "Track opening angle wrt beam axis vs momentum;Track momentum [GeV/c];#theta",
		       160, 0, 80, 100, 0.0, 0.02));
    BookHisto("hTrackXYMUV3",
	      new TH2F("TrackXYMUV3", "Track (x,y) at MUV3 plane;x [mm];y [mm]",
		       65, -1300, 1300, 65, -1300, 1300));

    BookHisto("hPhysicsEventsPerBurst", new
	      TH1F("PhysicsEventsPerBurst", "Physics events per burst;Burst ID",
		   fMaxNBursts, -0.5, fMaxNBursts-0.5));
    BookHisto("hKmu2EventsPerBurst", new
	      TH1F("Kmu2EventsPerBurst", "Kmu2 candidates per burst;Burst ID",
		   fMaxNBursts, -0.5, fMaxNBursts-0.5));

    BookHisto("hNCHODHits", new
	      TH1F("NCHODHits", "Number of CHOD hits;Number of hits", 50, -0.5, 49.5));
    BookHisto("hNRICHHits", new
	      TH1F("NRICHHits", "Number of RICH hits;Number of hits", 100, -0.5, 99.5));
    BookHisto("hNLKrCells", new
	      TH1F("NLKrCells", "Number of LKr cells with any signal;Number of cells", 100, -0.5, 99.5));
    BookHisto("hNLKrGoodCells", new
	      TH1F("NLKrGoodCells", "Number of LKr cells with E>40MeV;Number of cells", 100, -0.5, 99.5));
    BookHisto("hNLKrGoodCellsIn", new
	      TH1F("NLKrGoodCellsIn", "Number of LKr cells with E>40MeV;Number of cells", 100, -0.5, 99.5));
    BookHisto("hNLKrGoodCellsOut", new
	      TH1F("NLKrGoodCellsOut", "Number of LKr cells with E>40MeV;Number of cells", 100, -0.5, 99.5));
    BookHisto("hCellTrackDistance", new
	      TH1F("CellTrackDistance", "Good cell - track distance;Distance [mm]", 200, 0, 2000));
    BookHisto("hNLKrClusters", new
	      TH1F("NLKrClusters", "Number of LKr clusters;Number of clusters", 10, -0.5, 9.5));
    BookHisto("hLKrClusterEnergy", new
	      TH1F("LKrClusterEnergy", "LKr cluster energy;Energy [GeV]", 100, 0, 50));
    BookHisto("hLKrCellTotalEnergy", new
	      TH1F("LKrCellTotalEnergy",
		   "LKr total cell (E>40MeV) energy;Total cell energy [GeV]", 70, 0, 70));
    BookHisto("hLKrCellTotalEnergyIn", new
	      TH1F("LKrCellTotalEnergyIn",
		   "LKr total cell (E>40MeV) energy near track;Total cell energy [GeV]", 70, 0, 70));
    BookHisto("hLKrCellTotalEnergyOut", new
	      TH1F("LKrCellTotalEnergyOut",
		   "LKr total cell (E>40MeV) energy far from track;Total cell energy [GeV]", 70, 0, 70));
    BookHisto("hLKrCellTotalEnergyInOut", new
	      TH2F("LKrCellTotalEnergyInOut",
		   "LKr total cell energies;Energy near track [GeV];Energy far from track [GeV]", 100, 0, 50, 100, 0, 50));

    BookHisto("hLKrCellClusterTotalEnergy", new
	      TH2F("LKrCellClusterTotalEnergy",
		   "LKr total cluster energy vs cell energy;Total cell (>40MeV) energy [GeV];Total cluster energy [GeV]",
		   70, 0, 70, 70, 0, 70));

    ////////////////////////////
    // Set up the online monitor

    CreateCanvas("Kmu2Canvas");
    PlacePlotOnCanvas("hNTracks",      "Kmu2Canvas");
    PlacePlotOnCanvas("hEoP",          "Kmu2Canvas");
    PlacePlotOnCanvas("hZvtx",         "Kmu2Canvas");
    PlacePlotOnCanvas("hCDA",          "Kmu2Canvas");
    PlacePlotOnCanvas("hdTTrackCedar", "Kmu2Canvas");
    PlacePlotOnCanvas("hMMiss2Mu",     "Kmu2Canvas");
    SetUpdateInterval(50000);
  }

  else {
    cout << user_normal() << "Reading my own output" << endl;
    fHPhysicsEventsPerBurst =
      static_cast<TH1F*>(RequestHistogram(fAnalyzerName, "PhysicsEventsPerBurst", true));
    fHKmu2EventsPerBurst =
      static_cast<TH1F*>(RequestHistogram(fAnalyzerName, "Kmu2EventsPerBurst", true));
    fHMass =
      static_cast<TH1F*>(RequestHistogram(fAnalyzerName, "MMiss2Mu", true));
    fHEOP =
      static_cast<TH1F*>(RequestHistogram(fAnalyzerName, "EoP", true));
    fHZvertex =
      static_cast<TH1F*>(RequestHistogram(fAnalyzerName, "Zvtx", true));
  }
}

void Kmu2Selection::InitOutput() {
  RegisterOutput("EventSelected", &fEventSelected);
  RegisterOutput("MuonMomentum",  &fMuonMomentum);
}

void Kmu2Selection::Process(Int_t) {

  SetOutputState("EventSelected", kOValid);
  fEventSelected = false;
  fMuonMomentum = 0.0;
  
  if (!fReadingData) return; // no action if reading its own output in --histo mode

  Int_t  L0DataType    = GetL0Data()->GetDataType();
  Int_t  L0TriggerWord = GetL0Data()->GetTriggerFlags();
  Bool_t PhysicsData   = L0DataType    & 0x1;
  Bool_t TriggerOK     = L0TriggerWord & fTriggerMask;

  Bool_t FillPlots = PhysicsData && TriggerOK;
  Int_t BurstID   = GetBurstID();

  if (FillPlots) FillHisto("hPhysicsEventsPerBurst", BurstID);

  TRecoCedarEvent*  CEDARevent  = GetEvent<TRecoCedarEvent>();;
  TRecoCHANTIEvent* CHANTIevent = GetEvent<TRecoCHANTIEvent>();;
  TRecoCHODEvent*   CHODevent   = GetEvent<TRecoCHODEvent>();;
  TRecoRICHEvent*   RICHevent   = GetEvent<TRecoRICHEvent>();;
  TRecoLKrEvent*    LKRevent    = GetEvent<TRecoLKrEvent>();;
  TRecoLAVEvent*    LAVevent    = GetEvent<TRecoLAVEvent>();;
  TRecoIRCEvent*    IRCevent    = GetEvent<TRecoIRCEvent>();;
  TRecoSACEvent*    SACevent    = GetEvent<TRecoSACEvent>();;
  
  ////////////////////////////////////////////////
  // Require a good Cedar candidates (>=5 sectors)

  Int_t NCEDARcand = CEDARevent->GetNCandidates();
  Int_t NGoodCedarCand = 0;
  for (Int_t i=0; i<NCEDARcand; i++) {
    TRecoCedarCandidate* Ccand = static_cast<TRecoCedarCandidate*>(CEDARevent->GetCandidate(i));
    if (Ccand->GetNSectors()>4) NGoodCedarCand++;
  }
  if (!NGoodCedarCand) return;

  //////////////////////////////////////////
  // Require exactly one track in acceptance

  std::vector<DownstreamTrack> Tracks =
    *GetOutput<std::vector<DownstreamTrack>>("DownstreamTrackBuilder.Output");
  if (FillPlots) FillHisto("hNTracks", Tracks.size());
  if (Tracks.size()!=1) return;

  TRecoSpectrometerCandidate* Scand = Tracks[0].GetSpectrometerCandidate();
  Int_t    Q            = Tracks[0].GetCharge();
  Double_t Ptrack       = Tracks[0].GetMomentum(); // spectrometer calibration included
  Double_t Ptrackbefore = Tracks[0].GetMomentumBeforeFit();
  Double_t Ttrack       =
    (Tracks[0].CHODAssociationExists()) ? Tracks[0].GetCHODTime() : Tracks[0].GetTrackTime();
  Double_t Chi2track    = Tracks[0].GetChi2();

  if (Q!=1) return;
  if (Chi2track>20.0) return;
  if (Scand->GetNChambers()!=4) return;
  if (fabs(Ptrack-Ptrackbefore)>20000.0) return; // 20 GeV
  if (Ptrack<5000 || Ptrack>70000) return;

  if (!GeometricAcceptance::GetInstance()->InAcceptance(Scand, kSpectrometer, 0)) return;
  if (!GeometricAcceptance::GetInstance()->InAcceptance(Scand, kSpectrometer, 1)) return;
  if (!GeometricAcceptance::GetInstance()->InAcceptance(Scand, kSpectrometer, 2)) return;
  if (!GeometricAcceptance::GetInstance()->InAcceptance(Scand, kSpectrometer, 3)) return;
  if (!GeometricAcceptance::GetInstance()->InAcceptance(Scand, kMUV3))            return;

  Double_t eop = Tracks[0].GetLKrEoP();
  if (FillPlots) FillHisto("hEoP", eop);

  /////////////////////////////////////////
  // Zvertex & CDA: track wrt the beam axis

  Double_t cda  = Tracks[0].GetBeamAxisCDA();
  Double_t Zvtx = Tracks[0].GetBeamAxisVertex().Z();

  TVector3 KaonThreeMomentum = BeamParameters::GetInstance()->GetBeamThreeMomentum();
  TLorentzVector Kaon;
  Kaon.SetVectM(KaonThreeMomentum, MKCH);
  TLorentzVector Muon;
  Muon.SetVectM(Scand->GetThreeMomentumBeforeMagnet(), MMU);
  Double_t Mmiss2Mu = (Kaon-Muon).M2();
  Double_t Theta    = Kaon.Angle(Scand->GetThreeMomentumBeforeMagnet());

  if (FillPlots) {
    FillHisto("hZvtx", 0.001*Zvtx);
    FillHisto("hCDA", cda);
    FillHisto("hZvtxCDA", 0.001*Zvtx, cda);
  }
  Bool_t pas_zvtx = (Zvtx>130000 && Zvtx<180000);
  Bool_t pas_cda  = (cda<40);

  if (!pas_zvtx) return;
  if (!pas_cda) return;
  
  ////////////////////////////////////////////////////////////////////
  // Track-Cedar timing: look for the closest Cedar candidate to track

  Double_t dT_Track_Cedar =  999.999;
  Double_t CedarTime      = -999.999;
  for (Int_t i=0; i<NCEDARcand; i++) {
    TRecoCedarCandidate* Ccand = static_cast<TRecoCedarCandidate*>(CEDARevent->GetCandidate(i));
    if (Ccand->GetNSectors()<5) continue;
    Double_t dT = Ttrack - Ccand->GetTime(); // CHOD-Cedar, or Track-Cedar if no CHOD association
    if (fabs(dT) < fabs(dT_Track_Cedar)) {
      CedarTime = Ccand->GetTime();
      dT_Track_Cedar = dT;
    }
  }
  if (pas_zvtx && FillPlots) FillHisto("hdTTrackCedar", dT_Track_Cedar);
  if (fabs(dT_Track_Cedar)>10.0) return;
  
  ////////////////////////////////////////////////////////////////////
  // No LKr candidates not associated to the track.
  // There are no timing cuts: this works for low intensity data only.

  /*
  Int_t NLKRcand = LKRevent->GetNCandidates();
  if (!Tracks[0].LKrAssociationExists() && NLKRcand) return;
  if ( Tracks[0].LKrAssociationExists() && NLKRcand>1) return;
  */

  /////////////////////////
  // LAV veto (with timing)

  LAVMatching* pLAVMatching = *GetOutput<LAVMatching*>("PhotonVetoHandler.LAVMatching");
  pLAVMatching->SetReferenceTime(CedarTime);
  if (pLAVMatching->LAVHasTimeMatching(LAVevent)) return;

  /////////////////////////////////
  // IRC and SAC veto (with timing)

  SAVMatching* pSAVMatching = *GetOutput<SAVMatching*>("PhotonVetoHandler.SAVMatching");
  pSAVMatching->SetReferenceTime(CedarTime);
  pSAVMatching->SetIRCTimeCuts(10.0, 10.0); // half time window; default = 5ns
  pSAVMatching->SetSACTimeCuts(10.0, 10.0); // half time window; default = 5ns
  Bool_t SAVmatched = pSAVMatching->SAVHasTimeMatching(IRCevent, SACevent);
  if (SAVmatched) return;

  ////////////////////////////////////////////////
  // CHANTI veto: reject all candidates, no timing

  if (CHANTIevent->GetNCandidates()) return;

  /////////////////////////////////////////////////////////////////////////
  // MUV3 trigger performance check: all cuts except Track-MUV3 association

  if (FillPlots) {
    if (eop<0.2 && fabs(Mmiss2Mu*1e-6)<0.02) {
      FillHisto("hTrackXYMUV3", Tracks[0].xAtAfterMagnet(246800), Tracks[0].yAtAfterMagnet(246800));
    }
  }

  /////////////////////////////////////////////////////////////////////
  // Muon-Cedar timing: looking for the closest muon in time with Cedar

  if (!Tracks[0].MUV3AssociationExists()) return;

  Double_t dTmin = 9999;
  for (Int_t iMu=0; iMu<Tracks[0].GetNMUV3AssociationRecords(); iMu++) {
    Double_t dT = Tracks[0].GetMUV3Time(iMu) - CedarTime;
    if (fabs(dT)<fabs(dTmin)) dTmin = dT;
  }
  if (dTmin<-15.0 || dTmin>10.0) return;

  ///////////////////
  // Muon ID with E/p

  if (eop>0.2) return;

  if (FillPlots) {
    FillHisto("hMMiss2Mu",  Mmiss2Mu*1e-6); // [GeV^2]
    FillHisto("hPMMiss2Mu", 0.001*Ptrack, Mmiss2Mu*1e-6);
    FillHisto("hPTheta",    0.001*Ptrack, Theta);
    if (fabs(Mmiss2Mu*1e-6)<0.02) FillHisto("hKmu2EventsPerBurst", BurstID);
  }
  if (fabs(Mmiss2Mu*1e-6)>0.02) return;

  // CHOD response studies
  Int_t NCHODHits = CHODevent->GetNHits();
  if (FillPlots) FillHisto("hNCHODHits", NCHODHits);

  // RICH response studies
  Int_t NRICHHitsAll = RICHevent->GetNHits(); // including super-cells
  Int_t NRICHHits = 0;
  for (Int_t i=0; i<NRICHHitsAll; i++) {
    TRecoRICHHit* hit = static_cast<TRecoRICHHit*>(RICHevent->GetHit(i));
    if (hit->GetOrSuperCellID()==0) NRICHHits++; // no super-cells
  }
  if (FillPlots) FillHisto("hNRICHHits", NRICHHits);

  // LKr response studies
  if (GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[0], kLKr)) {
    Int_t NCells = LKRevent->GetNHits();
    Int_t NClusters = LKRevent->GetNCandidates();
    if (FillPlots) {
      FillHisto("hNLKrCells", NCells);
      FillHisto("hNLKrClusters", NClusters);
    }

    Double_t TotalCellEnergy = 0.0;
    Double_t TotalCellEnergy_in = 0.0;
    Double_t TotalCellEnergy_out = 0.0;
    Double_t TotalClusterEnergy = 0.0;
    Int_t NGoodCells = 0, NGoodCells_in = 0, NGoodCells_out = 0;
    for (Int_t i=0; i<NCells; i++) {
      TRecoLKrHit *hit = static_cast<TRecoLKrHit*>(LKRevent->GetHit(i));
      Double_t energy = hit->GetEnergy();
      if (energy<40.0) continue;
      NGoodCells++;
      TotalCellEnergy += energy;
      Double_t dx = hit->GetPosition().x() - Tracks[0].xAt(241093.0);
      Double_t dy = hit->GetPosition().y() - Tracks[0].yAt(241093.0);
      Double_t dist = sqrt(dx*dx+dy*dy);
      if (FillPlots) FillHisto("hCellTrackDistance", dist);
      if (dist<80.0) {
	NGoodCells_in++;
	TotalCellEnergy_in += energy;
      }
      else {
	NGoodCells_out++;
	TotalCellEnergy_out += energy;
      }
    }
    if (FillPlots) {
      FillHisto("hNLKrGoodCells", NGoodCells);
      FillHisto("hNLKrGoodCellsIn", NGoodCells_in);
      FillHisto("hNLKrGoodCellsOut", NGoodCells_out);
    }
    for (Int_t i=0; i<NClusters; i++) {
      TRecoLKrCandidate* Lcand = static_cast<TRecoLKrCandidate*>(LKRevent->GetCandidate(i));
      Double_t energy = Lcand->GetClusterEnergy();
      if (FillPlots) FillHisto("hLKrClusterEnergy", 0.001*energy);
      TotalClusterEnergy += energy;
    }
    if (FillPlots) {
      FillHisto("hLKrCellTotalEnergy", 0.001*TotalCellEnergy);
      FillHisto("hLKrCellTotalEnergyIn", 0.001*TotalCellEnergy_in);
      FillHisto("hLKrCellTotalEnergyOut", 0.001*TotalCellEnergy_out);
      FillHisto("hLKrCellTotalEnergyInOut",
		0.001*TotalCellEnergy_in, 0.001*TotalCellEnergy_out);
      FillHisto("hLKrCellClusterTotalEnergy",
		0.001*TotalCellEnergy, 0.001*TotalClusterEnergy);
    }
  }

  // Save the outputs
  fEventSelected = true;
  fMuonMomentum = Ptrack;
}

void Kmu2Selection::EndOfJobUser() {
  if (fReadingData) { // Data mode: save output
    SaveAllPlots();
    return;
  }
  if (!fHPhysicsEventsPerBurst) { // Histo mode required but no histograms found
    cout << user_normal() << "Asked to read my own output but cannot found it" << endl;
    return;
  }

  BuildPDFReport();
}

void Kmu2Selection::BuildPDFReport() {

  TString OutputPDFFileName = fAnalyzerName + ".pdf";
  gErrorIgnoreLevel = 5000; // suppress messages generated for each page printed
  gStyle->SetOptStat(11);

  TCanvas *Canvas0 = new TCanvas("Kmu2Canvas0");
  Canvas0->Print(Form(OutputPDFFileName + "["), "pdf"); // open file

  Canvas0->Divide(2,2);
  for (Int_t i=1; i<=4; i++) {
    Canvas0->GetPad(i)->SetLeftMargin(0.04);
    Canvas0->GetPad(i)->SetRightMargin(0.01);
    Canvas0->GetPad(i)->SetTopMargin(0.06);
    Canvas0->GetPad(i)->SetBottomMargin(0.10);
  }
  fHMass->SetLineColor(kBlue);
  fHMass->SetFillColor(kYellow);
  fHEOP->SetLineColor(kBlue);
  fHEOP->SetFillColor(kYellow);
  fHZvertex->SetLineColor(kBlue);
  fHZvertex->SetFillColor(kYellow);

  Canvas0->cd(1);
  fHZvertex->Draw();
  Canvas0->cd(2); gPad->SetLogy();
  fHEOP->Draw();
  Canvas0->cd(3);
  fHMass->Draw();

  Canvas0->Print(OutputPDFFileName, "pdf");
  Canvas0->Print(Form(OutputPDFFileName + "]"), "pdf"); // close file
  gErrorIgnoreLevel = -1; // restore the default 

  // PrintStatisticsPerBurst();
}

void Kmu2Selection::PrintStatisticsPerBurst() {
  for (Int_t i=1; i<=fHPhysicsEventsPerBurst->GetNbinsX(); i++) {
    Double_t N = fHPhysicsEventsPerBurst->GetBinContent(i);
    if (!N) continue;
    Double_t n = fHKmu2EventsPerBurst->GetBinContent(i);
    Double_t e = n/N;
    Double_t de = sqrt(e*(1.0-e)/N);
    cout << user_standard() << "@@Kmu2 "<<i-1<<" "<<n <<" "<<N<<" "<<e<<" "<<de<<endl;
  }
}
