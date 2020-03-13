#ifndef TRIGGERSTUDY_HH
#define TRIGGERSTUDY_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include <TCanvas.h>
#include "TriggerConditions.hh"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <stdio.h>     
#include <math.h>  
#include "TProfile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"

class TH1I;
class TH2F;
class TGraph;
class TTree;


class TriggerStudy : public NA62Analysis::Analyzer
{
	public:
		explicit TriggerStudy(NA62Analysis::Core::BaseAnalysis *ba);
		~TriggerStudy();
		void InitHist();
		void InitOutput();
		void DefineMCSimple();
		void ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType);
		void Process(int iEvent);
		void PostProcess();
		void StartOfBurstUser();
		void EndOfBurstUser();
		void StartOfRunUser();
		void EndOfRunUser();
        void EndOfJobUser();
		void DrawPlot();
		Double_t RICHlims(Double_t p, Double_t lim);
	protected:

Int_t fTriggerMaskPNN;
Double_t Frichval = 17000;
Double_t nrichval = 1.000061;
Double_t massofpion2 = 139.570*139.570;
Int_t highlim = 8;
Int_t lowlim = -2; 



};
#endif
