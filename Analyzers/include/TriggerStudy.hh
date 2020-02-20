#ifndef TRIGGERSTUDY_HH
#define TRIGGERSTUDY_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include <TCanvas.h>
#include "TriggerConditions.hh"

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
	protected:

Int_t fTriggerMaskPNN;



};
#endif
