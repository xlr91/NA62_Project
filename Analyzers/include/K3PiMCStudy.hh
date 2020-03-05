#ifndef K3PIMCSTUDY_HH
#define K3PIMCSTUDY_HH

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


class K3PiMCStudy : public NA62Analysis::Analyzer
{
	public:
		explicit K3PiMCStudy(NA62Analysis::Core::BaseAnalysis *ba);
		~K3PiMCStudy();
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
