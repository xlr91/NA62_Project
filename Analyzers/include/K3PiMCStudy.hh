#ifndef K3PIMCSTUDY_HH
#define K3PIMCSTUDY_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include <TCanvas.h>
#include "TriggerConditions.hh"

#include "L0PrimitiveHandler.hh"
#include "K3piSelection.hh" ///change this 
#include "GeometricAcceptance.hh"
#include "DownstreamTrack.hh"
#include "SpectrometerTrackVertex.hh"
#include "TProfile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "BeamParameters.hh"
#include "LAVMatching.hh"
#include "ConfigSettings.hh"

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
L0PrimitiveHandler* fPrimitiveHandler;
};
#endif
