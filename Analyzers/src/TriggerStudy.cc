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
using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

/// \class TriggerStudy
/// \Brief
/// Short description of your Analyzer
/// \EndBrief
///
/// \Detailed
/// Detailed description of your Analyzer\n\n
/// For examples of working Analyzer you can have a look at the examples in the Examples/ directory:\n
/// LKrPhotonMC \n
/// Pi0Reconstruction \n
/// SkimmingNoStrawMuonsTracks \n
/// Or the framework analyzers that can be found in the Analyzers/ directories: \n
/// CedarMCTester \n
/// VertexCDA \n
/// \n
/// All the following classes are available for you to use. Check their documentation for more information:\n
/// NA62Analysis::manip::TermManip \n
/// NA62Analysis::Analyzer \n
/// NA62Analysis::CounterHandler \n
/// NA62Analysis::EventFraction \n
/// NA62Analysis::MCSimple \n
/// NA62Analysis::NeuralNetwork \n
/// NA62Analysis::ParticleInterface \n
/// NA62Analysis::ParticleTree \n
/// NA62Analysis::StringBalancedTable \n
/// NA62Analysis::StringTable \n
/// NA62Analysis::UserMethods \n
/// NA62Analysis::Verbose \n
///
/// You might also be interested in checking the documentation for the following classes. However you should not
/// in general have to create instances of these. If necessary a pointer to the existing instance is usually
/// available or provided by specific methods.\n
/// NA62Analysis::Core::IOHisto \n
/// NA62Analysis::Core::IOTree \n
/// NA62Analysis::Core::IOHandler \n
/// NA62Analysis::Core::HistoHandler \n
/// NA62Analysis::Core::HistoHandler::Iterator \n
/// NA62Analysis::Core::PrimitiveReader \n
///
/// \EndDetailed

TriggerStudy::TriggerStudy(Core::BaseAnalysis *ba) : Analyzer(ba, "TriggerStudy")
{

RequestL0Data();
RequestL1Data();
fTriggerMaskPNN  = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-nQX-UTMC-nMUV-nLKr30");
	/// \MemberDescr
	/// \param ba : parent BaseAnalysis
	///
	/// Specify the trees you want to use and the event class corresponding\n
	/// Don't try to load MCTruth tree (RUN_0 or Event). Use the MCTruthEvent in Process function instead. Problems when opening twice the same tree.\n
	/// Example with RecoEvent\n
	///	\code
	/// 	RequestTree(new TRecoGigaTrackerEvent);
	///		RequestTree("GigaTracker", new TRecoGigaTrackerEvent);
	///		RequestTree("GigaTracker", new TRecoGigaTrackerEvent, "Reco");
	///		RequestTree("GigaTracker", new TRecoGigaTrackerEvent, "Digis");
	/// \endcode
	/// The first form can be used if the detector name is present in the class
	/// name.\n
	/// Example with MC Event\n
	///	\code
	/// 	RequestTree(new TGigaTrackerEvent);
	///		RequestTree("GigaTracker", new TGigaTrackerEvent);
	/// \endcode
	/// Example with generic tree\n
	///	\code
	///		RequestTree<MyClass>("MyTree", "BranchName", "MyClass", new MyClass);
	///		RequestTree("MyTree", "BranchName", "MyClass", new MyClass);
	/// 	RequestTree("MyTree", "BranchName", "int", new int);
	/// \endcode
	/// Requesting Trigger data\n
	///	\code
	///		RequestL0Data();
	///		RequestL1Data();
	///		RequestL2Data();
	/// \endcode
	/// \n
	/// Call one of: \n
	///	\code
	/// 	AddParam("paramName", &variableName, defaultValue);
	/// \endcode
	/// for each parameter of the analyzer. These parameters can be set when starting the FW from the command line with the -p option.\n
	/// paramName is the name of the parameter in the command line\n
	/// variableName is the name of the variable that should be declared in the definition of the class\n
	/// defaultValue is the default value if not specified in the command line\n
	/// The allowed types for parameters are the following: bool, int, long, float, double, char, string, TString\n
	/// \n
	/// A primitive file can be read in parallel to the event file. You can request to read primitives for a sub-detector
	/// with\n
	/// \code
	/// AddPrimitiveReader("detName", true/false);
	/// \endcode
	/// The boolean flag indicates if the primitives should be time sorted or kept in the original order they arrived from
	/// the detector.\n
	/// You can set the L0 matching window with \n
	/// \code
	/// SetL0MatchingWindowWidth("detName", val);
	/// SetL0MatchingWindowWidth("detName", timestamp, finetime);
	/// \endcode
	/// where val is a value in nanoseconds.
	/// \EndMemberDescr



}

void TriggerStudy::InitOutput(){
	/// \MemberDescr
	/// Register the output variables of the analyzer.\n
	/// Call: \n
	///	\code
	/// 	RegisterOutput("outputName", &variableName)
	/// \endcode
	/// for each variable that should be in the output of the Analyzer\n
	/// The name of the analyzer will be prepended to the outputName (to avoid collisions with other analyzers)\n
	/// variableName should be the name of a variable declared in the definition of the class\n
	/// \n
	/// To create a new TTree in the output file, call: \n
	///	\code
	/// 	void OpenNewTree("TTreeName", "TTreeTitle");
	/// \endcode
	/// TTreeName is the name of the TTree (will be used to refer to this TTree later)\n
	/// TTreeTitle is the title of the TTree\n
	/// \n
	/// To add a branch to the newly created TTree, call: \n
	///	\code
	/// 	void AddBranch<VariableType>("TTreeName", "BranchName", &pointer);
	/// \endcode
	/// VariableType is the type of the variable for this branch (fundamental data type or class)\n
	/// TTreeName is the name of the TTree to add this branch\n
	/// BranchName is the name of the branch\n
	/// pointer is a pointer to the variable (should be declared in the header file)\n
	/// \n
	/// To create a standard TTree containing KineParts (for candidates) use\n
	///	\code
	///     CreateStandardTree("TTreeName", "TTreeTitle");
	/// \endcode
	/// \EndMemberDescr
}

void TriggerStudy::InitHist(){
	/// \MemberDescr
	/// Book and Initialize histograms in this function.\n
	/// Same function to Book TH1, TH2, TGraph and TGraphAsymmErrors (anything derived from TH1 or TGraph)
	///	\code
	/// 	BookHisto(histogram*)
	/// \endcode

	BookHisto("hFullTrigStudy", new TH1I("TriggerStudy", "Kaon_Decay_Characteristics_afterL0", 10, 0, 10));
	BookHisto("hL0PNNChar", new TH1I("L0PNNChar", "L0PNN_Decay_Characteristics", 10, 0, 10));
	BookHisto("hLKREnergy", new TH1D("LKREnergyTest", "Energy_Distro_of_LKR", 30, 0, 70000));
	BookHisto("hPMom", new TH1D("PMomTest", "Momentum_Distro_from_Somewhere", 30, 0, 70000));
	BookHisto("hEOPCalc", new TH1D("EOPTestCalc", "E/p_thing", 30, 0, 1.5));
	
	BookHisto("hLKREoP",new TH1D("LKREOP", "E/p_thing", 30, 0, 1.5));
	BookHisto("hTotE", new TH1D("LKREnergyTot", "Energy_Distro_of_LKR", 30, 0, 70000));
	BookHisto("hTotEoP", new TH1D("LKRTotEoP", "E/p_thing", 30, 0, 1.5));
	BookHisto("hRICHMissingMass", new TH1D("RichMissingMass", "MissingMassReco", 100, -0.02, 0.04));
	BookHisto("hRICHMissingMass_cuts", new TH1D("RichMissingMass_cuts", "MissingMassReco", 100, -0.02, 0.04));
	

	BookHisto("hRICHring", new TH2D("RichRing", "Radius_of_ring_function_of_particle_momentum", 30, 0, 70000, 30, 0, 240));
	BookHisto("hLKREoP_cuts", new TH1D("EOP_cuts", "E/p_thing", 30, 0, 1.5));
	BookHisto("hRICHring_cuts", new TH2D("RichRing_cuts", "Radius_of_ring_function_of_particle_momentum", 30, 0, 70000, 30, 0, 240));


	///Comments
		/// If isAutotUpdate is true, this histogram will be drawn and updated regularly during the processing (default=false).\n
		/// The refresh interval can be set with (default=10):
		/// \code
		/// 	SetUpdateInterval(interval)
		/// \endcode
		/// Defining plots as AutoUpdate and setting the interval can also be done at runtime with a configuration file.\n
		/// \n
		/// Example of booking an histogram: \n
		///	\code
		/// 	BookHisto(new TH2I("PartEnergy", "Energy as a function of particle", 0, 0, 0, Bins, MinEnergy, MaxEnergy));
		/// \endcode
		/// THe histogram can be saved in a subdirectory of the output files in several ways:\n
		/// \code
		/// 	BookHisto(new TH2I("PartEnergy", "Energy as a function of particle", 0, 0, 0, Bins, MinEnergy, MaxEnergy), refresh, "dir1/dir2");
		/// 	BookHisto(new TH2I("dir1/dir2/PartEnergy", "Energy as a function of particle", 0, 0, 0, Bins, MinEnergy, MaxEnergy));
		/// 	BookHisto("dir1/dir2/PartEnergy", new TH2I("hPartEnergy", "Energy as a function of particle", 0, 0, 0, Bins, MinEnergy, MaxEnergy));
		/// \endcode
		/// In the first case, the ID of the histo in the framework (to use with FillHisto, ...) and its name in the output file is "PartEnergy".
		/// In the second and third case, the ID is "dir1/dir2/PartEnergy" and the name in the output file is respectively "PartEnergy" and "hPartEnergy".
		/// Booking of counters and creation of EventFraction can be done here with\n
		///	\code
		///		BookCounter(name)

		///		NewEventFraction(name)
		/// \endcode
		/// Example\n
		///	\code
		/// 	BookCounter("Total");
		/// 	BookCounter("PassCuts");
		/// 	NewEventFraction("Cuts");
		/// \endcode
		/// Add the counters to the EventFraction\n
		///	\code
		/// 	AddCounterToEventFraction(EventFractionName, CounterName)
		/// \endcode
		/// Example\n
		///	\code
		/// 	AddCounterToEventFraction("Cuts", "Total");
		/// 	AddCounterToEventFraction("Cuts", "PassCuts");
		/// \endcode
		/// Then define which counter represents the sample size\n
		///	\code
		/// 	DefineSampleSizeCounter(EventFractionName, CounterName)
		/// \endcode
		/// Example\n
		///	\code
		/// 	DefineSampleSizeCounter("Cuts", "Total");
		/// \endcode
		/// You can check the content of the input file (directory, histograms, generic key) in any directory with the following methods
		/// \code
		/// 	GetListOfKeys("CedarMonitoring"); // CedarMonitoring directory
		/// 	GetListOfHisto(""); // Top directory
		/// 	GetListOfTH1("CedarMonitoring");
		/// 	GetListOfTH2("CedarMonitoring");
		/// 	GetListOfTGraph("CedarMonitoring");
		/// 	GetListOfDirs("");
		/// \endcode
		/// The first one returns a vector of IOHandler::keyPair containing the name and class name of the object
		/// \code
		/// 	vector<IOHandler::keyPair> keys = GetListOfKeys("CedarMonitoring");
		/// 	cout << keys[0].name << " " << keys[0].className << endl;
		/// \endcode
		/// The others returns a vector of TString containing the name of the objects.\n
		/// You can retrieve histograms from the input ROOT file under the directory "dir1/dir2" (Anything derived from TH1) with\n
		///	\code
		/// 	RequestHistogram("dir1/dir2", "HistogramName", appendOnNewFile);
		/// 	RequestHistogram("dir1", "dir2/HistogramName", appendOnNewFile);
		/// \endcode
		/// In the first case the ID in the framework is "HistogramName", while it is "dir2/HistogramName" in the second case.
		/// appendOnNewFile is a boolean. If set to true, each time a new file is opened the content
		/// of the histogram will be appended to the content of the previous one. If set to false, the content
		/// of the histogram is replaced each time a new file is opened.
		/// \EndMemberDescr
	///


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
	
	BookCounter("ThreeTrack");
	BookCounter("TwoTrackwLepton");
	BookCounter("notAutopass");

	BookCounter("TotalEoP_counts");
	BookCounter("MuonExcluded");
	BookCounter("ElectronExcluded");
	BookCounter("Remains");

	BookCounter("TotalRich_counts");
	BookCounter("Rich_Included");
	BookCounter("Rich_Excluded");

	BookCounter("RecoMass_Included");
	BookCounter("RecoMass_Excluded");




	/// Tables
	NewEventFraction("TriggerAnalysis");
	NewEventFraction("L0PNNAnalysis");
	NewEventFraction("EoPCuts");
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
	DefineSampleSizeCounter("L0PNNAnalysis", "L0PNN");



	

	AddCounterToEventFraction("EoPCuts", "TotalEoP_counts");
	AddCounterToEventFraction("EoPCuts", "MuonExcluded");
	AddCounterToEventFraction("EoPCuts", "ElectronExcluded");
	AddCounterToEventFraction("EoPCuts", "Remains");
	DefineSampleSizeCounter("EoPCuts", "TotalEoP_counts");

	AddCounterToEventFraction("RichCuts", "TotalRich_counts");
	AddCounterToEventFraction("RichCuts", "Rich_Included");
	AddCounterToEventFraction("RichCuts", "Rich_Excluded");
	DefineSampleSizeCounter("RichCuts", "TotalRich_counts");
	
	AddCounterToEventFraction("RecoMassCuts", "TotalRich_counts");
	AddCounterToEventFraction("RecoMassCuts", "RecoMass_Included");
	AddCounterToEventFraction("RecoMassCuts", "RecoMass_Excluded");
	DefineSampleSizeCounter("RecoMassCuts", "TotalRich_counts");
	
	


}

void TriggerStudy::DefineMCSimple(){
	/// \MemberDescr
	/// Setup of fMCSimple. You must specify the generated MC particles you want.\n
	/// Add particles you want to recover from fMCSimple\n
	///	\code
	/// 	int particleID = fMCSimple.AddParticle(parentID, pdgCode)
	///	\endcode
	/// parentID : 	0=no parent (=beam particle)\n
	/// 	...\n
	/// Example : you want to retrieve the kaon from the beam, the pi0 an pi+ from the beam kaon and the 2 photons coming from the previous pi0 decay :\n
	///	\code
	/// 	int kaonID = fMCSimple.AddParticle(0, 321) //Ask beam kaon (sequence ID=1)
	/// 	fMCSimple.AddParticle(kaonID, 211) //Ask pi+ from previous kaon (sequence ID=2)
	/// 	int pi0ID = fMCSimple.AddParticle(kaonID, 111) //Ask pi0 from previous kaon (sequence ID=3)
	/// 	fMCSimple.AddParticle(pi0ID, 22) //Ask first gamma from previous pi0 (sequence ID=4)
	/// 	fMCSimple.AddParticle(pi0ID, 22) //Ask second gamma from previous pi0 (sequence ID=4)
	///	\endcode
	///
	/// @see ROOT TDatabasePDG for a list of PDG codes and particle naming convention
	/// \EndMemberDescr
}

void TriggerStudy::StartOfRunUser(){
	/// \MemberDescr
	/// This method is called at the beginning of the processing (corresponding to a start of run in the normal NA62 data taking)\n
	/// Do here your start of run processing if any
	/// \EndMemberDescr
}

void TriggerStudy::StartOfBurstUser(){
	/// \MemberDescr
	/// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the beginning of the first file\n
	/// Do here your start/end of burst processing if any
	/// \EndMemberDescr
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

	///
		/// \MemberDescr
		/// \param iEvent : Event number
		///
		/// Main process method. Called on each event. Write you analysis here.\n
		/// You can retrieve MC particles from the fMCSimple set with (returns a vector<KinePart*>)\n
		/// \code
		/// 	fMCSimple["particleName"]
		/// 	fMCSimple[pdgID]
		/// \endcode
		/// Example\n
		/// \code
		/// 	fMCSimple["K+"][index]; //for the kaon
		/// 	fMCSimple["pi+"][index]; //for the positive pion
		/// 	fMCSimple["gamma"][index]; //for the photon
		/// \endcode
		/// The number in the brackets is the index of the particle (if you asked for two photons in the set, you can ask fMCSimple["gamma"][0] for the first one and fMCSimple["gamma"][1] for the second)\n
		/// \n
		/// If you need a property of a particle, you can make a call to fParticleInterface (instance of the ParticleInterface class).\n
		///	This class has two methods FindParticle that will return a TParticlePDG with the required particle. You can search by pdgID or by name.\n
		///	This class also provide two methods to switch between particle name and pdgID if necessary.\n
		///	Example\n
		/// \code
		/// 	double kaonMass = fParticleInterface->FindParticle(321).Mass();
		/// 	double pi0Lifetime = fParticleInterface->FindParticle("pi0").Lifetime();
		/// \endcode
		/// You can retrieve the events from the trees with\n
		/// \code
		/// 	GetEvent<eventClass>();
		/// 	GetEvent<eventClass>("Digis");
		/// 	(eventClass*)GetEvent("detectorName");
		/// 	(eventClass*)GetEvent("detectorName", "Digis");
		/// \endcode
		/// The first two forms can be used only if the class is of type
		/// T<i>DetectorName</i>Event or TReco<i>DetectorName</i>Event.\n
		/// You can retrieve data from generic TTrees with\n
		/// \code
		/// 	GetObject<MyClass>("treeName");			// For class objects
		/// 	GetPrimitiveObject<int>("treeName");	// For primitive datatype (int, double, Int_t, ...)
		/// \endcode
		/// You can retrieve full MC events if available ( GetWithMC() ) with\n
		/// \code
		/// 	GetMCEvent();
		/// 	GetMCEvent("Digis");
		/// \endcode
		/// You can retrieve EventHeader if available ( GetWithEventHeader() ) with\n
		/// \code
		/// 	GetEventHeader();
		/// 	GetEventHeader("Digis");
		/// \endcode
		/// You can retrieve Trigger data if requested with\n
		/// \code
		/// 	GetL0Data();
		/// 	GetL1Data();
		/// 	GetL2Data();
		/// \endcode
		/// You can retrieve the histograms you booked (for drawing, changing, filling, ...) with\n
		/// \code
		/// 	fHisto.GetTH1("histoName");// for TH1
		/// 	fHisto.GetTH2("histoName");// for TH2
		/// 	fHisto.GetTGraph("graphName");// for TGraph and TGraphAsymmErrors
		/// 	fHisto.GetHisto("histoName");// for TH1 or TH2 (returns a TH1 pointer)
		/// \endcode
		/// To fill the histograms you can use\n
		/// \code
		/// 	FillHisto("histoName", values)
		/// \endcode
		/// where values are the same parameters as if you call histogram->Fill(values) (x,y,weight,...)\n
		/// If the histogram is not found, an error message is printed\n
		/// \n
		/// Modify a counter with one of the following methods\n
		/// \code
		/// 	IncrementCounter(name)
		/// 	IncrementCounter(name, delta)
		/// 	DecrementCounter(name)
		/// 	DecrementCounter(name, delta)
		/// 	SetCounterValue(name, value)
		/// \endcode
		/// \n
		/// To use the output of a different analyzer, use\n
		/// \code
		/// 	outputType *var = GetOutput<outputType>("analyzerName.outputName", state);
		/// \endcode
		/// Where outputType is the variable type and state is of type outputState\n
		/// State is set with the state of the variable (kOUninit, kOInvalid ,kOValid). The value of the output should only be trusted if state == kOValid\n
		/// example :
		/// \code
		/// 	TLorentzVector vertex = *(TLorentzVector*)GetOutput("simpleVertexAnalyzer.vertex", state);
		/// \endcode
		/// Before starting the processing of an event, the state flag of each output variable is reset to kOUninit\n
		/// When setting the value of an output variable, don't forget to set appropriately the state flag to either kOValid or kOInvalid\n
		/// to indicate if the value can/can't be used in other analyzer\n
		/// \code
		/// 	SetOutputState("outputName", kOValid);
		/// \endcode
		/// If you want to append a candidate in one of your standard output Tree, use\n
		/// \code
		/// 	KinePart *candidate = CreateStandardCandidate("treeName");
		/// \endcode
		/// and fill the properties of your candidate. It will be automatically written in the output tree.\n
		///	\n
		/// If you want to save this event in your custom and standard TTrees (not the input tree replication), call\n
		/// \code
		/// 	FillTrees();
		/// 	FillTrees("treeName");
		/// \endcode
		/// This will call the Fill method of every TTree created in this analyzer.\n
		///	\n
		/// If you want to replicate this event in the output file, call\n
		/// \code
		/// 	FilterAccept();
		/// \endcode
		/// The structure of all the trees that have been opened (by all Analyzer) will be copied in the output file\n
		/// and the events for which at least one analyzer called FilterAccept() will be replicated in the output trees.\n
		/// For more information, please refer to the \ref eventFiltering "Event Filtering" page
		/// \n
		/// The primitives associated to the current event can be retrieved in two ways:\n
		/// \code
		/// FindAllPrimitiveInMatchingWindow("detName");
		/// \endcode
		/// allows to get all the primitives in a certain window around the event timestamp (L0MatchingWindowWidth)
		/// while
		/// \code
		/// FindMatchingPrimitive("detName");
		/// \endcode
		/// allows to get the primitive closest the the event timestamp (but nevertheless withing the
		/// L0MatchingWindowWidth).
		/// @see ROOT TParticlePDG for the particle properties
		/// @see ROOT TDatabasePDG for a list of PDG codes and particle naming convention
		/// \EndMemberDescr
		///if(fMCSimple.fStatus == MCSimple::kMissing){printIncompleteMCWarning(iEvent);return;}
		///if(fMCSimple.fStatus == MCSimple::kEmpty){printNoMCWarning();return;}

	///

	

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
	///Bool_t TwoTrack = *(Bool_t*)GetOutput("FilterTwoTrackVertexWithLepton.EventSelected"); ///not good


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

	///for the MC study, use the l0trigger emulator hh and cc to find the way to the specific emulator, and then the definition is there
	/// EOP tests (Help from KMu2Selection.cc)
	/// do loop over all the vectors if you want
	std::vector<DownstreamTrack> Tracks = *GetOutput<std::vector<DownstreamTrack>>("DownstreamTrackBuilder.Output");
	
	///L0TriggerOnPNN = true;

    if(L0TriggerOnPNN) {


		IncrementCounter("L0PNN"); 
		FillHisto("hFullTrigStudy", 2); ///L0PNN
		FillHisto("hL0PNNChar", 0);

	

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
		

		if (Tracks.size() != 1) {return;}
		///next up: try to loop over the vectors

		Double_t Ptrack = Tracks[0].GetMomentum(); ///good sht
		Double_t LKREnergy = Tracks[0].GetLKrEnergy(); ///good sht
		Double_t LKREoP = Tracks[0].GetLKrEoP(); ///good sht
		Double_t TotEnergy = Tracks[0].GetLKrTotalEnergy();
		Double_t TotEoP = Tracks[0].GetLKrTotalEoP();

		Double_t RichRing = Tracks[0].GetRICHRingRadius();
		Double_t RichMass = Tracks[0].GetRICHSingleRingTrkCentredMass();
		Double_t RichMass2 = RichMass*RichMass/1000000;
		Double_t LowMassLim = 0.005;
		Double_t HighMassLim = 0.015;

		FillHisto("hLKREnergy", LKREnergy);
		FillHisto("hPMom", Ptrack);
		FillHisto("hEOPCalc", LKREnergy/Ptrack);

		FillHisto("hLKREoP", LKREoP);
		FillHisto("hTotE", TotEnergy);
		FillHisto("hTotEoP",TotEoP);
		FillHisto("hRICHring", Ptrack, RichRing);
		FillHisto("hRICHMissingMass", RichMass2);

		

		///EoP Cuts
		IncrementCounter("TotalEoP_counts");
		if (LKREoP < 0.2) {
			IncrementCounter("MuonExcluded");
		}
		else if (LKREoP >0.9){
			IncrementCounter("ElectronExcluded");
		}
		else
		{
			IncrementCounter("Remains");
			FillHisto("hLKREoP_cuts", LKREoP);
		}

		///Rich Cuts
		IncrementCounter("TotalRich_counts");
		if(Ptrack > 15000 && Ptrack < 35000 && RichRing < 180){
			IncrementCounter("Rich_Included");
			FillHisto("hRICHring_cuts", Ptrack, RichRing);
		}
		else {
			IncrementCounter("Rich_Excluded");
		}
	
		///RecoMass Cuts
		if(LowMassLim < RichMass2 && RichMass2 < HighMassLim ){
			IncrementCounter("RecoMass_Included");
			FillHisto("hRICHMissingMass_cuts", RichMass2);
		}
		else{
			IncrementCounter("RecoMass_Excluded");
		}

	}

	///labelling histograms
	Int_t ix;
	const Int_t nx1 = 9;
	const Int_t nx2 = 7;
	const char *labels1[nx1] = {"TotalEvent","PhysicsEvent", "L0PNN", "Kmu2Selection","K2piCounter",
      "K3piCounter","Ke3Selection","Kmu3SelectionNoSpectrometer","Pi0Selection"};
	const char *labels2[nx2] = {"L0PNN", "Main5", "Kmu2Selection","K2piCounter",
      "K3piCounter","Ke3Selection","Kmu3SelectionNoSpectrometer"};

	TH1 *MyHisto1 = fHisto.GetHisto("hFullTrigStudy");
	TH1 *MyHisto2 = fHisto.GetHisto("hL0PNNChar"); 

	for (ix=1;ix<=nx1;ix++) MyHisto1->GetXaxis()->SetBinLabel(ix,labels1[ix-1]); 
	for (ix=1;ix<=nx2;ix++) MyHisto2->GetXaxis()->SetBinLabel(ix,labels2[ix-1]); 

}

void TriggerStudy::PostProcess(){
	/// \MemberDescr
	/// This function is called after an event has been processed by all analyzers. It could be used to free some memory allocated
	/// during the Process.
	/// \EndMemberDescr

}

void TriggerStudy::EndOfBurstUser(){
	/// \MemberDescr
	/// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the end of the last file\n
	/// Do here your start/end of burst processing if any.
	/// Be careful: this is called after the event/file has changed.
	/// \EndMemberDescr
}

void TriggerStudy::EndOfRunUser(){
	/// \MemberDescr
	/// This method is called at the end of the processing (corresponding to a end of run in the normal NA62 data taking)\n
	/// Do here your end of run processing if any\n
	/// \EndMemberDescr

}

void TriggerStudy::EndOfJobUser(){
    /// \MemberDescr
    /// Called at the end of job, just before exiting NA62Analysis.\n
    /// Do here all output operation needed (write plots, objects)
    /// \n
	/// If you want to save all plots, just call\n
	/// \code
	/// 	SaveAllPlots();
	SaveAllPlots();

	///PDF_Files/TriggerStudy
	///Could have used iterator but just need this to work

	TCanvas *c = new TCanvas;
	
	fHisto.GetTH1("hFullTrigStudy")->Draw();
	c->SaveAs("PDF_Files/TriggerStudy/hFullTrigStudy.pdf");

	fHisto.GetTH1("hL0PNNChar")->Draw();
	c->SaveAs("PDF_Files/TriggerStudy/hL0PNNChar.pdf");
	
	fHisto.GetTH1("hLKREnergy")->Draw();
	///c->SaveAs("PDF_Files/TriggerStudy/hLKREnergy.pdf");

	fHisto.GetTH1("hPMom")->Draw();
	///c->SaveAs("PDF_Files/TriggerStudy/hPMom.pdf");

	fHisto.GetTH1("hEOPCalc")->Draw();
	///c->SaveAs("PDF_Files/TriggerStudy/hEOPCalc.pdf");


	fHisto.GetTH1("hLKREoP")->Draw();
	c->SaveAs("PDF_Files/TriggerStudy/hLKREoP.pdf");

	fHisto.GetTH1("hTotE")->Draw();
	///c->SaveAs("PDF_Files/TriggerStudy/hTotE.pdf");

	fHisto.GetTH1("hTotEoP")->Draw();
	///c->SaveAs("PDF_Files/TriggerStudy/hTotEoP.pdf");

	fHisto.GetTH1("hRICHMissingMass_cuts")->Draw();
	c->SaveAs("PDF_Files/TriggerStudy/hRICHMissingMass_cuts.pdf");
	
	fHisto.GetTH1("hRICHMissingMass")->Draw();
	c->SaveAs("PDF_Files/TriggerStudy/hRICHMissingMass.pdf");


	fHisto.GetTH2("hRICHring")->Draw();
	c->SaveAs("PDF_Files/TriggerStudy/hRICHring.pdf");

	fHisto.GetTH1("hLKREoP_cuts")->Draw();
	c->SaveAs("PDF_Files/TriggerStudy/hLKREoP_cuts.pdf");
	
	fHisto.GetTH2("hRICHring_cuts")->Draw();
	c->SaveAs("PDF_Files/TriggerStudy/hRICHring_cuts.pdf");
	delete c;
   

	/// \endcode
	/// Or you can just save the ones you want with\n
	/// \code
	/// 	histogram->Write()\n
	///		fHisto.Get...("histoname")->Write();
	/// \endcode
	/// \n
	/// To run over a set of histograms you can use Iterators (HistoHandler::IteratorTH1,
	/// HistoHandler::IteratorTH2, HistoHandler::IteratorTGraph). You can use it to run over
	/// all the histograms or only a subset of histogram by using one of the two forms of
	/// GetIterator...  (replace ... by TH1, TH2 or TGraph)\n
	/// \code
	/// 	GetIterator...()
	/// \endcode
	/// will get an Iterator running over all histograms of this type while
	/// \code
	/// 	GetIterator...("baseName")
	/// \endcode
	/// will get an Iterator running only over the histograms of this type whose name starts
	/// with baseName.\n
	/// For more details and examples on how to use the Iterator after getting it, please refer
	/// to the HistoHandler::Iterator documentation.\n
	/// Although this is described here, Iterators can be used anywhere after the
	/// histograms have been booked.
    /// \EndMemberDescr
}

void TriggerStudy::DrawPlot(){
	/// \MemberDescr
	/// This method is called at the end of processing to draw plots when the -g option is used.\n
	/// If you want to draw all the plots, just call\n
	/// \code
	/// 	DrawAllPlots();
	/// \endcode
	/// Or get the pointer to the histogram with\n
	/// \code
	/// 	fHisto.GetTH1("histoName");// for TH1
	/// 	fHisto.GetTH2("histoName");// for TH2
	/// 	fHisto.GetGraph("graphName");// for TGraph and TGraphAsymmErrors
	///     fHisto.GetHisto("histoName");// for TH1 or TH2 (returns a TH1 pointer)
	/// \endcode
	/// and manipulate it as usual (TCanvas, Draw, ...)\n
	/// \EndMemberDescr

}

TriggerStudy::~TriggerStudy(){
	/// \MemberDescr
	/// Destructor of the Analyzer. If you allocated any memory for class
	/// members, delete them here.
	/// \EndMemberDescr
}
