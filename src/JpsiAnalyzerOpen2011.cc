//ACESS DATA WITH TRIGGER
// -*- C++ -*-
//
// Package:    JpsiAnalyzerOpen2011
// Class:      JpsiAnalyzerOpen2011
// 
/**\class JpsiAnalyzerOpen2011 JpsiAnalyzerOpen2011.cc CMSOpenData2011/JpsiAnalyzerOpen2011/src/JpsiAnalyzerOpen2011.cc
Description: [one line class summary]
Implementation:
[Notes on implementation]
*/
//
// Original Author:  Sandro Fonseca De Souza,32 4-C14,+41227674949,
//         Created:  Thu Apr 13 21:58:44 CEST 2017
// $Id$
//
//


// system include files
#include <memory>
#include <cmath>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TLorentzVector.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Muon Reco
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
/*
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
*/
//includes ROOT
#include <TROOT.h>
#include <TTree.h>

//
// class declaration
//

class JpsiAnalyzerOpen2011 : public edm::EDAnalyzer {
	public:
		explicit JpsiAnalyzerOpen2011(const edm::ParameterSet&);
		~JpsiAnalyzerOpen2011();


	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		bool triggerfired(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname);
		bool triggerfound(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname);

		unsigned int triggerIndex(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname);

		// ----------member data ---------------------------

		//Declaring tree
		TTree* AnalysisTree;
		TTree* PlotControl;
		
		bool verbose_;
		bool triggerflag_;
		edm::InputTag primaryVertexProducer_;
		edm::InputTag recoMuons_;
		// Reco configs
		double minMuPt_;
		double maxMuEta_;
		double muonLeadPt_, muonTrailPt_;
		double minJPsiMass_ ;
		// Trigger
		edm::InputTag hlTriggerEvent_;      // Input tag for TriggerEvent
		//std::string triggerName_;
				
		std::vector<std::string> triggerName_;
		std::vector<std::string> NameTrigger; 
		edm::InputTag hlTriggerResults_;    // Input tag for TriggerResults
		edm::InputTag collectionName_;
		bool debug;
		//Counters
		int irun, ievt,lumiBlock;
		int Total_Events = 0;
		int CounterMuons = 0;
		int TrackerMuon	= 0;
		int nDimuon = 0;
		int TagAndProbes = 0;
		int countLooseMuons = 0; 

		//Lorentz Vector
		TLorentzVector mu_1;
		TLorentzVector mu_2;
		TLorentzVector mu1mu2;

		double M = 0.;
		double Pt = 0.;
		double Eta = 0.;
		double Rapidity = 0.;
		
		//Trigger
		int countInAccepted = 0;
		int countInTriggered = 0 ;

		double TagMuon_Pt;
		double TagMuon_Eta;
		double TagMuon_Phi;
		
		double ProbeMuon_Pt;
		double ProbeMuon_Eta;
		double ProbeMuon_Phi;
		
		int PassingProbeTrackingMuon;
		double MassPassingTagProbeTracking;
		double MassFailingTagProbeTracking;

		double PassingProbeStandAloneMuon;

		double PassingProbeIDMuon;

		//Creating Vectors



		



};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
JpsiAnalyzerOpen2011::JpsiAnalyzerOpen2011(const edm::ParameterSet& iConfig):
	verbose_ (iConfig.getParameter< bool > ("verbose")),
	triggerflag_ (iConfig.getParameter< bool > ("triggerflag")),
	primaryVertexProducer_ (iConfig.getParameter<edm::InputTag>("primaryVertexProducer")),
	recoMuons_(iConfig.getParameter<edm::InputTag>("recoMuonsLabel")),
	// Reco config with the trigger
	minMuPt_ (iConfig.getParameter<double>("minMuPt")),
	maxMuEta_ (iConfig.getParameter<double>("maxMuEta")),
	muonLeadPt_ (iConfig.getParameter<double>("minMuonLeadPt")),
	muonTrailPt_ (iConfig.getParameter<double>("minMuonTrailPt")),
	minJPsiMass_ (iConfig.getParameter<double>("minJPsiMass")),
	hlTriggerEvent_ (iConfig.getUntrackedParameter<edm::InputTag>("TriggerEventTag", edm::InputTag("hltTriggerSummaryAOD", "", "HLT"))),
//	triggerName_ (iConfig.getUntrackedParameter<std::string>("PathName","HLT_Dimuon0_Jpsi_v")),
        triggerName_ (iConfig.getUntrackedParameter<std::vector<std::string>>("PathName")),
	hlTriggerResults_  (iConfig.getUntrackedParameter<edm::InputTag>("TriggerResultsTag", edm::InputTag("TriggerResults", "", "HLT")))
{
	// Histos File
	edm::Service<TFileService> fs;

	// Define Histos
	TH1D::SetDefaultSumw2();
	//HistoMuon_Pt = fs->make<TH1F>("HistoMuon_Pt"  , "Muon_Pt"  ,   100,   0., 200.);
	

	//Define Trees
	AnalysisTree = fs->make<TTree>("AnalysisTree","Muon Analysis Tree");
	PlotControl = fs->make<TTree>("PlotControl","Muon Analysis Tree");



	

}

//------------------------------------------------------------------------------------------------------------
JpsiAnalyzerOpen2011::~JpsiAnalyzerOpen2011()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
    std::vector<std::string>().swap(NameTrigger);
	
	TagMuon_Pt = 0.;
	TagMuon_Eta = 0.;
	TagMuon_Phi = 0.;
	
	ProbeMuon_Pt = 0.;
	ProbeMuon_Eta = 0.;
	ProbeMuon_Phi = 0.;
	
	PassingProbeTrackingMuon = 0;
	MassPassingTagProbeTracking = 0;

	PassingProbeStandAloneMuon = 0;
	PassingProbeIDMuon = 0;

	//See if this happen between each event
	std::cout << "Cleaning Variables: " << std::endl;
		

}
//------------------------------------------------------------------------------------------------------------

//TRIGGERS
bool JpsiAnalyzerOpen2011::triggerfired(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname){
	const edm::TriggerNames TrigNames_ = ev.triggerNames(*triggerResultsHandle_);
	const unsigned int ntrigs = triggerResultsHandle_->size();
	for (unsigned int itr=0; itr<ntrigs; itr++){
		TString trigName = TrigNames_.triggerName(itr);
		if (!triggerResultsHandle_->accept(itr)) continue;
		if(trigName.Contains(trigname))      return true;
		if(verbose_) std::cout << "HLT trigger Name avaliable (fired) : " << trigName << std::endl;
	}
	return false;
}


bool JpsiAnalyzerOpen2011::triggerfound(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname){
	const edm::TriggerNames TrigNames_ = ev.triggerNames(*triggerResultsHandle_);
	const unsigned int ntrigs = triggerResultsHandle_->size();
	for (unsigned int itr=0; itr<ntrigs; itr++){
		TString trigName = TrigNames_.triggerName(itr);
		if(trigName.Contains(trigname) )     return true;
		if(verbose_) std::cout << "HLT trigger Name avaliable : " << trigName << std::endl;

	}
	return false;
}

unsigned int JpsiAnalyzerOpen2011::triggerIndex(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname){
	const edm::TriggerNames TrigNames_ = ev.triggerNames(*triggerResultsHandle_);
	const unsigned int ntrigs = triggerResultsHandle_->size();
	unsigned int itr;
	for (itr=0; itr< ntrigs; itr++){
		TString trigName=TrigNames_.triggerName(itr);
		if(trigName.Contains(trigname))      return itr;

	}
	return itr;
}

//
// member functions
//

// ------------ method called for each event  ------------
void JpsiAnalyzerOpen2011::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ //begin analize
    if (verbose_) std::cout<<"******************* STARTING ANALYSIS************"<<std::endl;
	ievt = iEvent.id().event();
	irun = iEvent.run();
	lumiBlock = iEvent.luminosityBlock();

	if (verbose_) std::cout<<" Run # "<<irun<<" Evt # "<<ievt<<" lumiBlock # "<<lumiBlock<<std::endl;
          
	using namespace edm;
	
	//Look for the Primary Vertex (and use the BeamSpot instead, if you can't find it):
  reco::Vertex::Point posVtx;
  reco::Vertex::Error errVtx;
  unsigned int theIndexOfThePrimaryVertex = 999.;

  edm::Handle<reco::VertexCollection> vertex;
  //iEvent.getByToken(theVertexLabel_, vertex);
	iEvent.getByLabel("offlinePrimaryVertices", vertex);
  if (vertex.isValid()) {
    for (unsigned int ind = 0; ind < vertex->size(); ++ind) {
      if ((*vertex)[ind].isValid() && !((*vertex)[ind].isFake())) {
        theIndexOfThePrimaryVertex = ind;
        break;
      }
    }
  }

  if (theIndexOfThePrimaryVertex < 100) {
    posVtx = ((*vertex)[theIndexOfThePrimaryVertex]).position();
    errVtx = ((*vertex)[theIndexOfThePrimaryVertex]).error();
  } else {
    LogInfo("RecoMuonValidator") << "reco::PrimaryVertex not found, use BeamSpot position instead\n";

    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
   // iEvent.getByToken(theBeamSpotLabel_, recoBeamSpotHandle);
	iEvent.getByLabel("offlineBeamSpot", recoBeamSpotHandle);
    reco::BeamSpot bs = *recoBeamSpotHandle;

    posVtx = bs.position();
    errVtx(0, 0) = bs.BeamWidthX();
    errVtx(1, 1) = bs.BeamWidthY();
    errVtx(2, 2) = bs.sigmaZ();
  }

	 const reco::Vertex vtx(posVtx, errVtx);

	//int nVertices = recoVertices->size();
	//if (verbose_) std::cout << "nVertices "<< nVertices << std::endl;
	//if(nVertices>0) vertex = & ((*recoVertices)[0]);


	//Accessing Muon Collection
	edm::Handle< reco::MuonCollection > recoMuons;
	iEvent.getByLabel(recoMuons_, recoMuons);
	//is valid method
       
   Total_Events++;

	//------------------------------------------------------------------
	//if(triggerflag_){

    //std::cout << "Using Trigger selection "<< std::endl;

	// *** get Handle to the TriggerEvent
	edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
	iEvent.getByLabel(hlTriggerEvent_,triggerEventHandle_);
	if (!triggerEventHandle_.isValid()) 
	{
		std::cout << "Error in getting TriggerEvent product from Event!" << std::endl;
		return;
	}

	// *** get Handle to the TriggerResults
	edm::Handle<TriggerResults> HLTR; 
	iEvent.getByLabel(hlTriggerResults_, HLTR);
	if (!HLTR.isValid()) 
	{
		std::cout << "HLT TriggerResults with label " << hlTriggerResults_ << " not found!";
		return;
	}
    
	// Only events in which the path actually fired had stored the filter results and products:	  
    for (unsigned int i=0; i<triggerName_.size(); i++) 
	{   
       
		bool triggerFound = triggerfound(iEvent,HLTR,triggerName_[i]);
		if (triggerFound) countInAccepted++; 
		if (triggerFound == true && verbose_) std::cout<< "HLT trigger is found: " << triggerFound <<" Trigger Path found : " << triggerName_[i] << std::endl;
		bool triggerFired = triggerfired(iEvent,HLTR,triggerName_[i]);
		if (triggerFired) countInTriggered++; 
		if (triggerFired == true && verbose_) std::cout<< "TriggerFired : "<< triggerFired <<" Trigger Path :"<< triggerName_[i] << std::endl;
		const unsigned int numberOfHltPaths = HLTR->size();
		//const unsigned int numberOfHltFilters = triggerEventHandle_->sizeFilters();
		if (triggerFired == true) NameTrigger.push_back(triggerName_[i]);
		unsigned int pathIndex = triggerIndex(iEvent,HLTR,triggerName_[i]);
		if (pathIndex>=numberOfHltPaths) {
		std::cout << " WARNING: path " << triggerName_[i] << " out of range in TriggerResults" << std::endl;
		return;
		}

		if (HLTR->wasrun(pathIndex)) {
		if (!triggerFound) std::cout << " WARNING: path exists in HLTR but not found in TriggerResults" << std::endl;
		}
		else {
		if (triggerFound) std::cout << " WARNING: path found in TriggerResults but it does not exist in HLTR" << std::endl;
		}
	
    //}// end trigger slection
	//------------------------------------------------------------------

	if (triggerFired == true)
	{ //triggerfired

	//Specific vector for muons (it has more information than normal C++ vector)
	std::vector<reco::Muon> myLeptons;
	
	//for (reco::MuonCollection::const_iterator recoMu1 = muons->begin(); recoMu1!=muons->end(); ++recoMu1) {
	for (reco::MuonCollection::const_iterator muon = recoMuons->begin(); muon != recoMuons->end(); muon++) 
	{
		CounterMuons++;
		//std::cout<< " muon->mass(): "<< muon->mass()  << "   muon->pt(): "<< muon->pt()  <<std::endl;

		//---------------------------------------------------------------------------
		//Loose Muon Criteria
		//---------------------------------------------------------------------------
		//As loose muon criteria we only select muons of the muon collection
	
		myLeptons.push_back(*muon); //fill

	}//End  Muon Loop

	//==========================================================================
	// sort the muons of highest pT to lowest 
	//==========================================================================
	std::sort(myLeptons.begin(),myLeptons.end(), [](const reco::Muon &a, const reco::Muon &b)
	{
		return a.pt() > b.pt();
	});

	// Dimuon Selection - Get the two muons with highest pT to make a pair
	//==========================================================================
	if (myLeptons.size() >= 2) 
	{ // Dimuon
			
		if(verbose_) std::cout<<"RECO  Muons Multiplicity:  " << myLeptons.size() << std::endl;
		reco::Muon leadingMuon = myLeptons[0];//first muon - with highest Pt
		reco::Muon trailingMuon = myLeptons[1];//second muon - with the second highest Pt

		//Loretz Vector of the Lepton
		mu_1.SetPtEtaPhiM(leadingMuon.pt(), leadingMuon.eta(), leadingMuon.phi(), leadingMuon.mass());
		mu_2.SetPtEtaPhiM(trailingMuon.pt(), trailingMuon.eta(), trailingMuon.phi(), trailingMuon.mass());
	
		//Unvariant Mass of two Particles
		M = (mu_1+mu_2).M();	

		//Oppossite Charge - The ressonance has charge 0, so the two bodies decay must have oppostite charges
		if ( leadingMuon.charge() == trailingMuon.charge() )	continue;
		nDimuon++;
		
		//Ressonance window
		if ( (M > 2.5) && (M < 3.6) )
		{
			std::cout<< "----------------------------------"<< std::endl;
			
			//Simple variables. Not vectors.
			//==========================================================================
			//To apply the tight muon criteria
			if( muon::isTightMuon(leadingMuon, vtx) )
			{
				TagAndProbes++;			
				//std::cout<< "TagMuon.pt(): "<< leadingMuon.pt() << " ProbeMuon.pt(): " << trailingMuon.pt()<< std::endl;

				//For tracking Muon Efficience
				if ( trailingMuon.isTrackerMuon() )
				{	
					TagMuon_Pt = leadingMuon.pt();
					TagMuon_Eta = leadingMuon.eta();
					TagMuon_Phi = leadingMuon.phi();
					
					ProbeMuon_Pt = trailingMuon.pt();
					ProbeMuon_Eta = trailingMuon.eta();
					ProbeMuon_Phi = trailingMuon.phi();
					
					MassPassingTagProbeTracking = M;
					PassingProbeTrackingMuon = 1;

					
				}
				else
				{
					//Probe
					MassFailingTagProbeTracking = M;
					PassingProbeTrackingMuon = 0;
					std::cout<< " MassFailingTagProbeTracking: " << M << std::endl;
					
				}

				//------------------------------------------
				//For StandAlone Muon Efficience
				//------------------------------------------
				if (  trailingMuon.isStandAloneMuon() )
				{
            	PassingProbeStandAloneMuon = 1;
				}
				else
				{
					PassingProbeStandAloneMuon = 0;
				}
				//------------------------------------------
				//For ID Muon Efficience
				//------------------------------------------
				if (  trailingMuon.isTrackerMuon() && trailingMuon.isStandAloneMuon() && trailingMuon.isGlobalMuon() )
				{
            	PassingProbeIDMuon = 1;
				}
				else
				{
					PassingProbeIDMuon = 0;
				}

				//std::cout<< "Filling Varaibles: "<< std::endl;
				AnalysisTree->Fill();
				PlotControl->Fill();
			
			}
			//Notice that we fill as opposite way that the IF above
			else if (muon::isTightMuon(trailingMuon, vtx) ) 
			{
				TagAndProbes++;
				//Fill Vectors				
				//std::cout<< "TagMuon.pt(): "<< trailingMuon.pt() << " ProbeMuon.pt(): " << leadingMuon.pt()<< std::endl;
				
				//------------------------------------------
				//For tracking Muon Efficience
				//------------------------------------------
				if ( leadingMuon.isTrackerMuon() )
				{
					TagMuon_Pt = trailingMuon.pt();
					TagMuon_Eta = trailingMuon.eta();
					TagMuon_Phi = trailingMuon.phi();
					
					ProbeMuon_Pt = leadingMuon.pt();
					ProbeMuon_Eta = leadingMuon.eta();
					ProbeMuon_Phi = leadingMuon.phi();
					
					MassPassingTagProbeTracking = M;
					PassingProbeTrackingMuon = 1;
					
				}
				else
				{
					MassFailingTagProbeTracking = M;
					PassingProbeTrackingMuon = 0;
					std::cout<< " MassFailingTagProbeTracking: " << M << std::endl;
					
				}
				//------------------------------------------
				//For StandAlone Muon Efficience
				//------------------------------------------
				if (  leadingMuon.isStandAloneMuon() )
				{
            	PassingProbeStandAloneMuon = 1;
				}
				else
				{
					PassingProbeStandAloneMuon = 0;
				}
				//------------------------------------------
				//For ID Muon Efficience
				//------------------------------------------
				if (  leadingMuon.isPFMuon() )
				{
            	PassingProbeIDMuon = 1;
				}
				else
				{
					PassingProbeIDMuon = 0;
				}
		
				//std::cout<< "Filling Varaibles: "<< M << std::endl;

				AnalysisTree->Fill();
				PlotControl->Fill();
		
			}
			
		}

	}//end Dimuon
	
	}//end //triggerfired

	}// loop trigger names

	//Cleaning Variables
	//std::cout << "Cleaning Variables: " << std::endl;

	/*	
	TagMuon_Pt = 0.;
	TagMuon_Eta = 0.;
	TagMuon_Phi = 0.;
	
	ProbeMuon_Pt = 0.;
	ProbeMuon_Eta = 0.;
	ProbeMuon_Phi = 0.;
	
	PassingProbeTrackingMuon = 0;
	MassPassingTagProbeTracking = 0;

	PassingProbeStandAloneMuon = 0;
	PassingProbeIDMuon = 0;
	*/		

}//end analizer

// ------------ method called once each job just before starting event loop  ------------
	void
JpsiAnalyzerOpen2011::beginJob()
{
	AnalysisTree->Branch("ProbeMuon_Pt", &ProbeMuon_Pt, "ProbeMuon_Pt/D");
	AnalysisTree->Branch("ProbeMuon_Eta", &ProbeMuon_Eta, "ProbeMuon_Eta/D");

	AnalysisTree->Branch("PassingProbeTrackingMuon", &PassingProbeTrackingMuon, "PassingProbeTrackingMuon/I");
	AnalysisTree->Branch("MassPassingTagProbeTracking", &MassPassingTagProbeTracking, "MassPassingTagProbeTracking/D");
	AnalysisTree->Branch("MassFailingTagProbeTracking", &MassFailingTagProbeTracking, "MassFailingTagProbeTracking/D");

	AnalysisTree->Branch("PassingProbeStandAloneMuon", &PassingProbeStandAloneMuon, "PassingProbeStandAloneMuon/I");

	AnalysisTree->Branch("PassingProbeIDMuon", &PassingProbeIDMuon, "PassingProbeIDMuon/I");
	
	//Create and Fill Branchs
	PlotControl->Branch("TagMuon_Pt", &TagMuon_Pt, "TagMuon_Pt/D");
	PlotControl->Branch("TagMuon_Eta", &TagMuon_Eta, "TagMuon_Eta/D");
	PlotControl->Branch("TagMuon_Phi", &TagMuon_Phi, "TagMuon_Phi/D");
	
	PlotControl->Branch("ProbeMuon_Pt", &ProbeMuon_Pt, "ProbeMuon_Pt/D");
	PlotControl->Branch("ProbeMuon_Eta", &ProbeMuon_Eta, "ProbeMuon_Eta/D");
	PlotControl->Branch("ProbeMuon_Phi", &ProbeMuon_Phi, "ProbeMuon_Phi/D");
	
	PlotControl->Branch("PassingProbeTrackingMuon", &PassingProbeTrackingMuon, "PassingProbeTrackingMuon/I");
	PlotControl->Branch("MassPassingTagProbeTracking", &MassPassingTagProbeTracking, "MassPassingTagProbeTracking/D");
	PlotControl->Branch("MassFailingTagProbeTracking", &MassFailingTagProbeTracking, "MassFailingTagProbeTracking/D");

   //Trigger Info Branches
	PlotControl->Branch("ievt", &ievt, "ievt/I");
	PlotControl->Branch("irun", &irun, "irun/I");
	PlotControl->Branch("lumiblock", &lumiBlock, "lumiblock/I" );
	PlotControl->Branch("Total_Events", &Total_Events, "Total_Events/I");
	PlotControl->Branch("HLTriggerName", &NameTrigger);
	PlotControl->Branch("TriggerAcceptedID", &countInAccepted, "countInAccepted/I");
	PlotControl->Branch("TriggeredFiredID", &countInTriggered, "countInTriggered/I");    
	PlotControl->Branch("Muons", &CounterMuons, "CounterMuons/I");
	PlotControl->Branch("nDimuon", &nDimuon, "nDimuon/I");

}

// ------------ method called once each job just after ending the event loop  ------------
//
	void
JpsiAnalyzerOpen2011::endJob()
{
	//Fill the Trees
	//AnalysisTree->Fill();
 	//PlotControl->Fill();

	std::cout << " " << std::endl;
	std::cout << "*********************************************************************** " << std::endl;
	std::cout << "Total_Events: " << Total_Events << std::endl;  
    if(triggerflag_)
	{ 
    	for (unsigned int i=0; i<triggerName_.size(); i++) 
		{	
        	std::cout<<"Paths: " <<  triggerName_[i] << std::endl;
        }
        std::cout<<"N of Evts using Trigger Fired :"<< countInTriggered << std::endl;
   }
	std::cout << "Muons Multiplicity: " << CounterMuons << std::endl;
	std::cout << "nDimuon: " << nDimuon << std::endl;
	std::cout << "TagAndProbes: " << TagAndProbes << std::endl;
	std::cout << "*********************************************************************** " << std::endl;
       //}

}

// ------------ method called when starting to processes a run  ------------
	void
JpsiAnalyzerOpen2011::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void
JpsiAnalyzerOpen2011::endRun(edm::Run const&, edm::EventSetup const&)
{
}
//define this as a plug-in
DEFINE_FWK_MODULE(JpsiAnalyzerOpen2011);
