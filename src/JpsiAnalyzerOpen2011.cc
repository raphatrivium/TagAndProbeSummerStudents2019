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

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Muon Reco
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
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

		TTree* AnalysisTree;
		//TBranch* branch;

		//Histogrms
		//TH1F* leadingMuon_Pt ;


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
		int GlobalMuon = 0;
		int TMOneStationTight = 0;
		int	NumberOfValidMuonHits = 0;
		int pixelLayersWithMeasurement = 0;
		int	normalizedChi2 = 0;
		int db_dz = 0;
		int PFMuon = 0; 
		int TrackerGlobalMuon = 0;
		int nDimuon = 0;

		//Trigger
		int countInAccepted = 0;
		int countInTriggered = 0 ;

		//Creating Vectors
		std::vector<int> VectorEvent;
		std::vector<int> VectorRun;
		std::vector<int> VectorlumiBlock;

		std::vector<double> VectorMuon_Pt;
		std::vector<double>	VectorMuon_Eta;
		std::vector<double>	VectorMuon_Phi;
		std::vector<int>	VectorMuon_Charge;
		std::vector<double>	VectorMuon_Mass;

		std::vector<double> TrackerMuonPt;
		std::vector<double> TrackerMuonEta;
		std::vector<double> TrackerMuonPhi;
		std::vector<int> TrackerMuonCharge;

		std::vector<double> GlobalMuonPt;
		std::vector<double> GlobalMuonEta;
		std::vector<double> GlobalMuonPhi;
		std::vector<int> GlobalMuonCharge;

		std::vector<double> VectorMuonTight_Pt;
		std::vector<double> VectorMuonTight_Eta;
		std::vector<double> VectorMuonTight_Phi;
		std::vector<int> VectorMuonTight_Charge;
		std::vector<double> VectorMuonTight_Mass;

		std::vector<double> VectorMuonTightValidHits_Pt;
		std::vector<double> VectorMuonTightValidHits_Eta;
		std::vector<double> VectorMuonTightValidHits_Phi;
		std::vector<int> VectorMuonTightValidHits_Charge;
		std::vector<double> VectorMuonTightValidHits_Mass;

		std::vector<double> VectorMuonTightValidHitsPixelLayer_Pt;
		std::vector<double> VectorMuonTightValidHitsPixelLayer_Eta;
		std::vector<double> VectorMuonTightValidHitsPixelLayer_Phi;
		std::vector<int> VectorMuonTightValidHitsPixelLayer_Charge;
		std::vector<double> VectorMuonTightValidHitsPixelLayer_Mass;

		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2_Pt;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2_Eta;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2_Phi;
		std::vector<int> VectorMuonTightValidHitsPixelLayerChi2_Charge;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2_Mass;

		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi;
		std::vector<int> VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass;

		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi;
		std::vector<int> VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass;

		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi;
		std::vector<int> VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge;
		std::vector<double> VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass;

		std::vector<double> VectorleadingMuon_Pt;
		std::vector<double> VectorleadingMuon_Eta;
		std::vector<double> VectorleadingMuon_Phi;
		std::vector<int> VectorleadingMuon_Charge;
		std::vector<double> VectorleadingMuon_Mass;

		std::vector<double> VectortrailingMuon_Pt;
		std::vector<double> VectortrailingMuon_Eta;
		std::vector<double> VectortrailingMuon_Phi;
		std::vector<int> VectortrailingMuon_Charge;
		std::vector<double> VectortrailingMuon_Mass;

		std::vector<double> VectorMll;
		std::vector<double> VectorMllpT;
		std::vector<double> VectorMlleta;
		std::vector<double> VectorMllphi;

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

}

//------------------------------------------------------------------------------------------------------------
JpsiAnalyzerOpen2011::~JpsiAnalyzerOpen2011()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
    std::vector<std::string>().swap(NameTrigger);
	std::vector<double>().swap(VectorMuon_Pt);
	std::vector<double>().swap(VectorMuon_Eta);
	std::vector<double>().swap(VectorMuon_Phi);
	std::vector<int>().swap(VectorMuon_Charge);
	std::vector<double>().swap(VectorMuon_Mass);

	std::vector<double>().swap(TrackerMuonPt);
	std::vector<double>().swap(TrackerMuonEta);
	std::vector<double>().swap(TrackerMuonPhi);
	std::vector<int>().swap(TrackerMuonCharge);

	std::vector<double>().swap(GlobalMuonPt);
	std::vector<double>().swap(GlobalMuonEta);
	std::vector<double>().swap(GlobalMuonPhi);
	std::vector<int>().swap(GlobalMuonCharge);

	std::vector<double>().swap(VectorMuonTight_Pt);
	std::vector<double>().swap(VectorMuonTight_Eta);
	std::vector<double>().swap(VectorMuonTight_Phi);
	std::vector<int>().swap(VectorMuonTight_Charge);
	std::vector<double>().swap(VectorMuonTight_Mass);

	std::vector<double>().swap(VectorMuonTightValidHits_Pt);
	std::vector<double>().swap(VectorMuonTightValidHits_Eta);
	std::vector<double>().swap(VectorMuonTightValidHits_Phi);
	std::vector<int>().swap(VectorMuonTightValidHits_Charge);
	std::vector<double>().swap(VectorMuonTightValidHits_Mass);

	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayer_Pt);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayer_Eta);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayer_Phi);
	std::vector<int>().swap(VectorMuonTightValidHitsPixelLayer_Charge);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayer_Mass);

	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2_Pt);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2_Eta);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2_Phi);
	std::vector<int>().swap(VectorMuonTightValidHitsPixelLayerChi2_Charge);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2_Mass);	

	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi);
	std::vector<int>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass);	

	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi);
	std::vector<int>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass);

	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi);
	std::vector<int>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge);
	std::vector<double>().swap(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass);

	std::vector<double>().swap(VectorleadingMuon_Pt);
	std::vector<double>().swap(VectorleadingMuon_Eta);
	std::vector<double>().swap(VectorleadingMuon_Phi);
	std::vector<int>().swap(VectorleadingMuon_Charge);
	std::vector<double>().swap(VectorleadingMuon_Mass);

	std::vector<double>().swap(VectortrailingMuon_Pt);
	std::vector<double>().swap(VectortrailingMuon_Eta);
	std::vector<double>().swap(VectortrailingMuon_Phi);
	std::vector<int>().swap(VectortrailingMuon_Charge);
	std::vector<double>().swap(VectortrailingMuon_Mass);

	std::vector<double>().swap(VectorMll);
	std::vector<double>().swap(VectorMllpT);
	std::vector<double>().swap(VectorMlleta);
	std::vector<double>().swap(VectorMllphi);

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
	const reco::Vertex* vertex = 0;
	edm::Handle< reco::VertexCollection > recoVertices;
	iEvent.getByLabel(primaryVertexProducer_, recoVertices);
	int nVertices = recoVertices->size();
	if (verbose_) std::cout << "nVertices "<< nVertices << std::endl;
	if(nVertices>0) vertex = & ((*recoVertices)[0]);

	edm::Handle< reco::MuonCollection > recoMuons;
	iEvent.getByLabel(recoMuons_, recoMuons);
        
    Total_Events++;

	//------------------------------------------------------------------
	//if(triggerflag_){

    std::cout << "Using Trigger selection "<< std::endl;

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
	//  reco::Vertex vertex;
	std::vector<reco::Muon> myLeptons;
	// Reco Muons
	for (reco::MuonCollection::const_iterator muon = recoMuons->begin(); muon != recoMuons->end(); muon++) 
	{
		CounterMuons++;
		VectorMuon_Pt.push_back(muon->pt());
		VectorMuon_Eta.push_back(muon->eta());
		VectorMuon_Phi.push_back(muon->phi());
		VectorMuon_Charge.push_back(muon->charge());
		VectorMuon_Mass.push_back(muon->mass());

		//HistoMuon_Pt->Fill(muon->pt()) ;           
			
		if (muon->isTrackerMuon())
		{
			TrackerMuon++;
			TrackerMuonPt.push_back(muon->pt());
			TrackerMuonEta.push_back(muon->eta());
			TrackerMuonPhi.push_back(muon->phi());
			TrackerMuonCharge.push_back(muon->charge());
		}

		if (muon->isGlobalMuon())
		{
			GlobalMuon++;
			GlobalMuonPt.push_back(muon->pt());
			GlobalMuonEta.push_back(muon->eta());
			GlobalMuonPhi.push_back(muon->phi());
			GlobalMuonCharge.push_back(muon->charge());	
		}

			bool isTMOneStationTight = muon::isGoodMuon(*muon,muon::TMOneStationTight);	
			if (isTMOneStationTight) 
			{
				TMOneStationTight++;
				VectorMuonTight_Pt.push_back(muon->pt());
				VectorMuonTight_Eta.push_back(muon->eta());
				VectorMuonTight_Phi.push_back(muon->phi());
				VectorMuonTight_Charge.push_back(muon->charge());
				VectorMuonTight_Mass.push_back(muon->mass());

				
	
				if (verbose_) std::cout<< " dxy: "<< fabs(muon->innerTrack()->dxy(vertex->position()))  << std::endl; 
				if (verbose_) std::cout<< " dz: "<< fabs(muon->innerTrack()->dz(vertex->position()))  << std::endl;
				// if (verbose_) std::cout << muon->charge() << std::endl;

				if(muon->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10)
				{
					NumberOfValidMuonHits++;
					VectorMuonTightValidHits_Pt.push_back(muon->pt());
					VectorMuonTightValidHits_Eta.push_back(muon->eta());
					VectorMuonTightValidHits_Phi.push_back(muon->phi());
					VectorMuonTightValidHits_Charge.push_back(muon->charge());
					VectorMuonTightValidHits_Mass.push_back(muon->mass());

					


					if(muon->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 1)
					{
						pixelLayersWithMeasurement++;
						VectorMuonTightValidHitsPixelLayer_Pt.push_back(muon->pt());
						VectorMuonTightValidHitsPixelLayer_Eta.push_back(muon->eta());
						VectorMuonTightValidHitsPixelLayer_Phi.push_back(muon->phi());
						VectorMuonTightValidHitsPixelLayer_Charge.push_back(muon->charge());
						VectorMuonTightValidHitsPixelLayer_Mass.push_back(muon->mass());

						

						if(muon->innerTrack()->normalizedChi2() < 1.8)
						{
							normalizedChi2++;
							VectorMuonTightValidHitsPixelLayerChi2_Pt.push_back(muon->pt());
							VectorMuonTightValidHitsPixelLayerChi2_Eta.push_back(muon->eta());
							VectorMuonTightValidHitsPixelLayerChi2_Phi.push_back(muon->phi());
							VectorMuonTightValidHitsPixelLayerChi2_Charge.push_back(muon->charge());
							VectorMuonTightValidHitsPixelLayerChi2_Mass.push_back(muon->mass());

							

							if( (fabs(muon->innerTrack()->dxy(vertex->position())) < 3./*cm*/) && (fabs(muon->innerTrack()->dz(vertex->position())) < 15./*cm*/) )
							{
								db_dz++;
								VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt.push_back(muon->pt());
								VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta.push_back(muon->eta());
								VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi.push_back(muon->phi());
								VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge.push_back(muon->charge());
								VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass.push_back(muon->mass());

								
							
								if (muon->isPFMuon())
								{ //Particle-Flow muon id - Can be complemented by muon quality cuts similar to those used in the Tight Muon selection.
									PFMuon++;
									VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt.push_back(muon->pt());
									VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta.push_back(muon->eta());
									VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi.push_back(muon->phi());
									VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge.push_back(muon->charge());
									VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass.push_back(muon->mass());

									
								
									if (muon->isTrackerMuon() || muon->isGlobalMuon())//muon type selection
									{
										TrackerGlobalMuon++;
										VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt.push_back(muon->pt());
										VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta.push_back(muon->eta());
										VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi.push_back(muon->phi());
										VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge.push_back(muon->charge());
										VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass.push_back(muon->mass());

										
									
										myLeptons.push_back(*muon); //fill

									}  //end muon type selection					
								} //PF muon
							} //end db_dz
						} //end normalizedChi2
					}//pixelLayersWithMeasurement
				} //end numberOfValidTrackerHits         
			}//end TMOneStationTight 
	} //Muon Loop
	//}// Trigger selection

	std::sort(myLeptons.begin(),myLeptons.end(), [](const reco::Muon &a, const reco::Muon &b)
	{
		return a.pt() > b.pt();
	});

	// Print Reco Muons
	for (std::vector<reco::Muon>::const_iterator muon = myLeptons.begin(); muon != myLeptons.end(); muon++) 
	{
		if(verbose_) std::cout<<"RECO myLeptons->pt():  " << muon->pt() << std::endl;
	}// Muon loop

	if(verbose_) std::cout<<"RECO myLeptons.size() all  " << myLeptons.size() << std::endl;


	// dimuon selection
	if (myLeptons.size() >= 2) 
	{ //loop ndimuon
		
		nDimuon++;
		if(verbose_) std::cout<<"RECO  Muons Multiplicity:  " << myLeptons.size() << std::endl;
		reco::Muon leadingMuon = myLeptons[0];
		reco::Muon trailingMuon = myLeptons[1];
		
		if(verbose_) std::cout<< "Leading Muon pt, eta, phi, charge = " << leadingMuon.pt() << " "<< leadingMuon.eta() << " "<< leadingMuon.phi() << " " << leadingMuon.charge() << std::endl;
		if(verbose_) std::cout<< "Trailing Muon  pt, eta, phi,charge = " << trailingMuon.pt() << " " << trailingMuon.eta() << " " << trailingMuon.phi() << " " << trailingMuon.charge()<< std::endl;

		//Invariant mass of dimuons
		double Mll = (leadingMuon.p4() + trailingMuon.p4()).mass();
		double MllpT = (leadingMuon.p4() + trailingMuon.p4()).pt();
		double Mlleta = (leadingMuon.p4() + trailingMuon.p4()).eta();
		double Mllphi = (leadingMuon.p4() + trailingMuon.p4()).phi();
		if(verbose_) std::cout<< "Dimuons Invariant Mass Mll, pT, eta, phi: " << Mll << " " << MllpT << " " << Mlleta << " " << Mllphi << std::endl;

		//Fill Histograms
		//leadingMuon_Pt->Fill(leadingMuon.pt());
		
		//trailingMuon_Pt->Fill(trailingMuon.pt());
		
		//Fill Vectors
		VectorleadingMuon_Pt.push_back(leadingMuon.pt());
		VectorleadingMuon_Eta.push_back(leadingMuon.eta());
		VectorleadingMuon_Phi.push_back(leadingMuon.phi());
		VectorleadingMuon_Charge.push_back(leadingMuon.charge());
		VectorleadingMuon_Mass.push_back(leadingMuon.mass());

		VectortrailingMuon_Pt.push_back(trailingMuon.pt());
		VectortrailingMuon_Eta.push_back(trailingMuon.eta());
		VectortrailingMuon_Phi.push_back(trailingMuon.phi());
		VectortrailingMuon_Charge.push_back(trailingMuon.charge());
		VectortrailingMuon_Mass.push_back(trailingMuon.mass());

		VectorMll.push_back(Mll);
		VectorMllpT.push_back(MllpT);
		VectorMlleta.push_back(Mlleta);
		VectorMllphi.push_back(Mllphi);

	}//end loop ndimuons
	
	}//end //triggerfired

	}// loop trigger names

}//end analize


// ------------ method called once each job just before starting event loop  ------------
	void
JpsiAnalyzerOpen2011::beginJob()
{
	//Create and Fill Branchs
	AnalysisTree->Branch("ievt", &ievt, "ievt/I");
	AnalysisTree->Branch("irun", &irun, "irun/I");
	AnalysisTree->Branch("lumiblock", &lumiBlock, "lumiblock/I" );
	AnalysisTree->Branch("Total_Events", &Total_Events, "Total_Events/I");
        
    // Trigger Info Branches
	AnalysisTree->Branch("HLTriggerName", &NameTrigger);
	AnalysisTree->Branch("TriggerAcceptedID", &countInAccepted, "countInAccepted/I");
	AnalysisTree->Branch("TriggeredFiredID", &countInTriggered, "countInTriggered/I");
     
	AnalysisTree->Branch("Muons", &CounterMuons, "CounterMuons/I");
	AnalysisTree->Branch("TrackerMuon", &TrackerMuon, "TrackerMuon/I");
	AnalysisTree->Branch("GlobalMuon", &GlobalMuon, "GlobalMuon/I");

	AnalysisTree->Branch("VectorMuon_Pt","std::vector<double>", &VectorMuon_Pt);
	AnalysisTree->Branch("VectorMuon_Eta","std::vector<double>", &VectorMuon_Eta);
	AnalysisTree->Branch("VectorMuon_Phi","std::vector<double>", &VectorMuon_Phi);
	AnalysisTree->Branch("VectorMuon_Charge","std::vector<int>", &VectorMuon_Charge);
	AnalysisTree->Branch("VectorMuon_Mass","std::vector<double>", &VectorMuon_Mass);

	AnalysisTree->Branch("TMOneStationTight", &TMOneStationTight, "TMOneStationTight/I");

	AnalysisTree->Branch("VectorMuonTight_Pt","std::vector<double>", &VectorMuonTight_Pt);
	AnalysisTree->Branch("VectorMuonTight_Eta","std::vector<double>", &VectorMuonTight_Eta);
	AnalysisTree->Branch("VectorMuonTight_Phi","std::vector<double>", &VectorMuonTight_Phi);
	AnalysisTree->Branch("VectorMuonTight_Charge","std::vector<int>", &VectorMuonTight_Charge);
	AnalysisTree->Branch("VectorMuonTight_Mass","std::vector<double>", &VectorMuonTight_Mass);

	AnalysisTree->Branch("NumberOfValidMuonHits", &NumberOfValidMuonHits, "NumberOfValidMuonHits/I");

	AnalysisTree->Branch("VectorMuonTightValidHits_Pt","std::vector<double>", &VectorMuonTightValidHits_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHits_Eta","std::vector<double>", &VectorMuonTightValidHits_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHits_Phi","std::vector<double>", &VectorMuonTightValidHits_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHits_Charge","std::vector<int>", &VectorMuonTightValidHits_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHits_Mass","std::vector<double>", &VectorMuonTightValidHits_Mass);

	AnalysisTree->Branch("pixelLayersWithMeasurement", &pixelLayersWithMeasurement, "pixelLayersWithMeasurement/I");

	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayer_Pt","std::vector<double>", &VectorMuonTightValidHitsPixelLayer_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayer_Eta","std::vector<double>", &VectorMuonTightValidHitsPixelLayer_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayer_Phi","std::vector<double>", &VectorMuonTightValidHitsPixelLayer_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayer_Charge","std::vector<int>", &VectorMuonTightValidHitsPixelLayer_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayer_Mass","std::vector<double>", &VectorMuonTightValidHitsPixelLayer_Mass);

	AnalysisTree->Branch("normalizedChi2", &normalizedChi2, "normalizedChi2/I");

	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2_Pt","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2_Eta","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2_Phi","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2_Charge","std::vector<int>", &VectorMuonTightValidHitsPixelLayerChi2_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2_Mass","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2_Mass);

	AnalysisTree->Branch("db_dz", &db_dz, "db_dz/I");

	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge","std::vector<int>", &VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass);

	AnalysisTree->Branch("PFMuon", &PFMuon, "PFMuon/I");

	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge","std::vector<int>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass);

	AnalysisTree->Branch("TrackerGlobalMuon", &TrackerGlobalMuon, "TrackerGlobalMuon/I");

	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge","std::vector<int>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge);
	AnalysisTree->Branch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass","std::vector<double>", &VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass);

	AnalysisTree->Branch("nDimuon", &nDimuon, "nDimuon/I");

	AnalysisTree->Branch("VectorEvent","std::vector<int>", &VectorEvent);
	AnalysisTree->Branch("VectorRun","std::vector<int>", &VectorRun);
	AnalysisTree->Branch("VectorlumiBlock","std::vector<int>", &VectorlumiBlock);


	AnalysisTree->Branch("leadingMuon_Pt","std::vector<double>", &VectorleadingMuon_Pt);
	AnalysisTree->Branch("leadingMuon_Eta","std::vector<double>", &VectorleadingMuon_Eta);
	AnalysisTree->Branch("leadingMuon_Phi","std::vector<double>", &VectorleadingMuon_Phi);
	AnalysisTree->Branch("leadingMuon_Charge","std::vector<int>", &VectorleadingMuon_Charge);
	AnalysisTree->Branch("leadingMuon_Mass","std::vector<double>", &VectorleadingMuon_Mass);

	AnalysisTree->Branch("trailingMuon_Pt","std::vector<double>", &VectortrailingMuon_Pt);
	AnalysisTree->Branch("trailingMuon_Eta","std::vector<double>", &VectortrailingMuon_Eta);
	AnalysisTree->Branch("trailingMuon_Phi","std::vector<double>", &VectortrailingMuon_Phi);
	AnalysisTree->Branch("trailingMuon_Charge","std::vector<int>", &VectortrailingMuon_Charge);
	AnalysisTree->Branch("trailingMuon_Mass","std::vector<double>", &VectortrailingMuon_Mass);

	AnalysisTree->Branch("Mll","std::vector<double>", &VectorMll);
	AnalysisTree->Branch("MllpT","std::vector<double>", &VectorMllpT);
	AnalysisTree->Branch("Mlleta","std::vector<double>", &VectorMlleta);
	AnalysisTree->Branch("Mllphi","std::vector<double>", &VectorMllphi);

}

// ------------ method called once each job just after ending the event loop  ------------nalysisTree->Branch("VectorMll","std::vector<double>", &VectorMll);
//
	void
JpsiAnalyzerOpen2011::endJob()
{
	//Fill the Trees
	AnalysisTree->Fill();
 
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
	std::cout << "TMOneStationTight: "<< TMOneStationTight<<std::endl;
	std::cout << "NumberOfValidMuonHits: " << NumberOfValidMuonHits << std::endl;
	std::cout << "pixelLayersWithMeasurement > 1: " << pixelLayersWithMeasurement << std::endl;
	std::cout << "normalizedChi2: "<< normalizedChi2<<std::endl;
	std::cout << "db < 3cm and dz < 15cm: "<< db_dz<<std::endl;
	std::cout << "PFMuon: " << PFMuon << std::endl;
	std::cout << "TrackerGlobalMuon: " << TrackerGlobalMuon << std::endl;
	std::cout << "nDimuon: " << nDimuon << std::endl;
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
