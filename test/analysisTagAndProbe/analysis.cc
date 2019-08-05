//Program to extract data (variables and vector format) of the root file and create histograms. This root file is made of others root files, therefore, has multiple entries which we access with a diferent method than a normal root file.

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atof */
#include <math.h>       /* sin */
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <vector>
#ifndef ROOT_TLatex
#define ROOT_TLatex
#ifndef ROOTiosfwd
#include "Riosfwd.h"
#endif
#ifndef ROOT_TText
#include "TText.h"
#endif
#ifndef ROOT_TAttLine
#include "TAttLine.h"
#endif

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int analysis()
{
	//call a file for a histogram style (Optional)
	gROOT->LoadMacro("styleTDR.C"); 
	setTDRStyle();

	//counters for file f1 (pythia)
	int count_Total_Events_pythia = 0;
	int counter_Muon_pythia = 0; 
	int counter_PtTrigger7_pythia = 0;
	int counter_TMOneStationTight_pythia = 0;
	int counter_NumberOfValidMuonHits_pythia = 0;
	int counter_pixelLayersWithMeasurement_pythia = 0;
	int counter_normalizedChi2_pythia = 0;
	int counter_db_dz_pythia = 0;
	int counter_PFMuon_pythia = 0;
	int counter_TrackerGlobalMuon_pythia = 0;
	int counter_nDimuon_pythia = 0;
	int count_OppositeCharge_pythia = 0; 
	int count_Eta_pythia = 0;
	int count_Jpsi_pythia = 0; 
	int count_region1_pythia = 0;
	int count_region2_pythia = 0; 
	int count_region3_pythia = 0;
	int count_region4_pythia = 0;
	int count_region5_pythia = 0;
	int count_pico_massa = 0;

	//counters for file f2 (data)
	int count_Total_Events_data = 0;
	int counter_Muon_data = 0;
	int counter_TMOneStationTight_data = 0;
	int counter_NumberOfValidMuonHits_data = 0;
	int counter_pixelLayersWithMeasurement_data = 0;
	int counter_normalizedChi2_data = 0;
	int counter_db_dz_data = 0;
	int counter_PFMuon_data = 0;
	int counter_TrackerGlobalMuon_data = 0;
	int counter_nDimuon_data = 0;
	int count_OppositeCharge_data = 0; 
	int count_Eta_data = 0;
	int count_Jpsi_data = 0; 
	int count_region1_data = 0;
	int count_region2_data = 0; 
	int count_region3_data = 0;
	int count_region4_data = 0;
	int count_region5_data = 0;    
		
	//Lorentz Vector
	TLorentzVector mu_1;
	TLorentzVector mu_2;	
	
	//Variaveis	and Vectors
	int Total_Events = 0;
	int Muons = 0;
	int PtTrigger7 = 0;
	int TMOneStationTight = 0;
	int NumberOfValidMuonHits = 0;
	int pixelLayersWithMeasurement = 0;
	int normalizedChi2 = 0;
	int db_dz = 0;
	int PFMuon = 0;
	int TrackerGlobalMuon = 0;
	int nDimuon = 0;

	std::vector<double>* VectorMuon_Pt = 0.;
	std::vector<double>* VectorMuon_Eta = 0.;
	std::vector<double>* VectorMuon_Phi = 0.;
	std::vector<int>* VectorMuon_Charge = 0.;
	std::vector<double>* VectorMuon_Mass = 0.;

	std::vector<double>* TrackerMuonPt = 0.;
	std::vector<double>* TrackerMuonEta = 0.;
	std::vector<double>* TrackerMuonPhi = 0.;
	std::vector<int>* TrackerMuonCharge = 0.;

	std::vector<double>* GlobalMuonPt = 0.;
	std::vector<double>* GlobalMuonEta = 0.;
	std::vector<double>* GlobalMuonPhi = 0.;
	std::vector<int>* GlobalMuonCharge = 0.;

	std::vector<double>* VectorMuonTight_Pt = 0.;
	std::vector<double>* VectorMuonTight_Eta = 0.;
	std::vector<double>* VectorMuonTight_Phi = 0.;
	std::vector<int>* VectorMuonTight_Charge = 0.;
	std::vector<double>* VectorMuonTight_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHits_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHits_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHits_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHits_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHits_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHitsPixelLayer_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayer_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayer_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHitsPixelLayer_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayer_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHitsPixelLayerChi2_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDz_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDz_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHitsPixelLayerChi2DbDz_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDz_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Mass = 0.;

	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi = 0.;
	std::vector<int>* VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge = 0.;
	std::vector<double>* VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Mass = 0.;
		
	std::vector<double>* leadingMuon_Pt = 0.;
	std::vector<double>* leadingMuon_Eta = 0.;
	std::vector<double>* leadingMuon_Phi = 0.;
	std::vector<int>* leadingMuon_Charge = 0.;
	std::vector<double>* leadingMuon_Mass = 0.;

	std::vector<double>* trailingMuon_Pt = 0.;
	std::vector<double>* trailingMuon_Eta = 0.;
	std::vector<double>* trailingMuon_Phi = 0.;
	std::vector<int>* trailingMuon_Charge = 0.;
	std::vector<double>* trailingMuon_Mass = 0.;

	int Total_Events2 = 0;
	int Muons2 = 0;
	int TMOneStationTight2 = 0;
	int NumberOfValidMuonHits2 = 0;
	int pixelLayersWithMeasurement2 = 0;
	int normalizedChi22 = 0;
	int db_dz2 = 0;
	int PFMuon2 = 0;
	int TrackerGlobalMuon2 = 0;
	int nDimuon2 = 0;

	std::vector<double>* leadingMuon_Pt2 = 0.;
	std::vector<double>* leadingMuon_Eta2 = 0.;
	std::vector<double>* leadingMuon_Phi2 = 0.;
	std::vector<int>* leadingMuon_Charge2 = 0.;
	std::vector<double>* leadingMuon_Mass2 = 0.;

	std::vector<double>* trailingMuon_Pt2 = 0.;
	std::vector<double>* trailingMuon_Eta2 = 0.;
	std::vector<double>* trailingMuon_Phi2 = 0.;
	std::vector<int>* trailingMuon_Charge2 = 0.;
	std::vector<double>* trailingMuon_Mass2 = 0.;
	
	std::vector<double>* VectorMll = 0.;
	std::vector<double>* VectorMllpT = 0.;
	std::vector<double>* VectorMlleta = 0.;
	std::vector<double>* VectorMllphi = 0.;
	
	double M = 0.;
	double Pt = 0.;
	double Eta = 0.;
	double Rapidity = 0.;
	
	//***********************************************************	
	//Criacao dos histogramas

	//Histogramas cinematics quantities of the muons
	TH1F *h_Muon_Pt = new TH1F("h_Muon_Pt","h_Muon_Pt",100,0,50);
	h_Muon_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_Muon_Pt->SetName("h_Muon_Pt");
	TH1F *h_Muon_Eta = new TH1F("h_Muon_Eta","h_Muon_Eta",100,-4,4);
	h_Muon_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_Muon_Eta->SetName("h_Muon_Eta");
	TH1F *h_Muon_Phi = new TH1F("h_Muon_Phi","h_Muon_Phi",100,-4,4);
	h_Muon_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_Muon_Phi->SetName("h_Muon_Phi");
	TH1F *h_Muon_Charge = new TH1F("h_Muon_Charge","h_Muon_Charge",100,-2,2);
	h_Muon_Charge->SetTitle("Distribuicao das cargas Muons; Charge  ; Eventos ");
	h_Muon_Charge->SetName("h_Muon_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTight_Pt = new TH1F("h_MuonTight_Pt","h_MuonTight_Pt",100,0,50);
	h_MuonTight_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTight_Pt->SetName("h_MuonTight_Pt");
	TH1F *h_MuonTight_Eta = new TH1F("h_MuonTight_Eta","h_MuonTight_Eta",100,-4,4);
	h_MuonTight_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTight_Eta->SetName("h_MuonTight_Eta");
	TH1F *h_MuonTight_Phi = new TH1F("h_MuonTight_Phi","h_MuonTight_Phi",100,-4,4);
	h_MuonTight_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTight_Phi->SetName("h_MuonTight_Phi");
	TH1F *h_MuonTight_Charge = new TH1F("h_MuonTight_Charge","h_MuonTight_Charge",100,-2,2);
	h_MuonTight_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTight_Charge->SetName("h_MuonTight_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHits_Pt = new TH1F("h_MuonTightValidHits_Pt","h_MuonTightValidHits_Pt",100,0,50);
	h_MuonTightValidHits_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHits_Pt->SetName("h_MuonTightValidHits_Pt");
	TH1F *h_MuonTightValidHits_Eta = new TH1F("h_MuonTightValidHits_Eta","h_MuonTightValidHits_Eta",100,-4,4);
	h_MuonTightValidHits_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHits_Eta->SetName("h_MuonTightValidHits_Eta");
	TH1F *h_MuonTightValidHits_Phi = new TH1F("h_MuonTightValidHits_Phi","h_MuonTightValidHits_Phi",100,-4,4);
	h_MuonTightValidHits_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHits_Phi->SetName("h_MuonTightValidHits_Phi");
	TH1F *h_MuonTightValidHits_Charge = new TH1F("h_MuonTightValidHits_Charge","h_MuonTightValidHits_Charge",100,-2,2);
	h_MuonTightValidHits_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTightValidHits_Charge->SetName("h_MuonTightValidHits_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHitsPixelLayer_Pt = new TH1F("h_MuonTightValidHitsPixelLayer_Pt","h_MuonTightValidHitsPixelLayer_Pt",100,0,50);
	h_MuonTightValidHitsPixelLayer_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHitsPixelLayer_Pt->SetName("h_MuonTightValidHitsPixelLayer_Pt");
	TH1F *h_MuonTightValidHitsPixelLayer_Eta = new TH1F("h_MuonTightValidHitsPixelLayer_Eta","h_MuonTightValidHitsPixelLayer_Eta",100,-4,4);
	h_MuonTightValidHitsPixelLayer_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHitsPixelLayer_Eta->SetName("h_MuonTightValidHitsPixelLayer_Eta");
	TH1F *h_MuonTightValidHitsPixelLayer_Phi = new TH1F("h_MuonTightValidHitsPixelLayer_Phi","h_MuonTightValidHitsPixelLayer_Phi",100,-4,4);
	h_MuonTightValidHitsPixelLayer_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHitsPixelLayer_Phi->SetName("h_MuonTightValidHitsPixelLayer_Phi");
	TH1F *h_MuonTightValidHitsPixelLayer_Charge = new TH1F("h_MuonTightValidHitsPixelLayer_Charge","h_MuonTightValidHitsPixelLayer_Charge",100,-2,2);
	h_MuonTightValidHitsPixelLayer_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTightValidHitsPixelLayer_Charge->SetName("h_MuonTightValidHitsPixelLayer_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHitsPixelLayerChi2_Pt = new TH1F("h_MuonTightValidHitsPixelLayerChi2_Pt","h_MuonTightValidHitsPixelLayerChi2_Pt",100,0,50);
	h_MuonTightValidHitsPixelLayerChi2_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2_Pt->SetName("h_MuonTightValidHitsPixelLayerChi2_Pt");
	TH1F *h_MuonTightValidHitsPixelLayerChi2_Eta = new TH1F("h_MuonTightValidHitsPixelLayerChi2_Eta","h_MuonTightValidHitsPixelLayerChi2_Eta",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2_Eta->SetName("h_MuonTightValidHitsPixelLayerChi2_Eta");
	TH1F *h_MuonTightValidHitsPixelLayerChi2_Phi = new TH1F("h_MuonTightValidHitsPixelLayerChi2_Phi","h_MuonTightValidHitsPixelLayerChi2_Phi",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2_Phi->SetName("h_MuonTightValidHitsPixelLayerChi2_Phi");
	TH1F *h_MuonTightValidHitsPixelLayerChi2_Charge = new TH1F("h_MuonTightValidHitsPixelLayerChi2_Charge","h_MuonTightValidHitsPixelLayerChi2_Charge",100,-2,2);
	h_MuonTightValidHitsPixelLayerChi2_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2_Charge->SetName("h_MuonTightValidHitsPixelLayerChi2_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDz_Pt = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDz_Pt","h_MuonTightValidHitsPixelLayerChi2DbDz_Pt",100,0,50);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDz_Pt->SetName("h_MuonTightValidHitsPixelLayerChi2DbDz_Pt");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDz_Eta = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDz_Eta","h_MuonTightValidHitsPixelLayerChi2DbDz_Eta",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDz_Eta->SetName("h_MuonTightValidHitsPixelLayerChi2DbDz_Eta");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDz_Phi = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDz_Phi","h_MuonTightValidHitsPixelLayerChi2DbDz_Phi",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDz_Phi->SetName("h_MuonTightValidHitsPixelLayerChi2DbDz_Phi");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDz_Charge = new TH1F("h_MuonTightValidHitsPixelLayerChi2_Charge","h_MuonTightValidHitsPixelLayerChi2DbDz_Charge",100,-2,2);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Charge->SetTitle("Distribuicao das cargas dos muons; Charge ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDz_Charge->SetName("h_MuonTightValidHitsPixelLayerChi2DbDz_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt","h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt",100,0,50);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta","h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Eta");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi","h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Phi");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge","h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge",100,-2,2);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPf_Charge");

//-----------------------------------------------------------------------------
	
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt","h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt",100,0,50);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->SetTitle("Distribuicao p_{T} dos Muons; pT [GeV] ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta","h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->SetTitle("Distribuicao Eta dos Muons; Eta ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Eta");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi","h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi",100,-4,4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->SetTitle("Distribuicao Phi dos Muons; Phi ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Phi");
	TH1F *h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge = new TH1F("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge","h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge",100,-2,2);
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->SetTitle("Distribuicao das cargas dos Muons; Charge ; Eventos ");
	h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge->SetName("h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Charge");

//-----------------------------------------------------------------------------
//Cinematics quantities of the Dimuons
	TH1F *h_Dimuons_M = new TH1F("h_Dimuons_M","h_Dimuons_M",100,0,10);
	h_Dimuons_M->SetTitle("Distribuicao Massa Invariante dos Dimuons; #mu#mu [GeV] ; Eventos ");
	h_Dimuons_M->SetName("h_Dimuons_M");

	TH1F *h_Dimuons_Pt = new TH1F("h_Dimuons_Pt","h_Dimuons_Pt",100,0,50);
	h_Dimuons_Pt->SetTitle("Distribuicao Pt dos Dimuons; #mu#mu p_{T} [GeV] ; Eventos ");
	h_Dimuons_Pt->SetName("h_Dimuons_Pt");

	TH1F *h_Dimuons_Eta = new TH1F("h_Dimuons_Eta","h_Dimuons_Eta",100,-4,4);
	h_Dimuons_Eta->SetTitle("Distribuicao Pseudo-Rapidez dos Dimuons; #eta ; Eventos ");
	h_Dimuons_Eta->SetName("h_Dimuons_Eta");

	TH1F *h_Dimuons_Rapidity = new TH1F("h_Dimuons_Rapidity","h_Dimuons_Rapidity",100,-4,4);
	h_Dimuons_Rapidity->SetTitle("Distribuicao Rapidez dos Muons; y ; Eventos ");
	h_Dimuons_Rapidity->SetName("h_Dimuons_Rapidity");

//-----------------------------------------------------------------------------

	TH1F *h_DimuonsOppositeCharge_M = new TH1F("h_DimuonsOppositeCharge_Pt","h_DimuonsOppositeCharge_M",100,0,10);
	h_DimuonsOppositeCharge_M->SetTitle("Distribuicao Massa Invariante dos muons; #mu#mu [GeV] ; Eventos ");
	h_DimuonsOppositeCharge_M->SetName("h_DimuonsOppositeCharge_M");

	TH1F *h_DimuonsOppositeCharge_Pt = new TH1F("h_DimuonsOppositeCharge_Pt","h_DimuonsOppositeCharge_Pt",100,0,50);
	h_DimuonsOppositeCharge_Pt->SetTitle("Distribuicao Pt dos muons; #mu#mu p_{T} [GeV] ; Eventos ");
	h_DimuonsOppositeCharge_Pt->SetName("h_DimuonsOppositeCharge_Pt");

	TH1F *h_DimuonsOppositeCharge_Eta = new TH1F("h_DimuonsOppositeCharge_Eta","h_DimuonsOppositeCharge_Eta",100,-4,4);
	h_DimuonsOppositeCharge_Eta->SetTitle("Distribuicao Pseudo-Rapidez dos Dimuons; #eta ; Eventos ");
	h_DimuonsOppositeCharge_Eta->SetName("h_DimuonsOppositeCharge_Eta");

	TH1F *h_DimuonsOppositeCharge_Rapidity = new TH1F("h_DimuonsOppositeCharge_Rapidity","h_DimuonsOppositeCharge_Rapidity",100,-4,4);
	h_DimuonsOppositeCharge_Rapidity->SetTitle("Distribuicao Rapidez dos muons; y ; Eventos ");
	h_DimuonsOppositeCharge_Rapidity->SetName("h_DimuonsOppositeCharge_Rapidity");

//-----------------------------------------------------------------------------

	TH1F *h_DimuonsOppositeChargeEta_M = new TH1F("h_DimuonsOppositeChargeEta_Pt","h_DimuonsOppositeChargeEta_M",100,0,10);
	h_DimuonsOppositeChargeEta_M->SetTitle("Distribuicao Massa Invariante dos muons; #mu#mu [GeV] ; Eventos ");
	h_DimuonsOppositeChargeEta_M->SetName("h_DimuonsOppositeChargeEta_M");

	TH1F *h_DimuonsOppositeChargeEta_Pt = new TH1F("h_DimuonsOppositeChargeEta_Pt","h_DimuonsOppositeChargeEta_Pt",100,0,50);
	h_DimuonsOppositeChargeEta_Pt->SetTitle("Distribuicao Pt dos muons; #mu#mu p_{T} [GeV] ; Eventos ");
	h_DimuonsOppositeChargeEta_Pt->SetName("h_DimuonsOppositeChargeEta_Pt");

	TH1F *h_DimuonsOppositeChargeEta_Eta = new TH1F("h_DimuonsOppositeChargeEta_Eta","h_DimuonsOppositeChargeEta_Eta",100,-4,4);
	h_DimuonsOppositeChargeEta_Eta->SetTitle("Distribuicao Pseudo-Rapidez dos Dimuons; #eta ; Eventos ");
	h_DimuonsOppositeChargeEta_Eta->SetName("h_DimuonsOppositeChargeEta_Eta");

	TH1F *h_DimuonsOppositeChargeEta_Rapidity = new TH1F("h_DimuonsOppositeChargeEta_Rapidity","h_DimuonsOppositeChargeEta_Rapidity",100,-4,4);
	h_DimuonsOppositeChargeEta_Rapidity->SetTitle("Distribuicao Rapidez dos muons; y ; Eventos ");
	h_DimuonsOppositeChargeEta_Rapidity->SetName("h_DimuonsOppositeChargeEta_Rapidity");

//-----------------------------------------------------------------------------
	
	TH1F *h_DimuonsOppositeChargeEtaJpsi_M = new TH1F("h_DimuonsOppositeChargeEtaJpsi_Pt","h_DimuonsOppositeChargeEtaJpsi_M",100,0,10);
	h_DimuonsOppositeChargeEtaJpsi_M->SetTitle("Distribuicao Massa Invariante dos muons; #mu#mu [GeV] ; Eventos ");
	h_DimuonsOppositeChargeEtaJpsi_M->SetName("h_DimuonsOppositeChargeEtaJpsi_M");

	TH1F *h_DimuonsOppositeChargeEtaJpsi_Pt = new TH1F("h_DimuonsOppositeChargeEtaJpsi_Pt","h_DimuonsOppositeChargeEtaJpsi_Pt",100,0,50);
	h_DimuonsOppositeChargeEtaJpsi_Pt->SetTitle("Distribuicao Pt dos muons; #mu#mu p_{T} [GeV] ; Eventos ");
	h_DimuonsOppositeChargeEtaJpsi_Pt->SetName("h_DimuonsOppositeChargeEtaJpsi_Pt");

	TH1F *h_DimuonsOppositeChargeEtaJpsi_Eta = new TH1F("h_DimuonsOppositeChargeEtaJpsi_Eta","h_DimuonsOppositeChargeEtaJpsi_Eta",100,-4,4);
	h_DimuonsOppositeChargeEtaJpsi_Eta->SetTitle("Distribuicao Pseudo-Rapidez dos Dimuons; #eta ; Eventos ");
	h_DimuonsOppositeChargeEtaJpsi_Eta->SetName("h_DimuonsOppositeChargeEtaJpsi_Eta");

	TH1F *h_DimuonsOppositeChargeEtaJpsi_Rapidity = new TH1F("h_DimuonsOppositeChargeEtaJpsi_Rapidity","h_DimuonsOppositeChargeEtaJpsi_Rapidity",100,-4,4);
	h_DimuonsOppositeChargeEtaJpsi_Rapidity->SetTitle("Distribuicao Rapidez dos muons; y  ; Eventos ");
	h_DimuonsOppositeChargeEtaJpsi_Rapidity->SetName("h_DimuonsOppositeChargeEtaJpsi_Rapidity");

//-----------------------------------------------------------------------------
	//TH2F (rapidity x  Transverse Momentum) for J/Psi candidates
	TH2F *h2_Jpsi = new TH2F("h2_Jpsi", "h2_Jpsi", 20, 0, 2.0, 20, 8, 30);
	h2_Jpsi->SetTitle("Distribuicao |y| x p_{T} na Janela J/#psi; |Y| ; p_{T} [GeV] ");
	h2_Jpsi->SetName("p_{T} #times y");

//-----------------------------------------------------------------------------
	//rapidity regions for pythia
	TH1F *h_dimuons_M_y1 = new TH1F("h_dimuons_M_y1","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y1->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y1->SetName("h_dimuons_M_y1");

	TH1F *h_dimuons_M_y2 = new TH1F("h_dimuons_M_y2","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y2->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y2->SetName("h_dimuons_M_y2");
	
	TH1F *h_dimuons_M_y3 = new TH1F("h_dimuons_M_y3","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y3->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y3->SetName("h_dimuons_M_y3");

	TH1F *h_dimuons_M_y4 = new TH1F("h_dimuons_M_y4","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y4->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y4->SetName("h_dimuons_M_y4");

	TH1F *h_dimuons_M_y5 = new TH1F("h_dimuons_M_y5","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y5->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y5->SetName("h_dimuons_M_y5");

	TH1F *h_pico_massa = new TH1F("h_pico_massa","h_dimuons_phi",100,2.8,3.4);
	h_pico_massa->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_pico_massa->SetName("h_pico_massa");

//-----------------------------------------------------------------------------
	//rapidity regions for data
	TH1F *h_dimuons_M_y1_trigger = new TH1F("h_dimuons_M_y1_trigger","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y1_trigger->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y1_trigger->SetName("h_dimuons_M_y1_trigger");

	TH1F *h_dimuons_M_y2_trigger = new TH1F("h_dimuons_M_y2_trigger","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y2_trigger->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y2_trigger->SetName("h_dimuons_M_y2_trigger");
	
	TH1F *h_dimuons_M_y3_trigger = new TH1F("h_dimuons_M_y3_trigger","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y3_trigger->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y3_trigger->SetName("h_dimuons_M_y3_trigger");

	TH1F *h_dimuons_M_y4_trigger = new TH1F("h_dimuons_M_y4_trigger","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y4_trigger->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y4_trigger->SetName("h_dimuons_M_y4_trigger");

	TH1F *h_dimuons_M_y5_trigger = new TH1F("h_dimuons_M_y5_trigger","h_dimuons_phi",100,2.8,3.4);
	h_dimuons_M_y5_trigger->SetTitle("; Massa Invariante J/#psi [GeV] ; Eventos ");
	h_dimuons_M_y5_trigger->SetName("h_dimuons_M_y5_trigger");

//-----------------------------------------------------------------------------
	//Behavior of the J/psi candidates
	TH1F *h_Jpsi_Pt = new TH1F("h_Jpsi_Pt","h_Jpsi_Pt",100,10,40);
	h_Jpsi_Pt->SetTitle("Distribuicao p_{T} dos Candidatos a J/#psi; p_{T} [GeV] ; Eventos ");
	h_Jpsi_Pt->SetName("h_Jpsi_Pt");

	TH1F *h_Jpsi_Eta = new TH1F("h_Jpsi_Eta","h_Jpsi_Eta",100,-3,3);
	h_Jpsi_Eta->SetTitle("Distribuicao #eta dos Candidatos a J/#psi; #eta ; Eventos ");
	h_Jpsi_Eta->SetName("h_Jpsi_Eta");
	
	TH1F *h_Jpsi_Rapidity = new TH1F("h_Jpsi_Rapidity","h_Jpsi_Rapidity",100,-3,3);
	h_Jpsi_Rapidity->SetTitle("Distribuicao y dos Candidatos a J/#psi; y ; Eventos ");
	h_Jpsi_Rapidity->SetName("h_Jpsi_Rapidity");

	//Dimuon Invariant Mass
	TH1F *h_dimuons_Mass = new TH1F("h_dimuons_Mass","h_dimuons_Mass",100,2.5,5);
	h_dimuons_Mass->SetTitle("Invariant Mass; Mass [GeV] ; Events ");
	h_dimuons_Mass->SetName("h_dimuons_Mass");

	TH1F *h_Jpsi_Mass = new TH1F("h_Jpsi_Mass","h_Jpsi_Mass",100,3,3.2);
	h_Jpsi_Mass->SetTitle("Invariant Mass; Mass [GeV] ; Events ");
	h_Jpsi_Mass->SetName("h_Jpsi_Mass");

	TH1F *h_Jpsi_peak = new TH1F("h_Jpsi_peak","h_Jpsi_peak",100,3,3.2);
	h_Jpsi_peak->SetTitle("Invariant Mass Jpsi; Mass [GeV] ; Events ");
	h_Jpsi_peak->SetName("h_Jpsi_peak");

	TH1F *h_Tagdimuons_Pt = new TH1F("h_Tagdimuons_Pt","h_Tagdimuons_Pt",100,0,60);
	h_Tagdimuons_Pt->SetTitle("Transverse Momentum ; pT [GeV] ; Events ");
	h_Tagdimuons_Pt->SetName("h_Tagdimuons_Pt");

	TH1F *h_Tagdimuons_Eta = new TH1F("h_Tagdimuons_Eta","h_Tagdimuons_Eta",100,-4,4);
	h_Tagdimuons_Eta->SetTitle("Pseudo rapidity ; Eta; Events ");
	h_Tagdimuons_Eta->SetName("h_Tagdimuons_Eta");

	TH1F *h_leadingmuon_pt = new TH1F("h_leadingmuon_pt","h_leadingmuon_pt",100,0,40);
	h_leadingmuon_pt->SetTitle("Tranverse momentum ; pT; Events ");
	h_leadingmuon_pt->SetName("h_leadingmuon_pt");

	TH1F *h_trailingMuon_Pt = new TH1F("h_trailingMuon_Pt","h_trailingMuon_Pt",100,0,40);
	h_trailingMuon_Pt->SetTitle("Tranverse momentum ; pT; Events ");
	h_trailingMuon_Pt->SetName("h_trailingMuon_Pt");

//end loop jpsi

	

	//3096.900±0.006
	

//End Histograms
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------	
	//Reading the root file and the tree
	//TFile *f1 = new TFile("histo_pythia_trigger.root");
	//TTree *t1 = (TTree*)f1->Get("demo/AnalysisTree");

	//for file2
	TFile *f2 = new TFile("data_histo_all.root");
	TTree *t2 = (TTree*)f2->Get("demo/AnalysisTree");

//---------------------------------------------------------------------------------
	

	//=======================================================================================
	//addressing the memory to vector and variables for file f2
	//For Variables
	TBranch *b_Total_Events2 = t2->GetBranch("Total_Events");
	b_Total_Events2->SetAddress(&Total_Events2);
	TBranch *b_Muons2 = t2->GetBranch("Muons");
	b_Muons2->SetAddress(&Muons2);
	TBranch *b_TMOneStationTight2 = t2->GetBranch("TMOneStationTight");
	b_TMOneStationTight2->SetAddress(&TMOneStationTight2);
	TBranch *b_NumberOfValidMuonHits2 = t2->GetBranch("NumberOfValidMuonHits");
	b_NumberOfValidMuonHits2->SetAddress(&NumberOfValidMuonHits2);
	TBranch *b_pixelLayersWithMeasurement2 = t2->GetBranch("pixelLayersWithMeasurement");
	b_pixelLayersWithMeasurement2->SetAddress(&pixelLayersWithMeasurement2);
	TBranch *b_normalizedChi22 = t2->GetBranch("normalizedChi2");
	b_normalizedChi22->SetAddress(&normalizedChi22);
	TBranch *b_db_dz2 = t2->GetBranch("db_dz");
	b_db_dz2->SetAddress(&db_dz2);
	TBranch *b_PFMuon2 = t2->GetBranch("PFMuon");
	b_PFMuon2->SetAddress(&PFMuon2);
	TBranch *b_TrackerGlobalMuon2 = t2->GetBranch("TrackerGlobalMuon");
	b_TrackerGlobalMuon2->SetAddress(&TrackerGlobalMuon2);
	TBranch *b_nDimuon2 = t2->GetBranch("nDimuon");
	b_nDimuon2->SetAddress(&nDimuon2);

	//For Vectors


	TBranch *b_VectorMuon_Pt = t2->GetBranch("VectorMuon_Pt");
	b_VectorMuon_Pt->SetAddress(&VectorMuon_Pt);
	
	TBranch *b_VectorMuonTight_Pt = t2->GetBranch("VectorMuonTight_Pt");
	b_VectorMuonTight_Pt->SetAddress(&VectorMuonTight_Pt);

	TBranch *b_VectorMuonTightValidHits_Pt = t2->GetBranch("VectorMuonTightValidHits_Pt");
	b_VectorMuonTightValidHits_Pt->SetAddress(&VectorMuonTightValidHits_Pt);

	TBranch *b_VectorMuonTightValidHitsPixelLayer_Pt = t2->GetBranch("VectorMuonTightValidHitsPixelLayer_Pt");
	b_VectorMuonTightValidHitsPixelLayer_Pt->SetAddress(&VectorMuonTightValidHitsPixelLayer_Pt);

	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2_Pt = t2->GetBranch("VectorMuonTightValidHitsPixelLayerChi2_Pt");
	b_VectorMuonTightValidHitsPixelLayerChi2_Pt->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2_Pt);

	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt = t2->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt);

	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt = t2->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt);

	TBranch *b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt = t2->GetBranch("VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt");
	b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->SetAddress(&VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt);

	TBranch *b_leadingMuon_Pt2 = t2->GetBranch("leadingMuon_Pt");
	b_leadingMuon_Pt2->SetAddress(&leadingMuon_Pt2);
	TBranch *b_leadingMuon_Eta2 = t2->GetBranch("leadingMuon_Eta");
	b_leadingMuon_Eta2->SetAddress(&leadingMuon_Eta2);
	TBranch *b_leadingMuon_Phi2 = t2->GetBranch("leadingMuon_Phi");
	b_leadingMuon_Phi2->SetAddress(&leadingMuon_Phi2);
	TBranch *b_leadingMuon_Charge2 = t2->GetBranch("leadingMuon_Charge");
	b_leadingMuon_Charge2->SetAddress(&leadingMuon_Charge2);
	TBranch *b_leadingMuon_Mass2 = t2->GetBranch("leadingMuon_Mass");
	b_leadingMuon_Mass2->SetAddress(&leadingMuon_Mass2);

	TBranch *b_trailingMuon_Pt2 = t2->GetBranch("trailingMuon_Pt");
	b_trailingMuon_Pt2->SetAddress(&trailingMuon_Pt2);
	TBranch *b_trailingMuon_Eta2 = t2->GetBranch("trailingMuon_Eta");
	b_trailingMuon_Eta2->SetAddress(&trailingMuon_Eta2);
	TBranch *b_trailingMuon_Phi2 = t2->GetBranch("trailingMuon_Phi");
	b_trailingMuon_Phi2->SetAddress(&trailingMuon_Phi2);
	TBranch *b_trailingMuon_Charge2 = t2->GetBranch("trailingMuon_Charge");
	b_trailingMuon_Charge2->SetAddress(&trailingMuon_Charge2);
	TBranch *b_trailingMuon_Mass2 = t2->GetBranch("trailingMuon_Mass");
	b_trailingMuon_Mass2->SetAddress(&trailingMuon_Mass2);
	
	//**********************************************************		
	
//**********************************************************		
	//Reading Number of tree entries for file f2
	Long64_t nentries2 = t2->GetEntries();
	cout<< "Numero de Entradas: "<< nentries2 <<std::endl;

	Long64_t nbytes = 0, nb = 0, i=0;
	for (Long64_t kentry=0; kentry<nentries2;kentry++) // loop tree entries for file f2
	{
		Long64_t ientry = t2->LoadTree(kentry);
		//std::cout << "nentries " << nentries2 << " ientry " << kentry << " kentry " <<kentry <<std::endl;
   
		if (ientry < 0) break;

		//For counters (Intergers)
		b_Total_Events2->GetEntry(ientry);
		b_Muons2->GetEntry(ientry);
		b_TMOneStationTight2->GetEntry(ientry);
		b_NumberOfValidMuonHits2->GetEntry(ientry);
		b_pixelLayersWithMeasurement2->GetEntry(ientry);
		b_normalizedChi22->GetEntry(ientry);
		b_db_dz2->GetEntry(ientry);
		b_PFMuon2->GetEntry(ientry);
		b_TrackerGlobalMuon2->GetEntry(ientry);
		b_nDimuon2->GetEntry(ientry);

		count_Total_Events_data += Total_Events2;
		counter_Muon_data += Muons2;
		counter_TMOneStationTight_data += TMOneStationTight2;
		counter_NumberOfValidMuonHits_data += NumberOfValidMuonHits2;
		counter_pixelLayersWithMeasurement_data += pixelLayersWithMeasurement2;
		counter_normalizedChi2_data += normalizedChi22;
		counter_db_dz_data += db_dz2;
		counter_PFMuon_data += PFMuon2;
		counter_TrackerGlobalMuon_data += TrackerGlobalMuon2;
		counter_nDimuon_data += nDimuon2;

		//For vectors
		b_VectorMuon_Pt->GetEntry(ientry);
		
		b_VectorMuonTight_Pt->GetEntry(ientry);

		b_VectorMuonTightValidHits_Pt->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayer_Pt->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayerChi2_Pt->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt->GetEntry(ientry);

		b_VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->GetEntry(ientry);


		b_leadingMuon_Pt2->GetEntry(ientry);
		b_leadingMuon_Eta2->GetEntry(ientry);
		b_leadingMuon_Phi2->GetEntry(ientry);
		b_leadingMuon_Charge2->GetEntry(ientry);
		b_leadingMuon_Mass2->GetEntry(ientry);

		b_trailingMuon_Pt2->GetEntry(ientry);
		b_trailingMuon_Eta2->GetEntry(ientry);
		b_trailingMuon_Phi2->GetEntry(ientry);
		b_trailingMuon_Charge2->GetEntry(ientry);
		b_trailingMuon_Mass2->GetEntry(ientry);

		for(Long64_t i=0; i<VectorMuon_Pt->size();i++)
		{   
			h_Muon_Pt->Fill(VectorMuon_Pt->at(i));			
  		} 

		for(Long64_t i=0; i < VectorMuonTight_Pt->size(); i++)
		{  
			h_MuonTight_Pt->Fill(VectorMuonTight_Pt->at(i));
  		} 

		for(Long64_t i=0; i<VectorMuonTightValidHits_Pt->size();i++)
		{   		
			h_MuonTightValidHits_Pt->Fill(VectorMuonTightValidHits_Pt->at(i));
  		} 

		for(Long64_t i=0; i<VectorMuonTightValidHitsPixelLayer_Pt->size();i++) //loop nDimuon
		{
			h_MuonTightValidHitsPixelLayer_Pt->Fill(VectorMuonTightValidHitsPixelLayer_Pt->at(i));
		}

		for(Long64_t i=0; i<VectorMuonTightValidHitsPixelLayerChi2_Pt->size();i++) //loop nDimuon
		{
			h_MuonTightValidHitsPixelLayerChi2_Pt->Fill(VectorMuonTightValidHitsPixelLayerChi2_Pt->at(i));
		}

		for(Long64_t i=0; i<VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt->size();i++) //loop nDimuon
		{
			h_MuonTightValidHitsPixelLayerChi2DbDz_Pt->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDz_Pt->at(i));
		}

		for(Long64_t i=0; i<VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt->size();i++) //loop nDimuon
		{
			h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDzPf_Pt->at(i));
		}

		for(Long64_t i=0; i<VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->size();i++) //loop nDimuon
		{
			h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->Fill(VectorMuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->at(i));
		}

		
		for(Long64_t i=0; i<trailingMuon_Pt2->size();i++) //loop nDimuon
		{  
			//T		
			mu_1.SetPtEtaPhiM(leadingMuon_Pt2->at(i), leadingMuon_Eta2->at(i), leadingMuon_Phi2->at(i), leadingMuon_Mass2->at(i));
			mu_2.SetPtEtaPhiM(trailingMuon_Pt2->at(i), trailingMuon_Eta2->at(i), trailingMuon_Phi2->at(i), trailingMuon_Mass2->at(i));
			//calculo das grandezas cinematicas
			M = (mu_1+mu_2).Mag();		//Massa Invariante #mu#mu of two Particles
			Pt = (mu_1+mu_2).Pt();      //transverse momentum muon pair
			Eta = (mu_1+mu_2).Eta();      //Pseudo-Rapidity muon pair
			Rapidity = (mu_1+mu_2).Rapidity(); //Rapidity muon pair

			h_dimuons_Mass->Fill(M);
			h_Jpsi_Mass->Fill(M);
			
			//std::cout<< "Invariant Mass  " << i << ": "<< M << std::endl;

			//3.096900-3*0.000006
			
			
			if ( leadingMuon_Charge2->at(i) != trailingMuon_Charge2->at(i) ) // loop charge
			{
				count_OppositeCharge_data++;
				h_leadingmuon_pt->Fill(leadingMuon_Pt2->at(i));
				h_trailingMuon_Pt->Fill(trailingMuon_Pt2->at(i));
											
					if ( fabs(leadingMuon_Eta2->at(i)) < 2.4 && fabs(trailingMuon_Eta2->at(i)) < 2.4 ) // loop eta
					{
						count_Eta_data++;

						if ( (M > 3.08) && (M < 3.11)) // loop Jpsi
						{							
							count_Jpsi_data++;
							h_Jpsi_peak->Fill(M);
							h_Tagdimuons_Pt->Fill(Pt);
							h_Tagdimuons_Eta->Fill(Eta);
						
						}//end loop jpsi
						
					}//end loop eta
			}//end loop charge
		} //end loop nDimuon

}//End loop tree entries for file f2

	cout << "        " << endl;	
	cout << "======================================================== " << endl;
	cout << "count_Total_Events_data: "<< count_Total_Events_data << endl;	
	cout << "counter_Muon_data: "<< counter_Muon_data << endl;	
	cout << "counter_TMOneStationTight_data: "<< counter_TMOneStationTight_data << endl;	
	cout << "counter_NumberOfValidMuonHits_data: "<< counter_NumberOfValidMuonHits_data << endl;	
	cout << "counter_pixelLayersWithMeasurement_data: "<< counter_pixelLayersWithMeasurement_data << endl;	
	cout << "counter_normalizedChi2_data: "<< counter_normalizedChi2_data << endl;	
	cout << "counter_db_dz_data: "<< counter_db_dz_data << endl;	
	cout << "counter_PFMuon_data: "<< counter_PFMuon_data << endl;	
	cout << "counter_TrackerGlobalMuon_data: "<< counter_TrackerGlobalMuon_data << endl;	
	cout << "counter_nDimuon_data: "<< counter_nDimuon_data << endl;
	cout << "count_OppositeCharge_data: "<< count_OppositeCharge_data << endl;
	cout << "count_Eta_data: "<< count_Eta_data << endl;
	cout << "count_Jpsi_data: "<< count_Jpsi_data << endl;
	cout << "count_region1_data: "<< count_region1_data << endl;
	cout << "count_region2_data: "<< count_region2_data << endl;
	cout << "count_region3_data: "<< count_region3_data << endl;
	cout << "count_region4_data: "<< count_region4_data << endl;
	cout << "count_region5_data: "<< count_region5_data << endl;
	cout << "======================================================== " << endl;
	cout << "        " << endl; 

   //=========================================================================		
	//Creating Canvas
	TCanvas* c0 = new TCanvas("c0","Canvas 0 - behavior of the dimuons after quality selection",1200,600);
	//c0->Divide(2,2);
	//c0->cd(1);

	h_Muon_Pt->SetLineColor(kRed);
	h_MuonTight_Pt->SetLineColor(kBlue);
	h_MuonTightValidHits_Pt->SetLineColor(kGreen);
	h_MuonTightValidHitsPixelLayer_Pt->SetLineColor(kBlack);
	h_MuonTightValidHitsPixelLayerChi2_Pt->SetLineColor(kYellow);
	h_MuonTightValidHitsPixelLayerChi2DbDz_Pt->SetLineColor(+4);
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt->SetLineColor(kViolet);
	//h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->SetLineColor(kBlue);

	TLegend* leg_dimuons_Pt = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_Pt->SetFillColor(kWhite);
	leg_dimuons_Pt->SetFillStyle(1001);
	leg_dimuons_Pt->AddEntry(h_Muon_Pt,"#mu","L");
	leg_dimuons_Pt->AddEntry(h_MuonTight_Pt,"Tight #mu","L");
	leg_dimuons_Pt->AddEntry(h_MuonTightValidHits_Pt,"#mu NumberOfValidTrackerHits > 10","L");
	leg_dimuons_Pt->AddEntry(h_MuonTightValidHitsPixelLayer_Pt,"#mu PixelLayersWithMeasurement > 1","L");
	leg_dimuons_Pt->AddEntry(h_MuonTightValidHitsPixelLayerChi2_Pt,"#mu Chi2","L");
	leg_dimuons_Pt->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDz_Pt,"#mu db < 3cm and dz < 15cm","L");
	leg_dimuons_Pt->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt,"Soft #mu","L");
	//leg_dimuons_Pt->AddEntry(h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt,"#mu TrackerGlobal","L");
		
	h_Muon_Pt->Draw();
	h_MuonTight_Pt->Draw("sames");
	h_MuonTightValidHits_Pt->Draw("sames");
	h_MuonTightValidHitsPixelLayer_Pt->Draw("sames");
	h_MuonTightValidHitsPixelLayerChi2_Pt->Draw("sames");
	h_MuonTightValidHitsPixelLayerChi2DbDz_Pt->Draw("sames");
	h_MuonTightValidHitsPixelLayerChi2DbDzPf_Pt->Draw("sames");
	//h_MuonTightValidHitsPixelLayerChi2DbDzPfTG_Pt->Draw("sames");
	leg_dimuons_Pt->Draw();

	//=========================================================================
	//Creating Canvas
	TCanvas* c1 = new TCanvas("c1","Canvas 2 - behavior of the dimuons after quality selection",1200,600);
	h_dimuons_Mass->SetLineColor(kBlack);
	TLegend* leg_dimuons_Mass = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_Mass->SetFillColor(kWhite);
	leg_dimuons_Mass->SetFillStyle(1001);
	leg_dimuons_Mass->AddEntry(h_dimuons_Mass,"#mu#mu","L");
	h_dimuons_Mass->Draw("ep1");
	leg_dimuons_Mass->Draw();

	//=========================================================================	
	//Creating Canvas
	TCanvas* c2 = new TCanvas("c2","Canvas 2 - behavior of the dimuons after quality selection",1200,600);
	c2->Divide(2,2);
	c2->cd(1);
	h_Jpsi_Mass->SetLineColor(kBlack);
	h_Jpsi_peak->SetLineColor(kRed);
	TLegend* leg_dimuons_Mass = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_Mass->SetFillColor(kWhite);
	leg_dimuons_Mass->SetFillStyle(1001);
	leg_dimuons_Mass->AddEntry(h_Jpsi_Mass,"#mu#mu","L");
	leg_dimuons_Mass->AddEntry(h_Jpsi_peak,"Around Jpsi peak","L");
	h_Jpsi_Mass->Draw();
	h_Jpsi_peak->Draw("same");
	leg_dimuons_Mass->Draw();
	
	c2->cd(2);
	h_Tagdimuons_Pt->SetLineColor(kRed);	
	TLegend* leg_Tagdimuons_Pt = new TLegend(0.75,0.81,0.97,0.97);
	leg_Tagdimuons_Pt->SetFillColor(kWhite);
	leg_Tagdimuons_Pt->SetFillStyle(1001);
	leg_Tagdimuons_Pt->AddEntry(h_Tagdimuons_Pt,"Around Jpsi peak","L");	
	h_Tagdimuons_Pt->Draw();
	leg_Tagdimuons_Pt->Draw();
	
	c2->cd(3);
	h_Tagdimuons_Eta->SetLineColor(kRed);
	TLegend* leg_Tagdimuons_Eta = new TLegend(0.75,0.81,0.97,0.97);
	leg_Tagdimuons_Eta->SetFillColor(kWhite);
	leg_Tagdimuons_Eta->SetFillStyle(1001);
	leg_Tagdimuons_Eta->AddEntry(h_Tagdimuons_Eta,"Around Jpsi peak","L");	
	h_Tagdimuons_Eta->Draw();
	leg_Tagdimuons_Eta->Draw();

	//=========================================================================	
	//Creating Canvas
	TCanvas* c3 = new TCanvas("c3","Canvas 3 - behavior of the leading muons and trailing muons",1200,600);
	c3->Divide(2,2);
	c3->cd(1);
	h_leadingmuon_pt->SetLineColor(kBlack);
	h_trailingMuon_Pt->SetLineColor(kRed);
	TLegend* leg_leadingmuon_pt = new TLegend(0.75,0.81,0.97,0.97);
	leg_leadingmuon_pt->SetFillColor(kWhite);
	leg_leadingmuon_pt->SetFillStyle(1001);
	leg_leadingmuon_pt->AddEntry(h_leadingmuon_pt,"leading #mu_{1}","L");
	leg_leadingmuon_pt->AddEntry(h_trailingMuon_Pt,"trailing #mu_{2}","L");
	h_trailingMuon_Pt->Draw();
	h_leadingmuon_pt->Draw("same");
	leg_leadingmuon_pt->Draw();
	
	c3->cd(2);
		
	c3->cd(3);
	

	/*	
	//Creating Canvas
	TCanvas* c2 = new TCanvas("c2","Canvas 2 - behavior of the dimuons after quality selection",1200,600);
	c2->Divide(2,2);
	c2->cd(1);

	h_Dimuons_M->SetLineColor(kRed);
	h_DimuonsOppositeCharge_M->SetLineColor(kBlue);
	h_DimuonsOppositeChargeEta_M->SetLineColor(kBlack);
	h_DimuonsOppositeChargeEtaJpsi_M->SetLineColor(kViolet);
	
	TLegend* leg_dimuons_M = new TLegend(0.75,0.81,0.97,0.97);
	leg_dimuons_M->SetFillColor(kWhite);
	leg_dimuons_M->SetFillStyle(1001);
	leg_dimuons_M->AddEntry(h_Dimuons_M,"#mu#mu","L");
	leg_dimuons_M->AddEntry(h_DimuonsOppositeCharge_M,"#mu^{+}#mu^{-}","L");
	leg_dimuons_M->AddEntry(h_DimuonsOppositeChargeEta_M,"#mu^{+}#mu^{-} Eta<2.4","L");
	leg_dimuons_M->AddEntry(h_DimuonsOppositeChargeEtaJpsi_M,"Candidatos a J/#psi","L");
	
	h_Dimuons_M->Draw();
	h_DimuonsOppositeCharge_M->Draw("sames");
	h_DimuonsOppositeChargeEta_M->Draw("sames");
	h_DimuonsOppositeChargeEtaJpsi_M->Draw("sames");
	leg_dimuons_M->Draw();*/
	
	
	
	//-------------------------------------------------------------------------

}//end program

