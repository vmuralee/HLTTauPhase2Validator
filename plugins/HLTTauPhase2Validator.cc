// -*- C++ -*-
//
// Package:    Phase2Monitoring/HLTTauPhase2Validator
// Class:      HLTTauPhase2Validator
//
/**\class HLTTauPhase2Validator HLTTauPhase2Validator.cc Phase2Monitoring/HLTTauPhase2Validator/plugins/HLTTauPhase2Validator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vinaya Krishnan Nair
//         Created:  Sat, 24 May 2025 13:01:56 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
#include "TMath.h"


//ROOT inclusion

#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TMath.h"
#include "TList.h"
#include "TString.h"


class HLTTauPhase2Validator : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit HLTTauPhase2Validator(const edm::ParameterSet&);
  ~HLTTauPhase2Validator() override;

  

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  double deltaPhi(double phi1, double phi2);
  double deltaR(double eta1, double phi1, double eta2, double phi2);

  
  // ----------member data ---------------------------

  
  TTree* RelValTree;

  TH1F* hTrig1PassPt;
  TH1F* hTrig1PassEta;
  TH1F* hTrig2PassPt;
  TH1F* hTrig2PassEta;

  TH1F* hTrig1TotalPt;
  TH1F* hTrig1TotalEta;
  TH1F* hTrig2TotalPt;
  TH1F* hTrig2TotalEta;
  
  
  edm::Service<TFileService> fs;
  edm::EventNumber_t eventN;

  std::string channel_;
  std::string triggerfilter1_;  // Read from config  {ditau, mutau and etau}
  std::string triggerfilter2_;
  
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticlesToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjectsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;

  const static int nMax = 40000000;
  
  int nTaus;
  int nMuons;
  int nElectrons;


  float TrigObj1Pt;
  float TrigObj1Eta;
  float TrigObj1Phi;

  float TrigObj2Pt;
  float TrigObj2Eta;
  float TrigObj2Phi;

  
  float GenTauPt;
  float GenTauEta;
  float GenTauPhi;
  float GenTauM;

  float GenMuonPt;
  float GenMuonEta;
  float GenMuonPhi;
  float GenMuonM;
  
  float GenElePt;
  float GenEleEta;
  float GenElePhi;
  float GenEleM;
  
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

HLTTauPhase2Validator::HLTTauPhase2Validator(const edm::ParameterSet& iConfig){
  usesResource("TFileService");
  RelValTree = fs->make<TTree>("RelValTree","RelValTree");

  hTrig1PassPt  = fs->make<TH1F>("hTrig1PassPt","",20,0,500);
  hTrig1PassEta = fs->make<TH1F>("hTrig1PassEta","",9,-2.3,2.3);
  hTrig2PassPt  = fs->make<TH1F>("hTrig2PassPt","",20,0,500);
  hTrig2PassEta = fs->make<TH1F>("hTrig2PassEta","",9,-2.3,2.3);
  
  hTrig1TotalPt  = fs->make<TH1F>("hTrig1TotalPt","",20,0,500);
  hTrig1TotalEta = fs->make<TH1F>("hTrig1TotalEta","",9,-2.3,2.3);
  hTrig2TotalPt  = fs->make<TH1F>("hTrig2TotalPt","",20,0,500);
  hTrig2TotalEta = fs->make<TH1F>("hTrig2TotalEta","",9,-2.3,2.3);
  
  RelValTree->Branch("event", &eventN, "event/i");

  RelValTree->Branch("nTaus",&nTaus,"nTaus/i");
  RelValTree->Branch("nMuons",&nMuons,"nMuons/i");
  RelValTree->Branch("nElectrons",&nElectrons,"nElectrons/i");
  
  RelValTree->Branch("TrigObj1Pt", &TrigObj1Pt, "TrigObj1Pt/F");
  RelValTree->Branch("TrigObj1Eta", &TrigObj1Eta, "TrigObj1Eta/F");
  RelValTree->Branch("TrigObj1Phi", &TrigObj1Phi, "TrigObj1Phi/F");

  RelValTree->Branch("TrigObj2Pt", &TrigObj2Pt, "TrigObj2Pt/F");
  RelValTree->Branch("TrigObj2Eta", &TrigObj2Eta, "TrigObj2Eta/F");
  RelValTree->Branch("TrigObj2Phi", &TrigObj2Phi, "TrigObj2Phi/F");
  
  RelValTree->Branch("GenTauPt", &GenTauPt, "GenTauPt/F");
  RelValTree->Branch("GenTauEta", &GenTauEta, "GenTauEta/F");
  RelValTree->Branch("GenTauPhi", &GenTauPhi, "GenTauPhi/F");
  RelValTree->Branch("GenTauM", &GenTauM, "GenTauM/F");

  RelValTree->Branch("GenMuonPt", &GenMuonPt, "GenMuonPt/F");
  RelValTree->Branch("GenMuonEta", &GenMuonEta, "GenMuonEta/F");
  RelValTree->Branch("GenMuonPhi", &GenMuonPhi, "GenMuonPhi/F");
  RelValTree->Branch("GenMuonM", &GenMuonM, "GenMuonM/F");

  RelValTree->Branch("GenElePt", &GenElePt, "GenElePt/F");
  RelValTree->Branch("GenEleEta", &GenEleEta, "GenEleEta/F");
  RelValTree->Branch("GenElePhi", &GenElePhi, "GenElePhi/F");
  RelValTree->Branch("GenEleM", &GenEleM, "GenEleM/F");

  genParticlesToken_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));
  triggerObjectsToken_ = consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"));
  triggerBitsToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"));

  triggerfilter1_ = iConfig.getParameter<std::string>("triggerfilter1");
  triggerfilter2_ = iConfig.getParameter<std::string>("triggerfilter2");


  channel_ = iConfig.getParameter<std::string>("channel");
  
  
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

HLTTauPhase2Validator::~HLTTauPhase2Validator() {
  
}

//
// member functions
//
double HLTTauPhase2Validator::deltaPhi(double phi1, double phi2) {
    double dphi = std::fmod(phi1 - phi2, 2 * M_PI);
    if (dphi > M_PI)
        dphi -= 2 * M_PI;
    if (dphi < -M_PI)
        dphi += 2 * M_PI;
    return dphi;
}

// Compute deltaR between two objects in eta-phi space
double HLTTauPhase2Validator::deltaR(double eta1, double phi1, double eta2, double phi2) {
    double dEta = eta1 - eta2;
    double dPhi = deltaPhi(phi1, phi2);
    return std::sqrt(dEta * dEta + dPhi * dPhi);
}
// ------------ method called for each event  ------------
void HLTTauPhase2Validator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

   // Load trigger objects                                                                                                                                                                                 
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjectsToken_, triggerObjects);

  // Load trigger results (needed for unpacking)                                                                                                                                                            
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsToken_, triggerBits);

  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  std::vector<const reco::GenParticle*> sortedGenParticles;

  // pT ordered Leptons
  for (const auto& p : *genParticles) {
    if(!p.isLastCopy())continue;
    if (std::abs(p.pdgId()) == 11 || std::abs(p.pdgId()) == 13 || std::abs(p.pdgId()) == 15) {
      sortedGenParticles.push_back(&p);
    }
  }

  // Sort by ascending pT
  std::sort(sortedGenParticles.begin(), sortedGenParticles.end(),
            [](const reco::GenParticle* a, const reco::GenParticle* b) {
              return a->pt() < b->pt();
            });

  nTaus = std::count_if(sortedGenParticles.begin(), sortedGenParticles.end(),[](const reco::GenParticle* x){return std::abs(x->pdgId()) == 15;});
  nMuons = std::count_if(sortedGenParticles.begin(), sortedGenParticles.end(),[](const reco::GenParticle* x){return std::abs(x->pdgId()) == 13;});
  nElectrons = std::count_if(sortedGenParticles.begin(), sortedGenParticles.end(),[](const reco::GenParticle* x){return std::abs(x->pdgId()) == 11;});

  std::vector<const reco::GenParticle*> taus;
  std::vector<const reco::GenParticle*> muons;
  std::vector<const reco::GenParticle*> electrons;
  
  std::copy_if(sortedGenParticles.begin(), sortedGenParticles.end(),std::back_inserter(taus),[](const reco::GenParticle* x){ return std::abs(x->pdgId()) == 15 && x->status() == 2;});
  std::copy_if(sortedGenParticles.begin(), sortedGenParticles.end(),std::back_inserter(muons),[](const reco::GenParticle* x){ return std::abs(x->pdgId()) == 13;});
  std::copy_if(sortedGenParticles.begin(), sortedGenParticles.end(),std::back_inserter(electrons),[](const reco::GenParticle* x){ return std::abs(x->pdgId()) == 11;});


  for (auto &tau : taus){
    bool hasLeptonicDecay = false;
    for (size_t i = 0; i < tau->numberOfDaughters(); ++i) {
      const reco::Candidate* dau = tau->daughter(i);
      int dauId = std::abs(dau->pdgId());
      
      // Leptonic decay if e or mu present
      if (dauId == 11 || dauId == 13){
	bool hasLeptonicDecay = true;
        break;}
    }
    if (tau->pt() < 30) continue;
    if (std::fabs(tau->eta()) > 2.1) continue;
    
    if(!hasLeptonicDecay){  // Hadronic decay only
      GenTauPt = tau->pt();
      GenTauEta = tau->eta();
      GenTauPhi = tau->phi();
      GenTauM = tau->energy();
    }
  }

  for (auto &mu : muons){
    if (mu->pt() < 21) continue;
    if (std::fabs(mu->eta()) > 2.1) continue;
    for (auto &tau :taus){
      if(deltaR(mu->eta(),mu->phi(),tau->eta(),tau->phi()) < 0.4)
	continue;
      GenMuonPt = mu->pt();
      GenMuonEta = mu->eta();
      GenMuonPhi = mu->phi();
      GenMuonM = mu->energy();
    }
  }

  for (auto &ele : electrons){
    if (ele->pt() < 31) continue;
    if (std::fabs(ele->eta()) > 2.1) continue;
    for (auto &tau :taus){
      if(deltaR(ele->eta(),ele->phi(),tau->eta(),tau->phi()) < 0.4)
	continue;
      GenElePt = ele->pt();
      GenEleEta = ele->eta();
      GenElePhi = ele->phi();
      GenEleM = ele->energy();
    }
  }


  TLorentzVector genTau; genTau.SetPtEtaPhiE(GenTauPt,GenTauEta,GenTauPhi,GenTauM);
  TLorentzVector genMuon; genMuon.SetPtEtaPhiE(GenMuonPt,GenMuonEta,GenMuonPhi,GenMuonM);
  TLorentzVector genEle; genEle.SetPtEtaPhiE(GenElePt,GenEleEta,GenElePhi,GenEleM);

  
  if(channel_ == "ditau"){
    hTrig1TotalPt->Fill(genTau.Pt());
    hTrig1TotalEta->Fill(genTau.Eta());
    hTrig2TotalPt->Fill(genTau.Pt());
    hTrig2TotalEta->Fill(genTau.Eta());
  }  
  else if(channel_ == "mutau"){
    hTrig1TotalPt->Fill(genMuon.Pt());
    hTrig1TotalEta->Fill(genMuon.Eta());
    hTrig2TotalPt->Fill(genTau.Pt());
    hTrig2TotalEta->Fill(genTau.Eta());
  }
  else{
    hTrig1TotalPt->Fill(genEle.Pt());
    hTrig1TotalEta->Fill(genEle.Eta());
    hTrig2TotalPt->Fill(genTau.Pt());
    hTrig2TotalEta->Fill(genTau.Eta());
  }

  for (const auto& obj : *triggerObjects) {
    pat::TriggerObjectStandAlone unpackedObj = obj;
    unpackedObj.unpackPathNames(iEvent.triggerNames(*triggerBits));
    unpackedObj.unpackFilterLabels(iEvent, *triggerBits);

    
    if (std::find(unpackedObj.filterLabels().begin(), unpackedObj.filterLabels().end(), triggerfilter1_) != unpackedObj.filterLabels().end()) {
      TLorentzVector trigObj1(unpackedObj.px(), unpackedObj.py(), unpackedObj.pz(), unpackedObj.energy());
      TrigObj1Pt  = trigObj1.Pt();
      TrigObj1Eta = trigObj1.Eta();
      TrigObj1Phi = trigObj1.Phi();

      if (channel_ == "ditau"){
	if (genTau.DeltaR(trigObj1) < 0.3){
	  hTrig1PassPt->Fill(genTau.Pt());
	  hTrig1PassEta->Fill(genTau.Eta());
	  break;
	}   
      }
      else if (channel_ == "mutau"){
	if (genMuon.DeltaR(trigObj1) < 0.3){
	  hTrig1PassPt->Fill(genMuon.Pt());
	  hTrig1PassEta->Fill(genMuon.Eta());
	  break;
	}
	
      }
      else {
	if (genEle.DeltaR(trigObj1) < 0.3){
	  hTrig1PassPt->Fill(genEle.Pt());
	  hTrig1PassEta->Fill(genEle.Eta());
	  break;
	}
      }
    }
  }

  for (const auto& obj : *triggerObjects) {
    pat::TriggerObjectStandAlone unpackedObj = obj;
    unpackedObj.unpackPathNames(iEvent.triggerNames(*triggerBits));
    unpackedObj.unpackFilterLabels(iEvent, *triggerBits);
  
    if (std::find(unpackedObj.filterLabels().begin(), unpackedObj.filterLabels().end(), triggerfilter2_) != unpackedObj.filterLabels().end()) {
      TLorentzVector trigObj2(unpackedObj.px(), unpackedObj.py(), unpackedObj.pz(), unpackedObj.energy());
      TrigObj2Pt  = trigObj2.Pt();
      TrigObj2Eta = trigObj2.Eta();
      TrigObj2Phi = trigObj2.Phi();

      if (genTau.DeltaR(trigObj2) < 0.3){
	hTrig2PassPt->Fill(genTau.Pt());
	hTrig2PassEta->Fill(genTau.Eta());
	break;
      }
    }    
  }

  
  RelValTree->Fill();
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}


//define this as a plug-in
DEFINE_FWK_MODULE(HLTTauPhase2Validator);
