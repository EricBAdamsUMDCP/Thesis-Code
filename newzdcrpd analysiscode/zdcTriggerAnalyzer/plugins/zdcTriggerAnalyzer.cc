// -*- C++ -*-
//
// Package:    ZDC/zdcTriggerAnalyzer
// Class:      zdcTriggerAnalyzer
//
/**\class zdcTriggerAnalyzer zdcTriggerAnalyzer.cc ZDC/zdcTriggerAnalyzer/plugins/zdcTriggerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Oliver Suranyi
//         Created:  Wed, 21 Nov 2018 14:36:08 GMT
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HcalDigi/interface/HFDataFrame.h"

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HFDataFrame.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"


#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class zdcTriggerAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit zdcTriggerAnalyzer(const edm::ParameterSet&);
    ~zdcTriggerAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
   
    edm::Service<TFileService> fs;
    //edm::EDGetTokenT<reco::Centrality> CentralityTag_;
    edm::EDGetTokenT<int> CentralityBinTag_;
    //edm::EDGetTokenT<QIE10DigiCollection> hfToken;
    edm::EDGetTokenT<QIE10DigiCollection> zdcToken;
    edm::EDGetTokenT<edm::TriggerResults> hltToken;

    TTree* zdcTriggerTree;
    int run, lumi, event, bxid;

    int hiBin;
    /*int hiNpix, hiNpixelTracks, hiNtracks, hiNtracksPtCut, hiNtracksEtaCut, hiNtracksEtaPtCut;
    int hiNpixPlus, hiNpixMinus, hiNpixelTracksPlus, hiNpixelTracksMinus;
    float hiHF, hiHFplus, hiHFminus, hiHFplusEta4, hiHFminusEta4, hiHFhit, hiHFhitPlus, hiHFhitMinus;
    float hiHFECut, hiHFECutPlus, hiHFECutMinus;
    float hiEB, hiET, hiEE, hiEEplus, hiEEminus;
    float hiZDC, hiZDCplus, hiZDCminus;*/

    int HF1pos, HF1neg, HF2pos, HF2neg;
    int ZDCpos, ZDCneg;
   
    bool firstEvent;
    int* trigflag;
};


//
// constructors and destructor
//

zdcTriggerAnalyzer::zdcTriggerAnalyzer(const edm::ParameterSet& iConfig) 
  : //CentralityTag_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("CentralitySrc"))),
    CentralityBinTag_(consumes<int>(iConfig.getParameter<edm::InputTag>("CentralityBinSrc"))),
    //hfToken(consumes<QIE10DigiCollection>(iConfig.getParameter<edm::InputTag>("hf"))),
    zdcToken(consumes<QIE10DigiCollection>(iConfig.getParameter<edm::InputTag>("zdc"))),
    hltToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltresults")))
    /*trgResultsProcess_(iConfig.getParameter<edm::InputTag>("hltresults").process())*/{}


zdcTriggerAnalyzer::~zdcTriggerAnalyzer(){}


//
// member functions
//

// ------------ method called for each event  ------------
void
zdcTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  edm::Handle<int> cbin_;
  iEvent.getByToken(CentralityBinTag_,cbin_);
  hiBin = *cbin_;

  /*edm::Handle<reco::Centrality> centrality;
  iEvent.getByToken(CentralityTag_, centrality);

  hiNpix = centrality->multiplicityPixel();
  hiNpixPlus = centrality->multiplicityPixelPlus();
  hiNpixMinus = centrality->multiplicityPixelMinus();
  hiNpixelTracks = centrality->NpixelTracks();
  hiNpixelTracksPlus = centrality->NpixelTracksPlus();
  hiNpixelTracksMinus = centrality->NpixelTracksMinus();
  hiNtracks = centrality->Ntracks();
  hiNtracksPtCut = centrality->NtracksPtCut();
  hiNtracksEtaCut = centrality->NtracksEtaCut();
  hiNtracksEtaPtCut = centrality->NtracksEtaPtCut();

  hiHF = centrality->EtHFtowerSum();
  hiHFplus = centrality->EtHFtowerSumPlus();
  hiHFminus = centrality->EtHFtowerSumMinus();
  hiHFECut = centrality->EtHFtowerSumECut();
  hiHFECutPlus = centrality->EtHFtowerSumECutPlus();
  hiHFECutMinus = centrality->EtHFtowerSumECutMinus();
  hiHFplusEta4 = centrality->EtHFtruncatedPlus();
  hiHFminusEta4 = centrality->EtHFtruncatedMinus();
  hiHFhit = centrality->EtHFhitSum();
  hiHFhitPlus = centrality->EtHFhitSumPlus();
  hiHFhitMinus = centrality->EtHFhitSumMinus();

  hiZDC = centrality->zdcSum();
  hiZDCplus = centrality->zdcSumPlus();
  hiZDCminus = centrality->zdcSumMinus();

  hiEEplus = centrality->EtEESumPlus();
  hiEEminus = centrality->EtEESumMinus();
  hiEE = centrality->EtEESum();
  hiEB = centrality->EtEBSum();
  hiET = centrality->EtMidRapiditySum();*/

  // Check HF results
  /*HF1pos = 0;
  HF1neg = 0;
  HF2pos = 0;
  HF2neg = 0; 

  edm::Handle<QIE10DigiCollection> hfDigiCollection;
  iEvent.getByToken(hfToken,hfDigiCollection); 

  for (QIE10DigiCollection::const_iterator it = hfDigiCollection->begin(); it != hfDigiCollection->end(); it++) {
    const QIE10DataFrame& frame(*it);
   
    HcalDetId cell = frame.id();

    if(cell.ieta() > 0){
      if(frame[1].adc() > 15)
        HF1pos = 1;
      if(frame[1].adc() > 19)
        HF2pos = 1;
    }
    else{
      if(frame[1].adc() > 15)
        HF1neg = 1;
      if(frame[1].adc() > 19)
        HF2neg = 1;  
    }  
  }*/

  // Processing ZDC QIE10 digis
  edm::Handle<QIE10DigiCollection> zdcDigiCollection;
  iEvent.getByToken(zdcToken, zdcDigiCollection);

  ZDCpos = 0;
  ZDCneg = 0;

  for(QIE10DigiCollection::const_iterator it = zdcDigiCollection->begin(); it != zdcDigiCollection->end(); it++) {
    const QIE10DataFrame& frame(*it);
  
    HcalZDCDetId cell = frame.id();

    int zside = cell.zside();
    int section = cell.section();
    int channel = cell.channel();

    if(zside == 1 && section == 2){
      if((channel == 1 || channel == 2) && frame[4].adc() > 100)
        ZDCpos = 1;
    }
    else if(zside == -1 && section == 2){
      if(channel == 1 && frame[4].adc() > 119)
        ZDCneg = 1;
      if(channel == 2 && frame[4].adc() > 100)
        ZDCneg = 1;
    }
  }

  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByToken(hltToken,hltresults);
  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*hltresults);

  int ntrigs = hltresults->size();

  if(firstEvent){
    for(int itrig = 0; itrig != ntrigs; ++itrig) {
      TString trigName = triggerNames.triggerName(itrig);
      zdcTriggerTree->Branch(trigName,trigflag+itrig, trigName+"/I");
    }
    firstEvent = false;
  }

  for(int itrig = 0; itrig != ntrigs; ++itrig){
    trigflag[itrig] = (hltresults->accept(itrig)) ? 1 : 0;
  //  prscflag[itrig] = (hltCfg.moduleType(hltCfg.moduleLabel(itrig,hltresults->index(itrig)))=="HLTPrescaler") ? 1 : 0;
      //if(!prscflag[n])
      //  save = true;+;
  }


  zdcTriggerTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
zdcTriggerAnalyzer::beginJob(){
  zdcTriggerTree = fs->make<TTree>("zdctrigger","v1");
  
  zdcTriggerTree->Branch("run",&run,"run/I");
  zdcTriggerTree->Branch("lumi",&lumi,"lumi/I");
  zdcTriggerTree->Branch("event",&event,"event/I");
  zdcTriggerTree->Branch("bxid",&bxid,"bxid/I");

  zdcTriggerTree->Branch("hiBin",&hiBin,"hiBin/I");

  /*zdcTriggerTree->Branch("HF1pos",&HF1pos,"HF1pos/I");
  zdcTriggerTree->Branch("HF1neg",&HF1neg,"HF1neg/I");
  zdcTriggerTree->Branch("HF2pos",&HF2pos,"HF2pos/I");
  zdcTriggerTree->Branch("HF2neg",&HF2neg,"HF2neg/I");
*/
  zdcTriggerTree->Branch("ZDCpos",&ZDCpos,"ZDCpos/I");
  zdcTriggerTree->Branch("ZDCneg",&ZDCneg,"ZDCneg/I");

  const int Max = 500;
  trigflag = new int[Max];

}

// ------------ method called once each job just after ending the event loop  ------------
void
zdcTriggerAnalyzer::endJob(){}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
zdcTriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(zdcTriggerAnalyzer);
