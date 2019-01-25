// -*- C++ -*-
//
// Package:    ZDC/newZDCAnalyzer
// Class:      newZDCAnalyzer
//
/**\class newZDCAnalyzer newZDCAnalyzer.cc ZDC/newZDCAnalyzer/plugins/newZDCAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Oliver Suranyi
//         Created:  Sat, 03 Nov 2018 14:24:57 GMT
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
//#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "TTree.h"

#include "nominal_fC.h"

//
// class declaration
//

const int MAX = 50000;
const int MAX_TS = 10;

class newZDCAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit newZDCAnalyzer(const edm::ParameterSet&);
    ~newZDCAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    //HLTPrescaleProvider hltPrescaleProvider_;
    string trgResultsProcess_;
    
    edm::Service<TFileService> fs;

    edm::EDGetTokenT<QIE10DigiCollection> zdcToken;
    edm::EDGetTokenT<CaloTowerCollection> caloTowerToken;
    edm::EDGetTokenT<edm::TriggerResults> hltToken;
    edm::EDGetTokenT<reco::TrackCollection> trackToken;
    //edm::EDGetTokenT<SiPixelRecHitCollection> pixelToken;

    TTree* zdcDigiTree;
    int run, lumi, event, bxid;
   
    int n;

    int zside[MAX];
    int section[MAX];
    int channel[MAX];
    int adc[MAX_TS][MAX];
    float nfC[MAX_TS][MAX];
    float rfC[MAX_TS][MAX];

    int rejected;

    int nPixel;
    int nTrack;
    int nHF_pos, nHF_neg;
    float eHF_pos, eHF_neg;

    bool firstEvent;
    int* trigflag;
    int* prscflag;
    //HLTPrescaleProvider hltPrescaleProvider;
};

newZDCAnalyzer::newZDCAnalyzer(const edm::ParameterSet& iConfig) /*:
  hltPrescaleProvider(iConfig, consumesCollector(), *this)*/{
  
  usesResource("TFileService");
  zdcToken = consumes<QIE10DigiCollection>(iConfig.getParameter<edm::InputTag>("zdc"));
  caloTowerToken = consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("tower"));
  hltToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltresults"));
  trgResultsProcess_ = iConfig.getParameter<edm::InputTag>("hltresults").process();
  //hltCfg = hltPrescaleProvider_.hltConfigProvider(iConfig, consumesCollector(), *this);

  trackToken = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("track"));
  //pixelToken = consumes<SiPixelRecHitCollection>(iConfig.getParameter<edm::InputTag>("pixel"));
}


newZDCAnalyzer::~newZDCAnalyzer(){
}


//
// member functions
//

// ------------ method called for each event  ------------
void newZDCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  event = iEvent.id().event();
  run = iEvent.id().run();
  lumi = iEvent.luminosityBlock();
  bxid = iEvent.bunchCrossing();

  // Processing ZDC QIE10 digis
  edm::Handle<QIE10DigiCollection> zdcDigiCollection;
  iEvent.getByToken(zdcToken, zdcDigiCollection);

  //edm::ESHandle<HcalDbService> conditions;
  //iSetup.get<HcalDbRecord>().get(conditions);

  // Add ZDC digi branches to event tree in the first cycle
  if(firstEvent){
    const QIE10DataFrame& frame = (*zdcDigiCollection)[0];
    std::cout << "Number of timeslices: " << frame.samples() << std::endl;

    for(int i=0; i < frame.samples();i++){
      TString adc_str("adc"), nfC_str("nfC"), rfC_str("rfC");
      adc_str+=i; nfC_str+=i, rfC_str+=i;

      zdcDigiTree->Branch(adc_str,adc[i],adc_str+"[n]/I");
      zdcDigiTree->Branch(nfC_str,nfC[i],nfC_str+"[n]/F");
      //zdcDigiTree->Branch(rfC_str,rfC[i],rfC_str+"[n]/F");
    }
  }

  n = 0;
  rejected = 0;

  bool valid = true;

  // Read ZDC QIE10 digis  
  for(QIE10DigiCollection::const_iterator it = zdcDigiCollection->begin(); it != zdcDigiCollection->end(); it++) {
    const QIE10DataFrame& frame(*it);
    CaloSamples caldigi;

    HcalZDCDetId cell = frame.id();

    /*const HcalQIECoder* qiecoder = conditions->getHcalCoder(cell);
    const HcalQIEShape* qieshape = conditions->getHcalShape(qiecoder);
   
    HcalCoderDb coder(*qiecoder,*qieshape);
    coder.adc2fC(frame,caldigi);*/

    zside[n] = cell.zside();
    section[n] = cell.section();
    channel[n] = cell.channel();

    int prevCapID;

    for(int isample = 0; isample < 10; isample++){
      if(isample != 0 && frame[isample].capid() != ((prevCapID+1)%4)){
        valid = false;
        rejected++;
        for(int jsample = 0; jsample < 10; jsample++){
          std::cout << frame[jsample].capid() << std::endl;
        }
        while(!getchar());
        break;
      }
      prevCapID = frame[isample].capid();

      adc[isample][n] = frame[isample].adc();
      nfC[isample][n] = QIE10_nominal_fC[frame[isample].adc()];
      //rfC[isample][n] = caldigi[isample];
    }

    n++;
  }

  // Processing HF towers
  edm::Handle<CaloTowerCollection> towerCollection;
  iEvent.getByToken(caloTowerToken,towerCollection);

  nHF_pos = 0;
  nHF_neg = 0;
  
  eHF_pos = 0.0;
  eHF_neg = 0.0;

  for(auto it = towerCollection->begin(); it != towerCollection->end(); it++){
    bool isHF = false;
    for(unsigned int i = 0; i < it->constituentsSize(); i++){
      DetId cell = it->constituent(i);
      if(cell.det() == DetId::Hcal && (HcalSubdetector)cell.subdetId() == HcalForward)
        isHF = true;
    }

    if(isHF && it->energy() > 4.0){
      if(it->eta() > 0){
        nHF_pos++;
        eHF_pos += it->energy();
      }
      else{
        nHF_neg++;
        eHF_neg += it->energy();    
      }
    }
  }

  // Processing tracks
  edm::Handle<reco::TrackCollection> trackCollection;
  iEvent.getByToken(trackToken, trackCollection);

  nTrack = trackCollection->size();

  // Processing pixels
  /*edm::Handle<SiPixelRecHitCollection> pixelCollection;
  iEvent.getByToken(pixelToken, pixelCollection);    

  nPixel = pixelCollection->size();
*/
  // Processing HLT results
  
  //HLTConfigProvider const&  hltCfg = hltPrescaleProvider.hltConfigProvider();

  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByToken(hltToken,hltresults);
  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*hltresults);

  int ntrigs = hltresults->size();

  if(firstEvent){
    for(int itrig = 0; itrig != ntrigs; ++itrig) {
      TString trigName = triggerNames.triggerName(itrig);
      zdcDigiTree->Branch(trigName+"_acc", trigflag+itrig, trigName+"_acc/I");
      zdcDigiTree->Branch(trigName+"_psc", prscflag+itrig, trigName+"_psc/I");
    }
    firstEvent = false;
  }

  for(int itrig = 0; itrig != ntrigs; ++itrig){
    trigflag[itrig] = (hltresults->accept(itrig)) ? 1 : 0;
  //  prscflag[itrig] = (hltCfg.moduleType(hltCfg.moduleLabel(itrig,hltresults->index(itrig)))=="HLTPrescaler") ? 1 : 0;
      //if(!prscflag[n])
      //  save = true;+;
  }

  // Fill event tree
  if(valid)
    zdcDigiTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void newZDCAnalyzer::beginJob(){
  zdcDigiTree = fs->make<TTree>("zdcdigi","v1");
  
  zdcDigiTree->Branch("run",&run,"run/I");
  zdcDigiTree->Branch("lumi",&lumi,"lumi/I");
  zdcDigiTree->Branch("event",&event,"event/I");
  zdcDigiTree->Branch("bxid",&bxid,"bxid/I");

  zdcDigiTree->Branch("n",&n,"n/I");
  zdcDigiTree->Branch("zside",zside,"zside[n]/I");
  zdcDigiTree->Branch("section",section,"section[n]/I");
  zdcDigiTree->Branch("channel",channel,"channel[n]/I");

  zdcDigiTree->Branch("rejected",rejected,"rejected[n]/I"); 

  zdcDigiTree->Branch("nHF_pos",&nHF_pos,"nHF_pos/I");
  zdcDigiTree->Branch("nHF_neg",&nHF_neg,"nHF_neg/I");

  zdcDigiTree->Branch("eHF_pos",&eHF_pos,"eHF_pos/F");
  zdcDigiTree->Branch("eHF_neg",&eHF_neg,"eHF_neg/F");

  zdcDigiTree->Branch("nTrack",&nTrack,"nTrack/I");
  zdcDigiTree->Branch("nPixel",&nPixel,"nPixel/I");

  firstEvent = true;

  const int Max = 500;
  trigflag = new int[Max];
  prscflag = new int[Max];
}

// ------------ method called once each job just after ending the event loop  ------------
void newZDCAnalyzer::endJob(){
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
newZDCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(newZDCAnalyzer);
