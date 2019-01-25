import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo",eras.Run2_2018_pp_on_AA)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.TFileService = cms.Service("TFileService",fileName=cms.string("zdc_digi_tree.root"))

process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring('/store/express/HIRun2018A/HIExpressPhysics/FEVT/Express-v1/000/326/381/00000/042E0A5B-1158-134C-BD04-F8B7D24F226C.root'))
#process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring('file:test_reco_py_RAW2DIGI_RECO.root'))

process.analyzer = cms.EDAnalyzer('newZDCAnalyzer',
                                zdc = cms.InputTag("hcalDigis","ZDC","RECO"),
                                tower = cms.InputTag("towerMaker"),
                                track = cms.InputTag("generalTracks"),
                                pixel = cms.InputTag("siPixelRecHits"),
                                hltresults = cms.InputTag("TriggerResults","","HLT")
                                )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v10'
#process.GlobalTag.globaltag = '103X_dataRun2_Express_v2'
#process.GlobalTag.globaltag = '103X_dataRun2_HLT_v1'
process.GlobalTag.globaltag = '103X_dataRun2_Prompt_v3'


#process.es_ascii = cms.ESSource(
#    'HcalTextCalibrations',
#    input = cms.VPSet(
#            cms.PSet(
#                object = cms.string('ElectronicsMap'),
#                file = cms.FileInPath("QWAna/QWNtrkOfflineProducer/run2018/HcalElectronicsMap_2018_v3.0_data_ext.txt")
#                )
#          )
#)
#process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')


#set digi and analyzer
#process.hcalDigis.InputLabel = "rawDataCollector"

#process.QWInfo = cms.EDProducer('QWEventInfoProducer')

# ZDC info
#process.load('QWZDC2018Producer_cfi')
#process.load('ZDC2018Pedestal_cfg')
#process.zdcdigi.SOI = cms.untracked.int32(4)

#process.digiPath = cms.Path(
#    process.hcalDigis *
#    process.zdcdigi
#)

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')

process.hltHighLevel.HLTPaths = cms.vstring("HLT_HIZeroBias_v1")
process.hltFilter = cms.Path(process.hltHighLevel)

#process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
#process.centrality = cms.Path(process.hiEvtAnalyzer)

process.analyze_step = cms.Path(process.analyzer)

process.schedule = cms.Schedule(
#                    process.digiPath,
                    process.hltFilter,
#                    process.centrality,
                    process.analyze_step)
