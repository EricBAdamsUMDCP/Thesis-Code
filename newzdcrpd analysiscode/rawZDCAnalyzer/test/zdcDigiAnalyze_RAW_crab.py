import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

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
#process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
#process.load('GeneratorInterface.Core.genFilterSummary_cff')
#process.load('Configuration.StandardSequences.SimIdeal_cff')
#process.load('Configuration.StandardSequences.Digi_cff')
#process.load('Configuration.StandardSequences.SimL1Emulator_cff')
#process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")


#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.TFileService = cms.Service("TFileService",fileName=cms.string("test.root"))

#mylist = FileUtils.loadListFromFile('file_326483.in')
#readFiles = cms.untracked.vstring(*mylist)

#process.source = cms.Source("PoolSource",
#                                 fileNames = cms.untracked.vstring(
#                                    *mylist
#                                    )
#                                 )

process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring('file:/eos/cms/store/hidata/HIRun2018/HIMinimumBias1/RAW/v1/000/326/237/00000/1B8BA9F6-3805-FA4B-95B2-6FE754244FA5.root'))
#process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring('file:test_reco_py_RAW2DIGI_RECO.root'))

process.analyzer = cms.EDAnalyzer('rawZDCAnalyzer',
                                zdc = cms.InputTag("hcalDigis","ZDC"),
                                hltresults = cms.InputTag("TriggerResults","","HLT")
                                )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v10'
#process.GlobalTag.globaltag = '103X_dataRun2_Express_v2'
process.GlobalTag.globaltag = '103X_dataRun2_HLT_v1'

process.hcalDigis.InputLabel = cms.InputTag("rawDataRepacker")

process.digiPath = cms.Path(
    process.hcalDigis
#    process.zdcdigi
)

process.analyze_step = cms.Path(process.analyzer)

process.schedule = cms.Schedule(
                    process.digiPath,
                    process.analyze_step)
