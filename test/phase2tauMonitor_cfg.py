import FWCore.ParameterSet.Config as cms


process = cms.Process("RelValTauMonitor")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))# import of standard configurations

process.source = cms.Source("PoolSource",
            fileNames = cms.untracked.vstring(
                'file:/eos/cms/store/relval/CMSSW_15_1_0_pre2/RelValTenTau_15_500/MINIAODSIM/PU_150X_mcRun4_realistic_v1_STD_Run4D110_PU-v1/2580000/4cd57bf8-69b9-4f88-91d9-7d00f39b8cf0.root'
                # 'file:/eos/cms/store/relval/CMSSW_15_1_0_pre1/RelValTenTau_15_500/MINIAODSIM/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v1/2580000/8743514b-c2b1-4892-be26-74835563d546.root',
                # 'file:/eos/cms/store/relval/CMSSW_15_1_0_pre1/RelValTenTau_15_500/MINIAODSIM/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v1/2580000/9c0b2456-f716-4c60-852e-37c554e875fb.root'
                # 'file:/eos/cms/store/relval/CMSSW_15_0_0_pre3/RelValTenTau_15_500/MINIAODSIM/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v2/2580000/0f6bdb79-9a9b-48a6-9f73-04af042b3af8.root',
                # 'file:/eos/cms/store/relval/CMSSW_15_0_0_pre3/RelValTenTau_15_500/MINIAODSIM/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v2/2580000/bca9d7bf-901f-4aa7-949c-32b96e7db6ea.root',
                # 'file:/eos/cms/store/relval/CMSSW_15_0_0_pre3/RelValTenTau_15_500/MINIAODSIM/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v2/2580000/f23c44ce-6e03-473a-8b6a-41359c0dabdc.root'
            ))

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'141X_mcRun4_realistic_v3' , '')#'140X_dataRun3_Prompt_v3'




process.ditau = cms.EDAnalyzer('HLTTauPhase2Validator',
                                       channel = cms.string("ditau"),
                                       genParticles = cms.InputTag('prunedGenParticles'),
                                       triggerObjects = cms.InputTag('slimmedPatTrigger'),
                                       triggerBits = cms.InputTag('TriggerResults', '', 'HLT'),
                                       triggerfilter1 = cms.string('hltHpsPFTauTrack'),
                                       triggerfilter2 = cms.string('hltHpsDoublePFTau35MediumDitauWPDeepTau')
                               )
process.mutau = cms.EDAnalyzer('HLTTauPhase2Validator',
                               channel = cms.string("mutau"),
                                       genParticles = cms.InputTag('prunedGenParticles'),
                                       triggerObjects = cms.InputTag('slimmedPatTrigger'),
                                       triggerBits = cms.InputTag('TriggerResults', '', 'HLT'),
                                       triggerfilter1 = cms.string('hltL3crIsoL1TkSingleMu22TrkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk'),
                                       triggerfilter2 = cms.string('hltHpsPFTau27LooseTauWPDeepTau')
                               )
process.eletau = cms.EDAnalyzer('HLTTauPhase2Validator',
                                channel = cms.string("eletau"),
                                       genParticles = cms.InputTag('prunedGenParticles'),
                                       triggerObjects = cms.InputTag('slimmedPatTrigger'),
                                       triggerBits = cms.InputTag('TriggerResults', '', 'HLT'),
                                       triggerfilter1 = cms.string('hltEle30WPTightGsfTrackIsoL1SeededFilter'),
                                       triggerfilter2 = cms.string('hltHpsPFTau30LooseTauWPDeepTau')
                                )

process.TFileService = cms.Service("TFileService", fileName=cms.string('Test.root'))

process.p = cms.Path(process.ditau + process.mutau + process.eletau)
