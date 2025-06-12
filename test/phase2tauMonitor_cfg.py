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
                '/store/relval/CMSSW_15_1_0_pre1/RelValTenTau_15_500/MINIAODSIM/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v1/2580000/8743514b-c2b1-4892-be26-74835563d546.root',
                '/store/relval/CMSSW_15_1_0_pre1/RelValTenTau_15_500/MINIAODSIM/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v1/2580000/9c0b2456-f716-4c60-852e-37c554e875fb.root'

            ))

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'141X_dataRun3_v6' , '')#'140X_dataRun3_Prompt_v3'



process.hltDoubleMediumChargedIsoPFTauHPS40 = cms.EDAnalyzer('HLTTauPhase2Validator',
                                       channel = cms.string("ditau"),
                                       genParticles = cms.InputTag('prunedGenParticles'),
                                       triggerObjects = cms.InputTag('slimmedPatTrigger'),
                                       triggerBits = cms.InputTag('TriggerResults', '', 'HLT'),
                                       triggerfilter1 = cms.string('hltHpsPFTauTrack'),
                                       triggerfilter2 = cms.string('hltHpsDoublePFTau40TrackPt1MediumChargedIsolation')
                               )
process.hltDoubleMediumDeepTauPFTauHPS35 = cms.EDAnalyzer('HLTTauPhase2Validator',
                                       channel = cms.string("ditau"),
                                       genParticles = cms.InputTag('prunedGenParticles'),
                                       triggerObjects = cms.InputTag('slimmedPatTrigger'),
                                       triggerBits = cms.InputTag('TriggerResults', '', 'HLT'),
                                       triggerfilter1 = cms.string('hltHpsPFTauTrack'),
                                       triggerfilter2 = cms.string('hltHpsDoublePFTau35MediumDitauWPDeepTau')
                               )
process.hltIsoMu20eta2p1LooseDeepTauPFTauHPS27CrossL1 = cms.EDAnalyzer('HLTTauPhase2Validator',
                               channel = cms.string("mutau"),
                                       genParticles = cms.InputTag('prunedGenParticles'),
                                       triggerObjects = cms.InputTag('slimmedPatTrigger'),
                                       triggerBits = cms.InputTag('TriggerResults', '', 'HLT'),
                                       triggerfilter1 = cms.string('hltL3crIsoL1TkSingleMu22TrkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk'),
                                       triggerfilter2 = cms.string('hltHpsPFTau27LooseTauWPDeepTau')
                               )
process.hltEle30WPTightL1SeededLooseDeepTauPFTauHPS30CrossL1 = cms.EDAnalyzer('HLTTauPhase2Validator',
                                channel = cms.string("eletau"),
                                       genParticles = cms.InputTag('prunedGenParticles'),
                                       triggerObjects = cms.InputTag('slimmedPatTrigger'),
                                       triggerBits = cms.InputTag('TriggerResults', '', 'HLT'),
                                       triggerfilter1 = cms.string('hltEle30WPTightGsfTrackIsoL1SeededFilter'),
                                       triggerfilter2 = cms.string('hltHpsPFTau30LooseTauWPDeepTau')
                                )

process.TFileService = cms.Service("TFileService", fileName=cms.string('TenTau15_1_0_pre1_PU_141X_mcRun4_realistic_v3_STD_Run4D110.root'))

process.p = cms.Path(process.hltDoubleMediumChargedIsoPFTauHPS40 + process.hltDoubleMediumDeepTauPFTauHPS35 + process.hltIsoMu20eta2p1LooseDeepTauPFTauHPS27CrossL1 + process.hltEle30WPTightL1SeededLooseDeepTauPFTauHPS30CrossL1)
