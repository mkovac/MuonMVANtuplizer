import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("MuonMVANtuplizer")

LEPTON_SETUP = 2016

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


outputFile = "muon_ntuple.root"

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         '/store/mc/RunIISpring18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/100X_upgrade2018_realistic_v10-v1/70000/AC5CFA46-C52E-E811-8D7D-FA163E7FD764.root'
    )
)

#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#         '/store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver2-v1/40000/9EEF7AA6-938C-E811-B8D0-B083FED42A6D.root'
#    )
#)


process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
                                           src = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                           cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
                                           filter = cms.bool(True)
                                           )

process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                   src = cms.InputTag("slimmedMuons"),
                                   preselection = cms.string("track.isNonnull"),
                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                   fractionOfSharedSegments = cms.double(0.499)
                                   )
                                   
process.bareSoftMuons = cms.EDFilter("PATMuonRefSelector",
                                     src = cms.InputTag("cleanedMu"),
                                     cut = cms.string("pt > 5 && abs(eta) < 2.4 && (isGlobalMuon || (isTrackerMuon && numberOfMatches > 0)) && muonBestTrackType != 2")
                                     )

process.ntuplizer = cms.EDAnalyzer('MuonMVANtuplizer',
#                                   srcMiniAOD           = cms.InputTag('slimmedMuons'),
                                   srcMiniAOD           = cms.InputTag('bareSoftMuons'),
#                                   verticesMiniAOD      = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   verticesMiniAOD      = cms.InputTag('goodPrimaryVertices'),
                                   pileupMiniAOD        = cms.InputTag('slimmedAddPileupInfo'),
#                                   rhoMiniAOD           = cms.InputTag('fixedGridRhoFastjetAll'),
                                   genParticlesMiniAOD  = cms.InputTag('prunedGenParticles'),
                                   isMC                 = cms.bool(True),
                                   deltaR               = cms.double(0.1),
                                   ptThreshold          = cms.double(5),
#                                   setup                = cms.int32(LEPTON_SETUP)
                                   )
                                   
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )
                                   
#process.PVfilter = cms.Path(process.goodPrimaryVertices)        
process.muons =  cms.Sequence(process.cleanedMu + process.bareSoftMuons)
process.proc = cms.Path(process.goodPrimaryVertices + process.muons + process.ntuplizer)
