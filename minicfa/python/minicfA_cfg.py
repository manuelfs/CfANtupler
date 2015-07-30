###########################################################
### Configuration file to run cfA in CMSSW 7+ on miniAOD
###########################################################

### General assignments
globalTags = {
# Global tag for early 2015 collisions data
    "PromptReco": "GR_P_V56::All", 
# GT for PHYS14 25ns (asymptotic alignment and calibration scenario)
    "PHYS14": "PHYS14_25_V3::All", 
# most 50 ns samples use asymptotic conditions, not startup conditions
#              "74X_MC_50ns_startup": "MCRUN2_74_V8::All", 
# standard set of conditions for 25 ns bunch spacing 74X MC
    "74X_MC-25ns": "MCRUN2_74_V9::All", 
# standard set of conditions for 50 ns bunch spacing 74X MC
    "74X_MC-50ns": "MCRUN2_74_V9A::All"}

######################################
# The following line must be changed #
######################################
datasetType = ""

import os
## Print out the cfA configuration information
if len(datasetType) == 0:
    print "Must set datasetType variable!"
    os._exit(1)
print "Using global tag " + globalTags[datasetType] + " selected from datasetType=" + datasetType

import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff import *
process = cms.Process("MinicfA")
process.load("FWCore.MessageService.MessageLogger_cfi")
# the following is done to avoid huge log files
process.MessageLogger.cerr.default.limit = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

##___________________________HCAL_Noise_Filter________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   taggingMode = cms.bool(True)
)


## process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
## process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True) )
## process.CPU = cms.Service("CPU")
## process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                            #'/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
                            #'/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/06843FC5-8370-E411-9B8C-0025905A60AA.root'
                            # '/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU40bx25_tsg_PHYS14_25_V1-v1/00000/06E41ADB-7870-E411-8850-0025905A605E.root'
                            # '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/004C6DA7-FB03-E411-96BD-0025905A497A.root'
                            #'file:/home/users/manuelf/cmssw/cfa/CMSSW_7_2_2_patch1/src/CfANtupler/minicfa/python/TT_Tune4C_13TeV-pythia8-tauola_MINIAODSIM.root'
                            #'/store/mc/Phys14DR/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/086903CE-2773-E411-A9B8-001E673967C5.root'
                            #'/store/mc/Phys14DR/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/682323F3-1774-E411-8F6A-002590A371D4.root'
                            '/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/06B5178E-F008-E511-A2CF-00261894390B.root'
                            #'/store/data/Run2015B/SingleMu/MINIAOD/PromptReco-v1/000/251/028/00000/705C6746-3C26-E511-92AC-02163E0139CF.root'
                            #'file:/home/users/manuelf/data/TT_Tune4C_13TeV-pythia8-tauola_MINIAODSIM.root'
                            #'file:/home/users/jbradmil/PHYS14/CMSSW_7_2_2_patch1/src/086903CE-2773-E411-A9B8-001E673967C5.root'
                                )
                                )



process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration/EventContent/EventContent_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process.GlobalTag.globaltag = globalTags[datasetType]   

### Loading branches
process.load("CfANtupler/minicfa.branchesminicfA_cfi")

process.InputTagDistributorService = cms.Service("InputTagDistributorService")
process.VariableHelperService = cms.Service("VariableHelperService")
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('cfA.root')
                                   )

## we don't really need this anymore since we save all pfcands
## process.trackIsolationMaker = cms.EDProducer("TrackIsolationMaker",
##                                              pfCandidatesTag = cms.InputTag("packedPFCandidates"),
##                                              vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
##                                              dR_ConeSize = cms.double(0.3),
##                                              dz_CutValue = cms.double(0.05),
##                                              minPt_PFCandidate = cms.double(5.0), #looser than the likely analysis selection
##                                              maxIso_PFCandidate = cms.double(0.25) #very loose
##                                              )

#Electron Identification for PHYS 14
## from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
## process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
## process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
## from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
## process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
## my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_cff']
## for idmod in my_id_modules:
##     setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


##     process.electronProducer = cms.EDProducer("ElectronProducer",
##                                               electronsInputTag = cms.InputTag("slimmedElectrons"),
##                                               electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-veto"),
##                                               electronLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-loose"),
##                                               electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-medium"),
##                                               electronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-tight")
##                                               )

process.jecCorL1Fast = cms.EDProducer("JECWrapper",
                                      jetsInputTag = cms.InputTag("slimmedJets"),
                                      jec_name = cms.string("ak4PFCHSL1Fastjet")
                                      )

process.jecCorL2L3 = cms.EDProducer("JECWrapper",
                                    jetsInputTag = cms.InputTag("slimmedJets"),
                                    jec_name = cms.string("ak4PFCHSL2L3")
                                    )

process.jecCorL1FastL2L3 = cms.EDProducer("JECWrapper",
                                          jetsInputTag = cms.InputTag("slimmedJets"),
                                          jec_name = cms.string("ak4PFCHSL1FastL2L3")
                                          )

process.photonProducer = cms.EDProducer("PhotonProducer",
                                        photonCollection =  cms.InputTag("slimmedPhotons"),
                                        electronCollection =  cms.InputTag("slimmedElectrons"),
                                        conversions = cms.InputTag("reducedEgamma", "reducedConversions"),
                                        beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),
                                        ecalRecHitsInputTag_EE = cms.InputTag("reducedEgamma","reducedEERecHits"),
                                        ecalRecHitsInputTag_EB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                        )

process.pdfProducer = cms.EDProducer("PDFProducer",
                                     genEventInfoInputTag = cms.string("generator"),
                                     hepmcHandle = cms.string("generator")
                                     )

from CfANtupler.Utils.trackIsolationMaker_cfi import trackIsolationFilter
from CfANtupler.Utils.trackIsolationMaker_cfi import trackIsolationCounter

process.IsolatedElectronTracksVeto = trackIsolationFilter.clone(
doTrkIsoVeto= False,
#vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
pfCandidatesTag = cms.InputTag("packedPFCandidates"),
dR_ConeSize = cms.double(0.3),
mini_ConeMax = cms.double(0.2),
mini_ConeMin = cms.double(0.05),
dz_CutValue = cms.double(0.1),
minPt_PFCandidate = cms.double(5.0),
isoCut = cms.double(0.2),
pdgId = cms.int32(11),
mTCut=cms.double(0)
)
process.IsolatedMuonTracksVeto = trackIsolationFilter.clone(
doTrkIsoVeto= False,
#vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
pfCandidatesTag = cms.InputTag("packedPFCandidates"),
dR_ConeSize = cms.double(0.3),
mini_ConeMax = cms.double(0.2),
mini_ConeMin = cms.double(0.05),
dz_CutValue = cms.double(0.1),
minPt_PFCandidate = cms.double(5.0),
isoCut = cms.double(0.2),
pdgId = cms.int32(13),
mTCut=cms.double(0)
)
process.IsolatedHadronicTracksVeto = trackIsolationFilter.clone(
doTrkIsoVeto= False,
#vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
pfCandidatesTag = cms.InputTag("packedPFCandidates"),
dR_ConeSize = cms.double(0.3),
mini_ConeMax = cms.double(0.2),
mini_ConeMin = cms.double(0.05),
dz_CutValue = cms.double(0.1),
minPt_PFCandidate = cms.double(10.0),
isoCut = cms.double(0.1),
pdgId = cms.int32(211),
mTCut=cms.double(0)
)

process.p = cms.Path(## process.trackIsolationMaker *
                     ## process.egmGsfElectronIDSequence *
                     ## process.electronProducer *
                     process.HBHENoiseFilterResultProducer*
                     process.jecCorL1Fast *
                     process.jecCorL2L3 *
                     process.jecCorL1FastL2L3 *
                     process.photonProducer *
                     process.pdfProducer *
                     process.IsolatedElectronTracksVeto *
                     process.IsolatedMuonTracksVeto *
                     process.IsolatedHadronicTracksVeto
                     )

process.outpath = cms.EndPath(cms.ignore(process.cfA))
