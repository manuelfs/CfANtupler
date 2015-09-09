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

JECs = {
# Global tag for early 2015 collisions data
    "PromptReco": "Summer15_50nsV4_DATA", 
# GT for PHYS14 25ns (asymptotic alignment and calibration scenario)
    "PHYS14": "PHYS14_25_V2", 
# most 50 ns samples use asymptotic conditions, not startup conditions
#              "74X_MC_50ns_startup": "MCRUN2_74_V8::All", 
# standard set of conditions for 25 ns bunch spacing 74X MC
    "74X_MC-25ns": "Summer15_25nsV2_MC", 
# standard set of conditions for 50 ns bunch spacing 74X MC
    "74X_MC-50ns": "Summer15_50nsV2_MC"}

######################################
# The following line must be changed #
######################################
datasetType = "74X_MC-25ns"

import os
## Print out the cfA configuration information
if len(datasetType) == 0:
    print "Must set datasetType variable!"
    os._exit(1)
print "Using global tag " + globalTags[datasetType] + " selected from datasetType=" + datasetType
print "Using JECs: " + JECs[datasetType]

collisionData = False
applyResidual = False
if datasetType=='PromptReco':
    collisionData = True
    applyResidual = True
    print "Applying residual JECs"

    
import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff import *
process = cms.Process("MinicfA")
process.load("FWCore.MessageService.MessageLogger_cfi")
# the following is done to avoid huge log files
process.MessageLogger.cerr.default.limit = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#process.dump=cms.EDAnalyzer('EventContentAnalyzer')


## process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
## process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True) )
## process.CPU = cms.Service("CPU")
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    #wantSummary = cms.untracked.bool(True)
)
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                             'file:pickevents.root'
        # '/store/mc/RunIISpring15DR74/TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/30000/0C1CDF5E-3B3B-E511-84DE-0CC47A4D99E6.root'        
                            #'/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
                            #'/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/06843FC5-8370-E411-9B8C-0025905A60AA.root'
                            # '/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU40bx25_tsg_PHYS14_25_V1-v1/00000/06E41ADB-7870-E411-8850-0025905A605E.root'
                            # '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/004C6DA7-FB03-E411-96BD-0025905A497A.root'
                            #'file:/home/users/manuelf/cmssw/cfa/CMSSW_7_2_2_patch1/src/CfANtupler/minicfa/python/TT_Tune4C_13TeV-pythia8-tauola_MINIAODSIM.root'
                            #'/store/mc/Phys14DR/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/086903CE-2773-E411-A9B8-001E673967C5.root'
                            #'/store/mc/Phys14DR/GJets_HT-600toInf_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/682323F3-1774-E411-8F6A-002590A371D4.root'
                            #'/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/06B5178E-F008-E511-A2CF-00261894390B.root'
                            #'/store/data/Run2015B/SingleMu/MINIAOD/PromptReco-v1/000/251/028/00000/705C6746-3C26-E511-92AC-02163E0139CF.root'
                            #'file:HTMHT_2A828FCA-182C-E511-9AF1-02163E01299A.root'
                            #'file:04412314-B92E-E511-97A6-002618943981.root'
                            #'file:/home/users/jbradmil/commissioningDPS/CMSSW_7_4_6_patch6/src/086903CE-2773-E411-A9B8-001E673967C5.root'
                                )
                                )



process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration/EventContent/EventContent_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')


from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process.GlobalTag.globaltag = globalTags[datasetType]


 ## ----------------------------------------------------------------------------------------------
## JECs
## ----------------------------------------------------------------------------------------------

# get the JECs (disabled by default)
# this requires the user to download the .db file from this twiki
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
JetTag = cms.InputTag('slimmedJets')
#if len(jecfile)>0:
JECPatch = cms.string('sqlite_file:'+JECs[datasetType]+'.db')
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
        connect = JECPatch,
        toGet   = cms.VPSet(
                cms.PSet(
                        record = cms.string("JetCorrectionsRecord"),
                        tag    = cms.string("JetCorrectorParametersCollection_"+JECs[datasetType]+"_AK4PFchs"),
                        label  = cms.untracked.string("AK4PFchs")
                ),
                cms.PSet(
                         record = cms.string("JetCorrectionsRecord"),
                         tag    = cms.string("JetCorrectorParametersCollection_"+JECs[datasetType]+"_AK4PF"),
                         label  = cms.untracked.string("AK4PF")
                )
        )
)
process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
    src     = cms.InputTag("slimmedJets"),
    levels  = ['L1FastJet',
               'L2Relative',
               'L3Absolute'],
               payload = 'AK4PFchs' # Make sure to choose the appropriate levels and payload here!
)
if applyResidual: process.patJetCorrFactorsReapplyJEC.levels.append('L2L3Residual')
        
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = patJetsUpdated.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
)
        
        #process.Baseline += process.patJetCorrFactorsReapplyJEC
        #process.Baseline += process.patJetsReapplyJEC
        
JetTag = cms.InputTag('patJetsReapplyJEC')

## Correct MET 
## ref: https://github.com/cms-met/cmssw/blob/METCorUnc74X/PhysicsTools/PatAlgos/test/corMETFromMiniAOD.py
METTag = cms.InputTag('slimmedMETs')
#    if is74X:
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=collisionData, # controls gen met
                           pfCandColl=cms.InputTag("packedPFCandidates"),
                           postfix="TypeICorr"
)
if collisionData:
    if not applyResidual: #skip residuals for data if not used
        process.patPFMetT1T2CorrTypeICorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT1T2SmearCorrTypeICorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT2CorrTypeICorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT2SmearCorrTypeICorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.shiftedPatJetEnDownTypeICorr.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
        process.shiftedPatJetEnUpTypeICorr.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
corrMETTag = cms.InputTag('slimmedMETsTypeICorr')


# use only PF cands with |eta|<3 (no HF) -- only for 74X
process.noHFCands = cms.EDFilter("CandPtrSelector",
                                 src=cms.InputTag("packedPFCandidates"),
                                 cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
                                 )
runMetCorAndUncFromMiniAOD(process,
                           isData=collisionData, # controls gen met
                           pfCandColl=cms.InputTag("noHFCands"),
                           postfix="TypeICorrNoHF"
)
if collisionData:
    if not applyResidual: #skip residuals for data if not used
        process.patPFMetT1T2CorrTypeICorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT1T2SmearCorrTypeICorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT2CorrTypeICorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.patPFMetT2SmearCorrTypeICorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
        process.shiftedPatJetEnDownTypeICorrNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
        process.shiftedPatJetEnUpTypeICorrNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
            

##___________________________HCAL_Noise_Filter________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

## process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
##    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
##    taggingMode = cms.bool(True)
## )




### Loading branches
process.load("CfANtupler/minicfa.branchesminicfA_cfi")

process.InputTagDistributorService = cms.Service("InputTagDistributorService")
process.VariableHelperService = cms.Service("VariableHelperService")
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('cfA.root')
                                   )




process.jecCorL1Fast = cms.EDProducer("JECWrapper",
                                      jetsInputTag = cms.InputTag("patJetsReapplyJEC"),
                                      jec_name = cms.string("ak4PFCHSL1Fastjet")
                                      )

process.jecCorL2L3 = cms.EDProducer("JECWrapper",
                                    jetsInputTag = cms.InputTag("patJetsReapplyJEC"),
                                    jec_name = cms.string("ak4PFCHSL2L3")
                                    )

process.jecCorL1FastL2L3 = cms.EDProducer("JECWrapper",
                                          jetsInputTag = cms.InputTag("patJetsReapplyJEC"),
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

from CfANtupler.Utils.metProducer_cfi import metProducer
process.pfTypeIMETLatestJEC = cms.EDProducer("METProducer",
                                             METTag=cms.InputTag("slimmedMETsTypeICorr")
)
process.pfTypeIMETLatestJECNoHF = cms.EDProducer("METProducer",
                                             METTag=cms.InputTag("slimmedMETsTypeICorrNoHF")
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
mTCut=cms.double(0),
METTag=corrMETTag
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
mTCut=cms.double(0),
METTag=corrMETTag
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
mTCut=cms.double(0),
METTag=corrMETTag
)

process.p = cms.Path(process.patJetCorrFactorsReapplyJEC*
                     process.patJetsReapplyJEC*
                     process.HBHENoiseFilterResultProducer*
                     #process.ApplyBaselineHBHENoiseFilter*
                     process.jecCorL1Fast *
                     process.jecCorL2L3 *
                     process.jecCorL1FastL2L3 *
                     process.photonProducer *
                     process.pdfProducer *
                     process.pfTypeIMETLatestJEC*
                     process.pfTypeIMETLatestJECNoHF*
                     process.IsolatedElectronTracksVeto *
                     process.IsolatedMuonTracksVeto *
                     process.IsolatedHadronicTracksVeto
#                     process.dump
                     )

process.outpath = cms.EndPath(cms.ignore(process.cfA))
