###########################################################
### Configuration file to run cfA in CMSSW 7+ on miniAOD
###########################################################

### General assignments
import FWCore.ParameterSet.Config as cms
process = cms.Process("MinicfA")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
        #'/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/06843FC5-8370-E411-9B8C-0025905A60AA.root'
        # '/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU40bx25_tsg_PHYS14_25_V1-v1/00000/06E41ADB-7870-E411-8850-0025905A605E.root'
        # '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/004C6DA7-FB03-E411-96BD-0025905A497A.root'
        #'file:/home/users/manuelf/cmssw/cfa/CMSSW_7_2_2_patch1/src/CfANtupler/minicfa/python/TT_Tune4C_13TeV-pythia8-tauola_MINIAODSIM.root'
        'file:/home/users/manuelf/data/TT_Tune4C_13TeV-pythia8-tauola_MINIAODSIM.root'
    )
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'   # GT for 25ns (asymptotic alignment and calibration scenario)

### Loading branches
process.load("CfANtupler/minicfa.branchesminicfA_cfi")

process.InputTagDistributorService = cms.Service("InputTagDistributorService")
process.VariableHelperService = cms.Service("VariableHelperService")
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('cfA.root')
)


process.trackIsolationMaker = cms.EDProducer("TrackIsolationMaker",
                                             pfCandidatesTag = cms.InputTag("packedPFCandidates"),
                                             vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                             dR_ConeSize = cms.double(0.3),
                                             dz_CutValue = cms.double(0.05),
                                             minPt_PFCandidate = cms.double(5.0), #looser than the likely analysis selection
                                             maxIso_PFCandidate = cms.double(0.25) #very loose
)

process.p = cms.Path(process.trackIsolationMaker)
process.outpath = cms.EndPath(cms.ignore(process.cfA))
