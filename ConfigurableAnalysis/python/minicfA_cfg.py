###########################################################
### Configuration file to run cfA in CMSSW 7+ on miniAOD
###########################################################

### General assignments
import FWCore.ParameterSet.Config as cms
process = cms.Process("MinicfA")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
        #'/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/004C6DA7-FB03-E411-96BD-0025905A497A.root'
        #'file:/afs/cern.ch/user/m/manuelf/work/store_mc_Spring14dr_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8_AODSIM_castor_PU_S14_POSTLS170_V6-v1_00000_0022A01C-E4E8-E311-97DA-003048FFD75C.root'
    )
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'   # GT for 25ns (asymptotic alignment and calibration scenario)

### Loading branches
process.load("CfANtupler/ConfigurableAnalysis.branchesminicfA_cfi")

process.InputTagDistributorService = cms.Service("InputTagDistributorService")
process.VariableHelperService = cms.Service("VariableHelperService")
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('configurableAnalysis.root')
)


process.p = cms.Path(process.configurableAnalysis)
