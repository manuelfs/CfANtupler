#########################################################
### Configuration file to run cfA
#########################################################

### General assignments
import FWCore.ParameterSet.Config as cms
process = cms.Process("MinicfA")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
        'file:/afs/cern.ch/user/m/manuelf/work/store_mc_Spring14dr_QCD_Pt-10to20_EMEnriched_Tune4C_13TeV_pythia8_AODSIM_castor_PU_S14_POSTLS170_V6-v1_00000_0022A01C-E4E8-E311-97DA-003048FFD75C.root'
    )
)

### configurableAnalysis assignments
basicKinematicLeaves = cms.PSet(
    status = cms.string('status'),
    phi = cms.string('phi'),
    pt = cms.string('pt'),
    pz = cms.string('pz'),
    px = cms.string('px'),
    py = cms.string('py'),
    eta = cms.string('eta'),
    theta = cms.string('theta'),
    et = cms.string('et'),
    energy = cms.string('energy')
)

process.configurableAnalysis = cms.EDFilter("ConfigurableAnalysis",
    Selections = cms.PSet(
        filters = cms.PSet(),
        selections = cms.PSet(
            minSelection = cms.PSet(
                filterOrder = cms.vstring(''),
                makeFinalPlots = cms.bool(True),
                makeSummaryTable = cms.bool(True),
                makeContentPlots = cms.bool(True),
                makeAllButOnePlots = cms.bool(True),
                ntuplize = cms.bool(True),
                nMonitor = cms.uint32(1000),
                makeCumulativePlots = cms.bool(True)
            )
        )
    ),
    Plotter = cms.PSet(
        TH1s = cms.PSet(),
        ComponentName = cms.string('VariablePlotter'),
        TProfiles = cms.PSet(),
        TH2s = cms.PSet()
    ),
    Variables = cms.PSet(
        L1Bit = cms.PSet(
            src = cms.InputTag("gtDigis"),
            method = cms.string('ComputedVariable'),
            computer = cms.string('L1BitComputer')
        )
    ),
    workAsASelector = cms.bool(True),
    flows = cms.vstring('minSelection'),
    Ntupler = cms.PSet(
        branchesPSet = cms.PSet(
            treeName = cms.string('eventB'),
            mus = cms.PSet(
                src = cms.InputTag('muons'),
                leaves = cms.PSet(
                    basicKinematicLeaves,
                    #vars = cms.vstring(
                    #    'mus_pt:muons.pt',
                    #)
                ),
                Class = cms.string('reco::Muon')
            ),
        ),
        ComponentName = cms.string('CompleteNTupler'),
        useTFileService = cms.bool(True), ## false for EDM; true for non EDM
    )
)




process.InputTagDistributorService = cms.Service("InputTagDistributorService")
process.VariableHelperService = cms.Service("VariableHelperService")
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('configurableAnalysis.root')
)


process.p = cms.Path(process.configurableAnalysis)
