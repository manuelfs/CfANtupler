import FWCore.ParameterSet.Config as cms

metProducer = cms.EDProducer("METProducer",
                             METTag=cms.InputTag("slimmedMETs")
)


