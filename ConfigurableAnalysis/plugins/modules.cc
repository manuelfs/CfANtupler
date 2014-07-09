#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"

//DEFINE_SEAL_MODULE();


#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
//#include "PhysicsTools/UtilAlgos/interface/SortCollectionSelector.h"
#include "CommonTools/UtilAlgos/interface/SortCollectionSelector.h"
#include "CommonTools/Utils/interface/PtComparator.h"

typedef ObjectSelector<SortCollectionSelector<reco::GenParticleCollection,GreaterByPt<reco::GenParticle> > > GenParticleSorterByPt;

DEFINE_FWK_MODULE(GenParticleSorterByPt);

//#include "CfANtupler/ConfigurableAnalysis/interface/ProcessIdSplitter.h"
//DEFINE_EDM_PLUGIN(CachingVariableFactory, ProcessIdSplitter, "ProcessIdSplitter");


#include "CfANtupler/ConfigurableAnalysis/interface/StringBasedNTupler.h"
DEFINE_EDM_PLUGIN(NTuplerFactory, StringBasedNTupler, "StringBasedNTupler");
#include "CfANtupler/ConfigurableAnalysis/interface/VariableNTupler.h"
DEFINE_EDM_PLUGIN(NTuplerFactory, VariableNTupler, "VariableNTupler");
#include "CfANtupler/ConfigurableAnalysis/interface/CompleteNTupler.h"
DEFINE_EDM_PLUGIN(NTuplerFactory, CompleteNTupler, "CompleteNTupler");
