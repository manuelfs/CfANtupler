#include "CfANtupler/minicfa/interface/miniStringBasedNTupler.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Hemisphere.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/METReco/interface/MET.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/METReco/interface/HcalNoiseRBX.h"

#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include <DataFormats/CaloRecHit/interface/CaloCluster.h>

#include <DataFormats/PatCandidates/interface/TriggerPath.h>

#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>


//--------------------------------------------------------------------------------
//just define here a list of objects you would like to be able to have a branch of
//--------------------------------------------------------------------------------
#define MINIANOTHER_VECTOR_CLASS(C) if (class_==#C) return StringBranchHelper<C>(*this, iEvent)()
#define MINIANOTHER_CLASS(C) if (class_==#C) return StringLeaveHelper<C>(*this, iEvent)()

miniTreeBranch::value miniTreeBranch::branch(const edm::Event& iEvent){
  MINIANOTHER_VECTOR_CLASS(pat::Jet);
  else MINIANOTHER_VECTOR_CLASS(pat::Muon);
  else MINIANOTHER_VECTOR_CLASS(reco::GenParticle);
  else MINIANOTHER_VECTOR_CLASS(pat::Electron);
  else MINIANOTHER_VECTOR_CLASS(pat::MET);
  else MINIANOTHER_VECTOR_CLASS(pat::Tau);
  else MINIANOTHER_VECTOR_CLASS(pat::Hemisphere);
  else MINIANOTHER_VECTOR_CLASS(pat::Photon);
  else MINIANOTHER_VECTOR_CLASS(reco::CaloMET);
  else MINIANOTHER_VECTOR_CLASS(reco::Muon);
  else MINIANOTHER_VECTOR_CLASS(reco::Track);
  else MINIANOTHER_VECTOR_CLASS(reco::GsfElectron);
  else MINIANOTHER_VECTOR_CLASS(SimTrack);
  else MINIANOTHER_VECTOR_CLASS(l1extra::L1ParticleMap);
  else MINIANOTHER_VECTOR_CLASS(reco::Vertex);
  else MINIANOTHER_VECTOR_CLASS(pat::GenericParticle);
  else MINIANOTHER_VECTOR_CLASS(reco::MET);
  else MINIANOTHER_CLASS(edm::HepMCProduct);
  else MINIANOTHER_CLASS(reco::BeamSpot);
  else MINIANOTHER_CLASS(HcalNoiseSummary);
  else MINIANOTHER_CLASS(GenEventInfoProduct);
  else MINIANOTHER_VECTOR_CLASS(reco::HcalNoiseRBX);
  else MINIANOTHER_VECTOR_CLASS(reco::BasicJet);
  else MINIANOTHER_VECTOR_CLASS(reco::CaloJet);
  else MINIANOTHER_VECTOR_CLASS(reco::GenJet);
  else MINIANOTHER_VECTOR_CLASS(pat::TriggerPath);
  else MINIANOTHER_VECTOR_CLASS(reco::PFCandidate);
  else MINIANOTHER_VECTOR_CLASS(reco::CaloCluster);
  else MINIANOTHER_VECTOR_CLASS(reco::Photon);
  else MINIANOTHER_VECTOR_CLASS(pat::PackedCandidate);
  else MINIANOTHER_VECTOR_CLASS(pat::PackedGenParticle);
  else {
    edm::LogError("miniTreeBranch")<<branchName()<<" failed to recognize class type: "<<class_<<". Shucks";
    return miniTreeBranch::value(new std::vector<float>());
  }
}
#undef MINIANOTHER_CLASS
