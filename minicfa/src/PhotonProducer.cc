//-----------------------------------------------------------------------------------------
//
// Get photon branches that can't be obtained from string-based ntupler.
//
//-----------------------------------------------------------------------------------------
// system include files

#include <memory>
#include <vector>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"


#include "CfANtupler/minicfa/interface/PhotonProducer.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"
typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// class declaration
//
//
// constructors and destructor
//
PhotonProducer::PhotonProducer(const edm::ParameterSet& iConfig) {
  photonCollection_ = iConfig.getParameter<InputTag>("photonCollection");
  electronCollection_ = iConfig.getParameter<InputTag>("electronCollection");
  conversionsToken_ =  consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"));
  beamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  ecalRecHitsInputTag_EE_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE");
  ecalRecHitsInputTag_EB_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB");
  ecalRecHitsInputTag_EE_Token_ = consumes<EcalRecHitCollection>(ecalRecHitsInputTag_EE_);
  ecalRecHitsInputTag_EB_Token_ = consumes<EcalRecHitCollection>(ecalRecHitsInputTag_EB_);
  produces<vector<float> >("photonsfull5x5sigmaIEtaIEta").setBranchAlias("photons_full5x5sigmaIEtaIEta");
  produces<vector<bool> >("photonspasselveto").setBranchAlias("photons_pass_el_veto");
}
PhotonProducer::~PhotonProducer()
{
  if (clusterTools_) delete clusterTools_;
}
void PhotonProducer::beginRun(edm::Run&, const edm::EventSetup& es) {}
void PhotonProducer::beginJob() {}
void PhotonProducer::endJob() {}
// ------------ method called to produce the data ------------
void PhotonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto_ptr<vector<float> > photons_full5x5sigmaIEtaIEta(new vector<float>);
  auto_ptr<vector<bool> > photons_pass_el_veto(new vector<bool>);
  //---------------------------------
  // Full sigmaIEtaIEta
  // apparently have to do this by hand in 7_2--https://github.com/cmstas/NtupleMaker/blob/d2d7682cbecf57abe25b5ce00a13ed51a47aecd4/src/PhotonMaker.cc#L336 */
  //---------------------------------

  // Get photons, electrons, beam spot

edm::Handle<pat::PhotonCollection> photons;
 iEvent.getByLabel(photonCollection_, photons);
 edm::Handle<pat::ElectronCollection> electrons;
 iEvent.getByLabel(electronCollection_, electrons);
 edm::Handle<vector<reco::Conversion> > conversions;
 iEvent.getByToken(conversionsToken_,conversions);
 edm::Handle<reco::BeamSpot> beamSpot;
 iEvent.getByToken(beamSpotToken_,beamSpot);

  clusterTools_ = new EcalClusterLazyTools(iEvent, iSetup, ecalRecHitsInputTag_EB_Token_, ecalRecHitsInputTag_EE_Token_);

  for (unsigned int iphoton(0); iphoton < photons->size(); iphoton++) {
    const pat::Photon &photon = (*photons)[iphoton];
    std::vector<float> vCov = clusterTools_->localCovariances( *(photon.superCluster()->seed()) ); 
    const float sieie = (isnan(vCov[0]) ? 0. : sqrt(vCov[0])); 
    photons_full5x5sigmaIEtaIEta->push_back(sieie);
    // Electron Veto
    photons_pass_el_veto->push_back(!hasMatchedPromptElectron(photon.superCluster(),electrons, conversions, beamSpot->position()));
  }
  
  // put everything back into event
  iEvent.put(photons_full5x5sigmaIEtaIEta,"photonsfull5x5sigmaIEtaIEta");
  iEvent.put(photons_pass_el_veto,"photonspasselveto");
}

bool PhotonProducer::hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<std::vector<pat::Electron> > &eleCol,
					      const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,
					      float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax) {
  // copied from https://github.com/RazorCMS/SUSYBSMAnalysis-RazorTuplizer/blob/6072ffb43bbeb3f6b34cf8a96426c7f104c5b902/plugins/RazorAux.cc#L127
  //check if a given SuperCluster matches to at least one GsfElectron having zero expected inner hits
  //and not matching any conversion in the collection passing the quality cuts
  if (sc.isNull()) return false;
  for (std::vector<pat::Electron>::const_iterator it = eleCol->begin(); it!=eleCol->end(); ++it) {
    //match electron to supercluster
    if (it->superCluster()!=sc) continue;
    //check expected inner hits
    if (it->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
    //check if electron is matching to a conversion
    if (ConversionTools::hasMatchedConversion(*it,convCol,beamspot,lxyMin,probMin,nHitsBeforeVtxMax)) continue;
    return true;
  }
  return false;
}
//define this as a plug-in
DEFINE_FWK_MODULE(PhotonProducer);
