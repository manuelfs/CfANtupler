// -*- C++ -*-

#ifndef PHOTONPRODUCER_H
#define PHOTONPRODUCER_H

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "Math/VectorUtil.h"

#include <DataFormats/PatCandidates/interface/Photon.h>
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//
// class decleration
//

class PhotonProducer : public edm::EDProducer {
public:
     explicit PhotonProducer (const edm::ParameterSet&);
     ~PhotonProducer();

private:
  //  virtual void beginJob() ;
  virtual void beginJob() ;
  virtual void beginRun(edm::Run&, const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<std::vector<pat::Electron> > &eleCol,
				const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,
				float lxyMin=2.0, float probMin=1e-6, unsigned int nHitsBeforeVtxMax=0);

  
  // ----------member data ---------------------------
  edm::InputTag photonCollection_;
  edm::InputTag electronCollection_;
  edm::EDGetTokenT<std::vector<reco::Conversion> > conversionsToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::InputTag ecalRecHitsInputTag_EE_;
  edm::InputTag ecalRecHitsInputTag_EB_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalRecHitsInputTag_EE_Token_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalRecHitsInputTag_EB_Token_;

  // noZS::EcalClusterLazyTools* clusterTools_;

};

#endif

