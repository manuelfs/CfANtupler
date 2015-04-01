// -*- C++ -*-

#ifndef JECWRAPPER_H
#define JECWRAPPER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "Math/VectorUtil.h"

//
// class decleration
//

class JECWrapper : public edm::EDProducer {
public:
     explicit JECWrapper (const edm::ParameterSet&);
     ~JECWrapper();

private:
  //  virtual void beginJob() ;
  virtual void beginJob() ;
  virtual void beginRun(edm::Run&, const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  
  // ----------member data ---------------------------

  edm::InputTag jetsTag_;

  std::string jec_name_;
  const pat::JetCollection *jets;

};

#endif

