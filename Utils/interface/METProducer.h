// adapted from https://github.com/cmstas/NtupleMaker

//
#ifndef METPRODUCER_H
#define METPRODUCER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class decleration
//

class METProducer : public edm::EDProducer {
 public:
  explicit METProducer (const edm::ParameterSet&);
  ~METProducer();

 private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag MetInputTag_;
  //    edm::InputTag pfMetCorInputTag;
  //  bool isData_;
};


#endif
