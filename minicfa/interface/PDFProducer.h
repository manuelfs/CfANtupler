// -*- C++ -*-

#ifndef PHOTONPRODUCER_H
#define PDFPRODUCER_H

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


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//
// class decleration
//

class PDFProducer : public edm::EDProducer {
public:
     explicit PDFProducer (const edm::ParameterSet&);
     ~PDFProducer();

private:
  //  virtual void beginJob() ;
  virtual void beginJob() ;
  virtual void beginRun(edm::Run&, const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  std::string genEventInfoInputTag_;
  std::string hepmcHandle_;

};

#endif

