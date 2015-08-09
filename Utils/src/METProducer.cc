// adapted from https://github.com/cmstas/NtupleMaker

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CfANtupler/Utils/interface/METProducer.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/PatCandidates/interface/MET.h"

typedef math::XYZTLorentzVectorF LorentzVector;

//
// constructors and destructor
//

METProducer::METProducer(const edm::ParameterSet& iConfig) {

  printf("METProducer\n");
  produces<float> ("pfmet"          ).setBranchAlias("pf_met"          );
  produces<float> ("pfmetPhi"       ).setBranchAlias("pf_metPhi"       );
  //produces<float> ("pfmetSig"       ).setBranchAlias("pf_metSig"       ); //this is just MET/sqrt(sumET). Use _pfmetSignificance unless you really want this
  produces<float> ("pfsumet"        ).setBranchAlias("pf_sumet"        );
  // produces<float> ("genmet"          ).setBranchAlias("gen_met"              );
  // produces<float> ("genmetPhi"       ).setBranchAlias("gen_metPhi"           );
  // produces<float> ("calomet"          ).setBranchAlias("calo_met"      );
  // produces<float> ("calometPhi"       ).setBranchAlias("calo_metPhi"   );

  MetInputTag_      = iConfig.getParameter<edm::InputTag>("METTag");
  //  isData_       = iConfig.getParameter<bool>     ( "isData"              );

}


METProducer::~METProducer() {}

void  METProducer::beginJob() {
}

void METProducer::endJob() {
}


// ------------ method called to produce the data  ------------
void METProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::auto_ptr<float>   pfmet         (new float   );
  std::auto_ptr<float>   pfmetPhi      (new float   );
  //    std::auto_ptr<float>   pfmetSig      (new float   ); //this is just MET/sqrt(sumET). Use pfmetSignificance unless you really want this branch
  std::auto_ptr<float>   pfsumet       (new float   );
  std::auto_ptr<float>   gen_met         (new float   );
  std::auto_ptr<float>   gen_metPhi      (new float   );
  std::auto_ptr<float>   calomet         (new float   );
  std::auto_ptr<float>   calometPhi      (new float   );

  //  printf("Get MET\n");

  edm::Handle< edm::View<pat::MET> > met_h;
  iEvent.getByLabel(MetInputTag_, met_h);

  // edm::Handle< edm::View<pat::MET> > genmet_h;
  // iEvent.getByLabel(MetInputTag_, genmet_h);
    
  if( !met_h.isValid() ) {
    throw cms::Exception("METProducer::produce: error getting particle-flow MET collection from Event!");
  }

  // if( !isData_ && !genmet_h.isValid() ) {
  //   throw cms::Exception("METProducer::produce: error getting gen particle-flow MET collection from Event!");
  // }

  *pfmet    = ( met_h->front() ).pt();
  // printf("MET = %3.2f\n", *pfmet);
  *pfmetPhi = ( met_h->front() ).phi();
  //    *pfmetSig = ( met_h->front() ).mEtSig();
  *pfsumet  = ( met_h->front() ).sumEt();       

  // if ( !isData_ ) {
  //   printf("Get gen MET\n");
  //   *gen_met      = ( genmet_h->front()).genMET()->pt();
  //   printf("Gen MET = %3.2f\n", *gen_met);
  //   *gen_metPhi   = ( genmet_h->front()).genMET()->phi();
  // }  
  // else {
  //   *gen_met      = -9999.;
  //   *gen_metPhi   = -9999.;
  // }
    
  // try {
  //   printf("Get calo MET\n");
  //   *calomet    = ( met_h->front() ).caloMETPt();
  //    printf("Gen MET = %3.2f\n", *calomet);
  //  *calometPhi = ( met_h->front() ).caloMETPhi();
  // }
  // catch ( cms::Exception& ex ) {
  //   *calomet    = -9999.;
  //   *calometPhi = -9999.;
  // }

  //  printf("Put MET back into event\n");
  iEvent.put(pfmet    , "pfmet"      );
  iEvent.put(pfmetPhi , "pfmetPhi"   );
  //    iEvent.put(pfmetSig , "pfmetSig"   );
  iEvent.put(pfsumet  , "pfsumet"    );  
  //iEvent.put(pfmetSignificance , "pfmetSignificance" );  
  // iEvent.put(gen_met      , "genmet"      );
  // iEvent.put(gen_metPhi   , "genmetPhi"   );
  // iEvent.put(calomet    , "calomet"      );
  // iEvent.put(calometPhi , "calometPhi"   );

}

//define this as a plug-in
DEFINE_FWK_MODULE(METProducer);
