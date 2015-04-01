//-----------------------------------------------------------------------------------------
//
// Since I can't seem to add an edm::EventSetup to AdHocNTupler,
// I'm creating a separate module to produce the JEC.
//
//-----------------------------------------------------------------------------------------
// system include files

#include <memory>
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

//For jet corrections
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "CfANtupler/minicfa/interface/JECWrapper.h"
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
JECWrapper::JECWrapper(const edm::ParameterSet& iConfig) {
  jetsTag_ = iConfig.getParameter<InputTag> ("jetsInputTag");
  jec_name_ = iConfig.getParameter<std::string>( "jec_name" );
  produces<vector<float> >("JEC").setBranchAlias("jets_JEC");
}
JECWrapper::~JECWrapper()
{
}
void JECWrapper::beginRun(edm::Run&, const edm::EventSetup& es) {}
void JECWrapper::beginJob() {}
void JECWrapper::endJob() {}
// ------------ method called to produce the data ------------
void JECWrapper::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto_ptr<vector<float> > jets_JEC(new vector<float>);
  //---------------------------------
  // get jets
  //---------------------------------
  edm::Handle<pat::JetCollection> jetsHandle;
  iEvent.getByLabel(jetsTag_, jetsHandle);
  jets = jetsHandle.product();

  // cout << "getJetCorrector( " << jec_name_ << ", iSetup)" << endl;

  const JetCorrector* corrector = JetCorrector::getJetCorrector ( jec_name_ , iSetup );

  for (unsigned int ijet(0); ijet < jets->size(); ijet++) { 
    const pat::Jet &jet = (*jets)[ijet];
    float jets_correction = corrector -> correction( jet, iEvent, iSetup );
    jets_JEC->push_back(jets_correction);
  }
  iEvent.put(jets_JEC,"JEC");
}
//define this as a plug-in
DEFINE_FWK_MODULE(JECWrapper);
