#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/ConditionsInEdm.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

using namespace std;

class AdHocNTupler : public NTupler {
 public:

  void fill(edm::Event& iEvent){

    // Trigger decisions and names
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::InputTag labtriggerBits("TriggerResults","","HLT");
    iEvent.getByLabel(labtriggerBits,triggerBits);  

    // Trigger prescales
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::InputTag labtriggerPrescales("patTrigger");
    iEvent.getByLabel(labtriggerPrescales,triggerPrescales);  

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      (*trigger_decision).push_back(triggerBits->accept(i));
      (*trigger_name).push_back(names.triggerName(i));
      (*trigger_prescalevalue).push_back(triggerPrescales->getPrescaleForIndex(i));
    }
   

   //fill the tree    
    if (ownTheTree_){ tree_->Fill(); }
    (*trigger_name).clear();
    (*trigger_decision).clear();
    (*trigger_prescalevalue).clear();
  }

  uint registerleaves(edm::ProducerBase * producer){
    uint nLeaves=0;
    if (useTFileService_){
      edm::Service<TFileService> fs;      
      if (ownTheTree_){
	ownTheTree_=true;
	tree_=fs->make<TTree>(treeName_.c_str(),"StringBasedNTupler tree");
      }else{
	TObject * object = fs->file().Get(treeName_.c_str());
	if (!object){
	  ownTheTree_=true;
	  tree_=fs->make<TTree>(treeName_.c_str(),"StringBasedNTupler tree");
	}
	tree_=dynamic_cast<TTree*>(object);
	if (!tree_){
	  ownTheTree_=true;
	  tree_=fs->make<TTree>(treeName_.c_str(),"StringBasedNTupler tree");
	}
      }
      
      //register the leaves by hand
      tree_->Branch("trigger_decision",&trigger_decision);
      tree_->Branch("trigger_name",&trigger_name);
      tree_->Branch("trigger_prescalevalue",&trigger_prescalevalue);
    }

    else{
      //EDM COMPLIANT PART
      //      producer->produce<ACertainCollection>(ACertainInstanceName);
    }


    return nLeaves;
  }

  void callBack(){
    //clean up whatever memory was allocated
  }

  AdHocNTupler (const edm::ParameterSet& iConfig){
    edm::ParameterSet adHocPSet = iConfig.getParameter<edm::ParameterSet>("AdHocNPSet");

    if (adHocPSet.exists("useTFileService"))
      useTFileService_=adHocPSet.getParameter<bool>("useTFileService");         
    else
      useTFileService_=iConfig.getParameter<bool>("useTFileService");

    if (useTFileService_){
      if (adHocPSet.exists("treeName")){
	treeName_=adHocPSet.getParameter<std::string>("treeName");
	ownTheTree_=true;
      }
      else{
	treeName_=iConfig.getParameter<std::string>("treeName");
	ownTheTree_=false;
      }
    }
    

    trigger_decision = new std::vector<bool>;
    trigger_name = new std::vector<std::string>;
    trigger_prescalevalue = new std::vector<float>;

  }

  ~AdHocNTupler(){
    delete trigger_decision;
    delete trigger_name;
    delete trigger_prescalevalue;
  }

 private:
  bool ownTheTree_;
  std::string treeName_;
  bool useTFileService_;

  std::vector<bool> * trigger_decision;
  std::vector<std::string> * trigger_name;
  std::vector<float> * trigger_prescalevalue;

};
