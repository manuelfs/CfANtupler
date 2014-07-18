// ADHOCNTUPLER: Creates eventA in the cfA ntuples, the tree that requires
//               Ad hoc c++ code to be filled.

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
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

using namespace std;

class AdHocNTupler : public NTupler {
 public:

  void fill(edm::Event& iEvent){

    //////////////// Pile up information //////////////////
    if(!iEvent.isRealData()) { //Access PU info in MC
      edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
      iEvent.getByLabel("addPileupInfo", PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator PVI;

      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
	// cout << " PU Information: bunch crossing " << PVI->getBunchCrossing() 
	//      << ", NumInteractions " << PVI->getPU_NumInteractions() 
	//      << ", TrueNumInteractions " << PVI->getTrueNumInteractions() 
	//      <<", evend ID   "<< iEvent.id().event() << std::endl;
	(*PU_NumInteractions_).push_back(PVI->getPU_NumInteractions());
	(*PU_bunchCrossing_).push_back(PVI->getBunchCrossing());
	(*PU_TrueNumInteractions_).push_back(PVI->getTrueNumInteractions());
	(*PU_zpositions_).push_back(PVI->getPU_zpositions());
	(*PU_sumpT_lowpT_).push_back(PVI->getPU_sumpT_lowpT());
	(*PU_sumpT_highpT_).push_back(PVI->getPU_sumpT_highpT());
	(*PU_ntrks_lowpT_).push_back(PVI->getPU_ntrks_lowpT());
	(*PU_ntrks_highpT_).push_back(PVI->getPU_ntrks_highpT());
      }
    } // if it's not real data


    //////////////// Trigger decisions and names //////////////////
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::InputTag labtriggerBits("TriggerResults","","HLT");
    iEvent.getByLabel(labtriggerBits,triggerBits);  

    //////////////// Trigger prescales //////////////////
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::InputTag labtriggerPrescales("patTrigger");
    iEvent.getByLabel(labtriggerPrescales,triggerPrescales);  

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      (*trigger_decision).push_back(triggerBits->accept(i));
      (*trigger_name).push_back(names.triggerName(i));
      (*trigger_prescalevalue).push_back(triggerPrescales->getPrescaleForIndex(i));
    }
   
    //////////////// HLT trigger objects //////////////////
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::InputTag labtriggerObjects("selectedPatTrigger");
    iEvent.getByLabel(labtriggerObjects,triggerObjects);  

    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(names);
      (*standalone_triggerobject_collectionname).push_back(obj.collection()); 
      (*standalone_triggerobject_pt).push_back(obj.pt());
      (*standalone_triggerobject_px).push_back(obj.px());
      (*standalone_triggerobject_py).push_back(obj.py());
      (*standalone_triggerobject_pz).push_back(obj.pz());
      (*standalone_triggerobject_et).push_back(obj.et());
      (*standalone_triggerobject_energy).push_back(obj.energy());
      (*standalone_triggerobject_phi).push_back(obj.phi());
      (*standalone_triggerobject_eta).push_back(obj.eta());
    }

    //////////////// L1 trigger objects --- TO BE UNDERSTOOD ---
    edm::Handle<L1GlobalTriggerReadoutRecord> L1trigger_h;
    edm::InputTag labL1trigger("gtDigis","","HLT");
    iEvent.getByLabel(labL1trigger, L1trigger_h);  

    std::vector<bool> gtbits;
    int ngtbits = 128;
    gtbits.reserve(ngtbits); for(int i=0; i<ngtbits; i++) gtbits[i]=false;
    if(L1trigger_h.isValid()) 
      gtbits = L1trigger_h->decisionWord();
    for(int i=0; i<ngtbits; i++) if(gtbits[i]) cout<<"Bit "<<i<<" is true"<<endl;

    const L1GlobalTriggerReadoutRecord* L1trigger = L1trigger_h.failedToGet () ? 0 : &*L1trigger_h;
    if(L1trigger) cout<<"Level 1 decision: "<<L1trigger->decision()<<endl;

    //fill the tree    
    if (ownTheTree_){ tree_->Fill(); }
    (*trigger_name).clear();
    (*trigger_decision).clear();
    (*trigger_prescalevalue).clear();
    (*standalone_triggerobject_pt).clear();
    (*standalone_triggerobject_px).clear();
    (*standalone_triggerobject_py).clear();
    (*standalone_triggerobject_pz).clear();
    (*standalone_triggerobject_et).clear();
    (*standalone_triggerobject_energy).clear();
    (*standalone_triggerobject_phi).clear();
    (*standalone_triggerobject_eta).clear();
    (*standalone_triggerobject_collectionname).clear();

    (*PU_zpositions_).clear();
    (*PU_sumpT_lowpT_).clear();
    (*PU_sumpT_highpT_).clear();
    (*PU_ntrks_lowpT_).clear();
    (*PU_ntrks_highpT_).clear();
    (*PU_NumInteractions_).clear();
    (*PU_bunchCrossing_).clear();
    (*PU_TrueNumInteractions_).clear();
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
      tree_->Branch("standalone_triggerobject_pt",&standalone_triggerobject_pt);
      tree_->Branch("standalone_triggerobject_px",&standalone_triggerobject_px);
      tree_->Branch("standalone_triggerobject_py",&standalone_triggerobject_py);
      tree_->Branch("standalone_triggerobject_pz",&standalone_triggerobject_pz);
      tree_->Branch("standalone_triggerobject_et",&standalone_triggerobject_et);
      tree_->Branch("standalone_triggerobject_energy",&standalone_triggerobject_energy);
      tree_->Branch("standalone_triggerobject_phi",&standalone_triggerobject_phi);
      tree_->Branch("standalone_triggerobject_eta",&standalone_triggerobject_eta);
      tree_->Branch("standalone_triggerobject_collectionname",&standalone_triggerobject_collectionname);

      tree_->Branch("PU_zpositions",&PU_zpositions_);
      tree_->Branch("PU_sumpT_lowpT",&PU_sumpT_lowpT_);
      tree_->Branch("PU_sumpT_highpT",&PU_sumpT_highpT_);
      tree_->Branch("PU_ntrks_lowpT",&PU_ntrks_lowpT_);
      tree_->Branch("PU_ntrks_highpT",&PU_ntrks_highpT_);
      tree_->Branch("PU_NumInteractions",&PU_NumInteractions_);
      tree_->Branch("PU_bunchCrossing",&PU_bunchCrossing_);
      tree_->Branch("PU_TrueNumInteractions",&PU_TrueNumInteractions_);
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
    standalone_triggerobject_pt = new std::vector<float>;
    standalone_triggerobject_px = new std::vector<float>;
    standalone_triggerobject_py = new std::vector<float>;
    standalone_triggerobject_pz = new std::vector<float>;
    standalone_triggerobject_et = new std::vector<float>;
    standalone_triggerobject_energy = new std::vector<float>;
    standalone_triggerobject_phi = new std::vector<float>;
    standalone_triggerobject_eta = new std::vector<float>;
    standalone_triggerobject_collectionname = new std::vector<std::string>;

    PU_zpositions_ = new std::vector<std::vector<float> >;
    PU_sumpT_lowpT_ = new std::vector<std::vector<float> >;
    PU_sumpT_highpT_ = new std::vector<std::vector<float> >;
    PU_ntrks_lowpT_ = new std::vector<std::vector<int> >;
    PU_ntrks_highpT_ = new std::vector<std::vector<int> >;
    PU_NumInteractions_ = new std::vector<int>;
    PU_bunchCrossing_ = new std::vector<int>;
    PU_TrueNumInteractions_ = new std::vector<float>;
  }

  ~AdHocNTupler(){
    delete trigger_decision;
    delete trigger_name;
    delete trigger_prescalevalue;
    delete standalone_triggerobject_pt;
    delete standalone_triggerobject_px;
    delete standalone_triggerobject_py;
    delete standalone_triggerobject_pz;
    delete standalone_triggerobject_et;
    delete standalone_triggerobject_energy;
    delete standalone_triggerobject_phi;
    delete standalone_triggerobject_eta;
    delete standalone_triggerobject_collectionname;

    delete PU_zpositions_;
    delete PU_sumpT_lowpT_;
    delete PU_sumpT_highpT_;
    delete PU_ntrks_lowpT_;
    delete PU_ntrks_highpT_;
    delete PU_NumInteractions_;
    delete PU_bunchCrossing_;
    delete PU_TrueNumInteractions_;
  }

 private:
  bool ownTheTree_;
  std::string treeName_;
  bool useTFileService_;

  std::vector<bool> * trigger_decision;
  std::vector<std::string> * trigger_name;
  std::vector<float> * trigger_prescalevalue;
  std::vector<float> * standalone_triggerobject_pt;
  std::vector<float> * standalone_triggerobject_px;
  std::vector<float> * standalone_triggerobject_py;
  std::vector<float> * standalone_triggerobject_pz;
  std::vector<float> * standalone_triggerobject_et;
  std::vector<float> * standalone_triggerobject_energy;
  std::vector<float> * standalone_triggerobject_phi;
  std::vector<float> * standalone_triggerobject_eta;
  std::vector<std::string> * standalone_triggerobject_collectionname;

  std::vector<std::vector<float> > * PU_zpositions_;
  std::vector<std::vector<float> > * PU_sumpT_lowpT_;
  std::vector<std::vector<float> > * PU_sumpT_highpT_;
  std::vector<std::vector<int> > * PU_ntrks_lowpT_;
  std::vector<std::vector<int> > * PU_ntrks_highpT_;
  std::vector<int> * PU_NumInteractions_;
  std::vector<int> * PU_bunchCrossing_;
  std::vector<float> * PU_TrueNumInteractions_;
 
};
