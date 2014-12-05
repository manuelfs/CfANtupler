// ADHOCNTUPLER: Creates eventA in the cfA ntuples, the tree that requires
//               Ad hoc c++ code to be filled.

#include <cmath>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/ConditionsInEdm.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>

using namespace std;
using namespace fastjet;

class miniAdHocNTupler : public NTupler {
 public:

  void fill(edm::Event& iEvent){

    nevents++;

    //////////////// Fat jets //////////////////
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByLabel("slimmedJets", jets);

    JetDefinition jet_def_12(antikt_algorithm, 1.2);
    //    vector<vector<PseudoJet>> fjets_vvector(0);
    vector<PseudoJet> fjets_constituents(0), fjets(0);
    vector<float> ptThresholds;
    ptThresholds.push_back(30);
    //const float etaThreshold(5);

    for(unsigned int ipt(0); ipt < ptThresholds.size(); ipt++){
      fjets_constituents.resize(0);
      //cout<<endl<<"SKINNY JETS"<<endl;
      for (unsigned int ijet(0); ijet < jets->size(); ijet++) {
	const pat::Jet &jet = (*jets)[ijet];
	//if(fabs(jet.eta()) > etaThreshold) continue;
	if(jet.pt() < ptThresholds[ipt]) continue;
	fjets_constituents.push_back(PseudoJet(jet.px(),jet.py(),jet.pz(),jet.energy()));
	//cout<<"pt "<<jet.pt()<<", eta "<<jet.eta()<<", phi "<<jet.phi()<<endl;
      }
      ClusterSequence cs_fjets(fjets_constituents, jet_def_12);
      fjets = sorted_by_pt(cs_fjets.inclusive_jets());
      for (unsigned int ifjet(0); ifjet < fjets.size(); ifjet++) {
	fjets30_pt->push_back(fjets[ifjet].pt());
	fjets30_eta->push_back(fjets[ifjet].eta());
	fjets30_phi->push_back(fjets[ifjet].phi());
	fjets30_energy->push_back(fjets[ifjet].E());
	fjets30_m->push_back(fjets[ifjet].m());
      }

//      cout<<endl<<"FAT JETS"<<endl;
//      for (unsigned int ifjet(0); ifjet < fjets.size(); ifjet++) {
//	cout<<"pt "<<fjets[ifjet].pt()<<", eta "<<fjets[ifjet].eta()<<", phi "<<fjets[ifjet].phi()<<endl;
//      }
//      fjets_vvector.push_back(fjets);
    }


    //////////////// pfcands shenanigans //////////////////
    edm::Handle<pat::PackedCandidateCollection> pfcands;
    iEvent.getByLabel("packedPFCandidates", pfcands);
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByLabel("slimmedMuons", muons);
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByLabel("slimmedElectrons", electrons);
    edm::Handle<pat::TauCollection> taus;
    iEvent.getByLabel("slimmedTaus", taus);

    vector<const pat::PackedCandidate*> el_pfmatch, mu_pfmatch; 
    for (const pat::PackedCandidate &pfc : *pfcands) {
      for (unsigned int ilep(0); ilep < electrons->size(); ilep++) {
	const pat::Electron &lep = (*electrons)[ilep];
	if(el_pfmatch.size() <= ilep) el_pfmatch.push_back(&pfc);
	else if(lep.pdgId()==pfc.pdgId() && deltaR(pfc, lep) < deltaR(*(el_pfmatch[ilep]), lep)) el_pfmatch[ilep] = &pfc;
      }
      for (unsigned int ilep(0); ilep < muons->size(); ilep++) {
	const pat::Muon &lep = (*muons)[ilep];
	if(mu_pfmatch.size() <= ilep) mu_pfmatch.push_back(&pfc);
	else if(lep.pdgId()==pfc.pdgId() && deltaR(pfc, lep) < deltaR(*(mu_pfmatch[ilep]), lep)) mu_pfmatch[ilep] = &pfc;
      }
    } // Loop over pfcands

    // Finding electron PF match
    for (unsigned int ilep(0); ilep < electrons->size(); ilep++) {
      const pat::Electron &lep = (*electrons)[ilep];
      els_isPF->push_back(deltaR(lep, *el_pfmatch[ilep]) < 0.1 && abs(lep.p()-el_pfmatch[ilep]->p())/lep.p()<0.05 &&
			  lep.pdgId() == el_pfmatch[ilep]->pdgId());
      els_jet_ind->push_back(-1);
    }

    // Finding muon PF match
    for (unsigned int ilep(0); ilep < muons->size(); ilep++) {
      const pat::Muon &lep = (*muons)[ilep];
      mus_isPF->push_back(lep.numberOfSourceCandidatePtrs()==1 && lep.sourceCandidatePtr(0)->pdgId()==lep.pdgId());
      mus_jet_ind->push_back(-1);
    }

    // Finding leptons in jets
    for (unsigned int ijet(0); ijet < jets->size(); ijet++) {
      const pat::Jet &jet = (*jets)[ijet];
      jets_AK4_mu_ind->push_back(-1);
      jets_AK4_el_ind->push_back(-1);

      float maxp(-99.), maxp_mu(-99.), maxp_el(-99.);
      int maxid(0);
      for (unsigned int i = 0, n = jet.numberOfSourceCandidatePtrs(); i < n; ++i) {
	const pat::PackedCandidate &pfc = dynamic_cast<const pat::PackedCandidate &>(*jet.sourceCandidatePtr(i));
	int pf_id = pfc.pdgId();
	float pf_p = pfc.p();
	if(pf_p > maxp){
	  maxp = pf_p;
	  maxid = pf_id;
	}

	if(abs(pf_id) == 11){
	  for (unsigned int ilep(0); ilep < electrons->size(); ilep++) {
	    if(&pfc == el_pfmatch[ilep]){
	      els_jet_ind->at(ilep) = ijet;
	      if(pf_p > maxp_el){
		maxp_el = pf_p;
		jets_AK4_el_ind->at(ijet) = ilep; // Storing the index of the highest pt electron in jet
	      }
	      break;
	    }
	  } // Loop over electrons
	} // If pfc is an electron

	if(abs(pf_id) == 13){
	  for (unsigned int ilep(0); ilep < muons->size(); ilep++) {
	    if(&pfc == mu_pfmatch[ilep]){
	      mus_jet_ind->at(ilep) = ijet;
	      if(pf_p > maxp_mu){
		maxp_mu = pf_p;
		jets_AK4_mu_ind->at(ijet) = ilep; // Storing the index of the highest pt muon in jet
	      }
	      break;
	    }
	  } // Loop over muons
	} // If pfc is an muon

      } // Loop over jet constituents
      jets_AK4_maxpt_id->push_back(maxid);
    } // Loop over jets

    // Finding leptons in taus
    for (unsigned int itau(0); itau < taus->size(); itau++) {
      const pat::Tau &tau = (*taus)[itau];
      taus_mu_ind->push_back(-1);
      taus_el_ind->push_back(-1);

      if(tau.numberOfSourceCandidatePtrs() == 1){
	const pat::PackedCandidate &pfc = dynamic_cast<const pat::PackedCandidate &> (*tau.sourceCandidatePtr(0));
	if(abs(pfc.pdgId())==11){
	  for (unsigned int ilep(0); ilep < electrons->size(); ilep++) {
	    if(&pfc == el_pfmatch[ilep]){
	      taus_el_ind->at(itau) = ilep;
	      break;
	    } 
	  } // Loop over electrons
	}
	if(abs(pfc.pdgId())==13){
	  for (unsigned int ilep(0); ilep < muons->size(); ilep++) {
	    if(&pfc == mu_pfmatch[ilep]){
	      taus_mu_ind->at(itau) = ilep;
	      break;
	    } 
	  } // Loop over electrons
	}
      } // If tau has one constituent
    } // Loop over taus

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


    //////////////// Filter decisions and names //////////////////
    edm::Handle<edm::TriggerResults> filterBits;
    edm::InputTag labfilterBits("TriggerResults","","PAT");
    iEvent.getByLabel(labfilterBits,filterBits);  
    int trackingfailurefilterResult(1);			    
    int goodVerticesfilterResult(1);				    
    int cschalofilterResult(1);						    
    int trkPOGfilterResult(1);						    
    int trkPOG_logErrorTooManyClustersfilterResult(1);	
    int EcalDeadCellTriggerPrimitivefilterResult(1);	
    int ecallaserfilterResult(1);						    
    int trkPOG_manystripclus53XfilterResult(1);		    
    int eebadscfilterResult(1);						    
    int METFiltersfilterResult(1);						    
    int HBHENoisefilterResult(1);						    
    int trkPOG_toomanystripclus53XfilterResult(1);		    
    int hcallaserfilterResult(1);
			               
    const edm::TriggerNames &fnames = iEvent.triggerNames(*filterBits);
    for (unsigned int i = 0, n = filterBits->size(); i < n; ++i) {
      string filterName = fnames.triggerName(i);
      int filterdecision = filterBits->accept(i);
      if (filterName=="Flag_trackingFailureFilter")		 trackingfailurefilterResult = filterdecision;
      if (filterName=="Flag_goodVertices")			 goodVerticesfilterResult = filterdecision;
      if (filterName=="Flag_CSCTightHaloFilter")		 cschalofilterResult = filterdecision;
      if (filterName=="Flag_trkPOGFilters")			 trkPOGfilterResult = filterdecision;
      if (filterName=="Flag_trkPOG_logErrorTooManyClusters")	 trkPOG_logErrorTooManyClustersfilterResult = filterdecision;
      if (filterName=="Flag_EcalDeadCellTriggerPrimitiveFilter") EcalDeadCellTriggerPrimitivefilterResult = filterdecision;
      if (filterName=="Flag_ecalLaserCorrFilter")		 ecallaserfilterResult = filterdecision;
      if (filterName=="Flag_trkPOG_manystripclus53X")		 trkPOG_manystripclus53XfilterResult = filterdecision;
      if (filterName=="Flag_eeBadScFilter")			 eebadscfilterResult = filterdecision;
      if (filterName=="Flag_METFilters")			 METFiltersfilterResult = filterdecision;
      if (filterName=="Flag_HBHENoiseFilter")			 HBHENoisefilterResult = filterdecision;
      if (filterName=="Flag_trkPOG_toomanystripclus53X")	 trkPOG_toomanystripclus53XfilterResult = filterdecision;
      if (filterName=="Flag_hcalLaserEventFilter")		 hcallaserfilterResult = filterdecision;
    }

    *trackingfailurefilter_decision_			=                trackingfailurefilterResult;	   
    *goodVerticesfilter_decision_			=		    goodVerticesfilterResult;	   
    *cschalofilter_decision_				=			    cschalofilterResult;   	    
    *trkPOGfilter_decision_				=			    trkPOGfilterResult;	   
    *trkPOG_logErrorTooManyClustersfilter_decision_	=  trkPOG_logErrorTooManyClustersfilterResult;  
    *EcalDeadCellTriggerPrimitivefilter_decision_	=    EcalDeadCellTriggerPrimitivefilterResult;    
    *ecallaserfilter_decision_				=			    ecallaserfilterResult; 	    
    *trkPOG_manystripclus53Xfilter_decision_		=	    trkPOG_manystripclus53XfilterResult;   
    *eebadscfilter_decision_				=			    eebadscfilterResult;   	    
    *METFiltersfilter_decision_				=			    METFiltersfilterResult;	    
    *HBHENoisefilter_decision_				=			    HBHENoisefilterResult; 	    
    *trkPOG_toomanystripclus53Xfilter_decision_		=	    trkPOG_toomanystripclus53XfilterResult;
    *hcallaserfilter_decision_				=                       hcallaserfilterResult;     

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


   //isolated pf candidates as found by TrackIsolationMaker                                                                               
    edm::Handle< vector<float> > pfcand_dzpv;
    iEvent.getByLabel("trackIsolationMaker","pfcandsdzpv", pfcand_dzpv);
    edm::Handle< vector<float> > pfcand_pt;
    iEvent.getByLabel("trackIsolationMaker","pfcandspt", pfcand_pt);
    edm::Handle< vector<float> > pfcand_eta;
    iEvent.getByLabel("trackIsolationMaker","pfcandseta", pfcand_eta);
    edm::Handle< vector<float> > pfcand_phi;
    iEvent.getByLabel("trackIsolationMaker","pfcandsphi", pfcand_phi);
    edm::Handle< vector<float> > pfcand_iso;
    iEvent.getByLabel("trackIsolationMaker","pfcandstrkiso", pfcand_iso);
    edm::Handle< vector<int> > pfcand_charge;
    iEvent.getByLabel("trackIsolationMaker","pfcandschg", pfcand_charge);

   for (size_t it=0; it<pfcand_pt->size(); ++it ) {
     isotk_pt_->push_back( pfcand_pt->at(it));
     isotk_phi_ -> push_back( pfcand_phi->at(it));
     isotk_eta_ -> push_back( pfcand_eta->at(it));
     isotk_iso_ -> push_back( pfcand_iso->at(it));
     isotk_dzpv_ -> push_back( pfcand_dzpv->at(it));
     isotk_charge_ -> push_back( pfcand_charge->at(it));
   }

   // tauID
    for (unsigned int itau(0); itau < taus->size(); itau++) {
      const pat::Tau &tau = (*taus)[itau];
      taus_byCombinedIsolationDeltaBetaCorrRaw3Hits_->push_back( tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") );
      taus_byLooseCombinedIsolationDeltaBetaCorr3Hits_->push_back( tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") );
      taus_byMediumCombinedIsolationDeltaBetaCorr3Hits_->push_back( tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") );
      taus_byTightCombinedIsolationDeltaBetaCorr3Hits_->push_back( tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") );
      taus_n_pfcands_->push_back( tau.numberOfSourceCandidatePtrs() );
      taus_decayMode_->push_back( tau.pfEssential().decayMode_ );
    } // Loop over taus

   

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

    (*els_isPF).clear();
    (*mus_isPF).clear();

    (*jets_AK4_maxpt_id).clear();
    (*jets_AK4_mu_ind).clear();
    (*jets_AK4_el_ind).clear();
    (*taus_el_ind).clear();
    (*taus_mu_ind).clear();
    (*els_jet_ind).clear();
    (*mus_jet_ind).clear();

    (*isotk_pt_).clear();
    (*isotk_phi_).clear();
    (*isotk_eta_).clear();
    (*isotk_iso_).clear();
    (*isotk_dzpv_).clear();
    (*isotk_charge_).clear();
    
    (*taus_byCombinedIsolationDeltaBetaCorrRaw3Hits_).clear();
    (*taus_byLooseCombinedIsolationDeltaBetaCorr3Hits_).clear();
    (*taus_byMediumCombinedIsolationDeltaBetaCorr3Hits_).clear();
    (*taus_byTightCombinedIsolationDeltaBetaCorr3Hits_).clear();
    (*taus_n_pfcands_).clear();
    (*taus_decayMode_).clear();

  }

  uint registerleaves(edm::ProducerBase * producer){
    uint nLeaves=0;
    if (useTFileService_){
      edm::Service<TFileService> fs;      
      if (ownTheTree_){
	ownTheTree_=true;
	tree_=fs->make<TTree>(treeName_.c_str(),"miniStringBasedNTupler tree");
      }else{
	TObject * object = fs->file().Get(treeName_.c_str());
	if (!object){
	  ownTheTree_=true;
	  tree_=fs->make<TTree>(treeName_.c_str(),"miniStringBasedNTupler tree");
	}
	tree_=dynamic_cast<TTree*>(object);
	if (!tree_){
	  ownTheTree_=true;
	  tree_=fs->make<TTree>(treeName_.c_str(),"miniStringBasedNTupler tree");
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

      tree_->Branch("trackingfailurefilter_decision", trackingfailurefilter_decision_ ,"trackingfailurefilter_decision/I");    
      tree_->Branch("goodVerticesfilter_decision", goodVerticesfilter_decision_	 ,"goodVerticesfilter_decision/I");
      tree_->Branch("cschalofilter_decision", cschalofilter_decision_,"cschalofilter_decision/I");			  
      tree_->Branch("trkPOGfilter_decision",  trkPOGfilter_decision_ ,"trkPOGfilter_decision/I");			  	 
      tree_->Branch("trkPOG_logErrorTooManyClustersfilter_decision", trkPOG_logErrorTooManyClustersfilter_decision_ ,"trkPOG_logErrorTooManyClustersfilter_decision/I");  	 
      tree_->Branch("EcalDeadCellTriggerPrimitivefilter_decision",	  EcalDeadCellTriggerPrimitivefilter_decision_	 ,"ecalDeadCellTriggerPrimitivefilter_decision/I");	  
      tree_->Branch("ecallaserfilter_decision",	ecallaserfilter_decision_ ,"ecallaserfilter_decision/I");			  
      tree_->Branch("trkPOG_manystripclus53Xfilter_decision",	  trkPOG_manystripclus53Xfilter_decision_	 ,"trkPOG_manystripclus53Xfilter_decision/I");	  
      tree_->Branch("eebadscfilter_decision",  eebadscfilter_decision_		 ,"eebadscfilter_decision/I");			  
      tree_->Branch("METFiltersfilter_decision", METFiltersfilter_decision_	 ,"METFiltersfilter_decision/I");
      tree_->Branch("HBHENoisefilter_decision",	 HBHENoisefilter_decision_ ,"HBHENoisefilter_decision/I");			  
      tree_->Branch("trkPOG_toomanystripclus53Xfilter_decision",	  trkPOG_toomanystripclus53Xfilter_decision_	 ,"trkPOG_toomanystripclus53Xfilter_decision/I");	  
      tree_->Branch("hcallaserfilter_decision",    hcallaserfilter_decision_,"hcallaserfilter_decision/I");   

      tree_->Branch("els_isPF",&els_isPF);
      tree_->Branch("mus_isPF",&mus_isPF);

      tree_->Branch("jets_AK4_maxpt_id", &jets_AK4_maxpt_id);
      tree_->Branch("jets_AK4_mu_ind",	&jets_AK4_mu_ind);  
      tree_->Branch("jets_AK4_el_ind",	&jets_AK4_el_ind);  
      tree_->Branch("taus_el_ind",	&taus_el_ind);      
      tree_->Branch("taus_mu_ind",	&taus_mu_ind);      
      tree_->Branch("els_jet_ind",	&els_jet_ind);      
      tree_->Branch("mus_jet_ind",	&mus_jet_ind);      
      tree_->Branch("isotk_pt",&isotk_pt_);
      tree_->Branch("isotk_phi",&isotk_phi_);
      tree_->Branch("isotk_eta",&isotk_eta_);
      tree_->Branch("isotk_iso",&isotk_iso_);
      tree_->Branch("isotk_dzpv",&isotk_dzpv_);
      tree_->Branch("isotk_charge",&isotk_charge_);

      tree_->Branch("taus_byCombinedIsolationDeltaBetaCorrRaw3Hits",&taus_byCombinedIsolationDeltaBetaCorrRaw3Hits_);
      tree_->Branch("taus_byLooseCombinedIsolationDeltaBetaCorr3Hits",&taus_byLooseCombinedIsolationDeltaBetaCorr3Hits_);
      tree_->Branch("taus_byMediumCombinedIsolationDeltaBetaCorr3Hits",&taus_byMediumCombinedIsolationDeltaBetaCorr3Hits_);
      tree_->Branch("taus_byTightCombinedIsolationDeltaBetaCorr3Hits",&taus_byTightCombinedIsolationDeltaBetaCorr3Hits_);
      tree_->Branch("taus_n_pfcands",&taus_n_pfcands_);
      tree_->Branch("taus_decayMode",&taus_decayMode_);

      tree_->Branch("fjets30_pt", &fjets30_pt);  
      tree_->Branch("fjets30_eta", &fjets30_eta);
      tree_->Branch("fjets30_phi", &fjets30_phi);
      tree_->Branch("fjets30_energy",&fjets30_energy);
      tree_->Branch("fjets30_m",&fjets30_m);


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

  miniAdHocNTupler (const edm::ParameterSet& iConfig){
    edm::ParameterSet adHocPSet = iConfig.getParameter<edm::ParameterSet>("AdHocNPSet");
    nevents = 0;

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

    trackingfailurefilter_decision_			= new int;   
    goodVerticesfilter_decision_			= new int;
    cschalofilter_decision_				= new int;    
    trkPOGfilter_decision_				= new int;
    trkPOG_logErrorTooManyClustersfilter_decision_	= new int;
    EcalDeadCellTriggerPrimitivefilter_decision_	= new int;
    ecallaserfilter_decision_				= new int;    
    trkPOG_manystripclus53Xfilter_decision_		= new int;
    eebadscfilter_decision_				= new int;    
    METFiltersfilter_decision_				= new int;    
    HBHENoisefilter_decision_				= new int;    
    trkPOG_toomanystripclus53Xfilter_decision_		= new int;
    hcallaserfilter_decision_				= new int;

    els_isPF = new std::vector<bool>;
    mus_isPF = new std::vector<bool>;

    jets_AK4_maxpt_id = new std::vector<int>;
    jets_AK4_mu_ind = new std::vector<int>;
    jets_AK4_el_ind = new std::vector<int>;
    taus_el_ind = new std::vector<int>;
    taus_mu_ind = new std::vector<int>;
    els_jet_ind = new std::vector<int>;
    mus_jet_ind = new std::vector<int>;

    //isolated tracks (charged pf candidates)                                                                                                               
    isotk_pt_ = new std::vector<float>;
    isotk_phi_ = new std::vector<float>;
    isotk_eta_ = new std::vector<float>;
    isotk_iso_ = new std::vector<float>;
    isotk_dzpv_ = new std::vector<float>;
    isotk_charge_ = new std::vector<int>;

    taus_byCombinedIsolationDeltaBetaCorrRaw3Hits_ = new std::vector<bool>;
    taus_byLooseCombinedIsolationDeltaBetaCorr3Hits_ = new std::vector<bool>;
    taus_byMediumCombinedIsolationDeltaBetaCorr3Hits_ = new std::vector<bool>;
    taus_byTightCombinedIsolationDeltaBetaCorr3Hits_ = new std::vector<bool>;
    taus_n_pfcands_ = new std::vector<int>;
    taus_decayMode_ = new std::vector<int>;
  
    fjets30_pt =     new std::vector<float>;
    fjets30_eta =    new std::vector<float>;
    fjets30_phi =    new std::vector<float>;
    fjets30_energy = new std::vector<float>;
    fjets30_m =     new std::vector<float>;

  }

  ~miniAdHocNTupler(){
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

    delete trackingfailurefilter_decision_		     ;
    delete goodVerticesfilter_decision_		     ;
    delete cschalofilter_decision_			     ;
    delete trkPOGfilter_decision_			     ;
    delete trkPOG_logErrorTooManyClustersfilter_decision_;
    delete EcalDeadCellTriggerPrimitivefilter_decision_  ;
    delete ecallaserfilter_decision_		     ;    
    delete trkPOG_manystripclus53Xfilter_decision_	     ;
    delete eebadscfilter_decision_			     ;
    delete METFiltersfilter_decision_		     ;    
    delete HBHENoisefilter_decision_		     ;    
    delete trkPOG_toomanystripclus53Xfilter_decision_    ;
    delete hcallaserfilter_decision_		     ;

    delete els_isPF;
    delete mus_isPF;

    delete jets_AK4_maxpt_id;
    delete jets_AK4_mu_ind;
    delete jets_AK4_el_ind;
    delete taus_el_ind;
    delete taus_mu_ind;
    delete els_jet_ind;
    delete mus_jet_ind;

    delete isotk_pt_;
    delete isotk_phi_;
    delete isotk_eta_;
    delete isotk_iso_;
    delete isotk_dzpv_;
    delete isotk_charge_;

    delete taus_byCombinedIsolationDeltaBetaCorrRaw3Hits_;
    delete taus_byLooseCombinedIsolationDeltaBetaCorr3Hits_;
    delete taus_byMediumCombinedIsolationDeltaBetaCorr3Hits_;
    delete taus_byTightCombinedIsolationDeltaBetaCorr3Hits_;
    delete taus_n_pfcands_;
    delete taus_decayMode_;
  
    delete fjets30_pt;
    delete fjets30_eta;
    delete fjets30_phi;
    delete fjets30_energy;
    delete fjets30_m;

  }

 private:
  bool ownTheTree_;
  std::string treeName_;
  bool useTFileService_;
  long nevents;


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
 
  int *trackingfailurefilter_decision_		     ;
  int *goodVerticesfilter_decision_		     ;
  int *cschalofilter_decision_			     ;
  int *trkPOGfilter_decision_			     ;
  int *trkPOG_logErrorTooManyClustersfilter_decision_;
  int *EcalDeadCellTriggerPrimitivefilter_decision_  ;
  int *ecallaserfilter_decision_		     ;    
  int *trkPOG_manystripclus53Xfilter_decision_	     ;
  int *eebadscfilter_decision_			     ;
  int *METFiltersfilter_decision_		     ;    
  int *HBHENoisefilter_decision_		     ;    
  int *trkPOG_toomanystripclus53Xfilter_decision_    ;
  int *hcallaserfilter_decision_		     ;

  std::vector<bool> * els_isPF;
  std::vector<bool> * mus_isPF;

  std::vector<int> * jets_AK4_maxpt_id;
  std::vector<int> * jets_AK4_mu_ind;
  std::vector<int> * jets_AK4_el_ind;

  std::vector<int> * taus_el_ind;
  std::vector<int> * taus_mu_ind;
  std::vector<int> * els_jet_ind;
  std::vector<int> * mus_jet_ind;

  std::vector<float> * isotk_pt_;
  std::vector<float> * isotk_phi_;
  std::vector<float> * isotk_eta_;
  std::vector<float> * isotk_iso_;
  std::vector<float> * isotk_dzpv_;
  std::vector<int> *   isotk_charge_;

  std::vector<bool> *  taus_byCombinedIsolationDeltaBetaCorrRaw3Hits_;
  std::vector<bool> *  taus_byLooseCombinedIsolationDeltaBetaCorr3Hits_;
  std::vector<bool> *  taus_byMediumCombinedIsolationDeltaBetaCorr3Hits_;
  std::vector<bool> *  taus_byTightCombinedIsolationDeltaBetaCorr3Hits_;
  std::vector<int> *  taus_n_pfcands_;
  std::vector<int> *  taus_decayMode_;
 
  std::vector<float> * fjets30_pt;
  std::vector<float> * fjets30_eta;
  std::vector<float> * fjets30_phi;
  std::vector<float> * fjets30_energy;
  std::vector<float> * fjets30_m;


};
