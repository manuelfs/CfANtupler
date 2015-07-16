//-----------------------------------------------------------------------------------------
//
// Computation of the trackIsolation, for use with the isolated track veto 
// used for the stop quark search in the single lepton channel
// Author: Ben Hooberman
// Adapted by Arne-Rasmus Drager, Jack Bradmiller-Feld, & Kevin Pedro
//
//
//-----------------------------------------------------------------------------------------

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
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

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CfANtupler/Utils/interface/TrackIsolationFilter.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
// miniAOD
//#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"

//#include "CfANtupler/minicfa/interface/getPFIsolation.h"




using namespace reco;
using namespace edm;
using namespace std;

//
// class decleration
//

//
// constructors and destructor
//

TrackIsolationFilter::TrackIsolationFilter(const edm::ParameterSet& iConfig) {

  pfCandidatesTag_		= iConfig.getParameter<InputTag>("pfCandidatesTag");
  //	vertexInputTag_               = iConfig.getParameter<InputTag>("vertexInputTag");
  MetInputTag_      = iConfig.getParameter<InputTag>("METTag");
  dR_               = iConfig.getParameter<double>          ("dR_ConeSize");       // dR value used to define the isolation cone                (default 0.3 )
  miniIsoMax_               = iConfig.getParameter<double>          ("mini_ConeMax");       // max cone size for mini isolation                  (default 0.2 )
  miniIsoMin_               = iConfig.getParameter<double>          ("mini_ConeMin");       // min cone size for mini isolation                  (default 0.05 )
  dzcut_            = iConfig.getParameter<double>          ("dz_CutValue");       // cut value for dz(trk,vtx) for track to include in iso sum (default 0.05)
  minPt_            = iConfig.getParameter<double>          ("minPt_PFCandidate"); // store PFCandidates with pt above this cut                 (default 0   )
  isoCut_           = iConfig.getParameter<double>          ("isoCut"); // isolation cut value
  doTrkIsoVeto_     = iConfig.getParameter<bool>            ("doTrkIsoVeto");
  pdgId_     = iConfig.getParameter<int>            ("pdgId");
  mTCut_=   iConfig.getParameter<double>   ("mTCut");
  maxEta_= iConfig.getParameter<double>("etaCut");
  debug_= iConfig.getParameter<bool>("debug");
	
  produces<std::vector<pat::PackedCandidate> >(""); 
  //	produces<bool>("GoodVtx");
  produces<vector<TLorentzVector> >("pfcands");
  produces<vector<double> >("pfcandstrkiso").setBranchAlias("pfcands_trkiso");
  produces<vector<double> >("pfcandsminiso").setBranchAlias("pfcands_miniso_chg_only");
  produces<vector<double> >("pfcandsminisochgonly").setBranchAlias("pfcands_miniso_chg_only");
  produces<vector<double> >("pfcandsdzpv"  ).setBranchAlias("pfcands_dzpv");
  produces<vector<double> >("pfcandsmT"    ).setBranchAlias("pfcands_mT");
  produces<vector<int>   >("pfcandschg"    ).setBranchAlias("pfcands_chg");
  produces<vector<int>   >("pfcandsfromPV"    ).setBranchAlias("pfcands_fromPV");
  produces<vector<int>   >("pfcandsid"     ).setBranchAlias("pfcands_id");
  produces<int>("isoTracks").setBranchAlias("isoTracks");

}

TrackIsolationFilter::~TrackIsolationFilter() {

}

// ------------ method called to produce the data  ------------

bool TrackIsolationFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<vector<TLorentzVector> > pfcands(new vector<TLorentzVector>);
  auto_ptr<vector<double> >  pfcands_trkiso(new vector<double>);
  auto_ptr<vector<double> >  pfcands_miniso(new vector<double>);
  auto_ptr<vector<double> >  pfcands_miniso_chg_only(new vector<double>);
  auto_ptr<vector<double> >  pfcands_dzpv  (new vector<double>);
  auto_ptr<vector<double> >  pfcands_mT    (new vector<double>);
  auto_ptr<vector<int>   >  pfcands_chg    (new vector<int>  );
  auto_ptr<vector<int>   >  pfcands_fromPV    (new vector<int>  );
  auto_ptr<vector<int>   >  pfcands_id     (new vector<int>  );
	
  edm::Handle< edm::View<pat::MET> > MET;
  iEvent.getByLabel(MetInputTag_,MET);
  reco::MET::LorentzVector metLorentz(0,0,0,0);
  if(MET.isValid() ){
    metLorentz=MET->at(0).p4();
  }

  //---------------------------------
  // get PFCandidate collection
  //---------------------------------
  
  edm::Handle<PFCandidateCollection> pfCandidatesHandle;
  //	edm::Handle<edm::View< pat::PackedCandidate> >pfCandidatesHandle;
  iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
  
  edm::Handle<edm::View< pat::PackedCandidate> >pfCandidates;
  iEvent.getByLabel(pfCandidatesTag_, pfCandidates);

  // Do I really need another one of these--to be compatible with mini isolation fucntion?
  edm::Handle<pat::PackedCandidateCollection> packedcands;
  iEvent.getByLabel("packedPFCandidates", packedcands);
  //---------------------------------
  // get Vertex Collection
  //---------------------------------
	
  // edm::Handle<edm::View<reco::Vertex> > vertices;
  // iEvent.getByLabel(vertexInputTag_, vertices);
  // //reco::Vertex::Point vtxpos = (vertices->size() > 0 ? (*vertices)[0].position() : reco::Vertex::Point());
  // int firstGoodVertexIdx = -1;
  // vtxSize = vertices->size();
  // for(int v=0; v<vtxSize;++v){
  // 	if ( !(*vertices)[v].isFake() && (*vertices)[v].ndof()>4. && (*vertices)[v].position().Rho()<2.0 && fabs((*vertices)[v].position().Z())<24.0) {
  // 		firstGoodVertexIdx=v;
  // 		//std::cout<<"GoodVertexIdV found: "<<firstGoodVertexIdx<<std::endl;
  // 		break;
  // 	}
  // }
  // bool hasGoodVtx = false;
  // if(vertices->size() > 0 && firstGoodVertexIdx>=0) hasGoodVtx = true;
  // std::auto_ptr<bool> GoodVtx(new bool(hasGoodVtx));
	
  //-------------------------------------------------------------------------------------------------
  // loop over PFCandidates and calculate the trackIsolation and dz w.r.t. 1st good PV for each one
  // for neutral PFCandidates, store trkiso = 999 and dzpv = 999
  //-------------------------------------------------------------------------------------------------
	
  std::auto_ptr<std::vector<pat::PackedCandidate> > prodminiAOD(new std::vector<pat::PackedCandidate>());
   
  // miniAOD
  for(size_t i=0; i<pfCandidates->size();i++)
    {
      const pat::PackedCandidate pfCand = (*pfCandidates)[i];
		
      //calculated mT value
      double dphiMET = fabs(pfCand.phi()-metLorentz.phi());
      double mT = sqrt(2 *metLorentz.pt() * pfCand.pt() * (1 - cos(dphiMET)));
		
      //to keep track of cuts in debug case (when continues are not used)
      bool goodCand = true;
		
      //-------------------------------------------------------------------------------------
      // skip events with no good vertices
      //-------------------------------------------------------------------------------------
      // if(!hasGoodVtx) {
      // 	if(debug_) goodCand &= false;
      // 	else continue;
      // }
      //-------------------------------------------------------------------------------------
      // only store PFCandidate values if PFCandidate.pdgId() == pdgId_
      //-------------------------------------------------------------------------------------
      if( pdgId_ != 0 && abs( pfCand.pdgId() ) != pdgId_ ) {
	if(debug_) goodCand &= false;
	else continue;
      }
      //-------------------------------------------------------------------------------------
      // only store PFCandidate values if pt > minPt
      //-------------------------------------------------------------------------------------
      if(pfCand.pt() <minPt_) {
	if(debug_) goodCand &= false;
	else continue;
      }
      if(fabs(pfCand.eta()) >maxEta_) {
	if(debug_) goodCand &= false;
	else continue;
      }
      //-------------------------------------------------------------------------------------
      // cut on mT of track and MET
      //-------------------------------------------------------------------------------------
      if(mTCut_>0.01 && mT>mTCut_) {
	if(debug_) goodCand &= false;
	else continue;
      }
		
      //store candidate values
      //(all values stored in debug case, otherwise just good candidates are stored)
      TLorentzVector p4(pfCand.px(),pfCand.py(),pfCand.pz(),pfCand.energy());
      pfcands->push_back(p4);
      pfcands_chg->push_back(pfCand.charge());
      pfcands_fromPV->push_back(pfCand.fromPV());
      pfcands_id->push_back(pfCand.pdgId());
      pfcands_mT->push_back(mT);
      //      std::cout << "Computing track mini isolation" << std::endl;
      pfcands_miniso->push_back(GetTrackMiniIsolation(p4, packedcands, abs(pfCand.pdgId()), miniIsoMin_, miniIsoMax_, 10., false));
      pfcands_miniso_chg_only->push_back(GetTrackMiniIsolation(p4, packedcands, abs(pfCand.pdgId()), miniIsoMin_, miniIsoMax_, 10., true));
      //      printf("miniso %3.2f, chgonly %3.2f\n", pfcands_miniso->back(), pfcands_miniso_chg_only->back());
      //loop over other PFCandidates in the event to calculate track isolation
      float trkiso = 0.0;
      float dz_it = 100;
      if( pfCand.charge() != 0 )
	{
	  for(size_t ii=0; ii<pfCandidates->size();ii++)
	    {  
	      // don't count the PFCandidate in its own isolation sum
	      if(i==ii) continue;
	      // require the PFCandidate to be charged
	      const pat::PackedCandidate otherCand = (*pfCandidates)[ii];
	      if( otherCand.charge() == 0 ) continue;
	      // cut on dR between the PFCandidates
	      float dR = deltaR(pfCand.eta(), pfCand.phi(), otherCand.eta(), otherCand.phi());
	      if( dR > dR_ ) continue;
	      // cut on the PFCandidate dz
	      float dz_other = otherCand.dz();
	      if( fabs(dz_other) > dzcut_ ) continue;
	      trkiso += otherCand.pt();
	    }
	  dz_it = pfCand.dz();
	}
      else { //neutral particle, store default values
	trkiso = 9999;
	dz_it = 9999;
      }
      pfcands_trkiso->push_back(trkiso);
      pfcands_dzpv->push_back(dz_it);
		
		
      //----------------------------------------------------------------------------
      // now make cuts on isolation and dz (for charged candidates only)
      //----------------------------------------------------------------------------		
      // if( pfCand.charge() != 0 )
      // {
      // 	if( debug_ && !goodCand) continue;
      // 	if( trkiso/pfCand.pt() > isoCut_ ) continue;
      // 	if( std::abs(dz_it) > dzcut_ ) continue;
			
      // 	//only keep tracks that pass all cuts
      // 	prodminiAOD->push_back( pfCand );
      // }
      //else neutral particle, nothing to do
		
    }

  bool result = (doTrkIsoVeto_ ? (prodminiAOD->size() == 0) : true);
	
  int isoTracks=prodminiAOD->size();
  //  std::cout << "Found " << isoTracks << " isoTracks" << std::endl;
  std::auto_ptr<int> htp(new int(isoTracks));
  iEvent.put(htp,"isoTracks" );

  //  cout << "calculated mini iso for " << pfcands_miniso->size() << " tracks" << endl;
  // put candidate values back into event
  iEvent.put(pfcands       ,"pfcands"      );
  iEvent.put(pfcands_trkiso,"pfcandstrkiso");
  iEvent.put(pfcands_miniso,"pfcandsminiso");
  iEvent.put(pfcands_miniso_chg_only,"pfcandsminisochgonly");
  iEvent.put(pfcands_dzpv  ,"pfcandsdzpv"  );
  iEvent.put(pfcands_mT    ,"pfcandsmT"    );
  iEvent.put(pfcands_chg   ,"pfcandschg"   );
  iEvent.put(pfcands_fromPV,"pfcandsfromPV");
  iEvent.put(pfcands_id    ,"pfcandsid"    );
  //	iEvent.put(GoodVtx       ,"GoodVtx"      );
	
  iEvent.put(prodminiAOD); 
  //  std::cout << "Putting isoTracks back into the event" << std::endl;
  return result;
}

double TrackIsolationFilter::GetTrackMiniIsolation(const TLorentzVector pfCandp4, edm::Handle<pat::PackedCandidateCollection> pfcands,                        
						   const int pdgId, double r_iso_min, double r_iso_max, double kt_scale,
						   bool charged_only) {

  
  if (pfCandp4.Pt()<5.) return 99999.;

  double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  double ptThresh(0.5);

  switch(abs(pdgId)) {
  case 11:
    ptThresh = 0;
    if(fabs(pfCandp4.Eta())>1.479){
      deadcone_ch = 0.015;
      deadcone_pu = 0.015;
      deadcone_ph = 0.08;
    }
    break;
  case 13:
    deadcone_ch = 0.0001;
    deadcone_pu = 0.01;
    deadcone_ph = 0.01;
    deadcone_nh = 0.01;
    break;
  default:
    deadcone_ch = 0.0001;
    deadcone_pu = 0.01;
    deadcone_ph = 0.01;
    deadcone_nh = 0.01; // Using muon cones for hadronic tracks
    break;
  }
  
  double iso_nh(0.); double iso_ch(0.); 
  double iso_ph(0.); double iso_pu(0.);
    
  double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/pfCandp4.Pt()));
  //  printf("Track 1: pt %3.2f, eta %3.2f, phi %3.2f, E %3.2f, r_iso %3.2f\n", pfCandp4.Pt(), pfCandp4.Eta(), pfCandp4.Phi(), pfCandp4.E(), r_iso);

  for (const pat::PackedCandidate &pfc : *pfcands) {
    if (abs(pfc.pdgId())<7) continue;


    double dr = deltaR(pfc.p4().Eta(), pfc.p4().Phi(), pfCandp4.Eta(), pfCandp4.Phi());
    if (pfCandp4.Pt()-pfc.pt()==0.&&dr==0.) continue; // don't double-count tracks
    //    printf("Track 2: pt %3.2f, eta %3.2f, phi %3.2f, E %3.2f, dR %3.2f\n", pfc.pt(), pfc.eta(), pfc.phi(), pfc.energy(), dr);
    if (dr > r_iso) continue;

    //////////////////  NEUTRALS  /////////////////////////
    if (pfc.charge()==0){
      if (pfc.pt()>ptThresh) {
	/////////// PHOTONS ////////////
	if (abs(pfc.pdgId())==22) {
	  if(dr < deadcone_ph) continue;
	  iso_ph += pfc.pt();
	  /////////// NEUTRAL HADRONS ////////////
	} else if (abs(pfc.pdgId())==130) {
	  if(dr < deadcone_nh) continue;
	  iso_nh += pfc.pt();
	}
      }
      //////////////////  CHARGED from PV  /////////////////////////
    } else if (pfc.fromPV()>1){
      if (abs(pfc.pdgId())==211) {
	if(dr < deadcone_ch) continue;
	iso_ch += pfc.pt();
      }
      //////////////////  CHARGED from PU  /////////////////////////
    } else {
      if (pfc.pt()>ptThresh){
	if(dr < deadcone_pu) continue;
	iso_pu += pfc.pt();
      }
    }
  }
  double iso(0.);
  if (charged_only){
    iso = iso_ch;
  } else {
    iso = iso_ph + iso_nh;
    iso -= 0.5*iso_pu;
    if (iso>0) iso += iso_ch;
    else iso = iso_ch;
  }
  iso = iso/pfCandp4.Pt();

  //  printf("ch %3.2f, nh %3.2f, ph %3.2f, pu %3.2f\n", iso_ch, iso_nh, iso_ph, iso_pu);

  return iso;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackIsolationFilter);

