//-----------------------------------------------------------------------------------------
//
// Computation of the trackIsolation, for use with the isolated track veto 
// used for the stop quark search in the single lepton channel
// Original Author: Ben Hooberman
// Adapted for subsequent use in the EWKino->hh(bbbb) search by Josh Thompson
// Adapted for compatibility with miniAOD by Jack Bradmiller-Feld
//
// For each PFCandidate above threshold and passing an isolation cut, store:
// pT, Eta, Phi, Isolation value
// charge of PFCandidate
// dz of PFCandidate w.r.t. the 1st good vertex
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

#include "CfANtupler/IsoTrackFinder/interface/TrackIsolationMaker.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// class decleration
//

//
// constructors and destructor
//

TrackIsolationMaker::TrackIsolationMaker(const edm::ParameterSet& iConfig) {

  pfCandidatesTag_		= iConfig.getParameter<InputTag>	("pfCandidatesTag");
  vertexInputTag_               = iConfig.getParameter<InputTag>        ("vertexInputTag");
  
  dR_               = iConfig.getParameter<double>          ("dR_ConeSize");       // dR value used to define the isolation cone                (default 0.3 )
  dzcut_            = iConfig.getParameter<double>          ("dz_CutValue");       // cut value for dz(trk,vtx) for track to include in iso sum (default 0.05)
  minPt_            = iConfig.getParameter<double>          ("minPt_PFCandidate"); // store PFCandidates with pt above this cut                 (default 0   )
  maxIso_           = iConfig.getParameter<double>          ("maxIso_PFCandidate");// store PFCandidates with iso below this cut                (default 0   )
  
  produces<vector<float> >("pfcandstrkiso").setBranchAlias("pfcands_trkiso");
  produces<vector<float> >("pfcandsdzpv"  ).setBranchAlias("pfcands_dzpv");
  produces<vector<float> >("pfcandspt"    ).setBranchAlias("pfcands_pt");
  produces<vector<float> >("pfcandseta"   ).setBranchAlias("pfcands_eta");
  produces<vector<float> >("pfcandsphi"   ).setBranchAlias("pfcands_phi");
  produces<vector<int>   >("pfcandschg"   ).setBranchAlias("pfcands_chg");
    
}

TrackIsolationMaker::~TrackIsolationMaker() 
{

}

void  TrackIsolationMaker::beginRun(edm::Run&, const edm::EventSetup& es) {}
void  TrackIsolationMaker::beginJob() {}
void  TrackIsolationMaker::endJob()   {}

// ------------ method called to produce the data  ------------

void TrackIsolationMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<vector<float> >  pfcands_trkiso(new vector<float>);
  auto_ptr<vector<float> >  pfcands_dzpv  (new vector<float>);
  auto_ptr<vector<float> >  pfcands_pt    (new vector<float>);
  auto_ptr<vector<float> >  pfcands_eta   (new vector<float>);
  auto_ptr<vector<float> >  pfcands_phi   (new vector<float>);
  auto_ptr<vector<int>   >  pfcands_chg   (new vector<int>  );

  //---------------------------------
  // get PFCandidate collection
  //---------------------------------
  
  edm::Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
  pfCandidates  = pfCandidatesHandle.product();

  //---------------------------------
  // get Vertex Collection
  //---------------------------------
  
  Handle<reco::VertexCollection> vertex_h;
  iEvent.getByLabel(vertexInputTag_, vertex_h);
  const reco::VertexCollection *vertices = vertex_h.product();

  //-----------------------------------
  // Find 1st good vertex
  //-----------------------------------

  VertexCollection::const_iterator firstGoodVertex = vertices->end();

  int firstGoodVertexIdx = 0;

  for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  //  cout<<" [TrackIsolationMaker] "<<endl;

  //must have a good vertex to want to store anything
  if ( firstGoodVertex!=vertices->end() ) {

    //-------------------------------------------------------------------------------------------------
    // loop over PFCandidates and calculate the trackIsolation and dz w.r.t. 1st good PV for each one
    //-------------------------------------------------------------------------------------------------
    for( pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++ ) {
      
      //-------------------------------------------------------------------------------------
      // only store PFCandidate values if pt > minPt
      //-------------------------------------------------------------------------------------
      
      if( pf_it->pt() < minPt_ ) continue;
      
      
      //-------------------------------------------------------
      // require PFCandidate is charged
      //-------------------------------------------------------
      if ( pf_it->charge() == 0 ) continue;
	
      //----------------------------------------------------------------------------
      // now loop over other PFCandidates in the event to calculate trackIsolation
      //----------------------------------------------------------------------------

      float trkiso = 0.0;

      for( pat::PackedCandidateCollection::const_iterator pf_other = pfCandidates->begin(); pf_other != pfCandidates->end(); pf_other++ ) {

	// don't count the PFCandidate in its own isolation sum
	if( pf_it == pf_other       ) continue;

	// require the PFCandidate to be charged
	if( pf_other->charge() == 0 ) continue;

	// cut on dR between the PFCandidates
	float dR = deltaR(pf_it->eta(), pf_it->phi(), pf_other->eta(), pf_other->phi());
	if( dR > dR_ ) continue;

	// cut on the PFCandidate dz
	float dz_other = pf_other->dz(firstGoodVertex->position());

	// if ( pf_other->trackRef().isNonnull()) {
	//   dz_other = pf_other->trackRef()->dz( firstGoodVertex->position() );
	// }

	if( fabs(dz_other) > dzcut_ ) continue;

	trkiso += pf_other->pt();
      }

      // calculate the dz of this candidate
      float dz_it = pf_it->dz( firstGoodVertex->position() ); //

      // if ( pf_it->trackRef().isNonnull()) {
      // 	dz_it = pf_it->trackRef()->dz( firstGoodVertex->position() );
      // }

      // key change from Ben's version: want to cut on iso already
      if ( trkiso / pf_it->pt() < maxIso_) {
	//	cout<<"\t"<<pf_it->pt()<<" "<<pf_it->eta()<<" "<<pf_it->phi()<<" "<<pf_it->charge()<<" "<<dz_it<<" "<<trkiso / pf_it->pt()<<endl;
	pfcands_trkiso->push_back(trkiso);
	pfcands_dzpv->push_back(dz_it);
	pfcands_pt->push_back(pf_it->pt());
	pfcands_eta->push_back(pf_it->eta());
	pfcands_phi->push_back(pf_it->phi());
	pfcands_chg->push_back(pf_it->charge());
      }

    } //end of loop of pf cands
          
  } //end of if good vtx

  // put trkiso and dz values back into event
  iEvent.put(pfcands_trkiso,"pfcandstrkiso");
  iEvent.put(pfcands_dzpv  ,"pfcandsdzpv"  );
  iEvent.put(pfcands_pt    ,"pfcandspt"    );
  iEvent.put(pfcands_eta   ,"pfcandseta"   );
  iEvent.put(pfcands_phi   ,"pfcandsphi"   );
  iEvent.put(pfcands_chg   ,"pfcandschg"   );
 
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackIsolationMaker);

