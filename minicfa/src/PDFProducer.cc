//-----------------------------------------------------------------------------------------
//
// Get PDF info branches--adapted from MT2 (https://github.com/cmstas/NtupleMaker/blob/master/src/PDFInfoMaker.cc)
//
//-----------------------------------------------------------------------------------------
// system include files

#include <memory>
#include <vector>
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

#include "CfANtupler/minicfa/interface/PDFProducer.h"

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
PDFProducer::PDFProducer(const edm::ParameterSet& iConfig) {

  genEventInfoInputTag_ = iConfig.getParameter<std::string>("genEventInfoInputTag");
  hepmcHandle_ = iConfig.getParameter<std::string>("hepmcHandle");
  
  produces<float> ("pdfinfox1" ).setBranchAlias("pdf_info_x1" );
  produces<float> ("pdfinfox2" ).setBranchAlias("pdf_info_x2" );
  produces<float> ("pdfinfoscale" ).setBranchAlias("pdf_info_scale");
  produces<float> ("pdfinfopdf1" ).setBranchAlias("pdf_info_pdf1" );
  produces<float> ("pdfinfopdf2" ).setBranchAlias("pdf_info_pdf2" );
  produces<int> ("pdfinfoid1" ).setBranchAlias("pdf_info_id1" );
  produces<int> ("pdfinfoid2" ).setBranchAlias("pdf_info_id2" );
  
}
PDFProducer::~PDFProducer()
{
  // if (clusterTools_) delete clusterTools_;
}
void PDFProducer::beginRun(edm::Run&, const edm::EventSetup& es) {}
void PDFProducer::beginJob() {}
void PDFProducer::endJob() {}
// ------------ method called to produce the data ------------
void PDFProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<float> pdfinfo_x1 ( new float );
  auto_ptr<float> pdfinfo_x2 ( new float );
  auto_ptr<float> pdfinfo_scale( new float );
  auto_ptr<float> pdfinfo_pdf1 ( new float );
  auto_ptr<float> pdfinfo_pdf2 ( new float );
  auto_ptr<int> pdfinfo_id1 ( new int );
  auto_ptr<int> pdfinfo_id2 ( new int );

  edm::Handle<GenEventInfoProduct> hEvtInfo;
  iEvent.getByLabel(genEventInfoInputTag_, hEvtInfo);
  // get MC particle collection
  edm::Handle<edm::HepMCProduct> hepmcHandle;
  iEvent.getByLabel( hepmcHandle_, hepmcHandle );
  const HepMC::GenEvent* evt = 0;
  const HepMC::PdfInfo* pdfinfo = 0;
  if(!hepmcHandle.failedToGet() ) {
    evt = hepmcHandle->GetEvent();
    pdfinfo = evt->pdf_info();
  }

  //try to get using the GenEventInfoProduct
  if(!hEvtInfo.failedToGet() && hEvtInfo->hasPDF()) {
    const gen::PdfInfo *pdf = hEvtInfo->pdf();
    *pdfinfo_x1 = pdf->x.first;
    *pdfinfo_x2 = pdf->x.second;
    *pdfinfo_scale = pdf->scalePDF;
    *pdfinfo_pdf1 = pdf->xPDF.first;
    *pdfinfo_pdf2 = pdf->xPDF.second;
    *pdfinfo_id1 = pdf->id.first;
    *pdfinfo_id2 = pdf->id.second;
  } else if(pdfinfo != 0) {
    //assign
    *pdfinfo_x1 = pdfinfo->x1();
    *pdfinfo_x2 = pdfinfo->x2();
    *pdfinfo_scale = pdfinfo->scalePDF();
    *pdfinfo_pdf1 = pdfinfo->pdf1();
    *pdfinfo_pdf2 = pdfinfo->pdf2();
    *pdfinfo_id1 = pdfinfo->id1();
    *pdfinfo_id2 = pdfinfo->id2();
  } else {
    *pdfinfo_x1 = -9999;
    *pdfinfo_x2 = -9999;
    *pdfinfo_scale = -9999;
    *pdfinfo_pdf1 = -9999;
    *pdfinfo_pdf2 = -9999;
    *pdfinfo_id1 = -9999;
    *pdfinfo_id2 = -9999;
    return;
   }

  // put everything back into event
  iEvent.put(pdfinfo_x1,"pdfinfox1");
  iEvent.put(pdfinfo_x2,"pdfinfox2");
  iEvent.put(pdfinfo_scale,"pdfinfoscale");
  iEvent.put(pdfinfo_pdf1,"pdfinfopdf1");
  iEvent.put(pdfinfo_pdf2,"pdfinfopdf2");
  iEvent.put(pdfinfo_id1,"pdfinfoid1");
  iEvent.put(pdfinfo_id2,"pdfinfoid2");

}

  //define this as a plug-in
  DEFINE_FWK_MODULE(PDFProducer);
