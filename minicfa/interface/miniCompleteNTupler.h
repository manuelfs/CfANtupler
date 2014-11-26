#include "PhysicsTools/UtilAlgos/interface/NTupler.h"

#include "CfANtupler/minicfa/interface/miniStringBasedNTupler.h"
#include "CfANtupler/minicfa/interface/miniVariableNTupler.h"
#include "CfANtupler/minicfa/interface/miniAdHocNTupler.h"

class miniCompleteNTupler : public NTupler {
 public:
  miniCompleteNTupler(const edm::ParameterSet& iConfig){
    sN = new miniStringBasedNTupler(iConfig);
    if (iConfig.exists("variablesPSet"))
      if (!iConfig.getParameter<edm::ParameterSet>("variablesPSet").empty())
	vN = new miniVariableNTupler(iConfig);
      else vN=0;
    else
      vN=0;
    if (iConfig.exists("AdHocNPSet"))
      if (!iConfig.getParameter<edm::ParameterSet>("AdHocNPSet").empty())
	aN = new miniAdHocNTupler(iConfig);
      else aN=0;
    else
      aN=0;
  }
  
  uint registerleaves(edm::ProducerBase * producer){
    uint nLeaves=0;
    nLeaves+=sN->registerleaves(producer);
    if (vN)
      nLeaves+=vN->registerleaves(producer);
    if (aN)
      nLeaves+=aN->registerleaves(producer);
    return nLeaves;
  }
  void fill(edm::Event& iEvent){
    sN->fill(iEvent);
    if (vN)
      vN->fill(iEvent);
    if (aN)
      aN->fill(iEvent);

    sN->callBack();
    if (vN)
      vN->callBack();
    if (aN)
      aN->callBack();
  }

 private:
  miniStringBasedNTupler * sN;
  miniVariableNTupler * vN;  
  miniAdHocNTupler * aN;

};

