CfANtupler
==========

Code to make cfA ntuples based on CMSSW 7XX and miniAOD.
It's been tested on `CMSSW_7_2_2_patch1` and in 7_4_X releases. 
CMSSW_7_2_2_patch1 should be used to process PHYS14 MC
and a 7_4_X release should be used to process 13 TeV data
and 7_4_X MC.

#### Running the code
Change the variable `datasetType` in `minicfA_cfg.py` to the appropriate value
for the dataset that you are processing.

Issue the following commands on lxplus:

    cmsrel CMSSW_7_2_2_patch1
    cd CMSSW_7_2_2_patch1/src
    cmsenv
    git clone git@github.com:manuelfs/CfANtupler
    scram b -j 4
    cmsRun CfANtupler/minicfa/python/minicfA_cfg.py

This will create a flat ntuple named configurableAnalysis.root in the
current directory.

In case of not having configured SSH in git, you can also check out the 
CfANtupler package with the http protocol

    git clone http://github.com/manuelfs/CfANtupler

#### Adding/changing tree content
Most of the branches are defined in `CfANtupler/minicfa/python/branchesminicfA_cfi.py`. 
You can easily modify the leaves parameter to add a branch with the format
`'name:member'`, where `name` will be the name of the branch, and `member` must
be a member function of the collection you are using.

Branches that require C++ code (e.g. triggers) are defined in 
`CfANtupler/minicfa/interface/AdHocNTupler.h`.
