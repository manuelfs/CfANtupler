CfANtupler
==========

Code to make cfA ntuples based on CMSSW 7XX and miniAOD.
It's been tested on `CMSSW_7_0_6_patch1`. 

#### Running the code
Issue the following commands on lxplus:

    cmsrel CMSSW_7_0_6_patch1
    cd CMSSW_7_0_6_patch1/src
    cmsenv
    git clone git@github.com:manuelfs/CfANtupler
    scram b -j 4
    cmsRun CfANtupler/ConfigurableAnalysis/python/minicfA_cfg.py

This will create a flat ntuple named configurableAnalysis.root in the
current directory.

In case of not having configured SSH in git, you can also check out the 
CfANtupler package with the http protocol

    git clone http://github.com/manuelfs/CfANtupler

#### Adding/changing tree content
Most of the branches are defined in `CfANtupler/ConfigurableAnalysis/python/branchesminicfA_cfi.py`. 
You can easily modify the leaves parameter to add a branch with the format
`'name:member'`, where `name` will be the name of the branch, and `member` must
be a member function of the collection you are using.

Branches that require C++ code (e.g. triggers) are defined in 
`CfANtupler/ConfigurableAnalysis/interface/AdHocNTupler.h`.
