CfANtupler
==========

Code to make cfA ntuples based on CMSSW 7XX and miniAOD.
It's been tested on CMSSW_7_0_6_patch1. 

To run it, issue the following commands on lxplus:

 cmsrel CMSSW_7_0_6_patch1
 cd CMSSW_7_0_6_patch1/src
 cmsenv
 git clone git@github.com:manuelfs/CfANtupler
 scram b -j 4
 cmsRun CfANtupler/ConfigurableAnalysis/python/minicfA_cfg.py

This will create a flat ntuple named configurableAnalysis.root in the
current directory.

In case of not having configured SSH on git, you can also check out the 
CfANtupler package with the http protocol

 git clone http://github.com/manuelfs/CfANtupler