## Configuration file to CRAB3 cfA jobs
## Submit for src with: crab submit -c crabcfA.py


dataset = '/SMS-T2tt_2J_mStop-850_mLSP-100_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/MINIAODSIM'

jobname = dataset[1:].replace('/','__')
jobname = jobname.replace(':','___')


from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = jobname
config.General.workArea = 'out_crab'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'CfANtupler/minicfa/python/minicfA_cfg.py'

config.section_("Data")
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = True
config.Data.publishDBS = 'phys03'

config.section_("Site")
config.Site.storageSite = 'T2_US_UCSD'
