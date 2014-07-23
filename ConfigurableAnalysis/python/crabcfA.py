from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'T1bbbb'
config.General.workArea = 'out_crab'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'CfANtupler/ConfigurableAnalysis/python/minicfA_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/SMS-T1bbbb_2J_mGl-1000_mLSP-900_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
config.Data.dbsUrl = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.publication = True
config.Data.publishDbsUrl = 'phys03'
config.Data.publishDataName = 'T1bbbb_cfA'

config.section_("Site")
config.Site.storageSite = 'T2_US_UCSD'
