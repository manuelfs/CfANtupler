## Configuration file to CRAB3 cfA jobs
## Submit for src with: crab submit -c crabcfA.py


dataset = ''

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
config.JobType.inputFiles = ['CfANtupler/JEC/Summer15_25nsV2_MC.db','CfANtupler/JEC/Summer15_50nsV2_MC.db', 'CfANtupler/JEC/Summer15_50nsV4_DATA.db']

config.section_("Data")
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = True
config.Data.publishDBS = 'phys03'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'

config.section_("Site")
config.Site.storageSite = 'T3_US_UCSB'
#config.Site.whitelist = ['T2_US_Caltech','T2_US_Florida', 'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD', 'T2_US_Vanderbilt', 'T2_US_Wisconsin', 'T1_US_FNAL','T2_US_MIT', 'T3_US_UCSB']
config.Site.blacklist = ['T1_RU_JINR']
# you may want to uncomment this line and force jobs to run in the US
# only a few datasets (mostly very new ones) will not be accessible
