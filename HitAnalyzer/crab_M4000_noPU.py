from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'btagHits_M4000_noPU'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'config.py'

config.Data.inputDataset = '/ZprimeToBBbar_M_4000/thaarres-ZprimeToBBbar_M_4000_TuneCUETP8M1_13TeV_pythia8_GEN-SIM-DIGI-RECO-86331fa3e5403a58492e6b8f4cc1a457/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'btagHits_noPU'

config.Site.storageSite = 'T3_CH_PSI'
config.Data.ignoreLocality = True
config.Site.whitelist = ["T2_CH*","T3_CH*",]