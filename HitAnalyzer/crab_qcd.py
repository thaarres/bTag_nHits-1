from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'btagHits_QCDwPU2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.maxMemoryMB=6000
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'config.py'

config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISummer17DRPremix-NZS_92X_upgrade2017_realistic_v10-v3/GEN-SIM-RECO'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = config.Data.unitsPerJob * 500
config.Data.outLFNDirBase = '/store/user/thaarres/'
config.Data.publication = False
config.Data.outputDatasetTag = 'QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_wPU-v2'

config.Site.storageSite = 'T3_CH_PSI'
# config.Data.ignoreLocality = True
# config.Site.whitelist = ["T2_CH*","T3_CH*",]