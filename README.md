# HitAnalyzer

ED analyzer to produce reco jet collection, pixel clusters and gen particles for b-tagging with nHits.

##Setup:
First fork this repository, then do
```
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src
cmsenv
export GITUSER=`git config user.github`
echo "Your github username has been set to \"$GITUSER\""
git clone git@github.com:${GITUSER}/HitAnalyzer.git
cd HitAnalyzer
git remote add MAIN git@github.com:thaarres/HitAnalyzer.git
scram b -j 12
```
To run, change input file in config.py, then
```
cmsRun config.py
```

The main code block to edit is in plugins/HitAnalyzer.cc.

To submit with crab
```
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit -c crab_config.py
```
