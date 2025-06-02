## Setup the frame work

The support of the package from the **CMSSW_15_1_X** onwards.
```
cmsrel CMSSW_15_1_0_pre1
cd CMSSW_15_1_0_pre1/src
cmsenv
mkdir Phase2Monitoring
cd Phase2Monitoring
git clone https://github.com/vmuralee/HLTTauPhase2Validator.git
scram b -j 32

git checkout <new branch>
```

To produce the flat ntuples make the necessay changes in the config file `test/phase2tauMonitor_cfg.py`.

possible changes,
  - **fileNames** : The list of input files ( mostly RelValTenTau samples)
  - **GlobalTag** : The global tag prescribed in the validation notification
  - **process.TFileService** : to change output file name
  - **triggerfilter** (optional): change the filtename which you want to monitor.

to run the config file, and it produce the validation ntuple
```
cmsRun test/phase2tauMonitor_cfg.py
```

In he `python`folder two scripts `makeValidationPlot.py` and `compareValidationPlots.py` one for make the validation plots and other for compare two validation ntuple

```
python3 makeValidationPlot.py file1.root
python3 compareValidationPlots.py file1.root file2.root

```
