#!/bin/bash
cd /afs/cern.ch/work/m/meshoaib/PatAnalysis/CMSSW_5_3_21/src/MyAnalysis/TbZ/
eval `scramv1 runtime -sh`
cmsRun tbzanalysisSignalMC2_cfg.py 
