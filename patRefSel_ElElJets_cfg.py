import sys

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "decide if run on MC or data")
options.register('outputFile', 'patRefSel_ElElJets.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, "name of output file")
if( hasattr(sys, "argv") ):
  options.parseArguments()



print sys.argv

process = cms.Process( 'PAT' )


### ======================================================================== ###
###                                                                          ###
###                                 Constants                                ###
###                            (user job steering)                           ###
###                                                                          ###
### ======================================================================== ###

###################
### Data or MC? ###
###################
runOnMC = options.runOnMC
#####################################
### Switch on/off selection steps ###
#####################################
# Step 0a
useTrigger           = True
# Step 0b
useGoodVertex        = True
# Step 0d
use2GoodElectrons    = True
# Step 0e
use2ElectronsVeto    = True
# Step 1
useDiElMassLT20Veto  = True
# Step 2
useDiElMassZVeto     = True
# Step 3
use2Jets             = True
# Step 4
useMET               = True
# Step 5a CSV low
use2BTagJetsL        = True
# Step 5b CSV medium
use2BTagJetsM        = False
# Step 5c CSV tight
use2BTagJetsT        = False

### Trigger matching?
addTriggerMatching = True
###########################
### Reference selection ###
###########################
from TopQuarkAnalysis.Configuration.patRefSel_refElElJets import *

# Electrons
electronVertexMaxDZ = 0.5
#electronCut = 'pt > 20. && abs(eta) < 2.5'

# Jets
#jetCut          = ''
#veryLooseJetCut = 'pt > 20.' # transverse momentum (all jets)
#looseJetCut     = 'pt > 35.' # transverse momentum (3rd jet, optional for 'use3JetsLoose = True')
#tightJetCut     = 'pt > 45.' # transverse momentum (leading jets)

# Trigger and trigger object
#triggerSelectionData       = ''
#triggerObjectSelectionData = ''
#triggerSelectionMC       = ''
#triggerObjectSelectionMC = ''
triggerSelectionMC        = 'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17'
triggerObjectSelectionMC  = 'type("TriggerElectron") && ( path("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17"))'

#####################
### Particle flow ###
#####################
postfix = 'PF'
################################################################
# subtract charged hadronic pile-up particles (from wrong PVs) #
################################################################
# effects also JECs #
#####################
usePFnoPU       = True # before any top projection
usePfIsoLessCHS = True # switch to new PF isolation with L1Fastjet CHS
###############################################################
# other switches for PF top projections (default: all 'True') #
###############################################################
useNoMuon     = True  # before electron top projection
useNoElectron = True  # before jet top projection
useNoJet      = True  # before tau top projection
useNoTau      = False # before MET top projection
################################
# cuts used in top projections #
################################
# vertices #
############
pfVertices  = 'goodOfflinePrimaryVertices'
pfD0Cut     = 0.2
pfDzCut     = 0.5
#############
# electrons #
#############
pfElectronSelectionCut  = 'abs(eta)<2.5 && pt>20. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits < 2'
useElectronCutBasePF  = False # use minimal (veto) electron selection cut on top of 'pfElectronSelectionCut'
electronCutPF = 'pt > 20. && abs(eta) < 2.5 && electronID("mvaTrigV0") > 0.5 && abs(dB()) < 0.04 && passConversionVeto() && gsfTrack().trackerExpectedHitsInner().numberOfHits() <= 0 && (chargedHadronIso+max(0.,neutralHadronIso+photonIso-1.0*userIsolation("User1Iso")))/et < 0.15'
#electronCutPF += ' && (chargedHadronIso+max(0.,neutralHadronIso)+photonIso-0.5*puChargedHadronIso)/et < 0.15' # relative isolation with Delta beta corrections
pfElectronIsoConeR03 = True
pfElectronCombIsoCut = 0.15
##################
### JEC levels ###
##################
# levels to be accessible from the jets
# jets are corrected to L3Absolute (MC), L2L3Residual (data) automatically, if enabled here
# and remain uncorrected, if none of these levels is enabled here
useL1FastJet    = True  # needs useL1Offset being off, error otherwise
useL1Offset     = False # needs useL1FastJet being off, error otherwise; not available from current GT!!!
useL2Relative   = True
useL3Absolute   = True
useL2L3Residual = True  # takes effect only on data
useL5Flavor     = False
useL7Parton     = False

typeIMetCorrections = True

metCut = 'et > 40.'
#############
### Input ###
#############
#######################
# list of input files #
#######################
useRelVals = False # if 'False', "inputFiles" is used
inputFiles = ['file:../../../../../../MC/000C5D15-AB1A-E211-8BDE-00215E22053A.root',
     	      'file:../../../../../../MC/0010005A-421A-E211-9E3C-E41F13181DA4.root'
	 # DY Samples for Syncro
              #'file:/afs/cern.ch/work/t/tkhurshi/Analysis/2013_Analysis/MC/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/70300E2E-27D2-E111-92BD-001E67397AE4.root',
              #'file:/afs/cern.ch/work/t/tkhurshi/Analysis/2013_Analysis/MC/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/7041870B-D3D2-E111-8CFE-001E67397008.root', 
              #'file:/afs/cern.ch/work/t/tkhurshi/Analysis/2013_Analysis/MC/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/70B638CE-C6D2-E111-8430-003048673F0A.root',
              #'file:/afs/cern.ch/work/t/tkhurshi/Analysis/2013_Analysis/MC/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/70CC5B25-C3D2-E111-85A7-001E6739751C.root',
              #'file:/afs/cern.ch/work/t/tkhurshi/Analysis/2013_Analysis/MC/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/70EA5873-5AD2-E111-AE4F-003048D462C4.root'] # overwritten, if "useRelVals" is 'True'
############################
# maximum number of events #
############################
maxEvents = options.maxEvents
##################
### Conditions ###
##################
##############
# GlobalTags #
##############
globalTagData = 'FT53_V21A_AN6::All'
globalTagMC   = 'START53_V27::All'
##############
### Output ###
##############
###############
# output file #
###############
outputFile = options.outputFile
#################################
# event frequency of Fwk report #
#################################
fwkReportEvery = 1000

# switch for 'TrigReport'/'TimeReport' at job end
wantSummary = True
### ======================================================================== ###
###									     ###	
###                              End of constants                            ###
###                                                                          ###
### ======================================================================== ###


###########################
### Basic configuration ###
###########################

process.load( "TopQuarkAnalysis.Configuration.patRefSel_basics_cff" )
process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery
process.options.wantSummary = wantSummary
if runOnMC:
  process.GlobalTag.globaltag = globalTagMC
else:
  process.GlobalTag.globaltag = globalTagData


###########################
### Input configuration ###
###########################

if useRelVals:
  if runOnMC:
    from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValProdTTbarAODSIM
    inputFiles = filesRelValProdTTbarAODSIM
  else:
    from PhysicsTools.PatAlgos.patInputFiles_cff import filesSingleMuRECO
    inputFiles = filesSingleMuRECO
process.load( "TopQuarkAnalysis.Configuration.patRefSel_inputModule_cfi" )
process.source.fileNames = inputFiles
process.maxEvents.input  = maxEvents


############################
### Output configuration ###
############################

process.load( "TopQuarkAnalysis.Configuration.patRefSel_outputModule_cff" )
####################
# output file name #
####################
process.out.fileName = outputFile
#################
# event content #
#################
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out.outputCommands += patEventContent
#########################
# clear event selection #
#########################
process.out.SelectEvents.SelectEvents = []


####################################################
### Cleaning and trigger selection configuration ###
####################################################

#########################
### Trigger selection ###
#########################
if runOnMC:
  triggerSelection = triggerSelectionMC
else:
  triggerSelection = triggerSelectionData
from TopQuarkAnalysis.Configuration.patRefSel_triggerSelection_cff import triggerResults
process.step0a = triggerResults.clone(
  triggerConditions = [ triggerSelection ]
)
#############################
### Good vertex selection ###
#############################
process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )
process.step0b = process.goodOfflinePrimaryVertices.clone( filter = True )
######################
### Event cleaning ###
######################
process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )
process.trackingFailureFilter.VertexSource = cms.InputTag( pfVertices )
process.step0c = process.eventCleaning
if runOnMC:
  process.step0c += process.eventCleaningMC
else:
  process.step0c += process.eventCleaningData


################################
### PAT/PF2PAT configuration ###
################################
process.load( "PhysicsTools.PatAlgos.patSequences_cff" )
##################
### Check JECs ###
##################
###########
# JEC set #
###########
jecSet = 'AK5PF'
if usePFnoPU:
  jecSet += 'chs'
##############
# JEC levels #
##############
if useL1FastJet and useL1Offset:
  sys.exit( 'ERROR: switch off either "L1FastJet" or "L1Offset"' )
jecLevels = []
if useL1FastJet:
  jecLevels.append( 'L1FastJet' )
if useL1Offset:
  jecLevels.append( 'L1Offset' )
if useL2Relative:
  jecLevels.append( 'L2Relative' )
if useL3Absolute:
  jecLevels.append( 'L3Absolute' )
if useL2L3Residual and not runOnMC:
  jecLevels.append( 'L2L3Residual' )
if useL5Flavor:
  jecLevels.append( 'L5Flavor' )
if useL7Parton:
  jecLevels.append( 'L7Parton' )

############################
### Switch configuration ###
############################
from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
usePF2PAT( process
         , runPF2PAT           = True
         , runOnMC             = runOnMC
         , jetAlgo             = jetAlgo
         , postfix             = postfix
         , jetCorrections      = ( jecSet
                                 , jecLevels
                                 )
         , typeIMetCorrections = typeIMetCorrections
         , pvCollection        = cms.InputTag( pfVertices )
         )

if useElectronCutBasePF:
  #process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
  from TopQuarkAnalysis.Configuration.patRefSel_pfIdentifiedElectrons_cfi import pfIdentifiedElectrons
  process.pfIdentifiedElectronsPF     = pfIdentifiedElectrons
  process.pfIdentifiedElectronsPF.src = cms.InputTag( 'pfElectronsFromVertexPF' )
  process.pfSelectedElectronsPF.src   = cms.InputTag( 'pfIdentifiedElectronsPF' )
  process.patPF2PATSequencePF.replace( process.pfSelectedElectronsPF, process.pfIdentifiedElectronsPF + process.pfSelectedElectronsPF )
  pfElectronSelectionCut += ' && %s'%( electronCutPF )

process.pfNoPileUpPF.enable   = usePFnoPU	# True
process.pfNoMuonPF.enable     = useNoMuon	# True
process.pfNoElectronPF.enable = useNoElectron	# True
process.pfNoJetPF.enable      = useNoJet	# True
process.pfNoTauPF.enable      = useNoTau	# False

if useL1FastJet:
  process.pfPileUpIsoPF.checkClosestZVertex = usePfIsoLessCHS

process.pfElectronsFromVertexPF.d0Cut               = pfD0Cut
process.pfElectronsFromVertexPF.dzCut               = pfDzCut
process.pfSelectedElectronsPF.cut                   = pfElectronSelectionCut
process.pfIsolatedElectronsPF.doDeltaBetaCorrection = True # applies EA corrections here!
process.pfIsolatedElectronsPF.deltaBetaFactor       = -1.0
process.pfIsolatedElectronsPF.isolationCut          = pfElectronCombIsoCut

if pfElectronIsoConeR03:
  from EgammaAnalysis.ElectronTools.electronIsolatorFromEffectiveArea_cfi import elPFIsoValueEA03
  process.elPFIsoValueEA03PF = elPFIsoValueEA03
  process.elPFIsoValueEA03PF.pfElectrons = cms.InputTag( 'pfSelectedElectronsPF' )
  process.patPF2PATSequencePF.replace( process.pfSelectedElectronsPF, process.pfSelectedElectronsPF + process.elPFIsoValueEA03PF )
  process.pfIsolatedElectronsPF.isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFIdPF' ) )

  process.pfIsolatedElectronsPF.deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValueEA03PF' ) # EA corrections
  process.pfIsolatedElectronsPF.isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFIdPF' )                                                                                                  , cms.InputTag( 'elPFIsoValueGamma03PFIdPF' ) )
  process.pfElectronsPF.isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFIdPF' ))

  process.pfElectronsPF.deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValueEA03PF' ) # EA corrections
  process.pfElectronsPF.isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFIdPF' )
                                                                  , cms.InputTag( 'elPFIsoValueGamma03PFIdPF' ) )

  process.patElectronsPF.isolationValues.pfNeutralHadrons   = cms.InputTag( 'elPFIsoValueNeutral03PFIdPF' )
  process.patElectronsPF.isolationValues.pfChargedAll       = cms.InputTag( 'elPFIsoValueChargedAll03PFIdPF' )
  process.patElectronsPF.isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFIdPF' )
  process.patElectronsPF.isolationValues.pfPhotons          = cms.InputTag( 'elPFIsoValueGamma03PFIdPF' )
  process.patElectronsPF.isolationValues.pfChargedHadrons   = cms.InputTag( 'elPFIsoValueCharged03PFIdPF' )
  process.patElectronsPF.isolationValues.user               = cms.VInputTag( cms.InputTag( "elPFIsoValueEA03%s"%( postfix ) ) )
  process.patElectronsPF.electronIDSources = cms.PSet( mvaTrigV0 = cms.InputTag("mvaTrigV0")) # just added
else:
  from EgammaAnalysis.ElectronTools.electronIsolatorFromEffectiveArea_cfi import elPFIsoValueEA04
  process.elPFIsoValueEA04PF = elPFIsoValueEA04
  process.elPFIsoValueEA04PF.pfElectrons = cms.InputTag( 'pfSelectedElectrons' + postfix )
  process.patPF2PATSequencePF.replace( process.pfSelectedElectronsPF, process.pfSelectedElectronsPF + process.elPFIsoValueEA04PF )
  process.pfIsolatedElectronsPF.deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValueEA04PF' ) # EA corrections
  process.pfElectronsPF.deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValueEA04PF' ) # EA corrections
  process.patElectronsPF.isolationValues.user = cms.VInputTag( cms.InputTag( "elPFIsoValueEA04%s"%( postfix ) ) )


from PhysicsTools.PatAlgos.tools.coreTools import *

from TopQuarkAnalysis.Configuration.patRefSel_refElElJets_cfi import *
#####################################################
# remove MC matching, object cleaning, objects etc. #
#####################################################
if not runOnMC:
  runOnData( process
           , names = [ 'PFAll' ]
           , postfix = postfix
           )
removeSpecificPATObjects( process
                        , names = [ 'Photons', 'Taus' ]
                        , postfix = postfix
                        ) # includes 'removeCleaning'

# additional event content has to be (re-)added _after_ the call to 'removeCleaning()':
process.out.outputCommands += [ 'keep edmTriggerResults_*_*_*'
                              , 'keep *_hltTriggerSummaryAOD_*_*'
                              # vertices and beam spot
                              , 'keep *_offlineBeamSpot_*_*'
                              , 'keep *_offlinePrimaryVertices*_*_*'
                              , 'keep *_goodOfflinePrimaryVertices*_*_*'
                              ]
if runOnMC:
  process.out.outputCommands += [ 'keep GenEventInfoProduct_*_*_*'
                                , 'keep recoGenParticles_*_*_*'
                                , 'keep *_addPileupInfo_*_*'
                                ]


################################
### Additional configuration ###
################################
#################
### Electrons ###
#################
intermediatePatElectrons.src       = cms.InputTag( 'selectedPatElectronsPF')
process.intermediatePatElectronsPF = intermediatePatElectrons

goodPatElectrons.electronSource = cms.InputTag( 'intermediatePatElectronsPF' )
goodPatElectrons.vertexSource   = cms.InputTag( pfVertices )
process.goodPatElectronsPF      = goodPatElectrons

process.patElectronsPF.usePV             = electronsUsePV
process.patElectronsPF.embedTrack        = electronEmbedTrack
process.patElectronsPF.electronIDSources = electronIDSources

process.selectedPatElectronsPF.cut = '''abs(eta)<2.5 && pt>10. && 
                                     (chargedHadronIso+max(0.,neutralHadronIso+photonIso-1.0*userIsolation("User1Iso")))/et < 0.15 && 
                                     electronID("mvaTrigV0") > 0.00'''

process.intermediatePatElectronsPF.cut = electronCutPF
process.goodPatElectronsPF.maxDZ       = electronVertexMaxDZ

step0d.src = cms.InputTag( 'goodPatElectronsPF' )
process.step0dPF = step0d
# for electron veto
step0e.src = cms.InputTag( 'selectedPatElectronsPF' )
process.step0ePF = step0e

##################### 
## DiElectron Mass ##
#####################                                   
from TopQuarkAnalysis.Configuration.DiElectronFilter_cfi import *
#########################
# Less than 20 GeV Veto #
#########################
process.step1PF              = filterElecPair.clone()
process.step1PF.electrons    = cms.InputTag("goodPatElectronsPF")
process.step1PF.filterCharge = -1
process.step1PF.Cut          = cms.vdouble(0.,20.)
process.step1PF.isVeto       = True
##############
# Zmass Veto #
##############
process.step2PF              = filterElecPair.clone()
process.step2PF.electrons    = cms.InputTag("goodPatElectronsPF")
process.step2PF.filterCharge = -1
process.step2PF.Cut          = cms.vdouble(76.,106.)
process.step2PF.isVeto       = True

############
### Jets ###
############
loosePatJets.src       = cms.InputTag( 'selectedPatJetsPF' )
loosePatJets.cut       = looseJetCut
process.loosePatJetsPF = loosePatJets

tightPatJets.src       = cms.InputTag( 'loosePatJetsPF' )
tightPatJets.cut       = tightJetCut
process.tightPatJetsPF = tightPatJets

process.selectedPatJetsPF.cut  = jetCut
process.loosePatJetsPF.cut     = looseJetCut
process.tightPatJetsPF.cut     = tightJetCut

process.cleanPatJets.checkOverlaps.muons.src = cms.InputTag( 'goodPatMuonsPF' )
process.cleanPatJets.checkOverlaps.muons.deltaR = 0.5
process.cleanPatJets.checkOverlaps.muons.requireNoOverlaps = True

step3.src        = cms.InputTag( 'loosePatJetsPF' )
process.step3PF  = step3

##########
## METs ##
##########
mets.src         = cms.InputTag( 'patMETsPF' )
mets.cut         = metCut
process.metsPF   = mets
step4.src        = cms.InputTag( 'metsPF' )
process.step4PF  = step4

#############
## BTaging ##
#############
bJetsCSVL.src       = cms.InputTag( 'loosePatJetsPF' )
process.bJetsCSVLPF = bJetsCSVL
step5a.src          = cms.InputTag( 'bJetsCSVLPF' )
process.step5aPF    = step5a

bJetsCSVM.src       = cms.InputTag( 'loosePatJetsPF' )
process.bJetsCSVMPF = bJetsCSVM
step5b.src          = cms.InputTag( 'bJetsCSVMPF' )
process.step5bPF    = step5b

bJetsCSVT.src       = cms.InputTag( 'loosePatJetsPF' )
process.bJetsCSVTPF = bJetsCSVT
step5c.src          = cms.InputTag( 'bJetsCSVTPF' )
process.step5cPF    = step5c

process.out.outputCommands.append( 'keep *_goodPatElectrons*_*_*' )
process.out.outputCommands.append( 'keep *_veryLoosePatJets*_*_*' )
process.out.outputCommands.append( 'keep *_loosePatJets*_*_*' )
process.out.outputCommands.append( 'keep *_tightPatJets*_*_*' )

########################
### Trigger matching ###
########################

if addTriggerMatching:

  if runOnMC:
    triggerObjectSelection = triggerObjectSelectionMC
  else:
    triggerObjectSelection = triggerObjectSelectionData
  ### Trigger matching configuration
  from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger
  from TopQuarkAnalysis.Configuration.patRefSel_triggerMatching_cfi import patElectronTriggerMatch
  from PhysicsTools.PatAlgos.tools.trigTools import *
  triggerProducerPF = patTrigger.clone()
  process.patTriggerPF = triggerProducerPF 
  triggerMatchPF = patElectronTriggerMatch.clone( matchedCuts = triggerObjectSelection )
  process.triggerMatchPF = triggerMatchPF 
  switchOnTriggerMatchEmbedding( process
                               , triggerProducer = 'patTrigger' + postfix
                               , triggerMatchers = [ 'triggerMatch' + postfix ]
                               , sequence        = 'patPF2PATSequence' + postfix
                               , postfix         = postfix
                               )
  removeCleaningFromTriggerMatching( process
                                   , sequence = 'patPF2PATSequencePF'
                                   )
  process.intermediatePatElectronsPF.src = cms.InputTag( 'selectedPatElectronsPFTriggerMatch' )


##################
### Scheduling ###
##################

###################
# MVA electron ID #
###################
process.load( "EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi" )
process.eidMVASequence = cms.Sequence(
  process.mvaTrigV0
+ process.mvaNonTrigV0
)
###########################
# The additional sequence #
###########################

patAddOnSequence = cms.Sequence(
  process.intermediatePatElectronsPF
* process.goodPatElectronsPF
* process.loosePatJetsPF
* process.tightPatJetsPF
* process.metsPF
* process.bJetsCSVLPF
* process.bJetsCSVMPF
* process.bJetsCSVTPF
)
process.patAddOnSequencePF = patAddOnSequence 

#############
# The paths #
#############
process.p = cms.Path()
if useTrigger:
  process.p += process.step0a
process.p += process.goodOfflinePrimaryVertices
if useGoodVertex:
  process.p += process.step0b
process.p += process.step0c
process.p += process.eidMVASequence
process.p += process.patPF2PATSequencePF
process.p += process.patAddOnSequencePF
if use2GoodElectrons:
  process.p += process.step0dPF
if use2ElectronsVeto:
  process.p += process.step0ePF
if useDiElMassLT20Veto:
  process.p += process.step1PF
if useDiElMassZVeto:
  process.p += process.step2PF
if use2Jets:
  process.p += process.step3PF
if useMET:
  process.p += process.step4PF
if use2BTagJetsL:
  process.p += process.step5aPF
if use2BTagJetsM:
  process.p += process.step5bPF
if use2BTagJetsT:
  process.p += process.step5cPF

process.out.SelectEvents.SelectEvents.append( 'p' )
