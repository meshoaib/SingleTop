import FWCore.ParameterSet.Config as cms

#Process name:
process = cms.Process("SingleTop")

#MessageLogger options:
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    FailPath = cms.untracked.vstring('ProductNotFound','Type Mismatch')
    )

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

ChannelName = "TChannel"

#Data or MC:
isData=False

#Input file:
process.source = cms.Source (
    "PoolSource",
    fileNames = cms.untracked.vstring (
      #"file:/afs/cern.ch/work/o/oiorio/public/xFrancescoFab/T_t-channel_Synch.root"
 "file:/afs/cern.ch/work/m/meshoaib/NewAnalysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/A27919CE-486A-E311-8EDD-002590D0B0C8.root"
#      "file:/afs/cern.ch/work/o/oiorio/public/xFrancescoFab/DataReRecoA.root"
#      "file:8425E88A-AED7-E211-8067-002481E14F5C.root"
    ),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

if isData: process.source.fileNames= cms.untracked.vstring ( "file:/afs/cern.ch/work/o/oiorio/public/xFrancescoFab/DataReRecoA.root" )

#process.MessageLogger.cerr.FwkReport.reportEvery = 100

#Include the hermetic top projection: 
hermeticTopProjection=True

#Gsf electron or PF electron: 
doGsfElectrons=False

#Add nJ >= 2 cut: 
addJetsCut=True 

#Run JetMET uncertainties from JME tool
doRunMETUncertainties=True

#Enable PF taus (true = default: adds the PF tau objects;false = embed them in jets)
EnablePFTaus=False

#Geometry:
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff") ### real data

#Tag:
### new GT - 10Sept2013
if isData:  process.GlobalTag.globaltag = cms.string('FT53_V21A_AN6::All')
else: process.GlobalTag.globaltag = cms.string('START53_V27::All')
#process.GlobalTag.globaltag = cms.string('FT_53_V6C_AN3::All')


###Jet probability calibration used for b-tagging must be rerun in simulation:
### https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG?rev=169#2012_Data_and_MC_EPS13_prescript
###Relevant code snippet is copied from this reference. Note that the calibration in 22Jan2013 rereco of real data is fine.

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
             tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
             connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
             connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
    )




# dummy output: needed to avoid crash
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('dummy.root'),
    outputCommands = cms.untracked.vstring(""),
)

# ---> PAT + Single top sequences <---
process.load("PhysicsTools.PatAlgos.patSequences_cff") 

### Filter on good vertices commented for tH analysis
process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter" ,
                                                   filterParams = cms.PSet( minNdof = cms.double( 4. ) , maxZ = cms.double( 24. ) , maxRho = cms.double( 2. ) ) ,
#filter = cms.bool( True) , src = cms.InputTag( 'offlinePrimaryVertices' ) )
filter = cms.bool( False) , src = cms.InputTag( 'offlinePrimaryVertices' ) )



# Configure PAT to use PFBRECO instead of AOD sources
# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *
from PhysicsTools.PatUtils.tools.metUncertaintyTools import *

postfix = ""
runOnMC = not(isData)
jetAlgoName = "AK5"

if runOnMC:#MC
    jetCorrections=['L1FastJet','L2Relative','L3Absolute']
else:#Data
    jetCorrections=['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

############ PRINTOUT ###################
sep_line = "-" * 50
print sep_line
print 'running the following PFBRECO sequence:'
print jetAlgoName
print 'run on MC        : ', runOnMC
print sep_line
print 'postfix       : ', postfix
print sep_line
print 'JEC        : ', jetCorrections
print sep_line
#########################################

usePF2PAT(process, runPF2PAT=True, jetAlgo=jetAlgoName, runOnMC=runOnMC, postfix=postfix,
          jetCorrections=('AK5PFchs',jetCorrections), pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
          # jetCorrections=('AK5PFchs',jetCorrections), pvCollection=cms.InputTag('offlinePrimaryVertices'),
          typeIMetCorrections=(not doRunMETUncertainties))

if doRunMETUncertainties:
    from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties
    if isData:
        runMEtUncertainties(process,electronCollection = "selectedPatElectrons", doSmearJets= False, muonCollection = "selectedPatMuons", tauCollection="selectedPatTaus", jetCollection = "selectedPatJets",jetCorrLabel="L2L3Residual")
        process.patPFMet.addGenMET = False
        process.patPFMetJetEnUp.addGenMET = False
        process.patPFMetJetEnDown.addGenMET = False
        process.patPFMetElectronEnUp.addGenMET = False
        process.patPFMetElectronEnDown.addGenMET = False
        process.patPFMetMuonEnUp.addGenMET = False
        process.patPFMetMuonEnDown.addGenMET = False
        process.patPFMetTauEnUp.addGenMET = False
        process.patPFMetTauEnDown.addGenMET = False
        process.patPFMetTauEnUp.addGenMET = False
        process.patPFMetTauEnDown.addGenMET = False
    else:
        runMEtUncertainties(process,electronCollection = "selectedPatElectrons", doSmearJets= False, muonCollection = "selectedPatMuons", tauCollection="selectedPatTaus", jetCollection = "selectedPatJets",)

# CONFIGURE LEPTONS for the analysis


# Use DR = 0.3 for PF electrons:
process.pfIsolatedElectrons.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFId"))
#process.pfIsolatedElectrons.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFId")
process.pfIsolatedElectrons.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"), cms.InputTag("elPFIsoValueGamma03PFId"))


    
# Switch to GsfElectrons if required
if doGsfElectrons:
    # Gsf Electrons
    # Use DR = 0.3:
    useGsfElectrons(process,postfix,"03")
if not(doGsfElectrons):
    # Set Isolation for PAT Electrons
    process.patElectrons.isolationValues = cms.PSet( pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFId"), pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFId"), pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFId"), pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFId"), pfPhotons = cms.InputTag("elPFIsoValueGamma03PFId") )

#Apply PFnoPU:
process.pfPileUp.Enable = True
process.pfPileUp.checkClosestZVertex = cms.bool(False)

#Disable top projection for taus
process.pfNoTau.enable = EnablePFTaus

# Prepare MVA electronId
process.load("EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi")
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectronIDs = cms.Sequence( process.mvaID )
# Add MVA electronId
process.electronIDSources = cms.PSet(
    mvaTrigV0    = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0    = cms.InputTag("mvaNonTrigV0")
    )

process.patElectrons.electronIDSources = process.electronIDSources

### Removing request of compatibility with the primary vertex
### Because info on the impact parameter are used in the MVA ID
### https://hypernews.cern.ch/HyperNews/CMS/get/egamma-elecid/72/1.html
process.pfElectronsFromVertex.d0Cut = 9999.
process.pfElectronsFromVertex.d0SigCut = 9999.
process.pfElectronsFromVertex.dzCut = 9999.
process.pfElectronsFromVertex.dzSigCut = 9999.

process.pfMuonsFromVertex.d0Cut = 9999.
process.pfMuonsFromVertex.d0SigCut = 9999.
process.pfMuonsFromVertex.dzCut = 9999.
process.pfMuonsFromVertex.dzSigCut = 9999.



### ========================
### HERMETIC TOP PROJECTION
### ========================

if(hermeticTopProjection):

    ### ELECTRONS
    process.pfSelectedElectrons.cut = 'abs(eta)<2.5 && gsfElectronRef.pt>20. && gsfTrackRef.isNonnull '
    
    process.load("EGamma.EGammaAnalysisTools.electronIsolatorFromEffectiveArea_cfi")

    process.pfIsolatedElectrons.deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValueEA03')
    process.pfIsolatedElectrons.doDeltaBetaCorrection = True
    process.pfIsolatedElectrons.deltaBetaFactor = -1.0
    process.pfIsolatedElectrons.isolationCut = 0.15

    process.patPF2PATSequence.replace( process.pfSelectedElectrons, process.patElectronIDs  + process.pfSelectedElectrons + process.elPFIsoValueEA03)
    ###  same definition as for the veto electrons
    process.selectedPatElectrons.cut = ''' isPF && ecalDrivenMomentum.pt > 20 && abs(eta) < 2.5 && userFloat(\"RhoCorrectedIso\") <0.15'''



    ### MUONS
    process.pfSelectedMuons.cut = 'abs(eta)<2.5 && pt>10.'
#    process.pfSelectedMuons.cut = 'abs(eta)<2.5 && pt>10. && ( isGlobalMuon || isTrackerMuon )'
    process.pfIsolatedMuons.doDeltaBetaCorrection = True
    process.pfIsolatedMuons.deltaBetaFactor = -0.5
    process.pfIsolatedMuons.isolationCut = 0.2
    
    process.patMuons.embedTrack = True
    #process.patMuons.usePV = False
    ###  same definition as for the veto muons
    process.selectedPatMuons.cut = '''abs(eta)<2.5 && pt>10. &&
    (chargedHadronIso+max(0.,neutralHadronIso+photonIso-0.50*puChargedHadronIso))/pt < 0.20 &&
    (isPFMuon && (isGlobalMuon || isTrackerMuon) )'''




# Define the PAT sequence:
process.patseq = cms.Sequence(
    process.goodOfflinePrimaryVertices *
    process.patElectronIDs *
    getattr(process,"patPF2PATSequence"+postfix)
    )

# Add PUJetID
process.load("CMGTools.External.pujetidsequence_cff")

# Define new PAT muons/electrons with no isolation in pf reco (ZeroIso suffix)
# They will be used to get anti-isolated leptons:


# MuonsZeroIso
process.pfIsolatedMuonsZeroIso = process.pfIsolatedMuons.clone(combinedIsolationCut =  cms.double(float("inf")),
                                                               isolationCut =  cms.double(float("inf"))
                                                               )
from TopQuarkAnalysis.SingleTop.AdaptPFMuonsFix_cff import adaptPFMuonsAnd
process.patMuonsZeroIso = process.patMuons.clone(pfMuonSource = cms.InputTag("pfIsolatedMuonsZeroIso"))
# use pf isolation, but do not change matching:
tmp = process.muonMatch.src

adaptPFMuonsAnd(process, process.patMuonsZeroIso, "")
process.muonMatch.src = tmp
process.muonMatchZeroIso = process.muonMatch.clone(src = cms.InputTag("pfIsolatedMuonsZeroIso"))
process.patMuonsZeroIso.genParticleMatch = cms.InputTag("muonMatchZeroIso")
process.patMuonsZeroIso.pfMuonSource = cms.InputTag("pfIsolatedMuonsZeroIso")


# ElectronsZeroIso
process.pfIsolatedElectronsZeroIso = process.pfIsolatedElectrons.clone(combinedIsolationCut = cms.double(float("inf")),
                                                                       isolationCut =  cms.double(float("inf")),
                                                                       )
process.patElectronsZeroIso = process.patElectrons.clone(pfElectronSource = cms.InputTag("pfIsolatedElectronsZeroIso"))

#Define ZeroIso leptons sequence:
if isData:
    process.ZeroIsoLeptonSequence = cms.Sequence(
        process.pfIsolatedMuonsZeroIso +
        process.patMuonsZeroIso +
        process.pfIsolatedElectronsZeroIso +
        process.patElectronsZeroIso
        )
else:    
    process.ZeroIsoLeptonSequence = cms.Sequence(
        process.pfIsolatedMuonsZeroIso +
        process.muonMatchZeroIso +
        process.patMuonsZeroIso +
        process.pfIsolatedElectronsZeroIso +
        process.patElectronsZeroIso
        )
    
##### Define leptons collections useful in single top analysis

# Veto leptons
process.load("TopQuarkAnalysis.SingleTop.userDataLeptonProducers_cfi") 

process.vetoMuons = process.userDataMuons.clone(
    cut = cms.string(" (isGlobalMuon || isTrackerMuon) " +
                     "& pt > 10 & abs(eta) < 2.5 " +
                     "& userFloat(\"DeltaCorrectedIso\") <0.2 ")
)

process.vetoElectrons = process.userDataElectrons.clone(
    cut = cms.string("ecalDrivenMomentum.pt > 20 " +
                     "& abs(eta) < 2.5 " +
                     "& userFloat(\"RhoCorrectedIso\") <0.15")# +
                    # "& userFloat(\"PassesVetoID\") >0.0")
)

process.vetoElectronsMVA = process.userDataElectrons.clone(
    cut =  cms.string(" ecalDrivenMomentum.pt > 20" +
                      "& abs(eta) < 2.5 && userFloat(\"RhoCorrectedIso\") <0.15")# +
                    #  "& electronID('mvaTrigV0') >0.0")
)

# Tight leptons
process.tightMuons = process.userDataMuons.clone(
    cut = cms.string(" pt > 26 & isGlobalMuon && isPFMuon & abs(eta) < 2.1 && normChi2 < 10 && track.hitPattern.trackerLayersWithMeasurement>5 "+
                     "& numberOfMatchedStations() > 1 && innerTrack.hitPattern.numberOfValidPixelHits > 0 " +
                     "& globalTrack.hitPattern.numberOfValidMuonHits > 0" +
                     #                     "& userFloat('VertexDxy')<0.02" +
                     "& abs(dB) < 0.2"
                     "& userFloat('VertexDz')<0.5" +
                     "& userFloat(\"DeltaCorrectedIso\") <0.12 " )
)

process.tightElectrons = process.userDataElectrons.clone(
    cut =  cms.string(" ecalDrivenMomentum.pt > 30  && abs(eta)<2.5" +
                      "& ( abs(superCluster.eta)> 1.5660 || abs(superCluster.eta)<1.4442)" +
                      "& gsfTrack.trackerExpectedHitsInner.numberOfHits <=0" +
                      "& passConversionVeto" +
                      #"& userFloat('VertexDxy')<0.02" +
                      "& userFloat('RhoCorrectedIso')<0.1" )
)

# Tight leptons ZeroIso
process.tightMuonsZeroIso = process.userDataMuons.clone(
    src = cms.InputTag("patMuonsZeroIso"),
    cut = cms.string(" pt > 26 & isGlobalMuon && isPFMuon & abs(eta) < 2.1 && normChi2 < 10 && track.hitPattern.trackerLayersWithMeasurement>5 "+
                     "& numberOfMatchedStations() > 1 && innerTrack.hitPattern.numberOfValidPixelHits > 0 " +
                     "& globalTrack.hitPattern.numberOfValidMuonHits > 0")
)

process.tightElectronsZeroIso = process.userDataElectrons.clone(
    src = cms.InputTag("patElectronsZeroIso"),
    cut =  cms.string(" ecalDrivenMomentum.pt > 30  && abs(eta)<2.5" +
                      "& ( abs(superCluster.eta)> 1.5660 || abs(superCluster.eta)<1.4442)" +
                      "& passConversionVeto")
)

##### Filtering on leptons numbers
process.load("TopQuarkAnalysis.SingleTop.leptonCounterFilter_cfi") 
# Select events with at least 1 tight lepton OR at least one tight leptonNoIso
process.countLeptons.minNumberLoose = 0
process.countLeptons.maxNumberLoose = 99
process.countLeptons.minNumberTight = 0
process.countLeptons.maxNumberTight = 99
process.countLeptons.minNumberQCD = 0
process.countLeptons.maxNumberQCD = 99

# define Jets for single top analysis
process.load("TopQuarkAnalysis.SingleTop.userDataJetsProducer_cfi") 

process.load("TopQuarkAnalysis.SingleTop.userDataMETsProducer_cfi") 

#definition: Jets Loose
process.topJetsPF.cut = cms.string("numberOfDaughters()>1 & pt()> 20 && abs(eta())<4.7 " +
                                   " & ((abs(eta())>=2.4) || ( chargedHadronEnergyFraction() > 0 & chargedMultiplicity()>0 " +
                                   " & chargedEmEnergyFraction()<0.99))" +
                                   " & (neutralHadronEnergy + HFHadronEnergy) / energy < 0.99" )
                                   #neutralEmEnergyFraction() < 0.99 & neutralHadronEnergyFraction() < 0.99 ")

process.UnclusteredMETPF = cms.EDProducer("SingleTopUnclusteredMETProducer",
                                  metSource = cms.InputTag("patType1CorrectedPFMet"),
                                  jetsSource = cms.InputTag("selectedPatJets"),
                                  electronsSource = cms.InputTag("selectedPatElectrons"),
                                  muonsSource = cms.InputTag("selectedPatMuons"),
                                  )


process.topJetsPF.isData = isData
process.topMETsPF.isData = isData
process.topMETsPF.addExternalUnclusteredMET = doRunMETUncertainties
if (not doRunMETUncertainties): 
    process.topMETsPF.metsSrc = cms.InputTag("patMETs")

if isData:
    process.basePath = cms.Sequence(
        process.vetoMuons +
        process.vetoElectrons +
        process.vetoElectronsMVA +
        process.topJetsPF +
        process.topMETsPF +
        process.tightMuonsZeroIso +
        process.tightElectronsZeroIso +
        process.tightMuons +
        process.tightElectrons
    )

else:    
    process.basePath = cms.Sequence(
        process.vetoMuons +
        process.vetoElectrons +
        process.vetoElectronsMVA +
        process.topJetsPF +
        process.topMETsPF +
        process.UnclusteredMETPF +
        process.tightMuonsZeroIso +
        process.tightElectronsZeroIso +
        process.tightMuons +
        process.tightElectrons
    )

#Trigger filter to be eventually used:
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as triggerFilter

process.HLTFilterMu2012  = triggerFilter.triggerResultsFilter.clone(
    hltResults = cms.InputTag( "TriggerResults","","HLT" ),
    triggerConditions = ["HLT_*"],#All trigger paths are included in the skim
#    triggerConditions = ["HLT_IsoMu24_eta2p1_v*"],
#    triggerConditions = ["HLT_Ele27_WP80_v*"],
    l1tResults = '',
    throw = False
    )

from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *
process.jetsCut = countPatJets.clone(src = 'topJetsPF', minNumber = 2)

# Overall skim path
process.singleTopSkimPath = cms.Path(
    process.HLTFilterMu2012 *
    process.patseq +
    process.puJetIdSqeuence +
    process.puJetIdSqeuenceChs *
    process.ZeroIsoLeptonSequence *
    process.basePath #+
# moved preselection in a standalone path
#    process.preselection# + 
#    process.nTuplesSkim 
    )

# Load recommended event filters
process.load("TopQuarkAnalysis.SingleTop.SingleTopEventFilters_cff") 

# Define event filtering path
#process.preselection = cms.Sequence(
process.preselection = cms.Path(
    process.HLTFilterMu2012 *
    process.HBHENoiseFilter *
    process.scrapingVeto *
    process.CSCTightHaloFilter *
    process.hcalLaserEventFilter *
    process.EcalDeadCellTriggerPrimitiveFilter *
    process.EcalDeadCellBoundaryEnergyFilter *
    process.goodVertices *
    process.trackingFailureFilter *
    process.eeBadScFilter *
    process.ecalLaserCorrFilter *
    ~process.manystripclus53X *
    ~process.toomanystripclus53X *
    ~process.logErrorTooManyClusters *
    process.countLeptons
    )

if addJetsCut:
    process.preselection += process.jetsCut

process.fullPath = cms.Schedule(
    process.singleTopSkimPath,
    process.preselection
    )

#Objects included in the pat-tuples
savePatTupleSkimLoose = cms.untracked.vstring(
    'drop *',
    'keep *_*Muons_*_*',
    'keep *_pfMuonsFromVertex_*_*',
    'keep patMuons_selectedPatMuons_*_*',
    'keep patElectrons_selectedPatElectrons_*_*',
    'keep patJets_selectedPatJets_*_*',
    'keep *_selectedPatJets_genJets_*', # to get embedded genJets
    'keep patMETs_patMETs_*_*',
    'keep *_kt6PFJets_rho_*',
    'keep *_topJetsPF_*_*',
    'keep *_topMETsPF_*_*',
    'keep patMuons_vetoMuons_*_*',
    'keep *_vetoElectrons_*_*',
    'keep *_vetoElectronsMVA_*_*',
    'keep patMuons_tightMuons_*_*',
    'keep *_tightElectrons_*_*',
    'keep *_tightElectronsZeroIso_*_*',
    'keep *_tightMuonsZeroIso_*_*',
# vertex
    'keep *_offlineBeamSpot_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep *_goodOfflinePrimaryVertices_*_*', # needed by SingleTopVertexInfoDumper module
# Trigger results
    "keep *_TriggerResults_*_*",
# gen particles
    'keep *_genParticles_*_*',
# gen info
    'keep PileupSummaryInfos_*_*_*',
    'keep GenEventInfoProduct_*_*_*',
    'keep GenRunInfoProduct_*_*_*',
    'keep LHEEventProduct_*_*_*',
    'keep *_genEventScale_*_*',
    'keep *_PDFInfo_*_*',
)


process.singleTopPatTuple = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('singleTopSkim_'+ChannelName+'.root'),
#    SelectEvents   = cms.untracked.PSet(
#      SelectEvents = cms.vstring(
#        'preselection')
#        'pathSelection')
#      ),
    outputCommands = savePatTupleSkimLoose
    )
process.singleTopPatTuple.dropMetaData = cms.untracked.string("DROPPED")


#### Ntuplization step ###
########################################################
process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")
########################################################

######### GET generator info ##############
#genJets:
process.genJetsPF = cms.EDProducer(
    "SingleTopGenJetPtEtaProducer",
    jetsSource = cms.InputTag("topJetsPF"),
)
#genAllJets:
process.genAllJetsPF = cms.EDProducer(
    "SingleTopGenJetPtEtaProducer",
    jetsSource = cms.InputTag("selectedPatJets"),
)
#PU Info
process.NVertices = cms.EDProducer("SingleTopPileUpProducer")

#n gen particles Info
process.NGenParticles = cms.EDProducer("SingleTopNGenParticlesProducer")

#PDF Info
process.PDFInfo = cms.EDProducer( "PDFInfoDumper" )

#Part of MC Truth particles production
process.MCTruthParticles = cms.EDProducer(
    "SingleTopMCProducer",
    genParticlesSource = cms.InputTag("genParticles"),
    lhes = cms.InputTag("source"),
    usesPythiaStatusConvention = cms.untracked.bool(True)
)
#############################################

######### EdmNtuples production ##############
process.load("TopQuarkAnalysis.SingleTop.SingleTopNtuplizers_cff")

# Ntuple sequence
process.genPath = cms.Sequence(
    process.genJetsPF +
    process.genAllJetsPF +
    process.NVertices +
    process.NGenParticles +
    process.PDFInfo           
)

process.singleTopNtuplePath = cms.Sequence(
    process.nTuplesSkim 
)

process.singleTopSkimPath += process.singleTopNtuplePath
if not(isData): process.singleTopSkimPath += process.genPath

from TopQuarkAnalysis.SingleTop.SingleTopNtuplizers_cff import saveNTuplesSkimLoose

#Add MC Truth information:
doMCTruth = True
if isData:
    doMCTruth = False

if doMCTruth:
    process.MCTruth = cms.Sequence(
        process.MCTruthParticles +
        process.nTuplesSkimMCTruth
    )
    process.singleTopSkimPath += process.MCTruth

    saveNTuplesSkimLoose.append('keep  floats_MCTruthParticles_*_*')
    saveNTuplesSkimLoose.append('keep  ints_MCTruthParticles_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCLeptons_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCQuarks_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCNeutrinos_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCBQuarks_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCTops_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCTopsW_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCTopsBQuark_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCTopsLepton_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCTopsNeutrino_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCTopsQuark_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCTopsQuarkBar_*_*')
    #Higgs
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCHiggs_*_*')
    saveNTuplesSkimLoose.append('keep  floats_singleTopMCHiggsBQuark_*_*')

    
## Output module configuration
process.singleTopNTupleOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('singleTopEdmNtuple_'+ChannelName+'.root'),
    #SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('preselection')),
    outputCommands = saveNTuplesSkimLoose,
    )

process.singleTopNTupleOut.dropMetaData = cms.untracked.string("ALL")

process.outpath = cms.EndPath(
    process.singleTopPatTuple +
    process.singleTopNTupleOut
    )

