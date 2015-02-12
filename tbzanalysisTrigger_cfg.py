import FWCore.ParameterSet.Config as cms

#from SignalMC import *

process = cms.Process("DEMO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('RecoBTag/Configuration/RecoBTag_cff')
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')

#JetMETCorrections
process.load("JetMETCorrections.Type1MET.caloMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsCaloMet_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0RecoTrack_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")
process.load("JetMETCorrections.Type1MET.correctedMet_cff")

#---------12-05-14-------
#process.load("CommonTools.ParticleFlow.pfIsolation_cfg")
#from CommonTools.ParticleFlow import pfIsolation_cfg
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff') Configuration.Geometry.GeometryIdeal_cff
#process.load("Configuration.StandardSequences.Reconstruction_cff")
# HLT bit HLT_L1_BscMinBiasOR_BptxPlusORMinus

#/////////////////////////////////////////////////////////////////////
#process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
#process.hltHighLevel.HLTPaths = cms.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus')
#triggername = "HLT_IsoMu24_eta2p1_v

#DoubleElectron
#HLT Ele17 CaloIdT CaloIsoVL TrkIdVL TrkIsoVL Ele8 CaloIdT CaloIsoVL TrkIdVL TrkIsoVL

#DoubleMu
#HLT_Mu17_Mu8
#HLT_Mu17_TkMu8

#MuEG
#HLT Mu17 Ele8 CaloIdT CaloIsoVL TrkIdVL TrkIsoVL
#HLT Mu8 Ele17 CaloIdT CaloIsoVL TrkIdVL TrkIsoVL
#////////////////////////////////////////////////////////////////////////

########################################
# Summer12 with CMSSW_5_3_X and GT START53_V7A
# DoubleMuon trigger 
# HLT_Mu17_Mu8_v17 or HLT_Mu17_TkMu8_v10 
# DoubleElectron trigger 
# HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17
########################################

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

    #DATA2012A6
    #DATA2012A7
    #DATA2012A8
    #DATA2012A9
    #Data2101A

    #SignalMCPat	
    #SignalMC
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/A27919CE-486A-E311-8EDD-002590D0B0C8.root'
    #'file:/afs/cern.ch/work/m/meshoaib/NewAnalysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/A27919CE-486A-E311-8EDD-002590D0B0C8.root'
    #2012_EE
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/60CA2BDE-71DA-E111-9275-0024E86E8D4C.root' 
    #2012_DataSet
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/F4FBA42A-A6CF-E111-BAAB-003048678A7E.root'
    #DoublMu
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/F4FBA42A-A6CF-E111-BAAB-003048678A7E.root'
    #MuEG
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/0A47AE1F-BFD9-E111-8AB4-20CF3019DF0C.root'
    #WZ_MC
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/38576F66-61DD-E111-A13B-002481E0DC4E.root'	
    ##'file:/afs/cern.ch/user/m/mwaqas/TBZ14/TBZMET/CMSSW_5_3_15/src/MyAnalysis/TbZ/001C0B68-536A-E311-B25F-002590D0B066.root'
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/001C0B68-536A-E311-B25F-002590D0B066.root'
    #'file:/afs/cern.ch/work/m/meshoaib/Analysis_TBZ/CMSSW_5_3_7/src/FirstAnalysis/TBZAnalysis/A00B400A-B1CF-E111-BCEF-0026189438C2.root'
    #'file:/afs/cern.ch/work/m/meshoaib/Analysis_TBZ/CMSSW_5_3_7/src/FirstAnalysis/TBZAnalysis/441B933C-3846-E211-9928-00266CFFA768.root'
    #'file:/ehep/ehep_data/datasets/DoubleMu_2012RunA/13Jul2012-v1/AC00CD43-B9D2-E111-97FD-001A92971B20.root'    
    
	#'file:/ehep/ehep_data/datasets/MC_2012/TT/0C8518DA-42F1-E111-BD2F-003048D3CB3C.root'
	
	# patTuple

	#'file:/afs/cern.ch/work/m/meshoaib/NewAnalysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/singleTopSkim_TChannel.root'
	#'file:/afs/cern.ch/work/m/meshoaib/PatAnalysis/CMSSW_5_3_21/src/MyAnalysis/TbZ/singleTopSkim_TChannel_74_1_Ajb.root'
	#"file:/afs/cern.ch/work/m/meshoaib/PatAnalysis/CMSSW_5_3_21/src/MyAnalysis/TbZ/singleTopSkim_TChannel.root"
        
	# "file:/afs/cern.ch/work/m/meshoaib/TupleSingleTop/CMSSW_5_3_11/src/TopQuarkAnalysis/SingleTop/test/singleTopSkim_TChannel.root"


#	"file:/afs/cern.ch/work/m/meshoaib/PatAnalysis/CMSSW_5_3_21/src/MyAnalysis/TbZ/singleTopSkim_TChannel.root"	
	
	#SignalMC
	"file:/afs/cern.ch/work/m/meshoaib/PatAnalysis/CMSSW_5_3_21/src/MyAnalysis/TbZ/singleTopSkim_TChannel_77_1_EL1.root"
	
	#QCD 
#"file:/afs/cern.ch/work/m/meshoaib/TupleSingleTop/CMSSW_5_3_11/src/TopQuarkAnalysis/SingleTop/test/crab_QCD_TestTuple/res/singleTopSkim_TChannel_100_1_TQr.root"

	#"file:/afs/cern.ch/work/m/meshoaib/PatAnalysis/CMSSW_5_3_21/src/MyAnalysis/TbZ/NoMetJetCorresingleTopSkim_TChannel.root"
	#"file:/afs/cern.ch/work/m/meshoaib/PatAnalysis/CMSSW_5_3_21/src/MyAnalysis/TbZ/WithAppliedCorr_singleTopSkim_TChannel.root"
#	"file:/afs/cern.ch/work/m/meshoaib/PatAnalysis/CMSSW_5_3_21/src/MyAnalysis/TbZ/singleTopSkim_TChannel_Synchro.root"
	#input from EOS directly

	#For-Systematic Checking 
	#"file:/afs/cern.ch/work/m/meshoaib/TupleSingleTop/CMSSW_5_3_11/src/TopQuarkAnalysis/SingleTop/test/WithoutJetEn_singleTopSkim_TChannel.root"
	#"file:/afs/cern.ch/work/m/meshoaib/TupleSingleTop/CMSSW_5_3_11/src/TopQuarkAnalysis/SingleTop/test/JetEnUp_singleTopSkim_TChannel.root"
	#"file:/afs/cern.ch/work/m/meshoaib/TupleSingleTop/CMSSW_5_3_11/src/TopQuarkAnalysis/SingleTop/test/JetEnDown_singleTopSkim_TChannel.root"


    )
)

#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')## To delete duplicate events, this included checkAllFilesOpened
process.source.duplicateCheckMode = cms.untracked.string('checkAllFilesOpened')
process.TFileService = cms.Service("TFileService",
        fileName= cms.string("ttbarMC_Pat.root")
)

print "before" 
#process.myBJetsProducer = cms.EDProducer('SingleTopBJetsProducer'
#                                          )
#process.myProducerLabel = cms.EDProducer('TbzMuonProducer',
#                                         muonTag = cms.InputTag("muons"),
#                                         muonPtCut = cms.double(20.)    ,
#                                         muonEtaCut = cms.double(2.1)                                         
#)                                                                                  
#process.myLoseMuons = cms.EDProducer('TbzLoseMuonProducer',
#                                         muonTag = cms.InputTag("muons"),
#                                         muonPtCut = cms.double(10.)    ,
#                                         muonEtaCut = cms.double(2.1)                                         
#)                                                                 
#process.myJetProdLabel = cms.EDProducer('JetProducer',
#                                         JetsPtCut  = cms.double(30.), #as suggested by jets contact person
#                                         JetsEtaCut = cms.double(3.0)
#                                        )
#process.selectedMuonsGenParticlesMatch = cms.EDProducer( "MCTruthDeltaRMatcher",
 #                                             src = cms.InputTag("TbzMuonProducer"),
 #                                             matched = cms.InputTag("genParticleCandidates"),
 #                                             distMin = cms.double(0.15),
 #                                             matchPDGId = cms.vint32(13)
 #                                             )
# rho value for isolation FastjetJetProducer::kt6PFJetsForIsolation
#

#from RecoJets.JetProducers.kt4PFJets_cfi import *
#process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
#process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

# particle flow isolation
#

#from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
#process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
#process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)
#
#process.myElectronProdLabel = cms.EDProducer('ElectronProducer'  ,
#            electronsInputTag = cms.InputTag("gsfElectrons")     ,
#            muonsInputTag = cms.InputTag("myLoseMuons")      ,
#            conversionsInputTag = cms.InputTag("allConversions") ,
#            beamSpotInputTag  = cms.InputTag("offlineBeamSpot")  ,
#            rhoIsoInputTag   = cms.InputTag("kt6PFJetsForIsolation","rho"),
#            primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
#            isoValInputTags =  cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
#                                                                 cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
#                                                                 cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),	
#                                             electronPtCut = cms.double(20.),
#                                             electronEtaCut = cms.double(2.1),
#                                             printDebug              = cms.bool(True)
#)
print "after" 

#process.demo = cms.EDAnalyzer('TBZAnalysis',
#triggerSelection = cms.string('HLT_Mu17_Mu8_v* OR HLT_IsoMu24_eta2p1_v*'),
#triggerConfiguration =  cms.PSet(
  #hltResults = cms.InputTag('TriggerResults','','HLT'),
  #l1tResults = cms.InputTag(''),
  #daqPartitions = cms.uint32(1),
 # l1tIgnoreMask = cms.bool( False ),
#  l1techIgnorePrescales = cms.bool( False ),
#  throw  = cms.bool( True )
#)
#)

#process.Z2mumu = cms.EDAnalyzer('Z2mumu')
#process.Z2ep   = cms.EDAnalyzer('Z2ep')
#process.Acceptance  = cms.EDAnalyzer('tbz_Acceptance')
#process.topFromWb = cms.EDAnalyzer('t2Wb')
#process.tbZ = cms.EDAnalyzer('tbZ_Final')


print "Before top analyzer with loose & tight"


process.topAnaLoosetight = cms.EDAnalyzer('TbZTopAnalyzer_tightLoose',

MutriggerSelection   = cms.string('(HLT_Mu17_Mu8_* OR HLT_Mu17_TkMu8_*)'),
MuEGtriggerSelection = cms.string('(HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_* OR HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),
EtriggerSelection    = cms.string('(HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),

triggerConfiguration =  cms.PSet(
  hltResults = cms.InputTag('TriggerResults','','HLT'),
  l1tResults = cms.InputTag(''),
  daqPartitions = cms.uint32(1),
  l1tIgnoreMask = cms.bool( False ),
  l1techIgnorePrescales = cms.bool( False ),
  throw  = cms.bool( True )
),

 #    Producer     = cms.InputTag("myBJetsProducer")              ,
                               beamSpotInputTag  = cms.InputTag("offlineBeamSpot")             ,
                               primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices")  ,

                               # rhoIsoInputTag   = cms.InputTag("kt6PFJetsForIsolation","rho"),

                               # isoValInputTags =  cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                               # cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                               #  cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),


                                metPtCut      = cms.double(30.)                                    ,# before it was 30 GeV
                                vertexSrc     = cms.InputTag("offlinePrimaryVertices")             ,

                                BtagEtaCut    = cms.double(2.4)                                    ,
                                BtagPtCut     = cms.double(30)                                     ,

                                #BtagDiscrCut  = cms.double(0.679)                                  ,
                               # DPHiENue      = cms.double(1.0)                                    ,
                                DPHiMuNue     = cms.double(1.0)                                    ,
                                MaxZMass      = cms.double(102.0)                                  ,
                                MinZMAss      = cms.double(78.0)                                   ,
                                JetsPtCut     = cms.double(30.)                                    ,
                                JetsEtaCut    = cms.double(3.0)                                    ,
                                ElecPtCut     = cms.double(20.)                                    ,
                                ElecEtaCut    = cms.double(2.5)                                    ,
                                muonPtCut     = cms.double(20)                                     ,
                                muonEtaCut    = cms.double(2.4)                                    ,
                                doTruthMatch  = cms.bool(False)                                    ,
                                #realdata      = cms.bool(True)                                    ,
                                realdata      = cms.bool(False)                                    ,
                                doPileup      = cms.bool(False)                                    ,
                                printDebug    = cms.bool(True)

                              )
print "After top analyzer with Loose & tight Lept"


print "Beffore top analyzer with Loose Lept"

############# top analyzer with Loose #####

process.topAnaLooseLept = cms.EDAnalyzer('TbZTopAnalyzer_LooseLept',

MutriggerSelection   = cms.string('(HLT_Mu17_Mu8_* OR HLT_Mu17_TkMu8_*)'),
MuEGtriggerSelection = cms.string('(HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_* OR HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),
EtriggerSelection    = cms.string('(HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),

triggerConfiguration =  cms.PSet(
  hltResults = cms.InputTag('TriggerResults','','HLT'),
  l1tResults = cms.InputTag(''),
  daqPartitions = cms.uint32(1),
  l1tIgnoreMask = cms.bool( False ),
  l1techIgnorePrescales = cms.bool( False ),
  throw  = cms.bool( True )
),


			   #    Producer     = cms.InputTag("myBJetsProducer")              ,
                               beamSpotInputTag  = cms.InputTag("offlineBeamSpot")             ,
                               primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices")  ,

                               # rhoIsoInputTag   = cms.InputTag("kt6PFJetsForIsolation","rho"),

                               # isoValInputTags =  cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                               # cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                               #  cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),


                                metPtCut      = cms.double(30.)                                    , # before it was 30. GeV
                                vertexSrc     = cms.InputTag("offlinePrimaryVertices")             ,

                                BtagEtaCut    = cms.double(2.4)                                    ,
                                BtagPtCut     = cms.double(30)                                     ,

                                #BtagDiscrCut  = cms.double(0.679)                                  ,
                               # DPHiENue      = cms.double(1.0)                                    ,
                                DPHiMuNue     = cms.double(1.0)                                    ,
                                MaxZMass      = cms.double(102.0)                                  ,
                                MinZMAss      = cms.double(78.0)                                   ,
                                JetsPtCut     = cms.double(30.)                                    ,
                                JetsEtaCut    = cms.double(3.0)                                    ,
                                ElecPtCut     = cms.double(20.)                                    ,
                                ElecEtaCut    = cms.double(2.5)                                    ,
                                muonPtCut     = cms.double(20)                                     ,
                                muonEtaCut    = cms.double(2.4)                                    ,
                                doTruthMatch  = cms.bool(False)                                    ,
                                #realdata      = cms.bool(True)                                    ,
                                realdata      = cms.bool(False)                                    ,
                                doPileup      = cms.bool(False)                                    ,
                                printDebug    = cms.bool(True)

                              )
print "After top analyzer with Loose Lept"
###### top Analyzer using tight leptons with Zero Isolation ######################

process.topAna_lepZeroIso = cms.EDAnalyzer('TbZTopAnalyzer_leptonsZeroIso',



MutriggerSelection   = cms.string('(HLT_Mu17_Mu8_* OR HLT_Mu17_TkMu8_*)'),
MuEGtriggerSelection = cms.string('(HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_* OR HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),
EtriggerSelection    = cms.string('(HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),

triggerConfiguration =  cms.PSet(
  hltResults = cms.InputTag('TriggerResults','','HLT'),
  l1tResults = cms.InputTag(''),
  daqPartitions = cms.uint32(1),
  l1tIgnoreMask = cms.bool( False ),
  l1techIgnorePrescales = cms.bool( False ),
  throw  = cms.bool( True )
),

			       bJetProducer     = cms.InputTag("myBJetsProducer")              ,
                               beamSpotInputTag  = cms.InputTag("offlineBeamSpot")             ,
                               primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices")  ,


			       metPtCut      = cms.double(30.)                                    ,# before it was 30. GeV
                               vertexSrc     = cms.InputTag("offlinePrimaryVertices")             ,
                               BtagEtaCut    = cms.double(2.4)                                    ,
                               BtagPtCut     = cms.double(30)                                     ,
			        DPHiMuNue     = cms.double(1.0)                                    ,
                                MaxZMass      = cms.double(102.0)                                  ,
                                MinZMAss      = cms.double(78.0)                                   ,
                                JetsPtCut     = cms.double(30.)                                    ,
                                JetsEtaCut    = cms.double(3.0)                                    ,
                                ElecPtCut     = cms.double(20.)                                    ,
                                ElecEtaCut    = cms.double(2.5)                                    ,
                                muonPtCut     = cms.double(20)                                     ,
                                muonEtaCut    = cms.double(2.4)                                    ,
                                doTruthMatch  = cms.bool(False)                                    ,
                                #realdata      = cms.bool(True)                                    ,
                                realdata      = cms.bool(False)                                    ,
                                doPileup      = cms.bool(False)                                    ,
                                printDebug    = cms.bool(True)

                              )
print "after topAna_leptonsZeroIso"

###################################  Systematcis    ##################################################################################
process.topAna_Systematic = cms.EDAnalyzer('TbZTopAnalyzer_Systematic',

MutriggerSelection   = cms.string('(HLT_Mu17_Mu8_* OR HLT_Mu17_TkMu8_*)'),
MuEGtriggerSelection = cms.string('(HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_* OR HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),
EtriggerSelection    = cms.string('(HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),

triggerConfiguration =  cms.PSet(
  hltResults = cms.InputTag('TriggerResults','','HLT'),
  l1tResults = cms.InputTag(''),
  daqPartitions = cms.uint32(1),
  l1tIgnoreMask = cms.bool( False ),
  l1techIgnorePrescales = cms.bool( False ),
  throw  = cms.bool( True )
),


  bJetProducer     = cms.InputTag("myBJetsProducer")              ,
                               beamSpotInputTag  = cms.InputTag("offlineBeamSpot")             ,
                               primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices")  ,


                               metPtCut      = cms.double(30.)                                    ,# before it was 30. GeV
                               vertexSrc     = cms.InputTag("offlinePrimaryVertices")             ,
                               BtagEtaCut    = cms.double(2.4)                                    ,
                               BtagPtCut     = cms.double(30)                                     ,
                                DPHiMuNue     = cms.double(1.0)                                    ,
                                MaxZMass      = cms.double(102.0)                                  ,
                                MinZMAss      = cms.double(78.0)                                   ,
                                JetsPtCut     = cms.double(30.)                                    ,
                                JetsEtaCut    = cms.double(3.0)                                    ,
                                ElecPtCut     = cms.double(20.)                                    ,
                                ElecEtaCut    = cms.double(2.5)                                    ,
                                muonPtCut     = cms.double(20)                                     ,
                                muonEtaCut    = cms.double(2.4)                                    ,
                                doTruthMatch  = cms.bool(False)                                    ,
                                #realdata      = cms.bool(True)                                    ,
                                realdata      = cms.bool(False)                                    ,
                                doPileup      = cms.bool(False)                                    ,
                                printDebug    = cms.bool(True)

                              )
print "after topAna_Systematic"

####top analyzer#####
process.topAna = cms.EDAnalyzer('TbZTopAnalyzer',

#MutriggerSelection   = cms.string('HLT_Mu17_Mu8_* OR HLT_Mu17_TkMu8_* AND !( HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),
#MutriggerSelection   = cms.string('(HLT_Mu17_Mu8_*)'),
#MutriggerSelection   = cms.string('(HLT_Mu17_Mu8_* OR HLT_Mu17_TkMu8_* )'),
#MutriggerSelection   = cms.string('( HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_* )'),
#MutriggerSelection = cms.string('(HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_* OR HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),
#EtriggerSelection    = cms.string('(HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),
# HLT_Mu17_Mu8_v16
# HLT_Mu17_TkMu8_v9

MutriggerSelection   = cms.string('(HLT_Mu17_Mu8_* OR HLT_Mu17_TkMu8_*)'),
MuEGtriggerSelection = cms.string('(HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_* OR HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),
EtriggerSelection    = cms.string('(HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),

triggerConfiguration =  cms.PSet(
  hltResults = cms.InputTag('TriggerResults','','HLT'),
  l1tResults = cms.InputTag(''),
  daqPartitions = cms.uint32(1),
  l1tIgnoreMask = cms.bool( False ),
  l1techIgnorePrescales = cms.bool( False ),
  throw  = cms.bool( True )
),

                               # electronsInputTag = cms.InputTag("myElectronProdLabel"),
                               # electronsInputTag = cms.InputTag("tightElectrons"),
                               # muonsInputTag = cms.InputTag("myProducerLabel"), 
                               # muonsInputTag = cms.InputTag("muons"),
                               # jetsInputTag  = cms.InputTag("ak5PFJets"),#ak5PFJetsak5CaloJets
                               # bjetInputTag  = cms.InputTag("combinedSecondaryVertexBJetTags"),
                               #  metInputTag   = cms.InputTag("pfMet")               ,
                               # conversionsInputTag = cms.InputTag("allConversions") ,
                               
                               
                               bJetProducer     = cms.InputTag("myBJetsProducer")              ,
                               beamSpotInputTag  = cms.InputTag("offlineBeamSpot")             ,
                               primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices")  ,
                                
                               # rhoIsoInputTag   = cms.InputTag("kt6PFJetsForIsolation","rho"),
           
                               # isoValInputTags =  cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                               # cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                               #  cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),                                
                                
                                
                                metPtCut      = cms.double(30.)                                    ,# before it was 30. GeV
                                vertexSrc     = cms.InputTag("offlinePrimaryVertices")             ,
                                
				BtagEtaCut    = cms.double(2.4)                                    ,
                                BtagPtCut     = cms.double(30)                                     ,

                                #BtagDiscrCut  = cms.double(0.679)                                  ,
                               # DPHiENue      = cms.double(1.0)                                    ,
                                DPHiMuNue     = cms.double(1.0)                                    ,
                                MaxZMass      = cms.double(102.0)                                  ,
                                MinZMAss      = cms.double(78.0)                                   ,
                                JetsPtCut     = cms.double(30.)                                    ,
                                JetsEtaCut    = cms.double(3.0)                                    ,
                                ElecPtCut     = cms.double(20.)                                    ,
				ElecEtaCut    = cms.double(2.5)                                    ,
			        muonPtCut     = cms.double(20)                                     ,		
				muonEtaCut    = cms.double(2.4)                                    ,
                                doTruthMatch  = cms.bool(False)                                    ,			                             
                                #realdata      = cms.bool(True)                                    ,
                                realdata      = cms.bool(False)                                    ,
                                doPileup      = cms.bool(False)                                    ,
                                printDebug    = cms.bool(True)

                              )
print "after 1" 
process.dump=cms.EDAnalyzer('EventContentAnalyzer')
print "after 2" 
#process.p = cms.Path(process.myProducerLabel*process.demo*process.dump )
#process.p = cms.Path(process.myProducerLabel*process.myJetProdLabel*process.kt6PFJetsForIsolation*process.pfiso*process.myElectronProdLabel*process.Z2mumu*process.Z2ep*process.Acceptance*process.topAna)
#process.p = cms.Path(process.myProducerLabel*process.myLoseMuons*process.myJetProdLabel*process.kt6PFJetsForIsolation*process.pfiso*process.myElectronProdLabel*process.topAna)
process.p = cms.Path(process.topAnaLoosetight*process.topAnaLooseLept*process.topAna_lepZeroIso*process.topAna_Systematic*process.topAna)
#process.p = cms.Path(process.myProducerLabel*process.myJetProdLabel*process.kt6PFJetsForIsolation*process.pfiso*process.myElectronProdLabel*process.Z2mumu*process.Z2ep*process.topAna)
print "after 3" 
