############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1TrackNtuple")
 
 
############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
Source_Files = cms.untracked.vstring(
    ## ttbar PU=140
    #'/store/group/upgrade/Tracker/L1Tracking/Synchro/Input/TTbar/CMSSW_6_2_0_SLHC26-DES23_62_V1_LHCCRefPU140-v1/FC5CDAC1-4E2F-E511-8085-0026189438D9.root'
    # single muons PU=140, pt=100
    #'/store/group/upgrade/Tracker/L1Tracking/Synchro/Input/Muon/RelValSingleMuPt100_CMSSW_6_2_0_SLHC26-DES23_62_V1_LHCCRefPU140-v1/D8D712F1-3A2F-E511-9716-0025905A605E.root',
    #'/store/group/upgrade/Tracker/L1Tracking/Synchro/Input/Muon/RelValSingleMuPt100_CMSSW_6_2_0_SLHC26-DES23_62_V1_LHCCRefPU140-v1/BA72AD12-392F-E511-8E29-0025905A60B0.root',
    #'/store/group/upgrade/Tracker/L1Tracking/Synchro/Input/Muon/RelValSingleMuPt100_CMSSW_6_2_0_SLHC26-DES23_62_V1_LHCCRefPU140-v1/BE12E362-3C2F-E511-BEBD-0025905A6136.root',
    #'/store/group/upgrade/Tracker/L1Tracking/Synchro/Input/Muon/RelValSingleMuPt100_CMSSW_6_2_0_SLHC26-DES23_62_V1_LHCCRefPU140-v1/BEAF343E-3B2F-E511-BADF-002618943800.root',
    #'/store/group/upgrade/Tracker/L1Tracking/Synchro/Input/Muon/RelValSingleMuPt100_CMSSW_6_2_0_SLHC26-DES23_62_V1_LHCCRefPU140-v1/C2D56BE8-392F-E511-814C-0025905B8598.root',
    
    ## new ttbar sample
    #'file:TTbar_new_1.root',
    #'file:TTbar_new_2.root',
    #'file:TTbar_new_3.root',
    #'file:TTbar_new_4.root',
    #'file:TTbar_new_5.root',
    #'file:TTbar_new_6.root',
    "root://xrootd.unl.edu//store/mc/TTI2023Upg14D/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1-v1/70000/0040AD92-0493-E611-9E16-0025907B4EDC.root"
	)
process.source = cms.Source("PoolSource", fileNames = Source_Files)

process.TFileService = cms.Service("TFileService", fileName = cms.string('ntuple_TTbar_PU0.root'), closeFileFast = cms.untracked.bool(True))


############################################################
# Path definitions & schedule
############################################################

#run the tracking (example for tracklet method)
BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
#process.TTTracksFromPixelDigis.phiWindowSF = cms.untracked.double(2.0)  ## uncomment this to run with wider projection windows
process.TT_step = cms.Path(process.TrackTriggerTTTracks)
process.TTAssociator_step = cms.Path(process.TrackTriggerAssociatorTracks)


############################################################
# primary vertex producer 
############################################################

process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkPrimaryVertexProducer_cfi")
process.L1TkPrimaryVertex.L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis","Level1TTTracks")
process.pL1TkPrimaryVertex = cms.Path( process.L1TkPrimaryVertex )

process.L1TkPrimaryVertexMC = process.L1TkPrimaryVertex.clone()
process.L1TkPrimaryVertexMC.MonteCarloVertex = cms.bool( True )
process.pL1TkPrimaryVertexMC = cms.Path( process.L1TkPrimaryVertexMC )


############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# Valid options are:
#      single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
#      pions in jets from ttbar = 6
#      pions from taus = 15
#      inclusively, store all TPs (also those not from primary interaction, if available in samples you are running on) = 1
############################################################

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       Slim = cms.bool(True),            # only keep the branches we really need
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
                                       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
                                       TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
                                       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
                                       TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.4),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis", "Level1TTTracks"),               # TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), # MCTruth input 
                                       ## isolation stuff 
                                       TrackIsolation = cms.bool(True),
                                       # cuts on the central object (the truth muon & track matched to it)
                                       PTmin = cms.double(20.0),           # central object pt > X GeV, ptmin < 0 means no cut applied
                                       ETAmax = cms.double(2.4),           # central object |eta| < X
                                       TrackPTmin = cms.double(20.0),      # for track matched to central object
                                       TrackETAmax = cms.double(2.4),
                                       TrackChi2max = cms.double(1e10),
                                       TrackNStubmin = cms.int32(4),
                                       # cuts on the tracks (used to determine isolation variable)
                                       IsoTrackZmax = cms.double(25.),     # in cm
                                       IsoTrackChi2max = cms.double(1e10),
                                       IsoTrackNStubmin = cms.int32(4),    
                                       IsoTrackPTmin = cms.double(3.0),    # in GeV
                                       IsoDRmin = cms.double(0.0),
                                       IsoDRmax = cms.double(0.3),
                                       IsoDZmax = cms.double(0.5),         # in cm
                                       ## tracking in jets stuff (--> requires AK4 genjet collection present!)
                                       TrackingInJets = cms.bool(True),
                                       ## save primary vertex information? (--> requires that you ran that above)
                                       PrimaryVertex = cms.bool(True),
                                       )
process.ana = cms.Path(process.L1TrackNtuple)


############################################################
# output module
############################################################

#process.out = cms.OutputModule( "PoolOutputModule",
#                                fileName = cms.untracked.string("FileOut.root"),
#                                fastCloning = cms.untracked.bool( False ),
#                                outputCommands = cms.untracked.vstring('drop *',
#                                                                       'keep *_*_Level1TTTracks_*',
#                                                                       'keep *_*_StubAccepted_*',
#                                                                       'keep *_*_ClusterAccepted_*',
#                                                                       'keep *_*_MergedTrackTruth_*',
#                                                                       'keep *_L1TkPrimaryVertex_*_*',
#                                                                       'keep *_L1TkPrimaryVertexMC_*_*',
#                                                                       'keep *_genParticles_*_*'
#)
#process.FEVToutput_step = cms.EndPath(process.out)


############################################################
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
############################################################

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D
process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)


############################################################
# process schedule
############################################################

process.schedule = cms.Schedule(process.TT_step,process.TTAssociator_step,process.pL1TkPrimaryVertex,process.pL1TkPrimaryVertexMC,process.ana)
