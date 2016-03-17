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
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V3::All', '')


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
Source_Files = cms.untracked.vstring(
    ## single muons PU=0
    #'/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC12/Extended2023TTI/Muons/NoPU/SingleMuon_E2023TTI_NoPU.root',
    #'/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC12/Extended2023TTI/Muons/NoPU/SingleMuPlus_E2023TTI_NoPU.root'
    ## ttbar PU=140
    '/store/group/upgrade/Tracker/L1Tracking/Synchro/Input/TTbar/CMSSW_6_2_0_SLHC26-DES23_62_V1_LHCCRefPU140-v1/FC5CDAC1-4E2F-E511-8085-0026189438D9.root'
	)
process.source = cms.Source("PoolSource", fileNames = Source_Files)

process.TFileService = cms.Service("TFileService", fileName = cms.string('TTbar_PU140_TrkPerf.root'), closeFileFast = cms.untracked.bool(True))


############################################################
# Path definitions & schedule
############################################################

#run the tracking
BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
#process.TTTracksFromPixelDigis.phiWindowSF = cms.untracked.double(2.0)  ## uncomment this to run with wider projection windows (for electrons)
process.TT_step = cms.Path(process.TrackTriggerTTTracks)
process.TTAssociator_step = cms.Path(process.TrackTriggerAssociatorTracks)


############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
#      pions in jets = 6
#      taus = 15
#      all TPs = 1
############################################################

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save all L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
                                       TP_minNStub = cms.int32(4),       # require TP to have >= minNstub (efficiency denominator)
                                       TP_minPt = cms.double(1.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis", "Level1TTTracks"),               # TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), # MCTruth input 
                                       )
process.ana = cms.Path(process.L1TrackNtuple)

#output module (can use this to store edm-file instead of root-ntuple)
#process.out = cms.OutputModule( "PoolOutputModule",
#                                fileName = cms.untracked.string("FileOut.root"),
#                                fastCloning = cms.untracked.bool( False ),
#                                outputCommands = cms.untracked.vstring('drop *',
#                                                                       'keep *_*_Level1TTTracks_*',
#                                                                       'keep *_*_StubAccepted_*',
#                                                                       'keep *_*_ClusterAccepted_*',
#                                                                       'keep *_*_MergedTrackTruth_*')
#)
#process.FEVToutput_step = cms.EndPath(process.out)


# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023TTI
process = cust_2023TTI(process)

process.schedule = cms.Schedule(process.TT_step,process.TTAssociator_step,process.ana)
