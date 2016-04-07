import FWCore.ParameterSet.Config as cms

process = cms.Process('AMFITBASE')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('L1Trigger.TrackFindingAM.L1AMTrack_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load("Extractors.RecoExtractor.MIB_extractor_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(NEVTS)
)

# Input source
#
# You can use as input file the result of the script AMPR_test.py of part 5.2.2 of the tutorial
#
# Any other EDM file containing patterns and produced with CMSSW 620_SLHC13 should also work
#

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('INPUTFILENAME'),     
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# Additional output definition
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MYGLOBALTAG', '')

process.TFileService = cms.Service("TFileService", fileName = cms.string('OFFFILENAME'), closeFileFast = cms.untracked.bool(True))

# The name of the stub container over which the association is done, please note that the filtered cluster container is
# not associated due to the lack of simPixelDigis in official samples

process.TTStubAssociatorFromPixelDigis.TTStubs        = cms.VInputTag( cms.InputTag("MergeFITOutputI", "StubInTrackI"))
process.TTStubAssociatorFromPixelDigis.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusterAccepted"))
process.TTTrackAssociatorFromPixelDigis.TTTracks      = cms.VInputTag( cms.InputTag("MergeFITOutputI", "AML1TracksI"))


process.TTTracksTAMUFromTC.ConstantsDir                = cms.FileInPath("L1Trigger/TrackFindingAM/data/PreEstimate_Transverse/matrixVD_2016.txt")

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('OUTPUTFILENAME'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    )
)

# For the moment need to explicitely keep the following containers
# (not yet in the customizing scripts)

process.MIBextraction.doMatch          = True
process.MIBextraction.doMC             = True
process.MIBextraction.doSTUB           = True
process.MIBextraction.doL1TRK          = True

process.MIBextraction.L1pattern_tag    = cms.InputTag( "MergePROutput", "AML1Patterns")
process.MIBextraction.L1tc_tag         = cms.InputTag( "MergeTCOutputb", "AML1BinTCs")
process.MIBextraction.L1track_tag      = cms.InputTag( "MergeFITOutputI", "AML1TracksI")
process.MIBextraction.CLUS_container   = cms.string( "TTStubsFromPixelDigis")
process.MIBextraction.CLUS_name        = cms.string( "ClusterAccepted" )
process.MIBextraction.extractedRootFile= cms.string('EXTRFILENAME')

# Keep the PR and TC outputs
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMPR*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMTC*')

# Keep the FIT output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMFIT*')
process.RAWSIMoutput.outputCommands.append('drop *_TTTracks*From*_*_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save all L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(5),         # use 4 or 5-parameter L1 track fit ??
                                       TP_minNStub = cms.int32(4),       # require TP to have >= minNstub (efficiency denominator)
                                       TP_minNStubLayer = cms.int32(4),
                                       TP_minPt = cms.double(1.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag("MergeFITOutputI", "AML1TracksI"),               ## TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "AML1TracksI") ## MCTruth input 

                                       )

# Path and EndPath definitions
process.L1AMFIT_step         = cms.Path(process.TTTracksIFromTCswStubs)
process.p                    = cms.Path(process.MIBextraction)
process.ana                  = cms.Path(process.L1TrackNtuple)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)

process.schedule = cms.Schedule(process.L1AMFIT_step,process.p,process.ana,process.endjob_step,process.RAWSIMoutput_step)

#process.schedule = cms.Schedule(process.L1AMFIT_step,process.p,process.endjob_step,process.RAWSIMoutput_step)


# Automatic addition of the customisation function

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D

process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

# End of customisation functions
