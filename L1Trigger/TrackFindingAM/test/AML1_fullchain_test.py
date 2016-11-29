#########################
#
# Configuration file for L1 PCA fit
# using a file with AMTC output 
#
# This script works on any official production sample
# (assuming that this sample contains a container of TTStubs,
# a container of TTClusters, and a container of TrackingParticles)
#
# And of course, a container of TCs.... (TTTracks) 
#
#
# Author: S.Viret (viret@in2p3.fr)
# Date        : 04/03/2016
#
# Script tested with release CMSSW_6_2_0_SLHC27
#
#########################

import FWCore.ParameterSet.Config as cms

process = cms.Process('AMFULL')

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

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
#
# You can use as input file the result of the script AMTC_test.py of part 6.1.2 of the tutorial
#
# Any other EDM file containing TCs and produced with CMSSW 620_SLHC27 should also work
#

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('root://lyogrid06.in2p3.fr//dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC/PR/ReSynchro/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV_PRonly_PU200_181116/FEF8DD48-289F-E611-B71E-0CC47A4DEEF2_with_AMPR.root'),
                            fileNames = cms.untracked.vstring('root://lyogrid06.in2p3.fr//dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC/PR/ReSynchro/MuGunFlatPt8to100_PRonly_PU0_181116/FA9F7BF8-80AB-E611-A47C-02163E0130E8_with_AMPR.root'),
                            #fileNames = cms.untracked.vstring('file:AM_output.root'),
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# The name of the stub container over which the association is done, please note that the filtered cluster container is
# not associated due to the lack of simPixelDigis in official samples

process.TTTrackAssociatorFromPixelDigis.TTTracks      = cms.VInputTag( cms.InputTag("MergeFITOutput", "AMTCL1Tracks"), cms.InputTag("MergeFITACBOutput", "AMACBL1Tracks"),cms.InputTag("MergeFITSCBOutput", "AMSCBL1Tracks"), cms.InputTag("MergeTCOutput", "AML1TCs"),cms.InputTag("MergeACBOutput", "AML1ACBs"),cms.InputTag("MergeSCBOutput", "AML1SCBs"), cms.InputTag("MergeDRTCBOutput", "AMTCL1TracksDR"), cms.InputTag("MergeDRACBOutput", "AMACBL1TracksDR"),cms.InputTag("MergeDRSCBOutput", "AMSCBL1TracksDR"))

process.TTTracksTAMUFromTC.ConstantsDir   = cms.FileInPath("L1Trigger/TrackFindingAM/data/PreEstimate_Transverse/matrixVD_2016.txt")
process.TTTracksTAMUFromACB.ConstantsDir  = cms.FileInPath("L1Trigger/TrackFindingAM/data/PreEstimate_Transverse/matrixVD_2016.txt")
process.TTTracksTAMUFromSCB.ConstantsDir  = cms.FileInPath("L1Trigger/TrackFindingAM/data/PreEstimate_Transverse/matrixVD_2016.txt")

#process.TTTCsFromPattern.maxStubsperRoad = cms.int32(18)
#process.TTTCsFromPattern.maxRoadsperTower= cms.int32(250)
#process.TTTCsFromPattern.maxSeedsperTC= cms.int32(12)
#process.TTTCsFromPattern.maxTCsperTower= cms.int32(100)
#DuplicateRemovalFromTCB0.ParameterBased= cms.bool(True),


# Additional output definition
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('mu0.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    )
)

# For the moment need to explicitely keep the following containers
# (not yet in the customizing scripts)

# Keep the necessary output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMPR')
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMFULL')
process.RAWSIMoutput.outputCommands.append('drop *_TT*From*_*_*')
process.RAWSIMoutput.outputCommands.append('drop *_Dupl*From*_*_*')
process.RAWSIMoutput.outputCommands.append('drop *_TriggerResults_*_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')


# Path and EndPath definitions
process.L1AMFULL_step        = cms.Path(process.TTTCsFromPattern*process.TTACBsFromPattern*process.TTSCBsFromPattern
                                        *process.MergeTCOutput*process.MergeACBOutput*process.MergeSCBOutput
                                        *process.TTTracksTAMUFromTC*process.TTTracksTAMUFromACB*process.TTTracksTAMUFromSCB
                                        *process.MergeFITOutput*process.MergeFITACBOutput*process.MergeFITSCBOutput
#                                        *process.DuplicateRemovalswStubs)
                                        *process.DuplicateRemovalsx2wStubs)

process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)

# Uncomment this to use the parameter-based duplicate removal
#process.DuplicateRemovalFromTCB.ParameterBased = True
#process.DuplicateRemovalFromACB.ParameterBased = True
#process.DuplicateRemovalFromSCB.ParameterBased = True
#process.DuplicateRemovalFromTCB0.ParameterBased = True
#process.DuplicateRemovalFromACB0.ParameterBased = True
#process.DuplicateRemovalFromSCB0.ParameterBased = True

# Use this to change the maximum number of common stubs for the stub-based duplicate removal
# process.DuplicateRemovalFromTCB.MaxCommonStubs = 2
# process.DuplicateRemovalFromACB.MaxCommonStubs = 2
# process.DuplicateRemovalFromSCB.MaxCommonStubs = 2

#process.DuplicateRemovalFromTCB.TTTrackName = cms.InputTag("MergeFITOutput", "AMTCL1Tracks")
#process.DuplicateRemovalFromSCB.TTTrackName = cms.InputTag("MergeFITSCBOutput", "AMSCBL1Tracks")
#process.DuplicateRemovalFromACB.TTTrackName = cms.InputTag("MergeFITACBOutput", "AMACBL1Tracks")

process.schedule = cms.Schedule(process.L1AMFULL_step,process.endjob_step,process.RAWSIMoutput_step)


# Automatic addition of the customisation function

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D

process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

# End of customisation functions
