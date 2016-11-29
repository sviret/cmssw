#########################
#
# Configuration file for L1 TAMU combination builder
# using a file with AMPR output 
#
# This script works on any official production sample
# (assuming that this sample contains a container of TTStubs,
# a container of TTClusters, and a container of TrackingParticles)
#
# And of course, a container of matched roads.... (TTTracks) 
#
#
# Author: M. De Mattia
# Date        : 10/2016
#
# Script tested with release CMSSW_6_2_0_SLHC28_patch1
#
#########################

import FWCore.ParameterSet.Config as cms

process = cms.Process('AMCB')

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
    input = cms.untracked.int32(-1)
)

# Input source
#
# You can use as input file the result of the script AMPR_test.py of part 5.2.2 of the tutorial
#
# Any other EDM file containing patterns and produced with CMSSW 620_SLHC28_patch1 should also work
#

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:AM_output.root'),
                            #fileNames = cms.untracked.vstring('root://lyogrid06.in2p3.fr//dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC/PR/ReSynchro/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV_PRonly_PU200_181116/FEF8DD48-289F-E611-B71E-0CC47A4DEEF2_with_AMPR.root'),
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )


# The name of the stub container over which the association is done, please note that the filtered cluster container is
# not associated due to the lack of simPixelDigis in official samples

process.TTStubAssociatorFromPixelDigis.TTStubs        = cms.VInputTag( cms.InputTag("MergeCBOutput", "StubInCB"))
process.TTStubAssociatorFromPixelDigis.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusterAccepted"))
process.TTTrackAssociatorFromPixelDigis.TTTracks      = cms.VInputTag( cms.InputTag("MergeCBOutput", "AML1CBs"))

# Select which combination builder (default False to use the simple combination builder)
process.TTCBsFromPattern.AdvancedCombinationBuilder   = cms.bool(True)

# Additional output definition
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('AMCB_output.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    )
)

# For the moment need to explicitely keep the following containers
# (not yet in the customizing scripts)

# Keep the PR output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMPR')
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMTC')


# Keep the CB output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMCB')
process.RAWSIMoutput.outputCommands.append('drop *_TT*From*_*_*')
process.RAWSIMoutput.outputCommands.append('drop *_TriggerResults_*_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')

# Path and EndPath definitions
process.L1AMCB_step          = cms.Path(process.TTCBsFromPatternswStubs)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)

process.schedule = cms.Schedule(process.L1AMCB_step,process.endjob_step,process.RAWSIMoutput_step)

# Automatic addition of the customisation function

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D

process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

# End of customisation functions
