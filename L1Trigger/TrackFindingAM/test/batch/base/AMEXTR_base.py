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

process = cms.Process('AMEXTR')

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
    input = cms.untracked.int32(-1)
)

# Input source
#
# You can use as input file the result of the script AMTC_test.py of part 6.1.2 of the tutorial
#
# Any other EDM file containing TCs and produced with CMSSW 620_SLHC27 should also work
#

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('INPUTFILENAME'),     
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# Additional output definition
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MYGLOBALTAG', '')


# For the moment need to explicitely keep the following containers
# (not yet in the customizing scripts)

process.MIBextraction.doMatch          = True
process.MIBextraction.doMC             = True
process.MIBextraction.doSTUB           = True
process.MIBextraction.doL1TRK          = True

process.MIBextraction.L1pattern_tag    = cms.InputTag( "MergePROutput", "AML1Patterns")
process.MIBextraction.L1tc_tag         = cms.InputTag( "MergeTCOutput", "AML1TCs")
process.MIBextraction.L1track_tag      = cms.InputTag( "MergeDRTCBOutput", "AMTCL1TracksDR")
process.MIBextraction.CLUS_container   = cms.string( "TTStubsFromPixelDigis")
process.MIBextraction.CLUS_name        = cms.string( "ClusterAccepted" )
process.MIBextraction.extractedRootFile= cms.string('OUTPUTFILENAME')


# Path and EndPath definitions
process.p                    = cms.Path(process.MIBextraction)

process.schedule = cms.Schedule(process.p)

# Automatic addition of the customisation function

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D

process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

# End of customisation functions
