#########################
#
# Base file for stub production
#
# This script works on any official production sample
#
# Author: S.Viret (viret@in2p3.fr)
# Date        : 02/12/2015
#
# Script tested with release CMSSW_6_2_0_SLHC26
#
#########################

import FWCore.ParameterSet.Config as cms

process = cms.Process('STUBS')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('INPUTFILENAME'),
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# Additional output definition
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MYGLOBALTAG', '')

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('output.root'), ## ADAPT IT ##
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    )
)

#Stub window tuning adapted to FE extraction capability
#
#
# The breakdown is the following
#
# 2 GeV official cut : all the 2S modules except TOB1 layer, rings 4 to 9 of the disks
# 2 GeV tight cut    : TIB3/TOB1 layers
# 3 GeV cut          : Rings 1 to 3 in the disks, TIB1 and TIB2

process.TTStubAlgorithm_tab2013_PixelDigi_.BarrelCut = cms.vdouble( 0, 1.5, 1.5, 2.5, 4, 5.5, 6.5)
process.TTStubAlgorithm_tab2013_PixelDigi_.EndcapCutSet = cms.VPSet(
                                                                    cms.PSet( EndcapCut = cms.vdouble( 0 ) ),
                                                                    cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1.5, 2.0, 2.5, 2.5, 2.5, 3.0, 3.5, 4.5, 3.0, 3.5, 4.0, 4.5, 5.0 ) ),
                                                                    cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1.5, 2.0, 2.0, 2.5, 2.5, 2.5, 3.0, 4.0, 2.5, 3.0, 3.5, 4.0, 4.5 ) ),
                                                                    cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1, 2.0, 2.0, 2.0, 2.5, 2.5, 2.5, 3.5, 4.0, 2.5, 3.0, 3.5, 4.0 ) ),
                                                                    cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1, 2.0, 2.0, 2.0, 2.0, 2.5, 2.5, 3.0, 3.5, 2.5, 2.5, 3.0, 3.5 ) ),
                                                                    cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1, 1.5, 1.5, 2.0, 2.0, 2.0, 2.5, 2.5, 3.0, 3.5, 2.5, 2.5, 3.0 ) ),
                                                                    )

# For the moment need to explicitely keep the following containers
# (not yet in the customizing scripts)

process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')

# Path and EndPath definitions
process.L1TrackTrigger_step  = cms.Path(process.TrackTriggerClustersStubs)
process.L1TTAssociator_step  = cms.Path(process.TrackTriggerAssociatorClustersStubs)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)

process.schedule = cms.Schedule(process.L1TrackTrigger_step,process.L1TTAssociator_step,process.endjob_step,process.RAWSIMoutput_step)


# Automatic addition of the customisation function

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D

process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

# End of customisation functions	


