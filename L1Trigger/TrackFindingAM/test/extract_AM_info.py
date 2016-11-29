#########
#
# Example script to run the extractor on events containing
# the L1 tracking info produced with the full geometry
# and official software
#
# This script extract the AM pattern reco output just the pattern 
# info or the full stuff if available
#
# Usage: cmsRun extract_AM_info.py
#
# More info:
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto620
#
# Look at part 6.2.1 of the tutorial
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 27/02/2014
# Maj. upd : 17/03/2016
#
# Script tested with release CMSSW_6_2_0_SLHC27
#
#########


import FWCore.ParameterSet.Config as cms

process = cms.Process("MIBextractor")

process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10D_cff')


# Other statements

# Global tag for PromptReco
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# The file you want to extract
process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('file:AMfull_fake.root'),
                            fileNames = cms.untracked.vstring('root://lyogrid06.in2p3.fr//dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC/PR/ReSynchro/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV_PU0_PDDS2_CB_FIT_DRx2_reopt/FA388D72-0693-E611-A088-002590D0B008_with_AMPR_and_FIT.root'),
                            #fileNames = cms.untracked.vstring('file:EDM_SLHC_extr_PILEUP4T_140_8.root'),
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# Load the extracto
process.load("Extractors.RecoExtractor.MIB_extractor_cff")

# Tune some options (see MIB_extractor_cfi.py for details)

process.MIBextraction.doMatch          = True
process.MIBextraction.doMC             = True
process.MIBextraction.doSTUB           = True


process.MIBextraction.doMatch          = True
process.MIBextraction.doMC             = True
process.MIBextraction.doSTUB           = True
process.MIBextraction.doL1TRK          = True

process.MIBextraction.L1pattern_tag    = cms.InputTag( "MergePROutput", "AML1Patterns")
process.MIBextraction.L1tc_tag         = cms.InputTag( "MergeTCOutput", "AML1TCs")
process.MIBextraction.L1track_tag      = cms.InputTag( "MergeDRTCBOutput", "AMTCL1TracksDR")
process.MIBextraction.CLUS_container   = cms.string( "TTStubsFromPixelDigis")
process.MIBextraction.CLUS_name        = cms.string( "ClusterAccepted" )
process.MIBextraction.extractedRootFile= cms.string('extr_tt0_reopt.root')



# You can choose to extract the info from filtered stubs only
#process.MIBextraction.STUB_container   = cms.string( "MergePROutput" )
#process.MIBextraction.STUB_name        = cms.string( "StubInPattern" )
#process.MIBextraction.CLUS_container   = cms.string( "TTStubsFromPixelDigis")
#process.MIBextraction.CLUS_name        = cms.string( "ClusterAccepted" )

#process.MIBextraction.doL1TRK          = True
#process.MIBextraction.L1pattern_tag    = cms.InputTag( "MergePROutput", "AML1Patterns")

# Choose the appropriate lines
#process.MIBextraction.L1track_tag      = cms.InputTag( "", "") # No L1Tracks
#process.MIBextraction.L1tc_tag         = cms.InputTag( "", "") # No TCs
#process.MIBextraction.L1track_tag      = cms.InputTag( "MergeFITOutput", "AML1Tracks") # Floating point tracks
#process.MIBextraction.L1track_tag      = cms.InputTag( "MergeFITOutputb", "AML1BinTracks") # Bit-wise tracks
#process.MIBextraction.L1tc_tag         = cms.InputTag( "MergeTCOutputb", "AML1BinTCs") # Bit-wise TCs
#process.MIBextraction.L1tc_tag         = cms.InputTag( "MergeTCOutput", "AML1TCs") # Floating point TCs

process.p = cms.Path(process.MIBextraction)

