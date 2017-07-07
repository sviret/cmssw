import FWCore.ParameterSet.Config as cms

# First register all the hit matching algorithms, then specify preferred ones at end.

# The stub windows used here refer to the definition provided in this link:
#
# http://sviret.web.cern.ch/sviret/docs/CMS/SLHC/L1TT_XX1215.pdf
#
# Extension to the tilted geometry is discussed here:
#
# https://indico.cern.ch/event/536881/contributions/2219856/
#
# This script is adapted to the very last Tilted Tracker geometry to date (tracker T5)
#
# The tuning is not yet fully validated
# But should enable a good turn-on at 2GeV with reasonable rates at PU200

# Tab2013 hit matching algorithm
TTStubAlgorithm_official_Phase2TrackerDigi_ = cms.ESProducer("TTStubAlgorithm_official_Phase2TrackerDigi_",
   zMatchingPS = cms.bool(True),
   zMatching2S = cms.bool(True),
   BarrelCut = cms.vdouble( 0, 1.5, 2.0, 3.0, 4.5, 5.5, 6.5), #Use 0 as dummy to have direct access using DetId to the correct element 
   NTiltedRings = cms.vdouble( 0., 12., 12., 12., 0., 0., 0.), #Number of tilted rings per side in barrel layers (for tilted geom only)
   TiltedBarrelCutSet = cms.VPSet(
        cms.PSet( TiltedCut = cms.vdouble( 0 ) ),
        cms.PSet( TiltedCut = cms.vdouble( 0, 2, 2, 2, 2, 2, 2, 1.5, 1.5, 1.5, 1.5, 1., 1.) ),
        cms.PSet( TiltedCut = cms.vdouble( 0, 2.5, 2.5, 2.5, 2.5, 2, 2, 2.5, 2.5, 2, 2, 2, 2) ),
        cms.PSet( TiltedCut = cms.vdouble( 0, 3.5, 3.5, 3, 3, 3, 3, 2.5, 2.5, 2.5, 2, 2, 2) ),
	),
   EndcapCutSet = cms.VPSet(
        cms.PSet( EndcapCut = cms.vdouble( 0 ) ),
        cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1.5, 2, 2, 2.5, 3, 3, 3.5, 4, 2.5, 3, 3.5, 4.5, 5.5) ),
        cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1.5, 2, 2, 2, 2.5, 3, 3, 3, 2, 3, 4, 5, 5.5) ),
        cms.PSet( EndcapCut = cms.vdouble( 0, 1.5, 1.5, 2, 2, 2.5, 2.5, 2.5, 3.5, 2.5, 5, 5.5, 6) ),
        cms.PSet( EndcapCut = cms.vdouble( 0, 1.5, 1.5, 1.5, 2, 2, 2, 2, 3, 3, 6, 6, 6.5) ),
        cms.PSet( EndcapCut = cms.vdouble( 0, 1.5, 1.5, 1.5, 1.5, 2, 2, 2, 3, 3, 6, 6, 6.5) ),
        )
)

# CBC3 hit matching algorithm
TTStubAlgorithm_cbc3_Phase2TrackerDigi_ = cms.ESProducer("TTStubAlgorithm_cbc3_Phase2TrackerDigi_",
   zMatchingPS = cms.bool(True),
   zMatching2S = cms.bool(True),
)


# Set the preferred hit matching algorithms.
# We prefer the global geometry algorithm for now in order not to break
# anything. Override with process.TTStubAlgorithm_PSimHit_ = ...,
# etc. in your configuration.
TTStubAlgorithm_Phase2TrackerDigi_ = cms.ESPrefer("TTStubAlgorithm_official_Phase2TrackerDigi_")

