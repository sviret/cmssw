import FWCore.ParameterSet.Config as cms

TTStubsFromPhase2TrackerDigis = cms.EDProducer("TTStubBuilder_Phase2TrackerDigi_",
    TTClusters = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterInclusive"),
    OnlyOnePerInputCluster = cms.bool(True),
    FEineffs      = cms.bool(True),
    CBClimit      = cms.uint32(3),
    MPAlimit      = cms.uint32(5),
    SS5GCIClimit  = cms.uint32(16),
    PS5GCIClimit  = cms.uint32(17),
    SS10GCIClimit = cms.uint32(32),
    PS10GCIClimit = cms.uint32(35)
)


