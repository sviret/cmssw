/*! \class TTStubBuilder
* \brief Plugin to load the Stub finding algorithm and produce the
* collection of Stubs that goes in the event content.
* \details After moving from SimDataFormats to DataFormats,
* the template structure of the class was maintained
* in order to accomodate any types other than PixelDigis
* in case there is such a need in the future.
*
* \author Andrew W. Rose
* \author Nicola Pozzobon
* \author Ivan Reid
* \date 2013, Jul 18
*
*/

#ifndef L1_TRACK_TRIGGER_STUB_BUILDER_H
#define L1_TRACK_TRIGGER_STUB_BUILDER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "L1Trigger/TrackTrigger/interface/TTStubAlgorithm.h"
#include "L1Trigger/TrackTrigger/interface/TTStubAlgorithmRecord.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include <memory>
#include <map>
#include <vector>

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

template< typename T >
class TTStubBuilder : public edm::EDProducer
{
  public:
    /// Constructor
    explicit TTStubBuilder( const edm::ParameterSet& iConfig );

    /// Destructor;
    ~TTStubBuilder();

  private:
    /// Data members
    edm::ESHandle< TTStubAlgorithm< T > > theStubFindingAlgoHandle;
    edm::EDGetTokenT< edmNew::DetSetVector< TTCluster< T > > > clustersToken;
    bool ForbidMultipleStubs;

    bool applyFE;
    unsigned int  maxStubs_2S; // Every BX
    unsigned int  maxStubs_PS; // Every 2BX
    unsigned int  maxStubs_2S_CIC_5; // Every 8BX
    unsigned int  maxStubs_PS_CIC_5; // Every 8BX
    unsigned int  maxStubs_2S_CIC_10; // Every 8BX
    unsigned int  maxStubs_PS_CIC_10; // Every 8BX

    /// Mandatory methods
    virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
    virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
    virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

    /// Sorting method for stubs
    /// NOTE: this must be static!
    static bool SortStubBendPairs( const std::pair< unsigned int, double >& left, const std::pair< unsigned int, double >& right );
    static bool SortStubsBend( const TTStub< T >& left, const TTStub< T >& right );

    int ievt;

    std::unordered_map< int, std::vector< TTStub< Ref_Phase2TrackerDigi_ > > > moduleStubs_CIC; /// Temporary storage for stubs before max check
    std::unordered_map< int, int > moduleStubs_MPA; /// Temporary storage for stubs before max check
    std::unordered_map< int, int > moduleStubs_CBC; /// Temporary storage for stubs before max check

}; /// Close class

/*! \brief Implementation of methods
* \details Here, in the header file, the methods which do not depend
* on the specific type <T> that can fit the template.
* Other methods, with type-specific features, are implemented
* in the source file.
*/

/// Constructors
template< typename T >
TTStubBuilder< T >::TTStubBuilder( const edm::ParameterSet& iConfig )
{
  clustersToken = consumes< edmNew::DetSetVector< TTCluster< T > > >(iConfig.getParameter< edm::InputTag >( "TTClusters" ));
  ForbidMultipleStubs = iConfig.getParameter< bool >( "OnlyOnePerInputCluster" );
  applyFE             = iConfig.getParameter< bool >( "FEineffs" );
  maxStubs_2S         = iConfig.getParameter< uint32_t >( "CBClimit" );
  maxStubs_PS         = iConfig.getParameter< uint32_t >( "MPAlimit" );
  maxStubs_2S_CIC_5   = iConfig.getParameter< uint32_t >( "SS5GCIClimit" );
  maxStubs_PS_CIC_5   = iConfig.getParameter< uint32_t >( "PS5GCIClimit" );
  maxStubs_2S_CIC_10  = iConfig.getParameter< uint32_t >( "SS10GCIClimit" );
  maxStubs_PS_CIC_10  = iConfig.getParameter< uint32_t >( "PS10GCIClimit" );
  produces< edmNew::DetSetVector< TTCluster< T > > >( "ClusterAccepted" );
  produces< edmNew::DetSetVector< TTStub< T > > >( "StubAccepted" );
  produces< edmNew::DetSetVector< TTStub< T > > >( "StubRejected" );
}

/// Destructor
template< typename T >
TTStubBuilder< T >::~TTStubBuilder(){}

/// Begin run
template< typename T >
void TTStubBuilder< T >::beginRun( const edm::Run& run, const edm::EventSetup& iSetup )
{
  /// Get the stub finding algorithm
  iSetup.get< TTStubAlgorithmRecord >().get( theStubFindingAlgoHandle );
  ievt=0;
  moduleStubs_CIC.clear();
  moduleStubs_MPA.clear();
  moduleStubs_CBC.clear();
}

/// End run
template< typename T >
void TTStubBuilder< T >::endRun( const edm::Run& run, const edm::EventSetup& iSetup ){}

/// Sort routine for stub ordering
template< typename T >
bool TTStubBuilder< T >::SortStubBendPairs( const std::pair< unsigned int, double >& left, const std::pair< unsigned int, double >& right )
{
  return fabs(left.second) < fabs(right.second);
}

/// Analogous sorting routine directly from stubs
template< typename T >
bool TTStubBuilder< T >::SortStubsBend( const TTStub< T >& left, const TTStub< T >& right )
{
  return fabs(left.getTriggerBend()) < fabs(right.getTriggerBend());
}

/// Implement the producer
template< >
void TTStubBuilder< Ref_Phase2TrackerDigi_ >::produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

#endif
