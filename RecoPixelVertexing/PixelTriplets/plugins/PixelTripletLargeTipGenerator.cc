#include "RecoPixelVertexing/PixelTriplets/plugins/PixelTripletLargeTipGenerator.h"
#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPair.h"

#include "RecoPixelVertexing/PixelTriplets/interface/ThirdHitPredictionFromCircle.h"
#include "ThirdHitRZPrediction.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoUtilities.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/ThirdHitCorrection.h"
#include "RecoTracker/TkHitPairs/interface/RecHitsSortedInPhi.h"

#include "MatchedHitRZCorrectionFromBending.h"
//#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerAlgo.h"
//#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerTools.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/KDTreeLinkerAlgo.h" //amend to point at your copy...
#include "RecoPixelVertexing/PixelTriplets/plugins/KDTreeLinkerTools.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>

#include "DataFormats/Math/interface/normalizedPhi.h"

#include "CommonTools/Utils/interface/DynArray.h"

using namespace std;

using Range=PixelRecoRange<float>;
using HelixRZ=ThirdHitPredictionFromCircle::HelixRZ;

namespace {
  struct LayerRZPredictions {
    ThirdHitRZPrediction<PixelRecoLineRZ> line;
    ThirdHitRZPrediction<HelixRZ> helix1, helix2;
    MatchedHitRZCorrectionFromBending rzPositionFixup;
    ThirdHitCorrection correction;
  };
}

constexpr double nSigmaRZ = 3.4641016151377544; // sqrt(12.)
constexpr double nSigmaPhi = 3.;
constexpr float fnSigmaRZ = nSigmaRZ;


PixelTripletLargeTipGenerator::PixelTripletLargeTipGenerator(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC)
  : HitTripletGeneratorFromPairAndLayers(cfg),
    useFixedPreFiltering(cfg.getParameter<bool>("useFixedPreFiltering")),
    extraHitRZtolerance(cfg.getParameter<double>("extraHitRZtolerance")),
    extraHitRPhitolerance(cfg.getParameter<double>("extraHitRPhitolerance")),
    useMScat(cfg.getParameter<bool>("useMultScattering")),
    useBend(cfg.getParameter<bool>("useBending")),
    dphi(useFixedPreFiltering ? cfg.getParameter<double>("phiPreFiltering") : 0) {}

PixelTripletLargeTipGenerator::~PixelTripletLargeTipGenerator() {}

namespace {
  inline
  bool intersect(Range &range, const Range &second)
  {
    if (range.first > second.max() || range.second < second.min())
      return false;
    if (range.first < second.min())
      range.first = second.min();
    if (range.second > second.max())
      range.second = second.max();
    return range.first < range.second;
  }
}

void PixelTripletLargeTipGenerator::hitTriplets(const TrackingRegion& region, 
						OrderedHitTriplets & result,
						const edm::Event & ev,
						const edm::EventSetup& es,
						SeedingLayerSetsHits::SeedingLayerSet pairLayers,
						const std::vector<SeedingLayerSetsHits::SeedingLayer>& thirdLayers)
{ 
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);

  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopoHand;
  es.get<TrackerTopologyRcd>().get(tTopoHand);
  const TrackerTopology *tTopo=tTopoHand.product();

  auto const & doublets = thePairGenerator->doublets(region,ev,es, pairLayers);
  
  if (doublets.empty()) return;
   
  auto outSeq =  doublets.detLayer(HitDoublets::outer)->seqNum();


  int size = thirdLayers.size();


  using NodeInfo = KDTreeNodeInfo<unsigned int>;
  std::vector<NodeInfo > layerTree; // re-used throughout
  std::vector<unsigned int> foundNodes; // re-used throughout
  foundNodes.reserve(100);

  declareDynArray(KDTreeLinkerAlgo<unsigned int>, size, hitTree);
  declareDynArray(LayerRZPredictions, size, mapPred);

  float rzError[size]; //save maximum errors

  const float maxDelphi = region.ptMin() < 0.3f ? float(M_PI)/4.f : float(M_PI)/8.f; // FIXME move to config??
  const float maxphi = M_PI+maxDelphi, minphi = -maxphi; // increase to cater for any range
  const float safePhi = M_PI-maxDelphi; // sideband


  const RecHitsSortedInPhi * thirdHitMap[size];

  for(int il = 0; il < size; il++) {
    thirdHitMap[il] = &(*theLayerCache)(thirdLayers[il], region, ev, es);
    auto const & hits = *thirdHitMap[il];
 
    const DetLayer *layer = thirdLayers[il].detLayer();
    LayerRZPredictions &predRZ = mapPred[il];
    predRZ.line.initLayer(layer);
    predRZ.helix1.initLayer(layer);
    predRZ.helix2.initLayer(layer);
    predRZ.line.initTolerance(extraHitRZtolerance);
    predRZ.helix1.initTolerance(extraHitRZtolerance);
    predRZ.helix2.initTolerance(extraHitRZtolerance);
    predRZ.rzPositionFixup = MatchedHitRZCorrectionFromBending(layer,tTopo);
    predRZ.correction.init(es, region.ptMin(), *doublets.detLayer(HitDoublets::inner), *doublets.detLayer(HitDoublets::outer), *thirdLayers[il].detLayer(), useMScat, false);


    layerTree.clear();
    float minv=999999.0; float maxv = -999999.0; // Initialise to extreme values in case no hits
    float maxErr=0.0f;
    for (unsigned int i=0; i!=hits.size(); ++i) {
      auto angle = hits.phi(i);
      auto v =  hits.gv(i);
      //use (phi,r) for endcaps rather than (phi,z)
      minv = std::min(minv,v);  maxv = std::max(maxv,v);
      float myerr = hits.dv[i];
      maxErr = std::max(maxErr,myerr);
      layerTree.emplace_back(i, angle, v); // save it
      // populate side-bands
      if (angle>safePhi) layerTree.emplace_back(i, angle-Geom::ftwoPi(), v);
      else if (angle<-safePhi) layerTree.emplace_back(i, angle+Geom::ftwoPi(), v);
    }
    KDTreeBox phiZ(minphi, maxphi, minv-0.01f, maxv+0.01f);  // declare our bounds
    //add fudge factors in case only one hit and also for floating-point inaccuracy
    hitTree[il].build(layerTree, phiZ); // make KDtree
    rzError[il] = maxErr; //save error
  }

  double curv = PixelRecoUtilities::curvature(1. / region.ptMin(), es);
  
  for (std::size_t ip =0;  ip!=doublets.size(); ip++) {
    auto xi = doublets.x(ip,HitDoublets::inner);
    auto yi = doublets.y(ip,HitDoublets::inner);
    auto zi = doublets.z(ip,HitDoublets::inner);
    // auto rvi = doublets.rv(ip,HitDoublets::inner);
    auto xo = doublets.x(ip,HitDoublets::outer);
    auto yo = doublets.y(ip,HitDoublets::outer);
    auto zo = doublets.z(ip,HitDoublets::outer);
    // auto rvo = doublets.rv(ip,HitDoublets::outer);
    GlobalPoint gp1(xi,yi,zi);
    GlobalPoint gp2(xo,yo,zo);

    auto toPos = std::signbit(zo-zi); 

    PixelRecoLineRZ line(gp1, gp2);
    PixelRecoPointRZ point2(gp2.perp(), zo);
    ThirdHitPredictionFromCircle predictionRPhi(gp1, gp2, extraHitRPhitolerance);

    Range generalCurvature = predictionRPhi.curvature(region.originRBound());
    if (!intersect(generalCurvature, Range(-curv, curv))) continue;

    for(int il = 0; il < size; il++) {
      if (hitTree[il].empty()) continue; // Don't bother if no hits
      const DetLayer *layer = thirdLayers[il].detLayer();
      bool barrelLayer = layer->isBarrel();

      if ( (!barrelLayer) & (toPos != std::signbit(layer->position().z())) ) continue;

      
      Range curvature = generalCurvature;

      LayerRZPredictions &predRZ = mapPred[il];
      predRZ.line.initPropagator(&line);

      auto & correction = predRZ.correction;
      correction.init(line, point2,  outSeq);


      Range rzRange;
      if (useBend) {
        // For the barrel region:
        // swiping the helix passing through the two points across from
        // negative to positive bending, can give us a sort of U-shaped
        // projection onto the phi-z (barrel) or r-z plane (forward)
        // so we checking minimum/maximum of all three possible extrema
        // 
        // For the endcap region:
        // Checking minimum/maximum radius of the helix projection
        // onto an endcap plane, here we have to guard against
        // looping tracks, when phi(delta z) gets out of control.
        // HelixRZ::rAtZ should not follow looping tracks, but clamp
        // to the minimum reachable r with the next-best lower |curvature|.
        // So same procedure as for the barrel region can be applied.
        //
        // In order to avoid looking for potential looping tracks at all
        // we also clamp the allowed curvature range for this layer,
        // and potentially fail the layer entirely
	
        if (!barrelLayer) {
          Range z3s = predRZ.line.detRange();
          double z3 = z3s.first < 0 ? std::max(z3s.first, z3s.second)
	    : std::min(z3s.first, z3s.second);
          double maxCurvature = HelixRZ::maxCurvature(&predictionRPhi,
                                                      gp1.z(), gp2.z(), z3);
          if (!intersect(curvature, Range(-maxCurvature, maxCurvature)))
            continue;
        }
	
        HelixRZ helix1(&predictionRPhi, gp1.z(), gp2.z(), curvature.first);
        HelixRZ helix2(&predictionRPhi, gp1.z(), gp2.z(), curvature.second);
	
        predRZ.helix1.initPropagator(&helix1);
        predRZ.helix2.initPropagator(&helix2);
	
        Range rzRanges[2] = { predRZ.helix1(), predRZ.helix2() };
        predRZ.helix1.initPropagator(nullptr);
        predRZ.helix2.initPropagator(nullptr);

        rzRange.first = std::min(rzRanges[0].first, rzRanges[1].first);
        rzRange.second = std::max(rzRanges[0].second, rzRanges[1].second);
	
        // if the allowed curvatures include a straight line,
        // this can give us another extremum for allowed r/z
        if (curvature.first * curvature.second < 0.0) {
          Range rzLineRange = predRZ.line();
          rzRange.first = std::min(rzRange.first, rzLineRange.first);
          rzRange.second = std::max(rzRange.second, rzLineRange.second);
        }
      } else {
        rzRange = predRZ.line();
      }

      if (rzRange.first >= rzRange.second)
        continue;

      correction.correctRZRange(rzRange);

      Range phiRange;
      if (useFixedPreFiltering) {
	float phi0 = doublets.phi(ip,HitDoublets::outer);
        phiRange = Range(phi0 - dphi, phi0 + dphi);
      } else {
        Range radius;
	
        if (barrelLayer) {
          radius = predRZ.line.detRange();
          if (!intersect(rzRange, predRZ.line.detSize()))
            continue;
        } else {
          radius = rzRange;
          if (!intersect(radius, predRZ.line.detSize()))
            continue;
        }
	
        auto rPhi1 = predictionRPhi(curvature, radius.first);
        bool ok1 = !rPhi1.empty();
        if (ok1) {
          correction.correctRPhiRange(rPhi1);
          rPhi1.first  /= radius.max();
          rPhi1.second /= radius.max();
        }
        auto rPhi2 = predictionRPhi(curvature, radius.second);
        bool ok2 = !rPhi2.empty();
        if (ok2) {
          correction.correctRPhiRange(rPhi2);
          rPhi2.first  /= radius.min();
          rPhi2.second /= radius.min();
        }

        if (ok1) {
          rPhi1.first = normalizedPhi(rPhi1.first);
          rPhi1.second = proxim(rPhi1.second,rPhi1.first);
          if(ok2) {
            rPhi2.first = proxim(rPhi2.first,rPhi1.first);
            rPhi2.second = proxim(rPhi2.second,rPhi1.first);
            phiRange = rPhi1.sum(rPhi2);
          } else phiRange=rPhi1;
        } else if(ok2) {
          rPhi2.first = normalizedPhi(rPhi2.first);
          rPhi2.second = proxim(rPhi2.second,rPhi2.first);
          phiRange=rPhi2;
        } else continue;

        // if (std::abs(phiRange.first)>float(M_PI)) std::cout << "bha1 " << phiRange.first << ' ' << phiRange.second << std::endl;
        if (std::abs(phiRange.first)>float(M_PI) && std::abs(phiRange.second)>float(M_PI)) std::cout << "bha2 " << phiRange.first << ' ' << phiRange.second << std::endl;

      }
      
      foundNodes.clear(); // Now recover hits in bounding box...
      float prmin=phiRange.min(), prmax=phiRange.max(); //get contiguous range

      if (prmax<prmin)  std::cout << "aarg " << phiRange.first << ' ' << phiRange.second << std::endl;
      if (prmax-prmin>maxDelphi) std::cout << "delphi " << ' ' << prmin << '/' << prmax << std::endl;


      if (barrelLayer) {
	Range regMax = predRZ.line.detRange();
	Range regMin = predRZ.line(regMax.min());
	regMax = predRZ.line(regMax.max());
	correction.correctRZRange(regMin);
	correction.correctRZRange(regMax);
	if (regMax.min() < regMin.min()) { swap(regMax, regMin);}
	KDTreeBox phiZ(prmin, prmax,
		       regMin.min()-fnSigmaRZ*rzError[il],
		       regMax.max()+fnSigmaRZ*rzError[il]);
	hitTree[il].search(phiZ, foundNodes);
      }
      else {
	KDTreeBox phiZ(prmin, prmax,
		       rzRange.min()-fnSigmaRZ*rzError[il],
		       rzRange.max()+fnSigmaRZ*rzError[il]);
	hitTree[il].search(phiZ, foundNodes);
      }
      
      MatchedHitRZCorrectionFromBending l2rzFixup(doublets.hit(ip,HitDoublets::outer)->det()->geographicalId(), tTopo);
      MatchedHitRZCorrectionFromBending l3rzFixup = predRZ.rzPositionFixup;

      thirdHitMap[il] = &(*theLayerCache)(thirdLayers[il], region, ev, es);
      auto const & hits = *thirdHitMap[il];
      for (auto KDdata : foundNodes) {
	GlobalPoint p3 = hits.gp(KDdata); 
	double p3_r = p3.perp();
	double p3_z = p3.z();
	float p3_phi =  hits.phi(KDdata); 

	Range rangeRPhi = predictionRPhi(curvature, p3_r);
	correction.correctRPhiRange(rangeRPhi);

	float ir = 1.f/p3_r;
	float phiErr = nSigmaPhi *  hits.drphi[KDdata]*ir;
	if (!checkPhiInRange(p3_phi, rangeRPhi.first*ir-phiErr, rangeRPhi.second*ir+phiErr))
	  continue;
	
	Basic2DVector<double> thc(p3.x(), p3.y());
	
	auto curv_ = predictionRPhi.curvature(thc);
	double p2_r = point2.r(); double p2_z = point2.z(); // they will be modified!
	
	l2rzFixup(predictionRPhi, curv_, *doublets.hit(ip,HitDoublets::outer), p2_r, p2_z, tTopo);
	l3rzFixup(predictionRPhi, curv_, *hits.theHits[KDdata].hit(), p3_r, p3_z, tTopo);
	
	Range rangeRZ;
	if (useBend) {
	  HelixRZ updatedHelix(&predictionRPhi, gp1.z(), p2_z, curv_);
	  rangeRZ = predRZ.helix1(barrelLayer ? p3_r : p3_z, updatedHelix);
	} else {
	  float tIP = predictionRPhi.transverseIP(thc);
	  PixelRecoPointRZ updatedPoint2(p2_r, p2_z);
	  PixelRecoLineRZ updatedLine(line.origin(), point2, tIP);
	  rangeRZ = predRZ.line(barrelLayer ? p3_r : p3_z, line);
	}
	correction.correctRZRange(rangeRZ);
	
	double err = nSigmaRZ * hits.dv[KDdata];
	
	rangeRZ.first -= err, rangeRZ.second += err;
	
	if (!rangeRZ.inside(barrelLayer ? p3_z : p3_r)) continue;

	if (theMaxElement!=0 && result.size() >= theMaxElement) {
	  result.clear();
	  edm::LogError("TooManyTriplets")<<" number of triples exceed maximum. no triplets produced.";
	  return;
	}
	result.emplace_back( doublets.hit(ip,HitDoublets::inner), doublets.hit(ip,HitDoublets::outer), hits.theHits[KDdata].hit()); 
      }
    }
  }
  // std::cout << "found triplets " << result.size() << std::endl;
}

