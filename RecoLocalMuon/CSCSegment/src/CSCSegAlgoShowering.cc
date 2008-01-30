/**
 * \file CSCSegAlgoShowering.cc
 *
 *  \author: D. Fortin - UC Riverside
 *
 * See header file for description.
 */

#include "RecoLocalMuon/CSCSegment/src/CSCSegAlgoShowering.h"

#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include <Geometry/CSCGeometry/interface/CSCChamber.h>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>


/* Constructor
 *
 */
CSCSegAlgoShowering::CSCSegAlgoShowering(const edm::ParameterSet& ps) {
  debug                  = ps.getUntrackedParameter<bool>("CSCSegmentDebug");
  dRPhiFineMax           = ps.getParameter<double>("dRPhiFineMax");
  dPhiFineMax            = ps.getParameter<double>("dPhiFineMax");
  tanThetaMax            = ps.getParameter<double>("tanThetaMax");
  tanPhiMax              = ps.getParameter<double>("tanPhiMax");	
  maxRatioResidual       = ps.getParameter<double>("maxRatioResidualPrune");
}


/* Destructor:
 *
 */
CSCSegAlgoShowering::~CSCSegAlgoShowering(){

}


/* showerSeg
 *
 */
CSCSegment CSCSegAlgoShowering::showerSeg( const CSCChamber* aChamber, ChamberHitContainer rechits ) {

  theChamber = aChamber;
  // Initialize parameters
  std::vector<float> x, y, gz;
  std::vector<int> n;
 
  for (int i = 0; i < 6; ++i) {
    x.push_back(0.);
    y.push_back(0.);
    gz.push_back(0.);
    n.push_back(0);
  }

  // Loop over hits to find center-of-mass position in each layer
  for (ChamberHitContainer::const_iterator it = rechits.begin(); it != rechits.end(); it++ ) {
    const CSCRecHit2D& hit = (**it);
    const CSCDetId id = hit.cscDetId();
    int l_id = id.layer();
    const CSCLayer* layer  = theChamber->layer(hit.cscDetId().layer());
    GlobalPoint gp         = layer->toGlobal(hit.localPosition());
    LocalPoint  lp         = theChamber->toLocal(gp);

    n[l_id -1]++;
    x[l_id -1] += lp.x();
    y[l_id -1] += lp.y();
    gz[l_id -1] += gp.z();
  }


  std::vector<LocalPoint> lpCOM;
  // Determine center of mass for each layer and average center of mass for chamber
  float avgChamberX = 0.;
  float avgChamberY = 0.;
  int n_lay = 0;

  for (unsigned i = 0; i < 6; ++i) {
    if (n[i] < 1 ) continue;
 
    x[i] = x[i]/n[i];
    y[i] = y[i]/n[i];
    avgChamberX += x[i];
    avgChamberY += y[i];
    n_lay++;

  }

  if ( n_lay > 0) {
    avgChamberX = avgChamberX / n_lay;
    avgChamberY = avgChamberY / n_lay;
  }
 

  // Find the 2 averages which are furthest away from average
  unsigned worse_lay_1 = 10;
  unsigned worse_lay_2 = 10;
  float worse_deltaR_1 = 0.;
  float worse_deltaR_2 = 0.;

  for (unsigned i = 0; i < 6; ++i) {
    if (n[i] < 1 ) continue;
 
    float deltaR = (avgChamberX-x[i])*(avgChamberX-x[i]) + (avgChamberY-y[i])*(avgChamberY-y[i]);

    if (deltaR > worse_deltaR_1 ) {
      worse_lay_1 = i;
      worse_deltaR_1 = deltaR;
    }
    else if (deltaR > worse_deltaR_2 ) {
      worse_lay_2 = i;
      worse_deltaR_2 = deltaR;
    }

  }


  // Recompute Global average without those 2 worse averages
  float avgChamberX2 = 0.;
  float avgChamberY2 = 0.;
  unsigned n_lay2 = 0;

  for (unsigned i = 0; i < 6; ++i) {
    if (n[i] < 1 ) continue;
 
    if ( n_lay > 3 && i == worse_lay_1 ) continue;
    if ( n_lay > 4 && i == worse_lay_2 ) continue;

    avgChamberX2 += x[i];
    avgChamberY2 += y[i];
    n_lay2++;
  }

  if ( n_lay2 > 0 ) {
    avgChamberX = avgChamberX2 / n_lay2;
    avgChamberY = avgChamberY2 / n_lay2;
  }



  // Store average value for each layer into a local point
  LocalPoint defaultLP(avgChamberX,avgChamberY,0.);

  for (unsigned i = 0; i < 6; ++i) {
    if (n[i] > 0) {

      if ( n_lay > 3 && i == worse_lay_1 ) {
        lpCOM.push_back(defaultLP);
      }
      else if ( n_lay > 4 && i == worse_lay_2 ) {;
        lpCOM.push_back(defaultLP);
      }
      else {
        LocalPoint lpt(x[i],y[i],0.);
        lpCOM.push_back(lpt);
      }  
    }
    else {
      lpCOM.push_back(defaultLP);
    }
  }


  std::vector<float> r_closest;
  std::vector<int> id;
  for (unsigned i = 0; i < 6; ++i ) {
    id.push_back(-1);
    r_closest.push_back(9999.);
  }

  int idx = 0;

  // Loop over all hits and find hit closest to com for that layer.
  for (ChamberHitContainer::const_iterator it = rechits.begin(); it != rechits.end(); it++ ) {    
    const CSCRecHit2D& hit = (**it);
    int layId = hit.cscDetId().layer();
    const CSCLayer* layer  = theChamber->layer(hit.cscDetId().layer());
    GlobalPoint gp         = layer->toGlobal(hit.localPosition());
    LocalPoint  lp         = theChamber->toLocal(gp);

    float d_x = lp.x() - lpCOM[layId-1].x();
    float d_y = lp.y() - lpCOM[layId-1].y();

    LocalPoint diff(d_x, d_y, 0.);
    
    if ( fabs(diff.mag() ) < r_closest[layId-1] ) {
       r_closest[layId-1] =  fabs(diff.mag());
       id[layId-1] = idx;
    }
    idx++;
  }

  // Now fill vector of rechits closest to center of mass:
  protoSegment.clear();
  idx = 0;

  // Loop over all hits and find hit closest to com for that layer.
  for (ChamberHitContainer::const_iterator it = rechits.begin(); it != rechits.end(); it++ ) {    
    const CSCRecHit2D& hit = (**it);
    int layId = hit.cscDetId().layer();

    if ( idx == id[layId-1] )protoSegment.push_back(*it);

    idx++;    
  }

  // Reorder hits in protosegment
  if ( gz[0] > 0. ) {
    if ( gz[0] > gz[5] ) { 
      reverse( protoSegment.begin(), protoSegment.end() );
    }    
  }
  else if ( gz[0] < 0. ) {
    if ( gz[0] < gz[5] ) {
      reverse( protoSegment.begin(), protoSegment.end() );
    }    
  }


  // Compute the segment properties
  updateParameters();

  // Clean up protosegment if there is one very bad hit on segment
  if (protoSegment.size() > 3) pruneFromResidual();

  // Look for better hits near segment  
  for (ChamberHitContainer::const_iterator it = rechits.begin(); it != rechits.end(); it++ ) {

    const CSCRecHit2D* h = *it;
    int layer = h->cscDetId().layer();

    if ( isHitNearSegment( h ) ) compareProtoSegment(h, layer);
  }


  // Prune worse hit if necessary
  if ( protoSegment.size() > 5 ) pruneFromResidual();

  // Update the parameters
  updateParameters();

  // Local direction
  double dz   = 1./sqrt(1. + protoSlope_u*protoSlope_u + protoSlope_v*protoSlope_v);
  double dx   = dz*protoSlope_u;
  double dy   = dz*protoSlope_v;
  LocalVector localDir(dx,dy,dz);
        
  // localDir may need sign flip to ensure it points outward from IP  
  double globalZpos    = ( theChamber->toGlobal( protoIntercept ) ).z();
  double globalZdir    = ( theChamber->toGlobal( localDir ) ).z();
  double directionSign = globalZpos * globalZdir;
  LocalVector protoDirection = (directionSign * localDir).unit();

  // Error matrix
  AlgebraicSymMatrix protoErrors = calculateError();     
        
  CSCSegment temp(protoSegment, protoIntercept, protoDirection, protoErrors, protoChi2); 

  return temp;

} 




/* isHitNearSegment
 *
 * Compare rechit with expected position from proto_segment
 */
bool CSCSegAlgoShowering::isHitNearSegment( const CSCRecHit2D* hit) const {

  const CSCLayer* layer = theChamber->layer(hit->cscDetId().layer());

  // hit phi position in global coordinates
  GlobalPoint Hgp = layer->toGlobal(hit->localPosition());
  double Hphi = Hgp.phi();                                
  if (Hphi < 0.) Hphi += 2.*M_PI;
  LocalPoint Hlp = theChamber->toLocal(Hgp);
  double z = Hlp.z();  

  double LocalX = protoIntercept.x() + protoSlope_u * z;
  double LocalY = protoIntercept.y() + protoSlope_v * z;
  LocalPoint Slp(LocalX, LocalY, z);
  GlobalPoint Sgp = theChamber->toGlobal(Slp); 
  double Sphi = Sgp.phi();
  if (Sphi < 0.) Sphi += 2.*M_PI;
  double R = sqrt(Sgp.x()*Sgp.x() + Sgp.y()*Sgp.y());
  
  double deltaPhi = Sphi - Hphi;
  if (deltaPhi >  2.*M_PI) deltaPhi -= 2.*M_PI;
  if (deltaPhi < -2.*M_PI) deltaPhi += 2.*M_PI;
  if (deltaPhi < 0.) deltaPhi = -deltaPhi; 

  double RdeltaPhi = R * deltaPhi;

  if (RdeltaPhi < dRPhiFineMax && deltaPhi < dPhiFineMax ) return true;

  return false;
}


/* Method addHit
 *
 * Test if can add hit to proto segment. If so, try to add it.
 *
 */
bool CSCSegAlgoShowering::addHit(const CSCRecHit2D* aHit, int layer) {
  
  // Return true if hit was added successfully and then parameters are updated.
  // Return false if there is already a hit on the same layer, or insert failed.
  
  bool ok = true;
  
  // Test that we are not trying to add the same hit again
  for ( ChamberHitContainer::const_iterator it = protoSegment.begin(); it != protoSegment.end(); it++ ) 
    if ( aHit == (*it)  ) return false;
  
  protoSegment.push_back(aHit);

  return ok;
}    


/* Method updateParameters
 *      
 * Perform a simple Least Square Fit on proto segment to determine slope and intercept
 *
 */   
void CSCSegAlgoShowering::updateParameters() {

  // Compute slope from Least Square Fit    
  HepMatrix M(4,4,0);
  HepVector B(4,0);

  ChamberHitContainer::const_iterator ih;
  
  for (ih = protoSegment.begin(); ih != protoSegment.end(); ++ih) {
    
    const CSCRecHit2D& hit = (**ih);
    const CSCLayer* layer  = theChamber->layer(hit.cscDetId().layer());
    GlobalPoint gp         = layer->toGlobal(hit.localPosition());
    LocalPoint  lp         = theChamber->toLocal(gp); 
    
    double u = lp.x();
    double v = lp.y();
    double z = lp.z();
    
    // ptc: Covariance matrix of local errors 
    HepMatrix IC(2,2);
    IC(1,1) = hit.localPositionError().xx();
    IC(1,2) = hit.localPositionError().xy();
    IC(2,2) = hit.localPositionError().yy();
    IC(2,1) = IC(1,2); // since Cov is symmetric
    
    // ptc: Invert covariance matrix (and trap if it fails!)
    int ierr = 0;
    IC.invert(ierr); // inverts in place
    if (ierr != 0) {
      LogDebug("CSC") << "CSCSegment::fitSlopes: failed to invert covariance matrix=\n" << IC << "\n";      
    }
    
    M(1,1) += IC(1,1);
    M(1,2) += IC(1,2);
    M(1,3) += IC(1,1) * z;
    M(1,4) += IC(1,2) * z;
    B(1)   += u * IC(1,1) + v * IC(1,2);
    
    M(2,1) += IC(2,1);
    M(2,2) += IC(2,2);
    M(2,3) += IC(2,1) * z;
    M(2,4) += IC(2,2) * z;
    B(2)   += u * IC(2,1) + v * IC(2,2);
    
    M(3,1) += IC(1,1) * z;
    M(3,2) += IC(1,2) * z;
    M(3,3) += IC(1,1) * z * z;
    M(3,4) += IC(1,2) * z * z;
    B(3)   += ( u * IC(1,1) + v * IC(1,2) ) * z;
    
    M(4,1) += IC(2,1) * z;
    M(4,2) += IC(2,2) * z;
    M(4,3) += IC(2,1) * z * z;
    M(4,4) += IC(2,2) * z * z;
    B(4)   += ( u * IC(2,1) + v * IC(2,2) ) * z;
  }
  
  HepVector p = solve(M, B);
  
  // Update member variables 
  // Note that origin has local z = 0

  protoIntercept = LocalPoint(p(1), p(2), 0.);
  protoSlope_u = p(3);
  protoSlope_v = p(4);

  // Determine Chi^2 for the proto wire segment
  
  double chsq = 0.;
  
  for (ih = protoSegment.begin(); ih != protoSegment.end(); ++ih) {
    
    const CSCRecHit2D& hit = (**ih);
    const CSCLayer* layer  = theChamber->layer(hit.cscDetId().layer());
    GlobalPoint gp         = layer->toGlobal(hit.localPosition());
    LocalPoint lp          = theChamber->toLocal(gp);
    
    double u = lp.x();
    double v = lp.y();
    double z = lp.z();
    
    double du = protoIntercept.x() + protoSlope_u * z - u;
    double dv = protoIntercept.y() + protoSlope_v * z - v;
    
    HepMatrix IC(2,2);
    IC(1,1) = hit.localPositionError().xx();
    IC(1,2) = hit.localPositionError().xy();
    IC(2,2) = hit.localPositionError().yy();
    IC(2,1) = IC(1,2);
    
    // Invert covariance matrix
    int ierr = 0;
    IC.invert(ierr);
    if (ierr != 0) {
      LogDebug("CSC") << "CSCSegment::fillChiSquared: failed to invert covariance matrix=\n" << IC << "\n";      
    }
    chsq += du*du*IC(1,1) + 2.*du*dv*IC(1,2) + dv*dv*IC(2,2);
  }
  protoChi2 = chsq;
}



/* Method compareProtoSegment
 *      
 * For hit coming from the same layer of an existing hit within the proto segment
 * test if achieve better chi^2 by using this hit than the other
 *
 */ 
void CSCSegAlgoShowering::compareProtoSegment(const CSCRecHit2D* h, int layer) {
  
  // Store old segment first
  double old_protoChi2                  = protoChi2;
  LocalPoint old_protoIntercept         = protoIntercept;
  float old_protoSlope_u                = protoSlope_u;
  float old_protoSlope_v                = protoSlope_v;
  LocalVector old_protoDirection        = protoDirection;
  ChamberHitContainer old_protoSegment  = protoSegment;
 

  // Try adding the hit to existing segment, and remove old one existing in same layer
  ChamberHitContainer::iterator it;
  for ( it = protoSegment.begin(); it != protoSegment.end(); ) {
    if ( (*it)->cscDetId().layer() == layer ) {
      it = protoSegment.erase(it);
    } else {
      ++it;
    }
  }
  bool ok = addHit(h, layer);

  if (ok) updateParameters();
  
  if ( (protoChi2 > old_protoChi2) || ( !ok ) ) {
    protoChi2       = old_protoChi2;
    protoIntercept  = old_protoIntercept;
    protoSlope_u    = old_protoSlope_u;
    protoSlope_v    = old_protoSlope_v;
    protoDirection  = old_protoDirection;
    protoSegment    = old_protoSegment;
  }
}


/* calculateError
 *
 */
AlgebraicSymMatrix CSCSegAlgoShowering::calculateError() const {

  // Blightly assume the following never fails
 
  std::vector<const CSCRecHit2D*>::const_iterator it;
  int nhits = protoSegment.size();
  int ierr; 

  AlgebraicSymMatrix weights(2*nhits, 0);
  AlgebraicMatrix A(2*nhits, 4);

  int row = 0;  
  for (it = protoSegment.begin(); it != protoSegment.end(); ++it) {
    const CSCRecHit2D& hit = (**it);
    const CSCLayer* layer = theChamber->layer(hit.cscDetId().layer());
    GlobalPoint gp = layer->toGlobal(hit.localPosition());      
    LocalPoint lp = theChamber->toLocal(gp); 
    float z = lp.z();
    ++row;
    weights(row, row)   = hit.localPositionError().xx();
    weights(row, row+1) = hit.localPositionError().xy();
    A(row, 1) = 1.;
    A(row, 3) = z;
    ++row;
    weights(row, row-1) = hit.localPositionError().xy();
    weights(row, row)   = hit.localPositionError().yy();
    A(row, 2) = 1.;
    A(row, 4) = z;
  }
  weights.invert(ierr);

  AlgebraicSymMatrix a = weights.similarityT(A);
  a.invert(ierr);
    
  // but reorder components to match what's required by TrackingRecHit interface 
  // i.e. slopes first, then positions 
    
  AlgebraicSymMatrix hold( a ); 
    
  // errors on slopes into upper left 
  a(1,1) = hold(3,3); 
  a(1,2) = hold(3,4); 
  a(2,1) = hold(4,3); 
  a(2,2) = hold(4,4); 
    
  // errors on positions into lower right 
  a(3,3) = hold(1,1); 
  a(3,4) = hold(1,2); 
  a(4,3) = hold(2,1); 
  a(4,4) = hold(2,2); 
    
  // off-diagonal elements remain unchanged 
  return a;    
} 



// Try to clean up segments by quickly looking at residuals
void CSCSegAlgoShowering::pruneFromResidual(){

  // Only prune if have at least 5 hits 
  if ( protoSegment.size() < 5 ) return ;


  // Now Study residuals
      
  float maxResidual = 0.;
  float sumResidual = 0.;
  int nHits = 0;
  int badIndex = -1;
  int j = 0;


  ChamberHitContainer::const_iterator ih;

  for ( ih = protoSegment.begin(); ih != protoSegment.end(); ++ih ) {
    const CSCRecHit2D& hit = (**ih);
    const CSCLayer* layer  = theChamber->layer(hit.cscDetId().layer());
    GlobalPoint gp         = layer->toGlobal(hit.localPosition());
    LocalPoint lp          = theChamber->toLocal(gp);

    double u = lp.x();
    double v = lp.y();
    double z = lp.z();

    double du = protoIntercept.x() + protoSlope_u * z - u;
    double dv = protoIntercept.y() + protoSlope_v * z - v;

    float residual = sqrt(du*du + dv*dv);

    sumResidual += residual;
    nHits++;
    if ( residual > maxResidual ) {
      maxResidual = residual;
      badIndex = j;
    }
    j++;
  }

  float corrAvgResidual = (sumResidual - maxResidual)/(nHits -1);

  // Keep all hits 
  if ( maxResidual/corrAvgResidual < maxRatioResidual ) return;


  // Drop worse hit and recompute segment properties + fill

  ChamberHitContainer newProtoSegment;

  j = 0;
  for ( ih = protoSegment.begin(); ih != protoSegment.end(); ++ih ) {
    if ( j != badIndex ) newProtoSegment.push_back(*ih);
    j++;
  }
  
  protoSegment.clear();

  for ( ih = newProtoSegment.begin(); ih != newProtoSegment.end(); ++ih ) {
    protoSegment.push_back(*ih);
  }

  // Update segment parameters
  updateParameters();

}



