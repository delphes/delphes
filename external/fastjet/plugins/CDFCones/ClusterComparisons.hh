#ifndef _CLUSTER_COMPARISONS_HH_
#define _CLUSTER_COMPARISONS_HH_

//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from
// http://www.pa.msu.edu/~huston/Les_Houches_2005/JetClu+Midpoint-StandAlone.tgz
//
// Permission to distribute it with FastJet has been granted by Joey
// Huston (see the COPYING file in the main FastJet directory for
// details).
// Changes from the original file are listed below.
//----------------------------------------------------------------------

// History of changes compared to the original ClusterComparison.hh file
// 
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * put the code in the fastjet::cdf namespace
// 
// 2007-03-10  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added support for the pttilde scale choice in the CDF midpoint code
// 
// 2007-02-21  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added option of choosing the scale used in the split-merge
//          procedure (pt [default], Et or mt)
// 
// 2006-09-24  Gavin Salam  <salam@lpthe.jussieu.fr>
// 
//        * added JetClu+MidPoint to FastJet

#include "Cluster.hh"

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace cdf{

class ClusterFourVectorEtGreater
{
 public:
  int operator()(const Cluster& c1, const Cluster& c2) const
  {
    return c1.fourVector.Et() > c2.fourVector.Et();
  }
};

class ClusterCentroidEtGreater
{
 public:
  int operator()(const Cluster& c1, const Cluster& c2) const
  {
    return c1.centroid.Et > c2.centroid.Et;
  }
};

class ClusterPtGreater
{
 public:
  int operator()(const Cluster& c1, const Cluster& c2) const
  {
    return c1.fourVector.pt() > c2.fourVector.pt();
  }
};

class ClusterMtGreater
{
 public:
  int operator()(const Cluster& c1, const Cluster& c2) const
  {
    return c1.fourVector.mt() > c2.fourVector.mt();
  }
};

class ClusterPtTildeGreater
{
 public:
  int operator()(const Cluster& c1, const Cluster& c2) const
  {
    return c1.pt_tilde > c2.pt_tilde;
  }
};

}  // namespace cdf

FASTJET_END_NAMESPACE

#endif
