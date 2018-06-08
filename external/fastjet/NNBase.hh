#ifndef __FASTJET_NNBASE_HH__
#define __FASTJET_NNBASE_HH__

//FJSTARTHEADER
// $Id: NNBase.hh 4355 2018-04-22 15:38:54Z salam $
//
// Copyright (c) 2016-2018, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER

#include<fastjet/ClusterSequence.hh>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup advanced_usage
/// \class _NoInfo
/// internal dummy class, used as a default template argument
class _NoInfo {};

/// @ingroup advanced_usage
/// \class NNInfo
///
/// internal helper template class to facilitate initialisation of a
/// BJ with a PseudoJet and extra information. Implementations of
/// NN-based clustering do not need to explicitly use or refer to
/// this class!
template<class I> class NNInfo {
public:
  NNInfo()         : _info(NULL) {}
  NNInfo(I * info) : _info(info) {}
  template<class BJ> void init_jet(BJ * briefjet, const fastjet::PseudoJet & jet, int index) { briefjet->init(jet, index, _info);}
private:
  I * _info;
};

/// @ingroup advanced_usage Internal helper specialisation of NNInfo
/// for cases where there is no extra info
template<> class NNInfo<_NoInfo>  {
public:
  NNInfo()           {}
  NNInfo(_NoInfo * ) {}
  template<class BJ> void init_jet(BJ * briefjet, const fastjet::PseudoJet & jet, int index) { briefjet->init(jet, index);}
};


//----------------------------------------------------------------------
/// @ingroup advanced_usage
/// \class NNBase
/// Helps solve closest pair problems with generic interparticle and
/// particle-beam distances.
///
/// \section Description Description and derived classes:
///
///   This is an abstract base class which defines the interface for
///   several classes that help carry out nearest-neighbour
///   clustering:
///  
///    - NNH          provides an implementation for generic measures,
///  
///    - NNFJN2Plain  provides an implementation for distances
///                   satisfying the FastJet lemma i.e. distances for
///                   which the minimum dij has the property that i is
///                   the geometrical nearest neighbour of j, or vice
///                   versa. I.e. the distance can be factorised in a
///                   momentum factor and a geometric piece. This is
///                   based on the fastjet N2Plain clustering strategy
///  
///    - NNFJN2Tiled  is a tiled version of NNFJN2Plain (based on the
///                   N2Tiled FastJet clustering strategy). Like
///                   NNPlain2 it applies to distance measures that
///                   satisfy the FastJet lemma, with the additional
///                   restriction that: (a) the underlying geometry
///                   should be cylindrical (e.g. rapidity--azimuth)
///                   and (b) the search for the geometric nearest
///                   neighbour of each particle can be limited to
///                   that particle's tile and its neighbouring tiles.
///
/// If you can use NNFJN2Plain it will usually be faster than
/// NNH. NNFJN2Tiled, where it can be used, will be faster for
/// multiplicities above a few tens of particles.
///
///   NOTE: IN ALL CASES, THE DISTANCE MUST BE SYMMETRIC (dij=dji)!!!
///
/// \section BJ Underlying BriefJet (BJ) class:
/// 
///   All derived classes must be templated with a BriefJet (BJ)
///   class --- BJ should basically cache the minimal amount of
///   information that is needed to efficiently calculate
///   interparticle distances and particle-beam distances.
///   
///   This class can be used with or without an extra "Information"
///   template, i.e. `NN*<BJ>` or `NN*<BJ,I>`. Accordingly BJ must provide
///   one of the two following init functions:
///
///   \code
///     void  BJ::init(const PseudoJet & jet);            // initialise with a PseudoJet
///     void  BJ::init(const PseudoJet & jet, I * info);  // initialise with a PseudoJet + info
///   \endcode
///
///   where info might be a pointer to a class that contains, e.g.,
///   information about R, or other parameters of the jet algorithm
///   
///   The BJ then provides information about interparticle and
///   particle-beam distances. The exact requirements depend on
///   whether you use NNH, NNFJN2Plain or NNFJN2Tiled. (See the
///   corresponding classes for details).
///
///
/// \section Workflow Workflow:
///
///   In all cases, the usage of NNBase classes works as follows:
///
///   First, from the list of particles, create an `NN*<BJ>`
///   object of the appropriate type with the appropriate BJ class
///   (and optional extra info).
///
///   Then, cluster using a loop like this (assuming a FastJet plugin)
///
///   \code
///     while (njets > 0) {
///       int i, j, k;
///       // get the i and j that minimize the distance
///       double dij = nn.dij_min(i, j);  
///     
///       // do the appropriate recombination and update the nn
///       if (j >= 0) {    // interparticle recombination
///         cs.plugin_record_ij_recombination(i, j, dij, k);
///         nn.merge_jets(i, j, cs.jets()[k], k); 
///       } else {         // bbeam recombination
///         double diB = cs.jets()[i].E()*cs.jets()[i].E(); // get new diB
///         cs.plugin_record_iB_recombination(i, diB);
///         nn.remove_jet(i);
///       }
///       njets--;
///     }
///   \endcode
///
/// For an example of how the NNH<BJ> class is used, see the
/// JadePlugin or EECambridgePlugin.
template<class I = _NoInfo> class NNBase : public NNInfo<I> {
public:
  /// Default constructor
  NNBase() {}
  /// Constuctor with additional Info 
  NNBase(I * info) : NNInfo<I>(info) {}

  /// initialisation from a given list of particles
  virtual void start(const std::vector<PseudoJet> & jets) = 0;
  
  /// returns the dij_min and indices iA, iB, for the corresponding jets.
  /// If iB < 0 then iA recombines with the beam
  virtual double dij_min(int & iA, int & iB) = 0;

  /// removes the jet pointed to by index iA
  virtual void remove_jet(int iA) = 0;

  /// merges the jets pointed to by indices A and B and replaces them with
  /// jet, assigning it an index jet_index.
  virtual void merge_jets(int iA, int iB, const PseudoJet & jet, int jet_index) =  0;

  virtual ~NNBase() {};
};


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh


#endif // __FASTJET_NNBASE_HH__
