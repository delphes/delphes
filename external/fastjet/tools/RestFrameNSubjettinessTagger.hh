#ifndef __FASTJET_RESTFRAMENSUBJETTINESS_TAGGER_HH__
#define __FASTJET_RESTFRAMENSUBJETTINESS_TAGGER_HH__

//FJSTARTHEADER
// $Id: RestFrameNSubjettinessTagger.hh 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2005-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/CompositeJetStructure.hh>
#include <fastjet/tools/Transformer.hh>

FASTJET_BEGIN_NAMESPACE

class RestFrameNSubjettinessTagger;
class RestFrameNSubjettinessTaggerStructure;

//----------------------------------------------------------------------
/// @ingroup tools_taggers
/// \class RestFrameNSubjettinessTagger
/// Class that helps perform 2-pronged boosted tagging using
/// a reclustering in the jet's rest frame, supplemented with a cut on N-subjettiness
/// (and a decay angle), as discussed by Ji-Hun Kim in arXiv:1011.1493.
///
/// To tag a fat jet, the tagger proceeds as follows:
///
///  - boost its constituents into the rest frame of the jet
///
///  - recluster them using another jet definition (the original
///    choice was SISCone in spherical coordinates with R=0.6 and
///    f=0.75.
///
///  - keep the 2 most energetic subjets (\f$q_{1,2}\f$) and compute
///    the 2-subjettiness
///    \f[ 
///      \tau_2^j = \frac{2}{m_{\rm jet}^2}\,
///                 \sum_{k\in {\rm jet}} {\rm min}(q_1.p_k,q_2.p_k)
///    \f]
///    where the sum runs over the constituents of the jet. 
///
///  - require \f$\tau_2^j < \tau_2^{\rm cut}\f$ [0.08 by default]
///
///  - impose that (in the rest frame of the fat jet), the angles
///    between the 2 most energetic subjets and the boost axis are
///    both large enough: \f$\cos(\theta_s)<c_\theta^{\rm cut}\f$ 
///    [0.8 by default]
///
/// Note that in the original version, the jets to be tagged were reconstructed
/// using SISCone with R=0.8 and f=0.75. Also, b-tagging was imposed
/// on the 2 subjets found in the rest-frame tagging procedure.
///
/// \section desc Options
/// 
/// The constructor has the following arguments:
///  - The first argument is the jet definition to be used to
///    recluster the constituents of the jet to be filtered (in the
///    rest frame of the tagged jet).
///  - The second argument is the cut on tau_2 [0.08 by default]
///  - The 3rd argument is the cut on cos(theta_s) [0.8 by default]
///  - If the 4th argument is true, 2 exclusive rest-frame jets will
///    be considered in place of the 2 most energetic inclusive jets
///
/// \section input Input conditions
/// 
///  - the original jet must have constituents
///
/// \section output Output/structure
/// 
///  - the 2 subjets are kept as pieces if some substructure is found,
///    otherwise a single 0-momentum piece
///  - the tau2 and maximal cos(theta_s) values computed during the
///    tagging can be obtained via the resulting jet's structure_of<...>() 
///    function
///
class RestFrameNSubjettinessTagger : public Transformer{
public:
  /// ctor with arguments (see the class description above)
  RestFrameNSubjettinessTagger(const JetDefinition subjet_def, 
                      const double tau2cut=0.08, 
                      const double costhetascut=0.8,
                      const bool use_exclusive = false)
    : _subjet_def(subjet_def), _t2cut(tau2cut), _costscut(costhetascut),
      _use_exclusive(use_exclusive){};

  /// returns a textual description of the tagger 
  virtual std::string description() const;

  /// runs the tagger on the given jet and
  /// returns the tagged PseudoJet if successful, a PseudoJet==0 otherwise
  /// (standard access is through operator()).
  virtual PseudoJet result(const PseudoJet & jet) const;

  /// the type of Structure returned
  typedef RestFrameNSubjettinessTaggerStructure StructureType;

protected:
  JetDefinition _subjet_def;
  double _t2cut, _costscut;
  bool _use_exclusive;
};


//------------------------------------------------------------------------
/// @ingroup tools_taggers
/// \class RestFrameNSubjettinessTaggerStructure
/// the structure returned by the RestFrameNSubjettinessTagger transformer.
///
/// See the RestFrameNSubjettinessTagger class description for the details of
/// what is inside this structure
///
class RestFrameNSubjettinessTaggerStructure : public CompositeJetStructure{
public:
  /// ctor with pieces initialisation
  RestFrameNSubjettinessTaggerStructure(const std::vector<PseudoJet> & pieces_in) :
    CompositeJetStructure(pieces_in), _tau2(0.0), _costhetas(1.0){}

  /// returns the associated N-subjettiness
  inline double tau2() const{return _tau2;}

  /// returns the associated angle with the boosted axis
  inline double costhetas() const {return _costhetas;}

//  /// returns the original jet (before tagging)
//  const PseudoJet & original() const {return _original_jet;}

protected:
  double _tau2;      ///< the value of the N-subjettiness
  double _costhetas; ///< the minimal angle between the dijets
                     ///< and the boost axis
//  PseudoJet _original_jet;  ///< the original jet (before tagging)

  // allow the tagger to set these
  friend class RestFrameNSubjettinessTagger;
};

FASTJET_END_NAMESPACE
#endif  //  __FASTJET_RESTFRAMENSUBJETTINESS_TAGGER_HH__

