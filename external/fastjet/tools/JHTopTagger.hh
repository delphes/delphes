#ifndef __FASTJET_JH_TOP_TAGGER_HH__
#define __FASTJET_JH_TOP_TAGGER_HH__

//FJSTARTHEADER
// $Id: JHTopTagger.hh 3433 2014-07-23 08:17:03Z salam $
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


#include <fastjet/tools/TopTaggerBase.hh>
#include <fastjet/CompositeJetStructure.hh>
#include <fastjet/LimitedWarning.hh>

FASTJET_BEGIN_NAMESPACE

class JHTopTagger;
class JHTopTaggerStructure;

//----------------------------------------------------------------------
/// @ingroup tools_taggers
/// \class JHTopTagger
/// Class that helps perform boosted top tagging using the "Johns Hopkins"
/// method from arXiv:0806.0848 (Kaplan, Rehermann, Schwartz
/// and Tweedie)
///
///The tagger proceeds as follows:
///  - start from a jet J obtained with the Cambridge/Aachen algorithm
///  - undo the last iteration j -> j_1,j_2 (with pt_1>pt_2) until the
///    two subjets satisfy pt_1 > delta_p pt_J (with pt_J the pt of
///    the original jet) and |y_1 - y_2| + |phi_1 - phi_2| > delta_r.
///  - if one of these criteria is not satisfied, carry on the
///    procedure with j_1 (discarding j_2)
///  - for each of the subjets found, repeat the procedure. If some
///    new substructure is found, keep these 2 new subjets, otherwise
///    keep the original subjet (found during the first iteration)
///  - at this stage, one has at most 4 subjets. If one has less than
///    3, the tagger has failed.
///  - reconstruct the W from the 2 subjets with a mass closest to the
///    W mass
///  - impose that the W helicity angle be less than a threshold
///    cos_theta_W_max.
///
/// \section input Input conditions
/// 
///  - the original jet must have an associated (and valid)
///    ClusterSequence
///  - the tagger is designed to work with jets formed by the
///    Cambridge/Aachen (C/A) algorithm; if a non-C/A jet is passed to
///    the tagger, a warning will be issued
///
/// \section Example
///
/// A  JHTopTagger can be used as follows:
///
/// \code
///    double delta_p = 0.10; // subjets must carry at least this fraction of the original jet's p_t
///    double delta_r = 0.19; // subjets must be separated by at least this Manhattan distance
///    double cos_theta_W_max = 0.7; // the maximal allowed value of the W helicity angle
///    JHTopTagger top_tagger(delta_p, delta_r, cos_theta_W_max);
///    // indicate the acceptable range of top, W masses (default: no limits)
///    top_tagger.set_top_selector(SelectorMassRange(150,200));
///    top_tagger.set_W_selector  (SelectorMassRange( 65, 95));
///    // now try and tag a jet
///    PseudoJet top_candidate = top_tagger(jet);  // jet should come from a Cambridge/Aachen clustering
///    if (top_candidate != 0) { // successful tagging
///      double top_mass = top_candidate.m();
///      double W_mass   = top_candidate.structure_of<JHTopTagger>().W().m();
///    }
/// \endcode
///
/// The full set of information available from the structure_of<JHTopTagger>() 
/// call is
///
/// - PseudoJet W()    : the W subjet of the top candidate
/// - PseudoJet non_W(): non-W subjet(s) of the top candidate (i.e. the b)
/// - double cos_theta_W(): the W helicity angle
/// - PseudoJet W1(): the harder of the two prongs of the W
/// - PseudoJet W2(): the softer of the two prongs of the W
///
/// The structure of the top_candidate can also be accessed through its
/// pieces() function:
///
/// - top_candidate.pieces()[0]: W
/// - top_candidate.pieces()[1]: non_W
///
/// The W itself has two pieces (corresponding to W1, W2). 
///
/// The existence of the first two of the structural calls (W(),
/// non_W()) and the fact that the top is made of two pieces (W,
/// non_W) are features that should be common to all taggers derived
/// from TopTaggerBase.
///
/// See also \subpage Example13 for a full usage example.
///
class JHTopTagger : public TopTaggerBase {
public:
  /// default ctor
  /// The parameters are the following:
  ///  \param delta_p          fractional pt cut imposed on the subjets
  ///                          (computed as a fraction of the original jet)
  ///  \param delta_r          minimal distance between 2 subjets
  ///                          (computed as |y1-y2|+|phi1-phi2|)
  ///  \param cos_theta_W_max  the maximal value for the polarisation 
  ///                          angle of the W
  ///  \param mW               the W mass
  ///
  /// The default values of all these parameters are taken from
  /// arXiv:0806:0848
  JHTopTagger(const double delta_p=0.10, const double delta_r=0.19, 
              double cos_theta_W_max=0.7, double mW=80.4)
    : _delta_p(delta_p), _delta_r(delta_r),
      _cos_theta_W_max(cos_theta_W_max), _mW(mW){};

  /// returns a textual description of the tagger
  virtual std::string description() const;

  /// runs the tagger on the given jet and
  /// returns the tagged PseudoJet if successful, or a PseudoJet==0 otherwise
  /// (standard access is through operator()).
  ///  \param jet   the PseudoJet to tag
  virtual PseudoJet result(const PseudoJet & jet) const;

  // the type of the associated structure
  typedef JHTopTaggerStructure StructureType;

protected:
  /// runs the Johns Hopkins decomposition procedure
  std::vector<PseudoJet> _split_once(const PseudoJet & jet_to_split,
                                     const PseudoJet & reference_jet) const;

  double _delta_p, _delta_r, _cos_theta_W_max, _mW;
  static LimitedWarning _warnings_nonca;
};


//------------------------------------------------------------------------
/// @ingroup tools_taggers
/// \class JHTopTaggerStructure
/// the structure returned by the JHTopTagger transformer.
///
/// See the JHTopTagger class description for the details of what
/// is inside this structure
///
class JHTopTaggerStructure : public CompositeJetStructure, public TopTaggerBaseStructure {
public:
  /// ctor with pieces initialisation
  JHTopTaggerStructure(std::vector<PseudoJet> pieces_in,
                 const JetDefinition::Recombiner *recombiner = 0) :
    CompositeJetStructure(pieces_in, recombiner), _cos_theta_w(0.0){}

  /// returns the W subjet
  inline const PseudoJet & W() const{ 
    return _pieces[0];
  }

  /// returns the first W subjet (the harder)
  inline PseudoJet W1() const{
    assert(W().pieces().size()>0);
    return W().pieces()[0];
  }
  
  /// returns the second W subjet
  inline PseudoJet W2() const{
    assert(W().pieces().size()>1);
    return W().pieces()[1];
  }

  /// returns the non-W subjet
  /// It will have 1 or 2 pieces depending on whether the tagger has
  /// found 3 or 4 pieces
  inline const PseudoJet & non_W() const{ 
    return _pieces[1];
  }

  /// returns the W helicity angle
  inline double cos_theta_W() const {return _cos_theta_w;}

//  /// returns the original jet (before tagging)
//  const PseudoJet & original() const {return _original_jet;}


protected:
  double _cos_theta_w;      ///< the W helicity angle
  //PseudoJet _W;             ///< the tagged W
  //PseudoJet _non_W;         ///< the remaining pieces
//  PseudoJet _original_jet;  ///< the original jet (before tagging)

  // allow the tagger to set these
  friend class JHTopTagger;
};



FASTJET_END_NAMESPACE

#endif  //  __FASTJET_JH_TOP_TAGGER_HH__

