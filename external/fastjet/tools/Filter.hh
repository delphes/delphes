#ifndef __FASTJET_TOOLS_FILTER_HH__
#define __FASTJET_TOOLS_FILTER_HH__

//STARTHEADER
// $Id: Filter.hh 2694 2011-11-14 22:27:51Z salam $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include <fastjet/CompositeJetStructure.hh> // to derive the FilterStructure from CompositeJetStructure
#include <fastjet/tools/Transformer.hh>     // to derive Filter from Transformer
#include <iostream>
#include <string>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// fwd declarations
class Filter;
class FilterStructure;

//----------------------------------------------------------------------
/// @ingroup tools_generic
/// \class Filter
/// Class that helps perform filtering (Butterworth, Davison, Rubin
/// and Salam, arXiv:0802.2470) and trimming (Krohn, Thaler and Wang,
/// arXiv:0912.1342) on jets, optionally in conjunction with
/// subtraction (Cacciari and Salam, arXiv:0707.1378).
///
/// For example, to apply filtering that reclusters a jet's
/// constituents with the Cambridge/Aachen jet algorithm with R=0.3
/// and then selects the 3 hardest subjets, one can use the following
/// code:
/// \code
///    Filter filter(JetDefinition(cambridge_algorithm, 0.3), SelectorNHardest(3));
///    PseudoJet filtered_jet = filter(original_jet);
/// \endcode
///
/// To obtain trimming, involving for example the selection of all
/// subjets carrying at least 3% of the original jet's pt, the
/// selector would be replaced by SelectorPtFractionMin(0.03).
///
/// To additionally perform subtraction on the subjets prior to
/// selection, either include a 3rd argument specifying the background
/// density rho, or call the set_subtractor(...) member function.  If
/// subtraction is requested, the original jet must be the result of a
/// clustering with active area with explicit ghosts support or a
/// merging of such pieces.
///
/// The information on the subjets that were kept and rejected can be
/// obtained using:
/// \code
///    vector<PseudoJet> kept_subjets = filtered_jet.pieces();
///    vector<PseudoJet> rejected_subjets = filtered_jet.structure_of<Filter>().rejected();
/// \endcode
///
/// \section impl Implementation Note
/// 
/// If the original jet was defined with the Cambridge/Aachen
/// algorithm (or is made of pieces each of which comes from the C/A
/// alg) and the filtering definition is C/A, then the filter does not
/// rerun the C/A algorithm on the constituents, but instead makes use
/// of the existent C/A cluster sequence in the original jet. This
/// increases the speed of the filter. 
///
/// See also \subpage Example11 for a further usage example.
///
/// Support for areas, reuse of C/A cluster sequences, etc.,
/// considerably complicates the implementation of Filter. For an
/// explanation of how a simpler filter might be coded, see the
/// "User-defined transformers" appendix of the manual.
class Filter : public Transformer{
public:
  /// trivial ctor
  /// Note: this is just for derived classes
  ///       a Filter initialised through this constructor will not work!
  Filter() : _Rfiltfunc(0){};

  /// define a filter that decomposes a jet into subjets using a
  /// generic JetDefinition and then keeps only a subset of these
  /// subjets according to a Selector. Optionally, each subjet may be
  /// internally bakground-subtracted prior to selection.
  ///
  ///  \param subjet_def   the jet definition applied to obtain the subjets
  ///  \param selector     the Selector applied to compute the kept subjets
  ///  \param rho          if non-zero, backgruond-subtract each subjet befor selection
  ///
  /// Note: internal subtraction only applies on jets that are
  /// obtained with a cluster sequence with area support and explicit
  /// ghosts
  Filter(JetDefinition subjet_def, Selector selector, double rho = 0.0) : 
    _subjet_def(subjet_def), _Rfiltfunc(0), _Rfilt(-1), _selector(selector), _rho(rho), _subtractor(0) {}

  /// Same as the full constructor (see above) but just specifying the radius
  /// By default, Cambridge-Aachen is used
  /// If the jet (or all its pieces) is obtained with a non-default
  /// recombiner, that one will be used
  ///  \param Rfilt   the filtering radius
  Filter(double Rfilt, Selector selector, double rho = 0.0) : 
    _Rfiltfunc(0), _Rfilt(Rfilt), _selector(selector), _rho(rho), _subtractor(0) { 
    if (_Rfilt<0)
      throw Error("Attempt to create a Filter with a negative filtering radius");
  }

  /// Same as the full constructor (see above) but just specifying a
  /// filtering radius that will depend on the jet being filtered
  /// As for the previous case, Cambridge-Aachen is used
  /// If the jet (or all its pieces) is obtained with a non-default
  /// recombiner, that one will be used
  ///  \param Rfilt_func   the filtering radius function of a PseudoJet
  Filter(FunctionOfPseudoJet<double> *Rfilt_func, Selector selector, double rho = 0.0) : 
    _Rfiltfunc(Rfilt_func), _Rfilt(-1), _selector(selector), _rho(rho), _subtractor(0) {}

  /// default dtor
  virtual ~Filter(){};

  /// Set a subtractor that is applied to all individual subjets before
  /// deciding which ones to keep. It takes precedence over a non-zero rho.
  void set_subtractor(const Transformer * subtractor) {_subtractor = subtractor;}

  /// runs the filtering and sets kept and rejected to be the jets of interest
  /// (with non-zero rho, they will have been subtracted).
  ///
  /// \param jet    the jet that gets filtered
  /// \return the filtered jet
  virtual PseudoJet result(const PseudoJet & jet) const;

  /// class description
  virtual std::string description() const;

  // the type of the associated structure
  typedef FilterStructure StructureType;

private:
  /// Sets filtered_elements to be all the subjets on which filtering will work.
  /// It also sets the subjet_def to be used in joining things (the bit of
  /// subjet def that is of interest for later is the recombiner).
  void _set_filtered_elements(const PseudoJet & jet,
                              std::vector<PseudoJet> & filtered_elements,
                              JetDefinition & subjet_def,
                              bool & discard_area) const;
  
  /// set the filtered elements in the simple case of C/A+C/A
  void _set_filtered_elements_cafilt(const PseudoJet & jet,
                                     std::vector<PseudoJet> & filtered_elements,
                                     double Rfilt) const;

  /// set the filtered elements in the generic re-clustering case
  void _set_filtered_elements_generic(const PseudoJet & jet, 
                                      std::vector<PseudoJet> & filtered_elements,
                                      const JetDefinition & subjet_def,
				      bool do_areas) const;

  /// gather the information about what is kept and rejected under the
  /// form of a PseudoJet with a special ClusterSequenceInfo
  PseudoJet _finalise(const PseudoJet & jet, 
                      std::vector<PseudoJet> & kept, 
                      std::vector<PseudoJet> & rejected,
                      const JetDefinition & subjet_def,
                      const bool discard_area) const;

  // a series of checks
  //--------------------------------------------------------------------
  /// get the pieces down to the fundamental pieces
  bool _get_all_pieces(const PseudoJet &jet, std::vector<PseudoJet> &all_pieces) const;

  /// get the common recombiner to all pieces (NULL if none)
  const JetDefinition::Recombiner* _get_common_recombiner(const std::vector<PseudoJet> &all_pieces) const;

  /// check if one can apply the simplified trick for C/A subjets
  bool _check_ca(const std::vector<PseudoJet> &all_pieces) const;

  /// check if the jet (or all its pieces) have explicit ghosts
  /// (assuming the jet has area support
  ///
  /// Note that if the jet has an associated cluster sequence that is no
  /// longer valid, an error will be thrown
  bool _check_explicit_ghosts(const std::vector<PseudoJet> &all_pieces) const;

  bool _uses_subtraction() const {return (_subtractor || _rho != 0);}

  JetDefinition _subjet_def;   ///< the jet definition to use to extract the subjets
  FunctionOfPseudoJet<double> *_Rfiltfunc; 
                               ///< a dynamic filtering radius function of the jet being filtered
  double _Rfilt;               ///< a constant specifying the subjet radius (with C/A)
  Selector _selector;  ///< the subjet selection criterium
  double _rho;                 ///< the background density (used for subtraction when possible)
  const Transformer * _subtractor; ///< for subtracting bkgd density from subjets
};



//----------------------------------------------------------------------
/// @ingroup tools_generic
/// \class FilterStructure
/// Class to contain structure information for a filtered jet.
class FilterStructure : public CompositeJetStructure {
public:
  /// constructor from an original ClusterSequenceInfo
  /// We just share the original ClusterSequenceWrapper and initialise
  /// the rest
  FilterStructure(const std::vector<PseudoJet> & pieces_in, 
                  const JetDefinition::Recombiner *rec = 0)
    : CompositeJetStructure(pieces_in, rec){}

  /// virtual dtor to allow further overloading  
  virtual ~FilterStructure(){}

  /// description
  virtual std::string description() const { return "Filtered PseudoJet"; }

  //------------------------------------------------------------------
  /// @name The filter-specific information
  //------------------------------------------------------------------

//  /// returns the original jet (the first of the original jets
//  /// if you filtered a collection of jets)
//  const PseudoJet & original() const {return _original_jet;}

  /// returns the subjets that were not kept during the filtering procedure
  /// (subtracted if the filter requests it, and valid in the original cs)
  const std::vector<PseudoJet> & rejected() const {return _rejected;}

  friend class Filter;  // allow the filter to change the protected/private members

protected:
//  PseudoJet _original_jet;           ///< the original jet
  std::vector<PseudoJet> _rejected;  ///< the subjets rejected by the filter
};


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif   // __FASTJET_TOOLS_FILTER_HH__
