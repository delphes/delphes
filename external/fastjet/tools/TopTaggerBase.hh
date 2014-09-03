#ifndef __FASTJET_TOP_TAGGER_BASE_HH__
#define __FASTJET_TOP_TAGGER_BASE_HH__

//FJSTARTHEADER
// $Id: TopTaggerBase.hh 3433 2014-07-23 08:17:03Z salam $
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

#include <fastjet/internal/base.hh>
#include <fastjet/tools/Transformer.hh>

FASTJET_BEGIN_NAMESPACE

class TopTaggerBase;
class TopTaggerBaseStructure;

//----------------------------------------------------------------------
/// @ingroup tools_taggers
/// \class TopTaggerBase
/// A base class that provides a common interface for top taggers
/// that are able to return a W (in addition to the top itself).
///
/// Top taggers that derive from this should satisfy the following 
/// criteria:
///
/// - their underlying structure should derive from TopTaggerBaseStructure
/// - tagged tops should have two pieces, the first of which is the W candidate
/// - they should apply the top and W selectors to decide if the top has been
///   tagged
class TopTaggerBase : public Transformer {
public:
  TopTaggerBase() : _top_selector(SelectorIdentity()),
                    _W_selector(SelectorIdentity()),
                    _top_selector_set(false),
                    _W_selector_set(false)              {}

  /// the type of the associated structure
  typedef TopTaggerBaseStructure StructureType;

  /// sets the selector that is applied to the top candidate
  void set_top_selector(const Selector & sel) {_top_selector = sel; _top_selector_set = true;}
  /// sets  the selector that is applied to the W candidate
  void set_W_selector  (const Selector & sel) {_W_selector   = sel; _W_selector_set = true;}
  
  /// returns a description of the top and W selectors
  virtual std::string description_of_selectors() const {
    std::string descr;
    if (_top_selector_set) descr = ", top selector: "+_top_selector.description();
    if (_W_selector_set) descr += ", W selector: "+_W_selector.description();
    return descr;
  }

protected:
  /// computes the W helicity angle
  double _cos_theta_W(const PseudoJet & result) const;

  Selector _top_selector, _W_selector;
  bool _top_selector_set, _W_selector_set;
};


//----------------------------------------------------------------------
/// @ingroup tools_taggers
/// \class TopTaggerBaseStructure
/// class that specifies the structure common to all top taggers
///
/// Note that this specifies only the W, non_W part of the
/// interface. An actual top tagger structure class will also need to
/// derive from a PseudoJetStructureBase type class
/// (e.g. CompositeJetStructure)
class TopTaggerBaseStructure {
public:
  virtual const PseudoJet & W() const = 0;
  virtual const PseudoJet & non_W() const = 0;
  virtual ~TopTaggerBaseStructure() {}
};


FASTJET_END_NAMESPACE

#endif  //  __FASTJET_TOP_TAGGER_BASE_HH__

