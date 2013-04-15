#ifndef D0RunIIconeJets_CONEJETINFO_HPP
#define D0RunIIconeJets_CONEJETINFO_HPP

// --------------------------------------------------------------------------
// ConeJetInfo.hpp
// Purpose: Hold informations about the cone jets that do not fit into
//   a CalTClusterChunk/IntEclusterChunk.
//
// Created: Laurent Duflot 31-JUL-2000
//
// Modified:
//  09-Aug-2000 Laurent Duflot
//   + add initial jet ET (i.e. before split/merge) 
//    1-May-2007 Lars Sonnenschein
//    extracted from D0 software framework and modified to remove subsequent dependencies
//
//
// This file is distributed with FastJet under the terms of the GNU
// General Public License (v2). Permission to do so has been granted
// by Lars Sonnenschein and the D0 collaboration (see COPYING for
// details)
//
// History of Changes in FastJet compared tothe original version of
// ConeJetInfo.hpp
//
// 2011-12-13  Gregory Soyez  <soyez@fastjet.fr>
// 
//        * added license information
//
// 2011-11-14  Gregory Soyez  <soyez@fastjet.fr>
//
//         * changed the name of a few parameters to avoid a gcc
//           -Wshadow warning
//
// 2009-01-17  Gregory Soyez  <soyez@fastjet.fr>
//
//        * put the code in the fastjet::d0 namespace
//
// --------------------------------------------------------------------------


//#define CONEJET_SPLITMERGE_MOD 100

#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace d0{

namespace D0RunIIconeJets_CONEJETINFO {

const int CONEJET_SPLITMERGE_MOD = 100;

class ConeJetInfo
{
public:
  ConeJetInfo(): _seedET(0.), _initial_jet_ET(0.), _nb_split_merge(0) {};
  ConeJetInfo( float seedET_in): _seedET(seedET_in),  _nb_split_merge(0) {}; 
  ConeJetInfo( float seedET_in, float initialET_in, int nb_split, int nb_merge): 
    _seedET(seedET_in), _initial_jet_ET(initialET_in), 
    _nb_split_merge(nb_merge + CONEJET_SPLITMERGE_MOD*nb_split) {};
  ~ConeJetInfo() {};

  float seedET() const {return _seedET;};
  float initialET() const { return _initial_jet_ET; };
  int nbSplit() const {return _nb_split_merge/CONEJET_SPLITMERGE_MOD;};
  int nbMerge() const {return _nb_split_merge%CONEJET_SPLITMERGE_MOD;};
  int SplitMergeWord() const {return _nb_split_merge;};

  void initialET(float ET) { _initial_jet_ET = ET;};
  void splitted() { _nb_split_merge += CONEJET_SPLITMERGE_MOD;};
  void merged() { _nb_split_merge += 1;};


private:
  float _seedET;
  float _initial_jet_ET;  // stable cone ET before split/merge
  int _nb_split_merge;
  
};

}

}  // namespace d0

FASTJET_END_NAMESPACE

#endif

