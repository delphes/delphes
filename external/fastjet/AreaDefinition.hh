//FJSTARTHEADER
// $Id: AreaDefinition.hh 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2006-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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


#ifndef __FASTJET_AREADEFINITION_HH__
#define __FASTJET_AREADEFINITION_HH__

#include "fastjet/GhostedAreaSpec.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
//
/// @ingroup area_classes
/// \class VoronoiAreaSpec
/// Specification for the computation of the Voronoi jet area
///
/// class for holding a "Voronoi area" specification; an area will be
/// assigned to each particle, which is the area of the intersection
/// of the particle's Voronoi cell with a circle of radius
/// R*effective_Rfact.
///
class VoronoiAreaSpec {
public:

  /// default constructor (effective_Rfact = 1);
  VoronoiAreaSpec() : _effective_Rfact(1.0) {};
  
  /// constructor that allows you to set effective_Rfact.
  VoronoiAreaSpec(double effective_Rfact_in) : 
    _effective_Rfact(effective_Rfact_in) {};

  /// return the value of effective_Rfact
  double effective_Rfact() const {return _effective_Rfact;}

  /// return a textual description of the area definition.
  std::string description() const;

private:
  double _effective_Rfact;
};


/// the different types of area that are supported
enum AreaType {invalid_area = -1, 
               active_area = 0, active_area_explicit_ghosts = 1,
               one_ghost_passive_area = 10, passive_area = 11, 
               voronoi_area=20};


//----------------------------------------------------------------------
/// @ingroup area_classes
/// \class AreaDefinition
/// class that holds a generic area definition
class AreaDefinition {
public:
  
  /// default constructor, which provides a ghosted active area, with
  /// sensible defaults for the ghosts.
  AreaDefinition() {
    _area_type  = active_area;
    _ghost_spec = GhostedAreaSpec();
  }

  /// constructor for an area definition based on an area type and a
  /// ghosted area specification
  AreaDefinition(AreaType type, const GhostedAreaSpec & spec) {
    _ghost_spec = spec;
    _area_type   = type;
    assert(type != voronoi_area);
  }

  /// constructor for an area definition based on an area type and a
  /// voronoi area specification (type must be voronoi_area)
  AreaDefinition(AreaType type, const VoronoiAreaSpec & spec) {
    _voronoi_spec = spec;
    _area_type   = type;
    assert(type == voronoi_area);
  }

  /// constructor for an area definition based on an area type and 
  /// which attempts to provide sensible defaults for everything else
  AreaDefinition(AreaType type) {
    _area_type   = type;
    if (type == voronoi_area) {
      _voronoi_spec = VoronoiAreaSpec();
    } else {
      _ghost_spec = GhostedAreaSpec();
    }
  }

  /// constructor for an area definition based on an ghosted area
  /// specification, and an option to select which ghosted area you want
  AreaDefinition(const GhostedAreaSpec & spec, AreaType type = active_area) {
    _ghost_spec = spec;
    _area_type   = type;
    assert(type != voronoi_area);
  }

  /// constructor for an area definition based on a voronoi area
  /// specification
  AreaDefinition(const VoronoiAreaSpec & spec) {
    _voronoi_spec = spec;
    _area_type    = voronoi_area;
  }

  /// return a description of the current area definition
  std::string description() const;

  /// return info about the type of area being used by this defn
  AreaType area_type() const {return _area_type;}

  /// return a reference to the active area spec
  const GhostedAreaSpec  & ghost_spec()  const {return _ghost_spec;}
  GhostedAreaSpec & ghost_spec()  {return _ghost_spec;}

  /// return a reference to the voronoi area spec
  const VoronoiAreaSpec & voronoi_spec() const {return _voronoi_spec;}
  
private:

  AreaType        _area_type;
  GhostedAreaSpec  _ghost_spec;
  VoronoiAreaSpec _voronoi_spec;
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh


#endif // __FASTJET_AREADEFINITION_HH__
