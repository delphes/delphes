// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: defines.h                                                           //
// Description: header file for generic parameters definitions               //
// This file is part of the SISCone project.                                 //
// For more details, see http://projects.hepforge.org/siscone                //
//                                                                           //
// Copyright (c) 2006 Gavin Salam and Gregory Soyez                          //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA //
//                                                                           //
// $Revision:: 225                                                          $//
// $Date:: 2008-05-20 16:59:47 +0200 (Tue, 20 May 2008)                     $//
///////////////////////////////////////////////////////////////////////////////

//! \file defines.h
#ifndef __DEFINES_H__
#define __DEFINES_H__

/// program name
// we get "SISCone" by calling
//  siscone::siscone_package_name
// defined in siscone.h
// Otherwise, config.h
// It is also defined as "PACKAGE_NAME" in config.h but this method 
// might lead to conflicts
//#define PROGRAM   PACKAGE_NAME 

// program version
// we get it from
//   siscone::siscone_version
// defined in siscone.h
// It is also defined as "VERSION" in config.h but this method 
// might lead to conflicts

/// perform final stability check using the quadtree
/// With the following define enabled, the final check for stability
/// is performed using the quadtree rather than the explicit list of 
/// particles (see Cstable_cone::proceed_with_stability())
//#define USE_QUADTREE_FOR_STABILITY_TEST


/// threshold for recomoutation of the cone (see Cstable_cones::update_cone())
/// When traversing cone candidates along the angular ordering,
/// the momentum of the protojet candidate is computed incrementally
/// from the particles that enter and leave the cone.
/// When the cumulative change in "|px|+|py|" exceeds the cone "|px|+|py|"
/// we explicitely recompute the cone contents
#define PT_TSHOLD 1000.0


/// The following parameter controls collinear safety. For the set of 
/// particles used in the search of stable cones, we gather particles
/// if their distance in eta and phi is smaller than EPSILON_COLLINEAR.
///
/// NB: for things to behave sensibly one requires 
///        1e-15 << EPSILON_COCIRCULAR << EPSILON_COLLINEAR << 1
///
/// among the scales that appear in practice (e.g. in deciding to use
/// special strategies), we have EPSILON_COLLINEAR, EPSILON_COCIRCULAR,
/// sqrt(EPSILON_COCIRCULAR) and EPSILON_COLLINEAR / EPSILON_COCIRCULAR.
///
#define EPSILON_COLLINEAR 1e-8


/// The following parameter controls cocircular situations.
/// When consecutive particles in the ordered vicinity list are separated
/// (in angle) by less that that limit, we consider that we face a situation
/// of cocircularity.
#define EPSILON_COCIRCULAR 1e-12


/// The following define enables you to allow for identical protocones
/// to be merged automatically after each split-merge step before
/// anything else happens. Whether this happens depends on teh value
/// of the merge_identical_protocones flag in Csplit_merge.
///
/// It we allow such a merging and define allow 
/// MERGE_IDENTICAL_PROTOCONES_DEFAULT_TRUE then the
/// 'merge_identical_protocones' flag in Csplit_merge to be set to
/// 'true'. It may be manually reset to false in which case the
/// merging of identical protocones (protojets) will be turned off.
///
/// Note that this merging identical protocones makes the algorithm
/// infrared-unsafe, so it should be disabled except for testing
/// purposes.
//#define ALLOW_MERGE_IDENTICAL_PROTOCONES
//#define MERGE_IDENTICAL_PROTOCONES_DEFAULT_TRUE


/// if EPSILON_SPLITMERGE is defined then, during the split-merge
/// step, when two jets are found with PTs that are identical to
/// within a relative difference of EPSILON_SPLITMERGE they are
/// compared element-by-element to see where the differences are, and
/// one then uses pt1^2-pt2^2 = (pt1-pt2).(pt1+pt2) as an estimator of
/// which is harder. NB: in unfortunate cases, this can take the split
/// merge step up to N n * ln N time, though on normal events there
/// don't seem to have been any major problems yet.
#define EPSILON_SPLITMERGE 1e-12

/// definition of 2*M_PI which is useful a bit everyhere!
const double twopi = 6.283185307179586476925286766559005768394;

/// debugging information
//#define DEBUG_STABLE_CONES   ///< debug messages in stable cones search
//#define DEBUG_SPLIT_MERGE    ///< debug messages in split-merge
//#define DEBUG                ///< all debug messages !

// in case all debug massages allowed, allow them in practice !
#ifdef DEBUG
#define DEBUG_STABLE_CONES
#define DEBUG_SPLIT_MERGE
#endif

#endif  //  __DEFINES_H__

