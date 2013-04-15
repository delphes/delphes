#ifndef __FASTJET_LIMITEDWARNING_HH__
#define __FASTJET_LIMITEDWARNING_HH__

//STARTHEADER
// $Id: LimitedWarning.hh 2577 2011-09-13 15:11:38Z salam $
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


#include "fastjet/internal/base.hh"
#include <iostream>
#include <string>
#include <list>

FASTJET_BEGIN_NAMESPACE

/// @ingroup error_handling
/// \class LimitedWarning
/// class to provide facilities for giving warnings up to some maximum
/// number of times and to provide global summaries of warnings that have
/// been issued.
class LimitedWarning {
public:
  
  /// constructor that provides a default maximum number of warnings
  LimitedWarning() : _max_warn(_max_warn_default), _n_warn_so_far(0), _this_warning_summary(0) {}

  /// constructor that provides a user-set max number of warnings
  LimitedWarning(int max_warn) : _max_warn(max_warn), _n_warn_so_far(0), _this_warning_summary(0) {}

  /// outputs a warning to standard error (or the user's default
  /// warning stream if set)
  void warn(const std::string & warning);

  /// outputs a warning to the specified stream
  void warn(const std::string & warning, std::ostream * ostr);

  /// sets the default output stream for all warnings (by default
  /// cerr; passing a null pointer prevents warnings from being output)
  static void set_default_stream(std::ostream * ostr) {
    _default_ostr = ostr;
  }

  /// sets the default maximum number of warnings of a given kind
  /// before warning messages are silenced.
  static void set_default_max_warn(int max_warn) {
    _max_warn_default = max_warn;
  }

  /// returns a summary of all the warnings that came through the
  /// LimiteWarning class
  static std::string summary();

private:
  int _max_warn, _n_warn_so_far;
  static int _max_warn_default;
  static std::ostream * _default_ostr;
  typedef std::pair<std::string, unsigned int> Summary;
  static std::list< Summary > _global_warnings_summary;
  Summary * _this_warning_summary;
  
};

FASTJET_END_NAMESPACE

#endif // __FASTJET_LIMITEDWARNING_HH__
