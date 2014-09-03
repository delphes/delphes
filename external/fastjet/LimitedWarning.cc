//FJSTARTHEADER
// $Id: LimitedWarning.cc 3619 2014-08-13 14:17:19Z salam $
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

#include "fastjet/LimitedWarning.hh"
#include <sstream>
#include <limits>

using namespace std;

FASTJET_BEGIN_NAMESPACE

ostream * LimitedWarning::_default_ostr = &cerr;
std::list< LimitedWarning::Summary > LimitedWarning::_global_warnings_summary;
int LimitedWarning::_max_warn_default = 5;


// /// output a warning to ostr
// void LimitedWarning::warn(const std::string & warning) {
//   warn(warning, _default_ostr);
// }

void LimitedWarning::warn(const char * warning, std::ostream * ostr) {
  if (_this_warning_summary == 0) {
    // prepare the information for the summary
    _global_warnings_summary.push_back(Summary(warning, 0));
    _this_warning_summary = & (_global_warnings_summary.back());
  }
  if (_n_warn_so_far < _max_warn) {
    // prepare the warning within a string stream
    ostringstream warnstr;
    warnstr << "WARNING from FastJet: ";
    warnstr << warning;
    _n_warn_so_far++;
    if (_n_warn_so_far == _max_warn) warnstr << " (LAST SUCH WARNING)";
    warnstr << std::endl;
    // arrange for the whole warning to be output in one go (that way
    // user can easily insert their own printout, e.g. event number
    // before the warning string).
    if (ostr) {
      (*ostr) << warnstr.str();
      ostr->flush(); // get something written to file even if the program aborts
    }
  }

  // maintain the count, but do not allow overflow
  if (_this_warning_summary->second < numeric_limits<unsigned>::max()) {
    _this_warning_summary->second++;
  }
}

//----------------------------------------------------------------------
string LimitedWarning::summary() {
  ostringstream str;
  for (list<Summary>::const_iterator it = _global_warnings_summary.begin();
       it != _global_warnings_summary.end(); it++) {
    str << it->second << " times: " << it->first << endl;
  }
  return str.str();
}

FASTJET_END_NAMESPACE
