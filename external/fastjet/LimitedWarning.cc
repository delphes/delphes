//FJSTARTHEADER
// $Id$
//
// Copyright (c) 2005-2025, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
atomic<ostream *> LimitedWarning::_default_ostr{&cerr};
atomic<mutex *> LimitedWarning::_stream_mutex{nullptr};
atomic<int> LimitedWarning::_max_warn_default{5};
std::mutex LimitedWarning::_global_warnings_summary_mutex;
#else
ostream * LimitedWarning::_default_ostr = &cerr;
int LimitedWarning::_max_warn_default = 5;
#endif // FASTJET_HAVE_LIMITED_THREAD_SAFETY

std::list< LimitedWarning::Summary > LimitedWarning::_global_warnings_summary;

// /// output a warning to ostr
// void LimitedWarning::warn(const std::string & warning) {
//   warn(warning, _default_ostr);
// }

/// the number of times so far that a warning has been registered
/// with this instance of the class.
int LimitedWarning::n_warn_so_far() const{
  // explicitly cast to the pointer type (useless wo thread-safety
  // features but works around an issue with the intel compiler
  // (v13.1.3) with thread-safety features
  if (((LimitedWarning::Summary *)_this_warning_summary) == 0) return 0;
  return (*_this_warning_summary).second;
}


void LimitedWarning::warn(const char * warning, std::ostream * ostr) {
  // update the summary
  if (((LimitedWarning::Summary *)_this_warning_summary) == 0){
#ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
    // Threadsafety note: 
    //   here we need to lock _this_warning_summary to be sure that we
    //   can quietly initialising things without taking the risk that
    //   another thread "regenerates" the Summary entry at the same time
    // 
    // first acquire the mutex! 
    // See e.g. http://baptiste-wicht.com/posts/2012/03/cp11-concurrency-tutorial-part-2-protect-shared-data.html
    // for a quick intro.
    std::lock_guard<std::mutex> guard(_global_warnings_summary_mutex);
    // then make sure we still need to create things (in case another
    // thread got us beaten (it's better to use the mutex as little as
    // possible, hence the repetition of the text which might
    // otherwise look stupid)
    if (((LimitedWarning::Summary *)_this_warning_summary) == 0){
      // prepare the information for the summary
      _global_warnings_summary.push_back(Summary(warning, 0));
      _this_warning_summary = & (_global_warnings_summary.back());
    }
    // the lock will automatically be released here
#else  
    // prepare the information for the summary
    _global_warnings_summary.push_back(Summary(warning, 0));
    _this_warning_summary = & (_global_warnings_summary.back());
#endif // FASTJET_HAVE_LIMITED_THREAD_SAFETY
  }


  // maintain the count, but do not allow overflow
  unsigned int count = (*_this_warning_summary).second.step();

  // print the warning if we have not done it enough already
  if ((_max_warn<0) || (count < (unsigned int)_max_warn)) {
    // prepare the warning within a string stream
    ostringstream warnstr;
    warnstr << "WARNING from FastJet: ";
    warnstr << warning;
    if ((_max_warn>0) && (count+1 == (unsigned int)_max_warn))
      warnstr << " (LAST SUCH WARNING)";
    warnstr << std::endl;
    // arrange for the whole warning to be output in one go (that way
    // user can easily insert their own printout, e.g. event number
    // before the warning string).
    if (ostr) {
      // if there is a mutex, use it to lock
#ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
      if (_stream_mutex){
        std::lock_guard<std::mutex> guard(*_stream_mutex);
        (*ostr) << warnstr.str();
        ostr->flush(); // get something written to file even if the program aborts
      } else 
#endif // FASTJET_HAVE_LIMITED_THREAD_SAFETY      
      {
        (*ostr) << warnstr.str();
        ostr->flush(); // get something written to file even if the program aborts
      }
    }
  }

}

//----------------------------------------------------------------------
string LimitedWarning::summary() {
  ostringstream str;
#ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
  {
    // is the lock here really necessary. The only potential issue I
    // see is if another thread adds a warning when the loop below
    // calls it++ on the previous last element. Is there a simpler way
    // to handle this?
    std::lock_guard<std::mutex> guard(_global_warnings_summary_mutex);
#endif
  for (list<Summary>::const_iterator it = _global_warnings_summary.begin();
       it != _global_warnings_summary.end(); it++) {
    str << it->second << " times: " << it->first << endl;
  }
#ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
  }
#endif
  return str.str();
}

FASTJET_END_NAMESPACE
