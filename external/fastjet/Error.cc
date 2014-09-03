//FJSTARTHEADER
// $Id: Error.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/Error.hh"
#include "fastjet/config.h"
#include <sstream>

// printing the stack would need execinfo
#ifdef FASTJET_HAVE_EXECINFO_H
#include <execinfo.h>
#include <cstdlib>
#endif

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

bool Error::_print_errors = true;
bool Error::_print_backtrace = false;
ostream * Error::_default_ostr = & cerr;

Error::Error(const std::string & message_in) {
  _message = message_in; 
  if (_print_errors && _default_ostr){
    ostringstream oss;
    oss << "fastjet::Error:  "<< message_in << endl;

    // only print the stack if execinfo is available and stack enabled
#ifdef FASTJET_HAVE_EXECINFO_H
    if (_print_backtrace){
      void * array[10];
      char ** messages;
 
      int size = backtrace(array, 10);
      messages = backtrace_symbols(array, size);
      
      oss << "stack:" << endl;
      for (int i = 1; i < size && messages != NULL; ++i){
	oss << "  #" << i << ": " << messages[i] << endl;
      }
      free(messages);
    }
#endif

    *_default_ostr << oss.str();
    // get something written to file even 
    // if the program aborts
    _default_ostr->flush(); 

    // // output error message either to cerr or to the user-set stream
    // if (_default_ostr) { *_default_ostr << oss.str();
    //                       // get something written to file even 
    // 			  // if the program aborts
    //                       _default_ostr->flush(); }
    // else               { std::cerr << oss.str(); }
    
  }
}

FASTJET_END_NAMESPACE

