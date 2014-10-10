//FJSTARTHEADER
// $Id: Error.cc 3695 2014-09-18 13:57:56Z cacciari $
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

#ifndef __FJCORE__
// printing the stack would need execinfo
#ifdef FASTJET_HAVE_EXECINFO_H
#include <execinfo.h>
#include <cstdlib>
#ifdef FASTJET_HAVE_DEMANGLING_SUPPORT
#include <cstdio>
#include <cxxabi.h>
#endif // FASTJET_HAVE_DEMANGLING_SUPPORT
#endif // FASTJET_HAVE_EXECINFO_H
#endif  // __FJCORE__

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

bool Error::_print_errors = true;
bool Error::_print_backtrace = false;
ostream * Error::_default_ostr = & cerr;
#if (!defined(FASTJET_HAVE_EXECINFO_H)) || defined(__FJCORE__)
  LimitedWarning Error::_execinfo_undefined;
#endif

//----------------------------------------------------------------------
#ifndef __FJCORE__ 

// demangling only is included, i.e. --enable-demangling is specified
// at configure time, execinfo.h is present and the GNU C++ ABI is
// supported
#if defined(FASTJET_HAVE_EXECINFO_H) && defined(FASTJET_HAVE_DEMANGLING_SUPPORT)
// demangle a given backtrace symbol
//
// Notes:
//  - at the moment, only the symbol is parsed.
//  - one can get the offset by using 
//      "%*[^(]%*[^_]%127[^+)]%64[+x0123456789abcdef]", symbol, offset
//    and checking if sscanf returns 0, 1 or 2
//    (offset includes the leading +)
//  - Similarly one could exctract the address and try to convert it
//    into a filename+line number like addr2line does but this seems
//    to require exteral dependencies. If we want to go down that
//    route, one could look into the inplementation o faddr2line(.c)
//    and/or dladdr.
string Error::_demangle(const char* symbol) {
  size_t size;
  int status;
  char temp[128];
  char* demangled;
  // first, try to demangle a c++ name
  // decryption:
  //   %*[^(]  matches any number of characters different from "("
  //           the * tells not to store in input var 
  //   %*[^_]  matches any number of characters different from "_"
  //           the * tells not to store in input var 
  //   %127[^)+]  matches at most 127 characters different from "+"
  //              match is stored
  if (1 == sscanf(symbol, "%*[^(]%*[^_]%127[^)+]", temp)) {
    //cout << symbol << " -> " << temp << endl;
    if (NULL != (demangled = abi::__cxa_demangle(temp, NULL, &size, &status))) {
      string result(demangled);
      free(demangled);
      return result;
    }
  }
  //if that didn't work, try to get a regular c symbol
  if (1 == sscanf(symbol, "%127s", temp)) {
    return temp;
  }
 
  //if all else fails, just return the symbol
  return symbol;
}
#endif  // FASTJET_HAVE_DEMANGLING_SUPPORT && FASTJET_HAVE_EXECINFO_H
#endif  // __FJCORE__


//----------------------------------------------------------------------
Error::Error(const std::string & message_in) {
  _message = message_in; 

  if (_print_errors && _default_ostr){
    ostringstream oss;
    oss << "fastjet::Error:  "<< message_in << endl;

#ifndef __FJCORE__
    // only print the stack if execinfo is available and stack enabled
#ifdef FASTJET_HAVE_EXECINFO_H
    if (_print_backtrace){
      void * array[10];
      char ** messages;
 
      int size = backtrace(array, 10);
      messages = backtrace_symbols(array, size);
      
      oss << "stack:" << endl;
      for (int i = 1; i < size && messages != NULL; ++i){
#ifdef FASTJET_HAVE_DEMANGLING_SUPPORT
	oss << "  #" << i << ": " << _demangle(messages[i])
	    << " [" << messages[i]  << "]" << endl;
#else
	oss << "  #" << i << ": " << messages[i] << endl;
#endif
      }
      free(messages);
    }
#endif  // FASTJET_HAVE_EXECINFO_H
#endif  // __FJCORE__

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

//----------------------------------------------------------------------
void Error::set_print_backtrace(bool enabled) {
#if (!defined(FASTJET_HAVE_EXECINFO_H)) || defined(__FJCORE__)
  if (enabled) {
    _execinfo_undefined.warn("Error::set_print_backtrace(true) will not work with this build of FastJet");
  }
#endif    
  _print_backtrace = enabled;
}

FASTJET_END_NAMESPACE

