#ifndef __FASTJET_ERROR_HH__
#define __FASTJET_ERROR_HH__

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

#include<iostream>
#include<string>
#include "fastjet/internal/base.hh"
#include "fastjet/config.h"
//#include <exception>
#if (!defined(FASTJET_HAVE_EXECINFO_H)) || defined(__FASTJET_ONLY_CORE__)
#include "fastjet/LimitedWarning.hh"
#endif
#ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
#include <atomic>
#include <mutex>
#endif // FASTJET_HAVE_LIMITED_THREAD_SAFETY

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup error_handling
/// \class Error
/// base class corresponding to errors that can be thrown by FastJet
class Error{
public:
  /// default constructors
  Error() {}

  /// ctor from an error message
  ///   \param message to be printed
  /// Note: in addition to the error message, one can choose to print the
  /// backtrace (showing the last few calls before the error) by 
  /// using set_print_backtrace(true). The default is "false".
  Error(const std::string & message);

  /// virtual dummy dtor
  virtual ~Error() {}

  /// the error message
  std::string message() const {return _message;}

  /// an alternative access to the error message (more standard)
  std::string description() const {return message();}
  
  /// controls whether the error message (and the backtrace, if its printing is enabled) 
  /// is printed out or not
  static void set_print_errors(bool print_errors) {_print_errors = print_errors;}

  /// controls whether the backtrace is printed out with the error message or not.
  /// The default is "false".
  static void set_print_backtrace(bool enabled);

  /// sets the default output stream for all errors; by default
  /// cerr; if it's null then error output is suppressed.
  static void set_default_stream(std::ostream * ostr) {
    _default_ostr = ostr;
  }

#ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
  /// sets the default output stream for all errors (by default
  /// cerr; passing a null pointer prevents errors from being output)
  /// The second argument is a mutex that would be used to guarantee 
  /// that only a single thread writes to the stream at a time
  static void set_default_stream_and_mutex(std::ostream * ostr, std::mutex * stream_mutex) {
    _default_ostr = ostr;
    _stream_mutex = stream_mutex;
  }
#endif // FASTJET_HAVE_LIMITED_THREAD_SAFETY

private:

#ifndef __FASTJET_ONLY_CORE__
#if defined(FASTJET_HAVE_EXECINFO_H) && defined(FASTJET_HAVE_DEMANGLING_SUPPORT)
  /// demangle a given backtrace symbol
  std::string _demangle(const char* symbol);
#endif
#endif

  std::string _message;                ///< error message

#ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY
  FASTJET_WINDLL static std::atomic<bool> _print_errors;           ///< do we print anything?
  FASTJET_WINDLL static std::atomic<bool> _print_backtrace;        ///< do we print the backtrace?
  FASTJET_WINDLL static std::atomic<std::ostream *> _default_ostr; ///< the output stream (cerr if not set)
  FASTJET_WINDLL static std::atomic<std::mutex *> _stream_mutex; ///< the mutex for the output stream (nullptr if not set)
#else
  FASTJET_WINDLL static bool _print_errors;           ///< do we print anything?
  FASTJET_WINDLL static bool _print_backtrace;        ///< do we print the backtrace?
  FASTJET_WINDLL static std::ostream * _default_ostr; ///< the output stream (cerr if not set)
#endif // FASTJET_HAVE_LIMITED_THREAD_SAFETY


#if (!defined(FASTJET_HAVE_EXECINFO_H)) || defined(__FASTJET_ONLY_CORE__)
  FASTJET_WINDLL static LimitedWarning _execinfo_undefined;
#endif
};


/// @ingroup error_handling
/// \class InternalError
/// class corresponding to critical internal errors
/// 
/// This is an error class (derived from Error) meant for serious,
/// critical, internal errors that we still want to be catchable by an
/// end-user [e.g. a serious issue in clustering where the end-user
/// can catch it and retry with a different strategy]
///
/// Please directly contact the FastJet authors if you see such an
/// error.
class InternalError : public Error{
public:
  /// ctor with error message:
  /// just add a bit of info to the message and pass it to the base class
  InternalError(const std::string & message_in) : Error(std::string("*** CRITICAL INTERNAL FASTJET ERROR *** CONTACT THE AUTHORS *** ") + message_in){ }
};

FASTJET_END_NAMESPACE

#endif // __FASTJET_ERROR_HH__
