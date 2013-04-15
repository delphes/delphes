#ifndef __FASTJET_ERROR_HH__
#define __FASTJET_ERROR_HH__

//STARTHEADER
// $Id: Error.hh 2577 2011-09-13 15:11:38Z salam $
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

#include<iostream>
#include<string>
#include "fastjet/internal/base.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup error_handling
/// \class Error
/// base class corresponding to errors that can be thrown by FastJet
class Error {
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

  /// controls whether the error message (and the backtrace, if its printing is enabled) 
  /// is printed out or not
  static void set_print_errors(bool print_errors) {_print_errors = print_errors;}

  /// controls whether the backtrace is printed out with the error message or not.
  /// The default is "false".
  static void set_print_backtrace(bool enabled) {_print_backtrace = enabled;}

  /// sets the default output stream for all errors; by default
  /// cerr; if it's null then error output is suppressed.
  static void set_default_stream(std::ostream * ostr) {
    _default_ostr = ostr;
  }

private:
  std::string _message;                ///< error message
  static bool _print_errors;           ///< do we print anything?
  static bool _print_backtrace;        ///< do we print the backtrace?
  static std::ostream * _default_ostr; ///< the output stream (cerr if not set)
};


FASTJET_END_NAMESPACE

#endif // __FASTJET_ERROR_HH__
