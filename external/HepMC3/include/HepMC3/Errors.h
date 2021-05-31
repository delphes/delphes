// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
/**
 *  @file Errors.h
 *  @brief Implementation of  error and HEPMC3_HEPMC3_WARNING macros
 *
 */
#ifndef HEPMC3_ERRORS_H
#define HEPMC3_ERRORS_H

#include <iostream>
#include <stdexcept>

namespace HepMC3 {


/// @name Printing macros
//@{

/** @brief Macro for printing error messages */
#define HEPMC3_ERROR(MESSAGE)   if ( Setup::print_errors() )   { std::cerr << "ERROR::" << MESSAGE << std::endl; }

/** @brief Macro for printing HEPMC3_HEPMC3_WARNING messages */
#define HEPMC3_WARNING(MESSAGE) if ( Setup::print_warnings() ) { std::cout << "WARNING::" << MESSAGE << std::endl; }

// Debug messages and code that will not go to the release version
#ifndef HEPMC3_RELEASE_VERSION

/** @brief Macro for printing debug messages with appropriate debug level */
#define HEPMC3_DEBUG(LEVEL,MESSAGE) if( Setup::debug_level()>=(LEVEL) ) { std::cout << "DEBUG(" << LEVEL <<")::" << MESSAGE << std::endl; }
/** @brief Macro for storing code useful for debugging */
#define HEPMC3_DEBUG_CODE_BLOCK( x ) x

#else

#define HEPMC3_DEBUG( x,y )
#define HEPMC3_DEBUG_CODE_BLOCK( x )

#endif

//@}




} // namespace HepMC3

#endif
