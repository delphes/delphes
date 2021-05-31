// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file Setup.cc
 *  @brief Implementation of Setup class
 *
 */
#include "HepMC3/Setup.h"

namespace HepMC3 {

const unsigned int Setup::DEFAULT_DOUBLE_ALMOST_EQUAL_MAXULPS = 10;
const double       Setup::DOUBLE_EPSILON                      = 10e-20;
bool Setup::print_errors()                      { return m_is_printing_errors;    }
void Setup::set_print_errors(const bool flag)   { m_is_printing_errors   = flag;  }
bool Setup::print_warnings()                    { return m_is_printing_warnings;  }
void Setup::set_print_warnings(const bool flag) { m_is_printing_warnings = flag;  }
int  Setup::debug_level()                       { return m_debug_level;           }
void Setup::set_debug_level(const int level)    { m_debug_level          = level; }
bool Setup::m_is_printing_errors    = true;
bool Setup::m_is_printing_warnings  = true;
int  Setup::m_debug_level           = 5;

} // namespace HepMC3
