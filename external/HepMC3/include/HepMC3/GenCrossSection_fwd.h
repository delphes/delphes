// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file GenCrossSection_fwd.h
/// @brief Minimal forward declarations for GenCrossSection
///

#ifndef HEPMC3_GENCROSSSECTION_FWD_H
#define HEPMC3_GENCROSSSECTION_FWD_H

#include <memory>

namespace HepMC3 {

class GenCrossSection;

using GenCrossSectionPtr = std::shared_ptr<GenCrossSection>;
using ConstGenCrossSectionPtr = std::shared_ptr<const GenCrossSection>;

}

#endif
