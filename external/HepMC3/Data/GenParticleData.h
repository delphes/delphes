// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_DATA_GENPARTICLEDATA_H
#define HEPMC3_DATA_GENPARTICLEDATA_H
/**
 *  @file GenParticleData.h
 *  @brief Definition of \b class GenParticleData
 *
 *  @struct HepMC3::GenParticleData
 *  @brief Stores serializable particle information
 *
 *  @ingroup data
 *
 */
#include "HepMC3/FourVector.h"

namespace HepMC3 {

// NOTE: Keep in mind the data alignment
//       Currently it's 8b alignment = 56b total
struct GenParticleData {
    int        pid;               ///< PDG ID
    int        status;            ///< Status
    bool       is_mass_set;       ///< Check if generated mass is set
    double     mass;              ///< Generated mass (if set)
    FourVector momentum;          ///< Momentum
};

} // namespace HepMC

#endif
