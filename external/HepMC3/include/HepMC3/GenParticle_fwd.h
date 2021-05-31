// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file GenParticle_fwd.h
/// @brief Minimal forward declarations for GenParticle
///

#ifndef HEPMC3_GENPARTICLE_FWD_H
#define HEPMC3_GENPARTICLE_FWD_H

#include <memory>
#include <vector>

namespace HepMC3 {

class GenParticle;

using GenParticlePtr = std::shared_ptr<GenParticle>;
using ConstGenParticlePtr = std::shared_ptr<const GenParticle>;

using GenParticles = std::vector<GenParticlePtr>;
using ConstGenParticles = std::vector<ConstGenParticlePtr>;

/// An alias to a vector of GenParticle pointers whose constness depends on the
/// constness of the template shared_ptr param
/// This is convenient for declaring the return types based on the input type
template<typename T>
using GenParticles_type = typename std::conditional<std::is_const<typename T::element_type>::value, ConstGenParticles, GenParticles>::type;

}

#endif
