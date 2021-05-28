// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file GenVertex_fwd.h
/// @brief Minimal forward declarations for GenVertex
///

#ifndef HEPMC3_GENVERTEX_FWD_H
#define HEPMC3_GENVERTEX_FWD_H

#include <memory>

namespace HepMC3 {

class GenVertex;

using GenVertexPtr = std::shared_ptr<GenVertex>;
using ConstGenVertexPtr = std::shared_ptr<const GenVertex>;

template<typename T>
using GenVertex_type = typename std::conditional<std::is_const<typename T::element_type>::value, ConstGenVertexPtr, GenVertexPtr>::type;

}

#endif
