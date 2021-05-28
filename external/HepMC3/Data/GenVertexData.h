// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_DATA_GENVERTEXDATA_H
#define HEPMC3_DATA_GENVERTEXDATA_H
/**
 *  @file GenVertexData.h
 *  @brief Definition of \b class GenVertexData
 *
 *  @struct HepMC3::GenVertexData
 *  @brief Stores serializable vertex information
 *
 *  @ingroup data
 *
 */
#include "HepMC3/FourVector.h"

namespace HepMC3 {

struct GenVertexData {
    int        status;   ///< Vertex status
    FourVector position; ///< Position in time-space
    
    /// @brief Check if this struct fields are zero
    bool is_zero() const {
        if( status ) return false;
        
        return position.is_zero();
    }
};

} // namespace HepMC

#endif
