// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_DATA_GENEVENTDATA_H
#define HEPMC3_DATA_GENEVENTDATA_H
/**
 *  @file GenEventData.h
 *  @brief Definition of \b struct GenEventData
 *
 *  @struct HepMC3::GenEventData
 *  @brief Stores serializable event information
 *
 *  @ingroup data
 *
 */
#include <vector>
#include <string>
#include "HepMC3/Data/GenParticleData.h"
#include "HepMC3/Data/GenVertexData.h"
#include "HepMC3/Units.h"

namespace HepMC3 {

struct GenEventData {
    int                 event_number;  ///< Event number
    Units::MomentumUnit momentum_unit; ///< Momentum unit
    Units::LengthUnit   length_unit;   ///< Length unit

    std::vector<GenParticleData> particles; ///< Particles
    std::vector<GenVertexData>   vertices;  ///< Vertices
    std::vector<double>          weights;   ///< Weights

    FourVector event_pos;                   ///< Event position

    /** @brief First id of the vertex links
     *
     *  If this id is positive - it is the incoming particle id
     *  of a vertex which id is written in GenEventData::links2
     *
     *  If this id is negative - it's the id of a vertex which
     *  outgoing particle id is written in GenEventData::links2
     *
     *  The links1[i] points to links2[i].
     *  In case links1[i] is particle, links2[i] is end vertex.
     *  In case links2[i] is vertex, links2[i] is outgoing particle.
     *  An example of usage is given in documentation.
     *
     */
    std::vector<int> links1;
    std::vector<int> links2; ///< Second id of the vertex links

    std::vector<int>         attribute_id;     ///< Attribute owner id
    std::vector<std::string> attribute_name;   ///< Attribute name
    std::vector<std::string> attribute_string; ///< Attribute serialized as string
};

} // namespace HepMC

#endif
