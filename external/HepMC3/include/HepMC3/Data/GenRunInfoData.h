// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_DATA_GENRUNINFODATA_H
#define HEPMC3_DATA_GENRUNINFODATA_H
/**
 *  @file GenRunInfoData.h
 *  @brief Definition of \b struct GenRunInfoData
 *
 *  @struct HepMC3::GenRunInfoData
 *  @brief Stores serializable run information
 *
 *  @ingroup data
 *
 */
#include <vector>
#include <string>

namespace HepMC3 {

struct GenRunInfoData {
    std::vector<std::string> weight_names;     ///< Weight names

    std::vector<std::string> tool_name;        ///< Tool names
    std::vector<std::string> tool_version;     ///< Tool versions
    std::vector<std::string> tool_description; ///< Tool descriptions

    std::vector<std::string> attribute_name;   ///< Attribute name
    std::vector<std::string> attribute_string; ///< Attribute serialized as string
};

} // namespace HepMC

#endif
