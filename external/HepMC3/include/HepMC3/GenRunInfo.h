// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file GenRunInfo.h
/// @brief Definition of \b class GenRunInfo
///
#ifndef HEPMC3_GENRUNINFO_H
#define HEPMC3_GENRUNINFO_H

#if !defined(__CINT__)
#include "HepMC3/Units.h"
#include "HepMC3/Attribute.h"
#include <mutex>
#endif // __CINT__

#ifdef HEPMC3_ROOTIO
class TBuffer;
#endif

namespace HepMC3 {
/** Deprecated */
using namespace std;

struct GenRunInfoData;

/// @brief Stores run-related information
///
/// Manages run-related information.
/// Contains run-wide attributes
class GenRunInfo {

public:

    /// @brief Interrnal struct for keeping track of tools.
    struct ToolInfo {

        /// @brief The name of the tool.
        std::string name;

        /// @brief The version of the tool.
        std::string version;

        /// @brief Other information about how the tool was used in
        /// the run.
        std::string description;
    };

public:

    /// @brief Default constructor
    GenRunInfo() {}
    /// @brief Copy constructor
    GenRunInfo(const GenRunInfo& r);
    /// @brief Assignmet
    GenRunInfo& operator=(const GenRunInfo& r);

#if !defined(__CINT__)

    /// @brief The vector of tools used to produce this run.
    const std::vector<ToolInfo> & tools() const {
        return m_tools;
    }
    /// @brief The vector of tools used to produce this run.
    std::vector<ToolInfo> & tools() {
        return m_tools;
    }

    /// @brief Check if a weight name is present.
    bool has_weight(const std::string& name) const {
        return m_weight_indices.find(name) !=  m_weight_indices.end();
    }

    /// @brief Return the index corresponding to a weight name.
    /// @return -1 if name was not found
    int weight_index(const std::string& name) const {
        std::map<std::string, int>::const_iterator it = m_weight_indices.find(name);
        return it == m_weight_indices.end()? -1: it->second;
    }

    /// @brief Get the vector of weight names.
    const std::vector<std::string> & weight_names() const {
        return m_weight_names;
    }

    /// @brief Set the names of the weights in this run.
    ///
    /// For consistency, the length of the vector should be the same as
    /// the number of weights in the events in the run.
    void set_weight_names(const std::vector<std::string> & names);

    /// @brief add an attribute
    /// This will overwrite existing attribute if an attribute
    /// with the same name is present
    void add_attribute(const std::string &name,
                       const std::shared_ptr<Attribute> &att) {
        std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
        if ( att ) m_attributes[name] = att;
    }

    /// @brief Remove attribute
    void remove_attribute(const std::string &name) {
        std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
        m_attributes.erase(name);
    }

    /// @brief Get attribute of type T
    template<class T>
    std::shared_ptr<T> attribute(const std::string &name) const;

    /// @brief Get attribute of any type as string
    std::string attribute_as_string(const std::string &name) const;

    /// @brief Get list of attribute names
    std::vector<std::string> attribute_names() const;

    /// @brief Get a copy of the list of attributes
    /// @note To avoid thread issues, this is returns a copy. Better solution may be needed.
    std::map< std::string, std::shared_ptr<Attribute> > attributes() const {
        return m_attributes;
    }


#endif // __CINT__

    /// @name Methods to fill GenRunInfoData and to read it back
    //@{

    /// @brief Fill GenRunInfoData object
    void write_data(GenRunInfoData &data) const;

    /// @brief Fill GenRunInfo based on GenRunInfoData
    void read_data(const GenRunInfoData &data);

#ifdef HEPMC3_ROOTIO
    /// @brief ROOT I/O streamer
    void Streamer(TBuffer &b);
    //@}
#endif

private:

    /// @name Fields
    //@{

#if !defined(__CINT__)

    /// @brief The vector of tools used to produce this run.
    std::vector<ToolInfo> m_tools;

    /// @brief A map of weight names mapping to indices.
    std::map<std::string, int> m_weight_indices;

    /// @brief A vector of weight names.
    std::vector<std::string> m_weight_names;

    /// @brief Map of attributes
    mutable std::map< std::string, std::shared_ptr<Attribute> > m_attributes;

    /// @brief Mutex lock for the m_attibutes map.
    mutable std::recursive_mutex m_lock_attributes;
    //@}

#endif // __CINT__
};

#if !defined(__CINT__)

//
// Template methods
//

template<class T>
std::shared_ptr<T> GenRunInfo::attribute(const std::string &name) const {
    std::lock_guard<std::recursive_mutex> lock(m_lock_attributes);
    std::map< std::string, std::shared_ptr<Attribute> >::iterator i =
        m_attributes.find(name);
    if( i == m_attributes.end() ) return std::shared_ptr<T>();

    if( !i->second->is_parsed() ) {

        std::shared_ptr<T> att = std::make_shared<T>();
        if ( att->from_string(i->second->unparsed_string()) &&
                att->init(*this) ) {
            // update map with new pointer
            i->second = att;

            return att;
        }
        else
            return std::shared_ptr<T>();
    }
    else return std::dynamic_pointer_cast<T>(i->second);
}

#endif // __CINT__

} // namespace HepMC3

#endif
