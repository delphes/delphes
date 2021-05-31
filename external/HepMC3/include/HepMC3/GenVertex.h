// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
/// @file GenVertex.h
/// @brief Definition of \b class GenVertex
//
#ifndef HEPMC3_GENVERTEX_H
#define HEPMC3_GENVERTEX_H
#include <string>
#include "HepMC3/GenParticle_fwd.h"
#include "HepMC3/GenVertex_fwd.h"
#include "HepMC3/Data/GenVertexData.h"
#include "HepMC3/FourVector.h"

namespace HepMC3 {

/** Deprecated */
using namespace std;

class Attribute;
class GenEvent;

/// Stores vertex-related information
class GenVertex : public std::enable_shared_from_this<GenVertex> {

    friend class GenEvent;

public:

    /// @name Constructors
    //@{

    /// Default constructor
    GenVertex( const FourVector& position = FourVector::ZERO_VECTOR() );

    /// Constructor based on vertex data
    GenVertex( const GenVertexData& data );

    //@}

public:

    /// @name Accessors
    //@{

    /// Get parent event
    GenEvent* parent_event() { return m_event; }

    /// Get parent event
    const GenEvent* parent_event() const { return m_event; }

    /// Check if this vertex belongs to an event
    bool in_event() const { return parent_event() != nullptr; }

    /// Get the vertex unique identifier
    ///
    /// @note This is not the same as id() in HepMC v2, which is now @c status()
    int id() const { return m_id; }

    /// @brief set the vertex identifier
    void set_id(int id);

    /// Get vertex status code
    int status() const { return m_data.status; }
    /// Set vertex status code
    void set_status(int stat) { m_data.status = stat; }

    /// Get vertex data
    const GenVertexData& data() const { return m_data; }

    /// Add incoming particle
    void add_particle_in ( GenParticlePtr p);
    /// Add outgoing particle
    void add_particle_out( GenParticlePtr p);
    /// Remove incoming particle
    void remove_particle_in ( GenParticlePtr p);
    /// Remove outgoing particle
    void remove_particle_out( GenParticlePtr p);

    /// Get list of incoming particles
    const std::vector<GenParticlePtr>& particles_in() { return m_particles_in; }
    /// Get list of incoming particles (for const access)
    const std::vector<ConstGenParticlePtr>& particles_in() const;
    /// Get list of outgoing particles
    const std::vector<GenParticlePtr>& particles_out() { return m_particles_out; }
    /// Get list of outgoing particles (for const access)
    const std::vector<ConstGenParticlePtr>& particles_out() const;

    /// @brief Get vertex position
    ///
    /// Returns the position of this vertex. If a position is not set on _this_ vertex,
    /// the production vertices of ancestors are searched to find the inherited position.
    /// FourVector(0,0,0,0) is returned if no position information is found.
    ///
    const FourVector& position() const;
    /// @brief Check if position of this vertex is set
    bool has_set_position() const { return !(m_data.position.is_zero()); }

    /// Set vertex position
    void set_position(const FourVector& new_pos); //!<

    /// @brief Add event attribute to this vertex
    ///
    /// This will overwrite existing attribute if an attribute with
    /// the same name is present. The attribute will be stored in the
    /// parent_event(). @return false if there is no parent_event();
    bool add_attribute(const std::string& name, std::shared_ptr<Attribute> att);

    /// @brief Get list of names of attributes assigned to this particle
    std::vector<std::string> attribute_names() const;

    /// @brief Remove attribute
    void remove_attribute(const std::string& name);

    /// @brief Get attribute of type T
    template<class T>
    std::shared_ptr<T> attribute(const std::string& name) const;

    /// @brief Get attribute of any type as string
    std::string attribute_as_string(const std::string& name) const;

    /// @name Deprecated functionality
    //@{


    /// Add incoming particle by raw pointer
    /// @deprecated Use GenVertex::add_particle_in( const GenParticlePtr &p ) instead
    void add_particle_in ( GenParticle *p ) { add_particle_in( GenParticlePtr(p) ); }

    /// Add outgoing particle by raw pointer
    /// @deprecated Use GenVertex::add_particle_out( const GenParticlePtr &p ) instead
    void add_particle_out( GenParticle *p ) { add_particle_out( GenParticlePtr(p) ); }


    //@}


private:

    /// @name Fields
    //@{
    GenEvent       *m_event;  //!< Parent event
    int             m_id;     //!< Vertex id
    GenVertexData   m_data;   //!< Vertex data

    std::vector<GenParticlePtr>  m_particles_in;  //!< Incoming particle list

    std::vector<GenParticlePtr>  m_particles_out; //!< Outgoing particle list
    //@}

};


} // namespace HepMC3

#include "HepMC3/GenEvent.h"
namespace HepMC3 {
/// @brief Get attribute of type T
template<class T> std::shared_ptr<T> GenVertex::attribute(const std::string& name) const {
    return parent_event()?
           parent_event()->attribute<T>(name, id()): std::shared_ptr<T>();
}
}

#endif
