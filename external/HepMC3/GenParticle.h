// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_GENPARTICLE_H
#define HEPMC3_GENPARTICLE_H
/**
 *  @file GenParticle.h
 *  @brief Definition of \b class GenParticle
 *
 *  @class HepMC3::GenParticle
 *  @brief Stores particle-related information
 *
 */
#include <string>
#include "HepMC3/Data/GenParticleData.h"
#include "HepMC3/FourVector.h"

#include "HepMC3/GenParticle_fwd.h"
#include "HepMC3/GenVertex_fwd.h"

namespace HepMC3 {

/** Deprecated */
using namespace std;

class GenEvent;
class Attribute;

class GenParticle : public std::enable_shared_from_this<GenParticle> {

    friend class GenVertex;
    friend class GenEvent;

//
// Constructors
//
public:
    /** @brief Default constructor */
    GenParticle( const FourVector &momentum = FourVector::ZERO_VECTOR(), int pid = 0, int status = 0 );

    /** @brief Constructor based on particle data */
    GenParticle( const GenParticleData &data );

//
// Functions
//
public:
    /** @brief Check if this particle belongs to an event */
    bool in_event() const { return (bool)(m_event); }

//
// Accessors
//
public:

    GenEvent*              parent_event() { return m_event; } //!< Get parent event
    const GenEvent*              parent_event() const { return m_event; } //!< Get parent event
    int                    id()           const { return m_id;    } //!< Get particle id
    const GenParticleData& data()         const { return m_data;  } //!< Get particle data


    ConstGenVertexPtr production_vertex() const;        //!< Get production vertex (const version)
    ConstGenVertexPtr end_vertex() const;               //!< Get end vertex (const version)

    GenVertexPtr production_vertex();        //!< Get production vertex
    GenVertexPtr end_vertex();               //!< Get end vertex

    /// @brief Convenience access to immediate incoming particles via production vertex
    /// @note Less efficient than via the vertex since return must be by value (in case there is no vertex)
    std::vector<GenParticlePtr> parents();

    /// @brief Convenience access to immediate incoming particles via production vertex (const version)
    /// @note Less efficient than via the vertex since return must be by value (in case there is no vertex)
    std::vector<ConstGenParticlePtr> parents() const;

    /// @brief Convenience access to immediate outgoing particles via end vertex
    /// @note Less efficient than via the vertex since return must be by value (in case there is no vertex)
    std::vector<GenParticlePtr> children();

    /// @brief Convenience access to immediate outgoing particles via end vertex
    /// @note Less efficient than via the vertex since return must be by value (in case there is no vertex)
    std::vector<ConstGenParticlePtr> children() const;

    int   pid()                   const { return m_data.pid;            } //!< Get PDG ID
    int   abs_pid()               const { return abs(pid());            } //!< Get absolute value of PDG ID
    int   status()                const { return m_data.status;         } //!< Get status code
    const FourVector& momentum()  const { return m_data.momentum;       } //!< Get momentum
    bool  is_generated_mass_set() const { return m_data.is_mass_set;    } //!< Check if generated mass is set

    /// @brief Get generated mass
    ///
    /// This function will return mass as set by a generator/tool.
    /// If not set, it will return momentum().m()
    double generated_mass() const;


    void set_pid(int pid);                         //!< Set PDG ID
    void set_status(int status);                   //!< Set status code
    void set_momentum(const FourVector& momentum); //!< Set momentum
    void set_generated_mass(double m);             //!< Set generated mass
    void unset_generated_mass();                   //!< Declare that generated mass is not set

    /** @brief Add an attribute to this particle
     *
     *  This will overwrite existing attribute if an attribute with
     *  the same name is present. The attribute will be stored in the
     *  parent_event(). @return false if there is no parent_event();
     */
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

    /// @brief Get PDG ID
    /// @deprecated Use pid() instead
    int pdg_id() const { return pid(); }

    /// @brief Set PDG ID
    /// @deprecated Use set_pid() instead
    void set_pdg_id(const int& pidin) { set_pid(pidin); }


    //@}
//
// Fields
//
private:
    GenEvent        *m_event; //!< Parent event
    int              m_id;    //!< Index
    GenParticleData  m_data;  //!< Particle data

    std::weak_ptr<GenVertex>    m_production_vertex; //!< Production vertex
    std::weak_ptr<GenVertex>    m_end_vertex;        //!< End vertex
//    weak_ptr<GenParticle>  m_this;              //!< Pointer to shared pointer managing this particle
};

} // namespace HepMC3

#include "HepMC3/GenEvent.h"
namespace HepMC3 {
/// @brief Get attribute of type T
template<class T> std::shared_ptr<T> GenParticle::attribute(const std::string& name) const {
    return parent_event()?
           parent_event()->attribute<T>(name, id()): std::shared_ptr<T>();
}
}
#endif
