// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file GenParticle.cc
 *  @brief Implementation of \b class GenParticle
 *
 */
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Setup.h"
#include "HepMC3/Attribute.h"

namespace HepMC3 {

GenParticle::GenParticle( const FourVector &mom, int pidin, int stat):
    m_event(nullptr),
    m_id(0) {
    m_data.pid               = pidin;
    m_data.momentum          = mom;
    m_data.status            = stat;
    m_data.is_mass_set       = false;
    m_data.mass              = 0.0;
}

GenParticle::GenParticle( const GenParticleData &dat ):
    m_event(nullptr),
    m_id(0),
    m_data(dat) {
}

double GenParticle::generated_mass() const {
    if(m_data.is_mass_set) return m_data.mass;
    else                   return m_data.momentum.m();
}

void GenParticle::set_pid(int pidin) {
    m_data.pid = pidin;
}

void GenParticle::set_status(int stat) {
    m_data.status = stat;
}

void GenParticle::set_momentum(const FourVector& mom) {
    m_data.momentum = mom;
}

void GenParticle::set_generated_mass(double m) {
    m_data.mass        = m;
    m_data.is_mass_set = true;
}

void GenParticle::unset_generated_mass() {
    m_data.mass        = 0.;
    m_data.is_mass_set = false;
}

GenVertexPtr GenParticle::production_vertex() {
    return m_production_vertex.lock();
}

ConstGenVertexPtr GenParticle::production_vertex() const {
    return std::const_pointer_cast<const GenVertex>(m_production_vertex.lock());
}

GenVertexPtr GenParticle::end_vertex() {
    return m_end_vertex.lock();
}

ConstGenVertexPtr GenParticle::end_vertex() const {
    return std::const_pointer_cast<const GenVertex>(m_end_vertex.lock());
}

std::vector<GenParticlePtr> GenParticle::parents() {
    return (m_production_vertex.expired())? std::vector<GenParticlePtr>() : production_vertex()->particles_in();
}

std::vector<ConstGenParticlePtr> GenParticle::parents() const {
    return (m_production_vertex.expired()) ? std::vector<ConstGenParticlePtr>() : production_vertex()->particles_in();
}

std::vector<GenParticlePtr> GenParticle::children() {
    return (m_end_vertex.expired())? std::vector<GenParticlePtr>() : end_vertex()->particles_out();
}

std::vector<ConstGenParticlePtr> GenParticle::children() const {
    return (m_end_vertex.expired()) ? std::vector<ConstGenParticlePtr>() : end_vertex()->particles_out();
}

bool GenParticle::add_attribute(const std::string& name, std::shared_ptr<Attribute> att) {
    if ( !parent_event() ) return false;
    parent_event()->add_attribute(name, att, id());
    return true;
}

std::vector<std::string> GenParticle::attribute_names() const {
    if ( parent_event() ) return parent_event()->attribute_names(id());

    return std::vector<std::string>();
}

void GenParticle::remove_attribute(const std::string& name) {
    if ( parent_event() ) parent_event()->remove_attribute(name, id());
}

std::string GenParticle::attribute_as_string(const std::string& name) const {
    return parent_event() ? parent_event()->attribute_as_string(name, id()) : std::string();
}

} // namespace HepMC3
