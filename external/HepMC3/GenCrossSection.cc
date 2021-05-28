// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file GenCrossSection.cc
 *  @brief Implementation of \b class GenCrossSection
 *
 */
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenEvent.h"
#include <cstring> // memcmp
#include <cstdlib> // atoi
#include <sstream>
#include <iomanip>

namespace HepMC3 {


int GenCrossSection::windx(std::string wName) const {
    if ( !event() || !event()->run_info() ) return 0;
    return event()->run_info()->weight_index(wName);
}

void GenCrossSection::set_cross_section(const double& xs, const double& xs_err,const long& n_acc, const long& n_att ) {
    double cross_section       = xs;
    double cross_section_error = xs_err;
    accepted_events     = n_acc;
    attempted_events    = n_att;
    size_t N=1;
    if ( event() ) N=std::max(event()->weights().size(),N);
    cross_sections = std::vector<double>(N, cross_section);
    cross_section_errors = std::vector<double>(N, cross_section_error);
}


bool GenCrossSection::from_string(const std::string &att) {
    const char *cursor = att.data();
    cross_sections.clear();
    cross_section_errors.clear();


    double cross_section = atof(cursor);
    cross_sections.push_back(cross_section);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    double cross_section_error = atof(cursor);
    cross_section_errors.push_back(cross_section_error);

    if( !(cursor = strchr(cursor+1,' ')) ) {accepted_events = -1; attempted_events = -1;}
    else
    {
        accepted_events = atoi(cursor);
        if( !(cursor = strchr(cursor+1,' ')) ) attempted_events = -1;
        else attempted_events = atoi(cursor);
    }
    size_t N=1;
    if ( event() ) N=std::max(event()->weights().size(),N);
    const size_t max_n_cross_sections=1000;
    while (cross_sections.size()<max_n_cross_sections) {
        if( !(cursor = strchr(cursor+1,' ')) ) break;
        cross_sections.push_back(atof(cursor));
        if( !(cursor = strchr(cursor+1,' ')) ) break;
        cross_section_errors.push_back(atof(cursor));
    }
    if (cross_sections.size()>=max_n_cross_sections)
        HEPMC3_WARNING( "GenCrossSection::from_string: too many optional cross-sections  N="<<cross_sections.size()<<" or ill-formed input:"<<att )
//        if (cross_sections.size()!=N)
//  So far it is not clear if there should be a warning or not
//  Frank suggests no.             HEPMC3_WARNING( "GenCrossSection::from_string: optional cross-sections are available not for all weights")
        size_t oldsize=cross_sections.size();
    for (size_t i=oldsize; i<N; i++) {cross_sections.push_back(cross_section); cross_section_errors.push_back(cross_section_error);}

    return true;
}

bool GenCrossSection::to_string(std::string &att) const {
    std::ostringstream os;

    os << std::setprecision(8) << std::scientific
       << cross_sections.at(0) << " "
       << cross_section_errors.at(0) << " "
       << accepted_events << " "
       << attempted_events;

    for (size_t i = 1; i < cross_sections.size(); ++i )
        os << " " << cross_sections.at(i)
           << " " << cross_section_errors.at(i);

    att = os.str();

    return true;
}

bool GenCrossSection::operator==( const GenCrossSection& a ) const {
    return ( memcmp( (void*)this, (void*) &a, sizeof(class GenCrossSection) ) == 0 );
}

bool GenCrossSection::operator!=( const GenCrossSection& a ) const {
    return !( a == *this );
}

bool GenCrossSection::is_valid() const {
    if( cross_sections.size()       == 0 ) return false;
    if( cross_section_errors.size() == 0 ) return false;
    if( cross_section_errors.size()!=cross_sections.size() ) return false;
    if( cross_sections.at(0)       != 0 ) return true;
    if( cross_section_errors.at(0) != 0 ) return true;
    return false;
}

} // namespace HepMC3
