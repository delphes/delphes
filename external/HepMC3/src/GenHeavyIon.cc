// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file GenHeavyIon.cc
 *  @brief Implementation of \b class GenHeavyIon
 *
 */
#include "HepMC3/GenHeavyIon.h"
#include <cstring> // memcmp
#include <cstdlib> // atoi
#include <cstdio> // sprintf
#include <sstream> // sprintf

namespace HepMC3 {

bool GenHeavyIon::from_string(const std::string &att) {

#ifdef HEPMC3_NO_DEPRECATED
    double spectator_neutrons, spectator_protons, eccentricity;
#endif

    std::istringstream is(att);
    std::string version;

    if ( att[0]  != 'v' ) {
        is >> Ncoll_hard >> Npart_proj >> Npart_targ >> Ncoll
           >> spectator_neutrons >> spectator_protons
           >> N_Nwounded_collisions >> Nwounded_N_collisions
           >> Nwounded_Nwounded_collisions >> impact_parameter
           >> event_plane_angle >> eccentricity >> sigma_inel_NN
           >> centrality;
        return !is.fail();
    } else
        is >> version;



    is >> Ncoll_hard >> Npart_proj >> Npart_targ >> Ncoll;
    if ( version == "v0" ) is >> spectator_neutrons >> spectator_protons;
    is >> N_Nwounded_collisions >> Nwounded_N_collisions
       >> Nwounded_Nwounded_collisions >> impact_parameter
       >> event_plane_angle;
    if ( version == "v0" ) is >> eccentricity;
    is >> sigma_inel_NN >> centrality;
    if ( version != "v0" ) is >> user_cent_estimate;
    is >> Nspec_proj_n >> Nspec_targ_n >> Nspec_proj_p >> Nspec_targ_p;

    int N, ord;
    is >> N;
    for ( int i = 0; i < N; ++i ) {
        is >> ord;
        is >> participant_plane_angles[ord];
    }
    is >> N;
    for ( int i = 0; i < N; ++i ) {
        is >> ord;
        is >> eccentricities[ord];
    }

    return !is.fail();
}

bool GenHeavyIon::to_string(std::string &att) const {
    std::ostringstream os;

#ifndef HEPMC3_NO_DEPRECATED
    if ( !forceoldformat ) os << "v0 ";
#else
    os << "v1 ";
#endif

    os << std::setprecision(8)
       << Ncoll_hard << " " << Npart_proj << " "
       << Npart_targ << " " << Ncoll << " "
#ifndef HEPMC3_NO_DEPRECATED
       << spectator_neutrons << " " << spectator_protons << " "
#endif
       << N_Nwounded_collisions << " " << Nwounded_N_collisions << " "
       << Nwounded_Nwounded_collisions << " " << impact_parameter << " "
       << event_plane_angle << " "
#ifndef HEPMC3_NO_DEPRECATED
       << eccentricity << " "
#endif
       << sigma_inel_NN << " " << centrality << " " << user_cent_estimate << " "
       << Nspec_proj_n << " " << Nspec_targ_n << " "
       << Nspec_proj_p << " " << Nspec_targ_p << " ";

    os << participant_plane_angles.size();
    for ( std::map<int,double>::const_iterator it = participant_plane_angles.begin();
            it != participant_plane_angles.end(); ++it )
        os << " " << it->first << " " << it->second;

    os << " " << eccentricities.size();
    for ( std::map<int,double>::const_iterator it = eccentricities.begin();
            it != eccentricities.end(); ++it )
        os << " " << it->first << " " << it->second;

    att = os.str();

    return true;
}


#ifndef HEPMC3_NO_DEPRECATED

bool GenHeavyIon::operator==( const GenHeavyIon& a ) const {
    return ( memcmp( (void*) this, (void*) &a, sizeof(class GenHeavyIon) ) == 0 );
}

bool GenHeavyIon::operator!=( const GenHeavyIon& a ) const {
    return !( a == *this );
}

void GenHeavyIon::set( const int&nh, const int&np, const int&nt, const int&nc, const int&ns, const int&nsp,
                       const int&nnw, const int&nwn, const int&nwnw,
                       const double& im, const double& pl, const double& ec, const double& s, const double& cent, const double& usrcent ) {
    Ncoll_hard                   = nh;
    Npart_proj                   = np;
    Npart_targ                   = nt;
    Ncoll                        = nc;
    spectator_neutrons           = ns;
    spectator_protons            = nsp;
    N_Nwounded_collisions        = nnw;
    Nwounded_N_collisions        = nwn;
    Nwounded_Nwounded_collisions = nwnw;
    impact_parameter             = im;
    event_plane_angle            = pl;
    eccentricity                 = ec;
    sigma_inel_NN                = s;
    centrality                   = cent;
    user_cent_estimate           = usrcent;
}

bool GenHeavyIon::is_valid() const {
    if( Ncoll_hard                   != 0 ) return true;
    if( Npart_proj                   != 0 ) return true;
    if( Npart_targ                   != 0 ) return true;
    if( Ncoll                        != 0 ) return true;
    if( spectator_neutrons           != 0 ) return true;
    if( spectator_protons            != 0 ) return true;
    if( N_Nwounded_collisions        != 0 ) return true;
    if( Nwounded_N_collisions        != 0 ) return true;
    if( Nwounded_Nwounded_collisions != 0 ) return true;
    if( impact_parameter             != 0 ) return true;
    if( event_plane_angle            != 0 ) return true;
    if( eccentricity                 != 0 ) return true;
    if( sigma_inel_NN                != 0 ) return true;
    if( centrality                   != 0 ) return true;
    return false;
}

#endif

} // namespace HepMC3
