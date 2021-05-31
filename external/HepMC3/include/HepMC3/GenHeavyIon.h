// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_HEAVYION_H
#define HEPMC3_HEAVYION_H
/**
 *  @file GenHeavyIon.h
 *  @brief Definition of attribute \b class GenHeavyIon
 *
 *  @class HepMC3::GenHeavyIon
 *  @brief Stores additional information about Heavy Ion generator
 *
 *  This is an example of event attribute used to store Heavy Ion information
 *
 *  @ingroup attributes
 *
 */
#include <iostream>
#include <map>
#include "HepMC3/Attribute.h"

namespace HepMC3 {
/** Deprecated */
using namespace std;

class GenHeavyIon : public Attribute {

public:

    /// Empty default constructor.
    GenHeavyIon()
        : Ncoll_hard(-1), Npart_proj(-1), Npart_targ(-1), Ncoll(-1),
#ifndef HEPMC3_NO_DEPRECATED
          spectator_neutrons(-1), spectator_protons(-1),
#endif
          N_Nwounded_collisions(-1), Nwounded_N_collisions(-1),
          Nwounded_Nwounded_collisions(-1), impact_parameter(-1.0),
          event_plane_angle(-1.0),
#ifndef HEPMC3_NO_DEPRECATED
          eccentricity(-1.0),
#endif
          sigma_inel_NN(-1.0), centrality(-1.0), user_cent_estimate(-1.0),
          Nspec_proj_n(-1), Nspec_targ_n(-1),
          Nspec_proj_p(-1), Nspec_targ_p(-1), forceoldformat(false) {}

//
// Fields
//
public:

    ///
    /// @brief the number of hard nucleon-nucleon collisions.
    ///
    /// Model-dependent. Usually the number of nucleon-nucleon
    /// collisions containing a special signal process. A negative
    /// value means that the information is not available.
    int    Ncoll_hard;

    /// @brief the number of participating nucleons in the projectile.
    ///
    /// The number of nucleons in the projectile participating in an
    /// inelastic collision (see Ncoll). A negative value means that
    /// the information is not available.
    int    Npart_proj;

    /// @brief the number of participating nucleons in the target.
    ///
    /// The number of nucleons in the target participating in an
    /// inelastic collision (see Ncoll). A negative value means that
    /// the information is not available.
    int    Npart_targ;

    /// @brief the number of inelastic nucleon-nucleon collisions.
    ///
    /// Note that a one participating nucleon can be involved in many
    /// inelastic collisions, and that inelastic also includes
    /// diffractive excitation. A negative value means that the
    /// information is not available.
    ///
    int    Ncoll;

#ifndef HEPMC3_NO_DEPRECATED
    /// @brief Total number of spectator neutrons.
    ///
    /// HEPMC3_DEPRECATED("Use Nspec_proj_n and Nspec_targ_n instead.")
    int    spectator_neutrons;

    /// @brief Total number of spectator protons.
    ///
    /// HEPMC3_DEPRECATED("Use Nspec_proj_p and Nspec_targ_p instead.")
    int    spectator_protons;
#endif

    /// @brief Collisions with a diffractively excited target nucleon.
    ///
    /// The number of single diffractive nucleon-nucleon collisions
    /// where the target nucleon is excited. A negative value means
    /// that the information is not available.
    int    N_Nwounded_collisions;

    /// @brief Collisions with a diffractively excited projectile nucleon.
    ///
    /// The number of single diffractive nucleon-nucleon collisions
    /// where the projectile nucleon is excited. A negative value
    /// means that the information is not available.
    int    Nwounded_N_collisions;

    /// @brief Non-diffractive or doubly diffractive collisions.
    ///
    /// The number of nucleon-nucleon collisions where both projectile
    /// and target nucleons are wounded. A negative value means that
    /// the information is not available.
    int    Nwounded_Nwounded_collisions;

    /// @brief The impact parameter.
    ///
    /// The impact parameter given in units of femtometer. A negative
    /// value means that the information is not available.
    double impact_parameter;

    /// @brief The event plane angle.
    ///
    /// The angle wrt. the x-axix of the impact parameter vector
    /// (pointing frm the target to the projectile). A positive number
    /// between 0 and two pi. A negative value means that the
    /// information is not available.
    double event_plane_angle;

#ifndef HEPMC3_NO_DEPRECATED
    /// @brief The eccentricity.
    ///
    /// HEPMC3_DEPRECATED("Use eccentricities insted.")
    double eccentricity;
#endif

    /// @brief The assumed inelastic nucleon-nucleon cross section
    ///
    /// in units of millibarn. As used in a Glauber calculation to
    /// simulate the distribution in Ncoll. A negative value means
    /// that the information is not available.
    double sigma_inel_NN;

    /// @brief The centrality.
    ///
    /// The generated centrality in percentiles, where 0 is the
    /// maximally central and 100 is the minimally central. A negative
    /// value means that the information is not available.
    double centrality;

    /// @brief A user defined centrality estimator.
    ///
    /// This variable may contain anything a generator feels is
    /// reasonable for estimating centrality. The value should be
    /// non-negative, and a low value corresponds to a low
    /// centrality. A negative value indicatess that the information
    /// is not available.
    double user_cent_estimate;


    /// @brief The number of spectator neutrons in the projectile
    ///
    /// ie. those that thave not participated in any inelastic
    /// nucleon-nucleon collision. A negative value indicatess that
    /// the information is not available.
    int Nspec_proj_n;

    /// @brief The number of spectator neutrons in the target
    ///
    /// ie. those that thave not participated in any inelastic
    /// nucleon-nucleon collision. A negative value indicatess that
    /// the information is not available.
    int Nspec_targ_n;

    /// @brief The number of spectator protons in the projectile
    ///
    /// ie. those that thave not participated in any inelastic
    /// nucleon-nucleon collision. A negative value indicatess that
    /// the information is not available.
    int Nspec_proj_p;

    /// @brief The number of spectator protons in the target
    ///
    /// ie. those that thave not participated in any inelastic
    /// nucleon-nucleon collision. A negative value indicatess that
    /// the information is not available.
    int Nspec_targ_p;

    /// @brief Participant plane angles
    ///
    /// calculated to different orders. The key of the map specifies
    /// the order, and the value gives to the angle wrt. the
    /// event plane.
    std::map<int,double> participant_plane_angles;

    /// @brief Eccentricities
    ///
    /// Calculated to different orders. The key of the map specifies
    /// the order, and the value gives the corresponding eccentricity.
    std::map<int,double> eccentricities;

//
// Functions
//
public:

    /// @brief Implementation of Attribute::from_string.
    bool from_string(const std::string &att) override;

    /// @brief Implementation of Attribute::to_string.
    bool to_string(std::string &att) const  override;

#ifndef HEPMC3_NO_DEPRECATED

    /// @brief Operator ==
    ///
    bool operator==( const GenHeavyIon& ) const;
    /// @brief Operator !=
    ///
    bool operator!=( const GenHeavyIon& ) const;

    /// @brief Set all fields.
    ///
    /// HEPMC3_DEPRECATED("Set individual fields directly instead.")
    /** @brief Set all fields */
    void set( const int&nh, const int&np, const int&nt, const int&nc, const int&ns, const int&nsp,
              const int&nnw=0, const int&nwn=0, const int&nwnw=0,
              const double& im=0., const double& pl=0., const double& ec=0., const double& s=0., const double& cent=0., const double& ucent=0. );

    /// @brief Verify that the instance contains non-zero information.
    ///
    /// HEPMC3_DEPRECATED("Each filed now have default values meaning
    /// that they have not been set")
    bool is_valid() const;

    /// @brief force writing in old format for compatibility purposes.
    ///
    /// HEPMC3_DEPRECATED("This should really not be needed");
    bool forceoldformat;

#endif

};


#ifndef HEPMC3_NO_DEPRECATED
typedef GenHeavyIon HeavyIon; ///< Backward compatibility typedef
#endif


} // namespace HepMC3

#endif
