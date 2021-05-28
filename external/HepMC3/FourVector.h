// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_FOURVECTOR_H
#define HEPMC3_FOURVECTOR_H
/**
 *  @file FourVector.h
 *  @brief Definition of \b class FourVector
 */
#include <cmath>
#ifndef M_PI
/** @brief Definition of PI. Needed on some platforms */
#define M_PI 3.14159265358979323846264338327950288
#endif
namespace HepMC3 {


/**
 *  @brief Generic 4-vector
 *
 *  Interpretation of its content depends on accessors used: it's much simpler to do this
 *  than to distinguish between space and momentum vectors via the type system (especially
 *  given the need for backward compatibility with HepMC2). Be sensible and don't call
 *  energy functions on spatial vectors! To avoid duplication, most definitions are only
 *  implemented on the spatial function names, with the energy-momentum functions as aliases.
 *
 *  This is @a not intended to be a fully featured 4-vector, but does contain the majority
 *  of common non-boosting functionality, as well as a few support operations on
 *  4-vectors.
 *
 *  The implementations in this class are fully inlined.
 */
class FourVector {
public:

    /** @brief Default constructor */
    FourVector()
        : m_v1(0.0),   m_v2(0.0), m_v3(0.0),    m_v4(0.0)  {}
    /** @brief Sets all FourVector fields */
    FourVector(double xx, double yy, double zz, double ee)
        : m_v1(xx),     m_v2(yy),   m_v3(zz),      m_v4(ee)    {}
    /** @brief Copy constructor */
    FourVector(const FourVector & v)
        : m_v1(v.m_v1), m_v2(v.m_v2), m_v3(v.m_v3), m_v4(v.m_v4) {}


    /// @name Component accessors
    //@{

    /** @brief Set all FourVector fields, in order x,y,z,t */
    void set(double x1, double x2, double x3, double x4) {
        m_v1 = x1;
        m_v2 = x2;
        m_v3 = x3;
        m_v4 = x4;
    }

    /// set component of position/displacement
    void set_component(const int i, const double x)   
    {
    if (i==0) {m_v1=x; return; }
    if (i==1) {m_v2=x; return; }
    if (i==2) {m_v3=x; return; }
    if (i==3) {m_v4=x; return; }
    }
    /// get component of position/displacement
    double get_component(const int i)   const
    {
    if (i==0) return m_v1;
    if (i==1) return m_v2;
    if (i==2) return m_v3;
    if (i==3) return m_v4;
    return 0.0;
    }


    /// x-component of position/displacement
    double x()        const { return m_v1; }
    /// Set x-component of position/displacement
    void set_x(double xx)   { m_v1 = xx;    }
    /// @deprecated Prefer the HepMC-style set_x() function
    void setX(double xx)   { set_x(xx);    }

    /// y-component of position/displacement
    double y()        const { return m_v2; }
    /// Set y-component of position/displacement
    void   set_y(double yy)   { m_v2 = yy;    }
    /// @deprecated Prefer the HepMC-style set_y() function
    void setY(double yy)   { set_y(yy);    }

    /// z-component of position/displacement
    double z()        const { return m_v3; }
    /// Set z-component of position/displacement
    void set_z(double zz)   { m_v3 = zz;    }
    /// @deprecated Prefer the HepMC-style set_z() function
    void setZ(double zz)   { set_z(zz);    }

    /// Time component of position/displacement
    double t()        const { return m_v4; }
    /// Set time component of position/displacement
    void set_t(double tt)   { m_v4 = tt;    }
    /// @deprecated Prefer the HepMC-style set_t() function
    void setT(double tt)   { set_t(tt);    }


    /// x-component of momentum
    double px() const { return x(); }
    /// Set x-component of momentum
    void set_px(double pxx) { set_x(pxx);   }
    /// @deprecated Prefer the HepMC-style set_px() function
    void setPx(double pxx) { set_px(pxx);    }

    /// y-component of momentum
    double py() const { return y(); }
    /// Set y-component of momentum
    void set_py(double pyy) { set_y(pyy);   }
    /// @deprecated Prefer the HepMC-style set_py() function
    void setPy(double pyy) { set_py(pyy);    }

    /// z-component of momentum
    double pz() const { return z(); }
    /// Set z-component of momentum
    void set_pz(double pzz) { set_z(pzz);   }
    /// @deprecated Prefer the HepMC-style set_pz() function
    void setPz(double pzz) { set_pz(pzz);    }

    /// Energy component of momentum
    double e() const { return t(); }
    /// Set energy component of momentum
    void set_e(double ee ) { this->set_t(ee);    }
    /// @deprecated Prefer the HepMC-style set_y() function
    void setE(double ee) { set_e(ee);    }

    //@}


    /// @name Computed properties
    //@{

    /// Squared magnitude of (x, y, z) 3-vector
    double length2()  const { return x()*x() + y()*y() + z()*z(); }
    /// Magnitude of spatial (x, y, z) 3-vector
    double length()  const { return sqrt(length2()); }
    /// Squared magnitude of (x, y) vector
    double perp2()      const { return  x()*x() + y()*y(); }
    /// Magnitude of (x, y) vector
    double perp()      const { return sqrt(perp2()); }
    /// Spacetime invariant interval s^2 = t^2 - x^2 - y^2 - z^2
    double interval() const { return t()*t() - length2(); }

    /// Squared magnitude of p3 = (px, py, pz) vector
    double p3mod2()       const { return length2(); }
    /// Magnitude of p3 = (px, py, pz) vector
    double p3mod()       const { return length(); }
    /// Squared transverse momentum px^2 + py^2
    double pt2()      const { return perp2(); }
    /// Transverse momentum
    double pt()      const { return perp(); }
    /// Squared invariant mass m^2 = E^2 - px^2 - py^2 - pz^2
    double m2()       const { return interval(); }
    /// Invariant mass. Returns -sqrt(-m) if e^2 - P^2 is negative
    double m() const { return (m2() > 0.0) ? std::sqrt(m2()) : -std::sqrt(-m2()); }

    /// Azimuthal angle
    double phi() const { return atan2( y(), x() ); }
    /// Polar angle w.r.t. z direction
    double theta() const {  return atan2( perp(), z() ); }
    // /// Cosine of polar angle w.r.t. z direction
    // double costheta() const {  return z() / p3mod(); }
    /// Pseudorapidity
    double eta() const  { return 0.5*std::log( (p3mod() + pz()) / (p3mod() - pz()) ); }
    /// Rapidity
    double rap() const {   return 0.5*std::log( (e() + pz()) / (e() - pz()) ); }
    /// Absolute pseudorapidity
    double abs_eta() const { return std::abs( eta() ); }
    /// Absolute rapidity
    double abs_rap() const { return std::abs( rap() ); }

    /// Same as eta()
    /// @deprecated Prefer 'only one way to do it', and we don't have equivalent long names for e.g. pid, phi or eta
    double pseudoRapidity() const { return eta(); }

    //@}


    /// @name Comparisons to another FourVector
    //@{

    /// Check if the length of this vertex is zero
    bool is_zero() const { return x() == 0 && y() == 0 && z() == 0 && t() == 0; }

    /// Signed azimuthal angle separation in [-pi, pi]
    double delta_phi(const FourVector &v) const {
        double dphi = phi() - v.phi();
        if (dphi != dphi) return dphi;
        while (dphi >=  M_PI) dphi -= 2.*M_PI;
        while (dphi <  -M_PI) dphi += 2.*M_PI;
        return dphi;
    }

    /// Pseudorapidity separation
    double delta_eta(const FourVector &v) const { return eta() - v.eta(); }

    /// Rapidity separation
    double delta_rap(const FourVector &v) const { return rap() - v.rap(); }

    /// R_eta^2-distance separation dR^2 = dphi^2 + deta^2
    double delta_r2_eta(const FourVector &v) const {
        return delta_phi(v)*delta_phi(v) + delta_eta(v)*delta_eta(v);
    }

    /// R_eta-distance separation dR = sqrt(dphi^2 + deta^2)
    double delta_r_eta(const FourVector &v) const {
        return sqrt( delta_r2_eta(v) );
    }

    /// R_rap^2-distance separation dR^2 = dphi^2 + drap^2
    double delta_r2_rap(const FourVector &v) const {
        return delta_phi(v)*delta_phi(v) + delta_rap(v)*delta_rap(v);
    }

    /// R-rap-distance separation dR = sqrt(dphi^2 + drap^2)
    double delta_r_rap(const FourVector &v) const {
        return sqrt( delta_r2_rap(v) );
    }

    //@}


    /// @name Operators
    //@{

    /// Equality
    bool operator==(const FourVector& rhs) const {
        return x() == rhs.x() && y() == rhs.y() && z() == rhs.z() && t() == rhs.t();
    }
    /// Inequality
    bool operator!=(const FourVector& rhs) const { return !(*this == rhs); }

    /// Arithmetic operator +
    FourVector  operator+ (const FourVector& rhs) const {
        return FourVector( x() + rhs.x(), y() + rhs.y(), z() + rhs.z(), t() + rhs.t() );
    }
    /// Arithmetic operator -
    FourVector  operator- (const FourVector& rhs) const {
        return FourVector( x() - rhs.x(), y() - rhs.y(), z() - rhs.z(), t() - rhs.t() );
    }
    /// Arithmetic operator * by scalar
    FourVector  operator* (const double rhs) const {
        return FourVector( x()*rhs, y()*rhs, z()*rhs, t()*rhs );
    }
    /// Arithmetic operator / by scalar
    FourVector  operator/ (const double rhs) const {
        return FourVector( x()/rhs, y()/rhs, z()/rhs, t()/rhs );
    }

    /// Arithmetic operator +=
    void operator += (const FourVector& rhs) {
        setX(x() + rhs.x());
        setY(y() + rhs.y());
        setZ(z() + rhs.z());
        setT(t() + rhs.t());
    }
    /// Arithmetic operator -=
    void operator -= (const FourVector& rhs) {
        setX(x() - rhs.x());
        setY(y() - rhs.y());
        setZ(z() - rhs.z());
        setT(t() - rhs.t());
    }
    /// Arithmetic operator *= by scalar
    void operator *= (const double rhs) {
        setX(x()*rhs);
        setY(y()*rhs);
        setZ(z()*rhs);
        setT(t()*rhs);
    }
    /// Arithmetic operator /= by scalar
    void operator /= (const double rhs) {
        setX(x()/rhs);
        setY(y()/rhs);
        setZ(z()/rhs);
        setT(t()/rhs);
    }

    //@}


    /// Static null FourVector = (0,0,0,0)
    static const FourVector& ZERO_VECTOR() {
        static const FourVector v;
        return v;
    }


private:

    double m_v1; ///< px or x. Interpretation depends on accessors used
    double m_v2; ///< py or y. Interpretation depends on accessors used
    double m_v3; ///< pz or z. Interpretation depends on accessors used
    double m_v4; ///< e  or t. Interpretation depends on accessors used

};


/// @name Unbound vector comparison functions
//@{

/// Signed azimuthal angle separation in [-pi, pi] between vecs @c a and @c b
inline double delta_phi(const FourVector &a, const FourVector &b) { return b.delta_phi(a); }

/// Pseudorapidity separation between vecs @c a and @c b
inline double delta_eta(const FourVector &a, const FourVector &b) { return b.delta_eta(a); }

/// Rapidity separation between vecs @c a and @c b
inline double delta_rap(const FourVector &a, const FourVector &b) { return b.delta_rap(a); }

/// R_eta^2-distance separation dR^2 = dphi^2 + deta^2 between vecs @c a and @c b
inline double delta_r2_eta(const FourVector &a, const FourVector &b) { return b.delta_r2_eta(a); }

/// R_eta-distance separation dR = sqrt(dphi^2 + deta^2) between vecs @c a and @c b
inline double delta_r_eta(const FourVector &a, const FourVector &b) { return b.delta_r_eta(a); }

/// R_rap^2-distance separation dR^2 = dphi^2 + drap^2 between vecs @c a and @c b
inline double delta_r2_rap(const FourVector &a, const FourVector &b) { return b.delta_r2_rap(a); }

/// R_rap-distance separation dR = sqrt(dphi^2 + drap^2) between vecs @c a and @c b
inline double delta_r_rap(const FourVector &a, const FourVector &b) { return b.delta_r_rap(a); }

//@}


} // namespace HepMC3


#endif
