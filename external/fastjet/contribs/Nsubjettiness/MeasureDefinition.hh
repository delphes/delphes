//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: MeasureDefinition.hh 828 2015-07-20 14:52:06Z jthaler $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef __FASTJET_CONTRIB_MEASUREDEFINITION_HH__
#define __FASTJET_CONTRIB_MEASUREDEFINITION_HH__

#include "fastjet/PseudoJet.hh"
#include <cmath>
#include <vector>
#include <list>
#include <limits>

#include "TauComponents.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{



// The following Measures are available (and their relevant arguments):
// Recommended for usage as jet shapes
class DefaultMeasure;               // Default measure from which next classes derive (should not be called directly)
class NormalizedMeasure;            // (beta,R0)
class UnnormalizedMeasure;          // (beta)
class NormalizedCutoffMeasure;      // (beta,R0,Rcutoff)
class UnnormalizedCutoffMeasure;    // (beta,Rcutoff)

// New measures as of v2.2
// Recommended for usage as event shapes (or for jet finding)
class ConicalMeasure;               // (beta,R)
class OriginalGeometricMeasure;     // (R)
class ModifiedGeometricMeasure;     // (R)
class ConicalGeometricMeasure;      // (beta, gamma, R)
class XConeMeasure;                 // (beta, R)
   
// Formerly GeometricMeasure, now no longer recommended, kept commented out only for cross-check purposes
//class DeprecatedGeometricMeasure;         // (beta)
//class DeprecatedGeometricCutoffMeasure;   // (beta,Rcutoff)

   
///////
//
// MeasureDefinition
//
///////

//This is a helper class for the Minimum Axes Finders. It is defined later.
class LightLikeAxis;                                          
   
///------------------------------------------------------------------------
/// \class MeasureDefinition
/// \brief Base class for measure definitions
///
/// This is the base class for measure definitions.  Derived classes will calculate
/// the tau_N of a jet given a specific measure and a set of axes.  The measure is
/// determined by various jet and beam distances (and possible normalization factors).
///------------------------------------------------------------------------
class MeasureDefinition {
   
public:
   
   /// Description of measure and parameters
   virtual std::string description() const = 0;
   
   /// In derived classes, this should return a copy of the corresponding derived class
   virtual MeasureDefinition* create() const = 0;

   //The following five functions define the measure by which tau_N is calculated,
   //and are overloaded by the various measures below
   
   /// Distanes to jet axis.  This is called many times, so needs to be as fast as possible
   /// Unless overloaded, it just calls jet_numerator
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return jet_numerator(particle,axis);
   }

   /// Distanes to beam.  This is called many times, so needs to be as fast as possible
   /// Unless overloaded, it just calls beam_numerator
   virtual double beam_distance_squared(const fastjet::PseudoJet& particle) const {
      return beam_numerator(particle);
   }
   
   /// The jet measure used in N-(sub)jettiness
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const = 0;
   /// The beam measure used in N-(sub)jettiness
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const = 0;
   
   /// A possible normalization factor
   virtual double denominator(const fastjet::PseudoJet& particle) const = 0;
   
   /// Run a one-pass minimization routine.  There is a generic one-pass minimization that works for a wide variety of measures.
   /// This should be overloaded to create a measure-specific minimization scheme
   virtual std::vector<fastjet::PseudoJet> get_one_pass_axes(int n_jets,
                                                             const std::vector<fastjet::PseudoJet>& inputs,
                                                             const std::vector<fastjet::PseudoJet>& seedAxes,
                                                             int nAttempts = 1000,         // cap number of iterations
                                                             double accuracy = 0.0001      // cap distance of closest approach
   ) const;
   
public:
   
   /// Returns the tau value for a give set of particles and axes
   double result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const {
      return component_result(particles,axes).tau();
   }
   
   /// Short-hand for the result() function
   inline double operator() (const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const {
      return result(particles,axes);
   }
   
   /// Return all of the TauComponents for specific input particles and axes
   TauComponents component_result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const;
   
   /// Create the partitioning according the jet/beam distances and stores them a TauPartition
   TauPartition get_partition(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const;

   /// Calculate the tau result using an existing partition
   TauComponents component_result_from_partition(const TauPartition& partition, const std::vector<fastjet::PseudoJet>& axes) const;

   

   /// virtual destructor
   virtual ~MeasureDefinition(){}
   
protected:
   
   /// Flag set by derived classes to choose whether or not to use beam/denominator
   TauMode _tau_mode;
   
   /// Flag set by derived classes to say whether cheap get_one_pass_axes method can be used (true by default)
   bool _useAxisScaling;

   /// This is the only constructor, which requires _tau_mode and _useAxisScaling to be manually set by derived classes.
   MeasureDefinition() : _tau_mode(UNDEFINED_SHAPE), _useAxisScaling(true) {}

   
   /// Used by derived classes to set whether or not to use beam/denominator information
   void setTauMode(TauMode tau_mode) {
      _tau_mode = tau_mode;
   }

   /// Used by derived classes to say whether one can use cheap get_one_pass_axes
   void setAxisScaling(bool useAxisScaling) {
      _useAxisScaling = useAxisScaling;
   }

   /// Uses denominator information?
   bool has_denominator() const { return (_tau_mode == NORMALIZED_JET_SHAPE || _tau_mode == NORMALIZED_EVENT_SHAPE);}
   /// Uses beam information?
   bool has_beam() const {return (_tau_mode == UNNORMALIZED_EVENT_SHAPE || _tau_mode == NORMALIZED_EVENT_SHAPE);}

   /// Create light-like axis (used in default one-pass minimization routine)
   fastjet::PseudoJet lightFrom(const fastjet::PseudoJet& input) const {
      double length = sqrt(pow(input.px(),2) + pow(input.py(),2) + pow(input.pz(),2));
      return fastjet::PseudoJet(input.px()/length,input.py()/length,input.pz()/length,1.0);
   }
   
   /// Shorthand for squaring
   static inline double sq(double x) {return x*x;}

};
   

///////
//
// Default Measures
//
///////
   
   
///------------------------------------------------------------------------
/// \enum DefaultMeasureType
/// \brief Options for default measure
///
/// Can be used to switch between pp and ee measure types in the DefaultMeasure
///------------------------------------------------------------------------
enum DefaultMeasureType {
   pt_R,       ///  use transverse momenta and boost-invariant angles,
   E_theta,    ///  use energies and angles,
   lorentz_dot, ///  use dot product inspired measure
   perp_lorentz_dot /// use conical geometric inspired measures
};

///------------------------------------------------------------------------
/// \class DefaultMeasure
/// \brief Base class for default N-subjettiness measure definitions
///
/// This class is the default measure as defined in the original N-subjettiness papers.
/// Based on the conical measure, but with a normalization factor
/// This measure is defined as the pT of the particle multiplied by deltaR
/// to the power of beta. This class includes the normalization factor determined by R0
///------------------------------------------------------------------------
class DefaultMeasure : public MeasureDefinition {
   
public:

   /// Description
   virtual std::string description() const;
   /// To allow copying around of these objects
   virtual DefaultMeasure* create() const {return new DefaultMeasure(*this);}

   /// fast jet distance
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return angleSquared(particle, axis);
   }
   
   /// fast beam distance
   virtual double beam_distance_squared(const fastjet::PseudoJet& /*particle*/) const {
      return _RcutoffSq;
   }
   
   /// true jet distance (given by general definitions of energy and angle)
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const{
      double jet_dist = angleSquared(particle, axis);
      if (jet_dist > 0.0) {
         return energy(particle) * std::pow(jet_dist,_beta/2.0);
      } else {
         return 0.0;
      }
   }
   
   /// true beam distance
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      return energy(particle) * std::pow(_Rcutoff,_beta);
   }
   
   /// possible denominator for normalization
   virtual double denominator(const fastjet::PseudoJet& particle) const {
      return energy(particle) * std::pow(_R0,_beta);
   }

   /// Special minimization routine (from v1.0 of N-subjettiness)
   virtual std::vector<fastjet::PseudoJet> get_one_pass_axes(int n_jets,
                                                             const std::vector<fastjet::PseudoJet>& inputs,
                                                             const std::vector<fastjet::PseudoJet>& seedAxes,
                                                             int nAttempts,   // cap number of iterations
                                                             double accuracy  // cap distance of closest approach
                                                             ) const;
   
protected:
   double _beta;      ///< Angular exponent
   double _R0;        ///< Normalization factor
   double _Rcutoff;   ///< Cutoff radius
   double _RcutoffSq; ///< Cutoff radius squared
   DefaultMeasureType _measure_type;  ///< Type of measure used (i.e. pp style or ee style)
   
   
   /// Constructor is protected so that no one tries to call this directly.
   DefaultMeasure(double beta, double R0, double Rcutoff, DefaultMeasureType measure_type = pt_R)
   : MeasureDefinition(), _beta(beta), _R0(R0), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)), _measure_type(measure_type)
   {
      if (beta <= 0) throw Error("DefaultMeasure:  You must choose beta > 0.");
      if (R0 <= 0) throw Error("DefaultMeasure:  You must choose R0 > 0.");
      if (Rcutoff <= 0) throw Error("DefaultMeasure:  You must choose Rcutoff > 0.");
   }
   
   /// Added set measure method in case it becomes useful later
   void setDefaultMeasureType(DefaultMeasureType measure_type) {
      _measure_type = measure_type;
   }

   /// Generalized energy value (determined by _measure_type)
   double energy(const PseudoJet& jet) const;
   /// Generalized angle value (determined by _measure_type)
   double angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const;

   /// Name of _measure_type, so description will include the measure type
   std::string measure_type_name() const {
      if (_measure_type == pt_R) return "pt_R";
      else if (_measure_type == E_theta) return "E_theta";
      else if (_measure_type == lorentz_dot) return "lorentz_dot";
      else if (_measure_type == perp_lorentz_dot) return "perp_lorentz_dot";
      else return "Measure Type Undefined";
   }      

   /// templated for speed (TODO: probably should remove, since not clear that there is a speed gain)
   template <int N> std::vector<LightLikeAxis> UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes,
                                                              const std::vector <fastjet::PseudoJet> & inputJets,
                                                              double accuracy) const;
   
   /// called by get_one_pass_axes to update axes iteratively
   std::vector<LightLikeAxis> UpdateAxes(const std::vector <LightLikeAxis> & old_axes,
                                         const std::vector <fastjet::PseudoJet> & inputJets,
                                         double accuracy) const;
};
   

///------------------------------------------------------------------------
/// \class NormalizedCutoffMeasure
/// \brief Dimensionless default measure, with radius cutoff
///
/// This measure is just a wrapper for DefaultMeasure
///------------------------------------------------------------------------
class NormalizedCutoffMeasure : public DefaultMeasure {

public:

   /// Constructor
   NormalizedCutoffMeasure(double beta, double R0, double Rcutoff, DefaultMeasureType measure_type = pt_R) 
   : DefaultMeasure(beta, R0, Rcutoff, measure_type) {
      setTauMode(NORMALIZED_JET_SHAPE);
   }

   /// Description
   virtual std::string description() const;

   /// For copying purposes
   virtual NormalizedCutoffMeasure* create() const {return new NormalizedCutoffMeasure(*this);}

};

///------------------------------------------------------------------------
/// \class NormalizedMeasure
/// \brief Dimensionless default measure, with no cutoff
///
/// This measure is the same as NormalizedCutoffMeasure, with Rcutoff taken to infinity.
///------------------------------------------------------------------------
class NormalizedMeasure : public NormalizedCutoffMeasure {

public:

   /// Constructor
   NormalizedMeasure(double beta, double R0, DefaultMeasureType measure_type = pt_R)
   : NormalizedCutoffMeasure(beta, R0, std::numeric_limits<double>::max(), measure_type) {
      _RcutoffSq = std::numeric_limits<double>::max();
      setTauMode(NORMALIZED_JET_SHAPE);
   }
   
   /// Description
   virtual std::string description() const;
   /// For copying purposes
   virtual NormalizedMeasure* create() const {return new NormalizedMeasure(*this);}

};
   

///------------------------------------------------------------------------
/// \class UnnormalizedCutoffMeasure
/// \brief Dimensionful default measure, with radius cutoff
///
/// This class is the unnormalized conical measure. The only difference from NormalizedCutoffMeasure
/// is that the denominator is defined to be 1.0 by setting _has_denominator to false.
/// class UnnormalizedCutoffMeasure : public NormalizedCutoffMeasure {
///------------------------------------------------------------------------
class UnnormalizedCutoffMeasure : public DefaultMeasure {
   
public:
   
   /// Since all methods are identical, UnnormalizedMeasure inherits directly
   /// from NormalizedMeasure. R0 is a dummy value since the value of R0 is unecessary for this class,
   /// and the "false" flag sets _has_denominator in MeasureDefinition to false so no denominator is used.
   UnnormalizedCutoffMeasure(double beta, double Rcutoff, DefaultMeasureType measure_type = pt_R)
   : DefaultMeasure(beta, std::numeric_limits<double>::quiet_NaN(), Rcutoff, measure_type) {
      setTauMode(UNNORMALIZED_EVENT_SHAPE);
   }

   /// Description
   virtual std::string description() const;
   /// For copying purposes
   virtual UnnormalizedCutoffMeasure* create() const {return new UnnormalizedCutoffMeasure(*this);}

};

   
///------------------------------------------------------------------------
/// \class UnnormalizedMeasure
/// \brief Dimensionless default measure, with no cutoff
///
/// This measure is the same as UnnormalizedCutoffMeasure, with Rcutoff taken to infinity.
///------------------------------------------------------------------------
class UnnormalizedMeasure : public UnnormalizedCutoffMeasure {
   
public:
   /// Since all methods are identical, UnnormalizedMeasure inherits directly
   /// from NormalizedMeasure. R0 is a dummy value since the value of R0 is unecessary for this class,
   /// and the "false" flag sets _has_denominator in MeasureDefinition to false so no denominator is used.
   UnnormalizedMeasure(double beta, DefaultMeasureType measure_type = pt_R)
   : UnnormalizedCutoffMeasure(beta, std::numeric_limits<double>::max(), measure_type) {
      _RcutoffSq = std::numeric_limits<double>::max();
      setTauMode(UNNORMALIZED_JET_SHAPE);
   }

   /// Description
   virtual std::string description() const;
   
   /// For copying purposes
   virtual UnnormalizedMeasure* create() const {return new UnnormalizedMeasure(*this);}
   
};


///------------------------------------------------------------------------
/// \class ConicalMeasure
/// \brief Dimensionful event-shape measure, with radius cutoff
///
/// Very similar to UnnormalizedCutoffMeasure, but with different normalization convention
/// and using the new default one-pass minimization algorithm.
/// Axes are also made to be light-like to ensure sensible behavior
/// Intended to be used as an event shape.
///------------------------------------------------------------------------
class ConicalMeasure : public MeasureDefinition {
   
public:
   
   /// Constructor
   ConicalMeasure(double beta, double Rcutoff)
   :   MeasureDefinition(), _beta(beta), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
      if (beta <= 0) throw Error("ConicalMeasure:  You must choose beta > 0.");
      if (Rcutoff <= 0) throw Error("ConicalMeasure:  You must choose Rcutoff > 0.");
      setTauMode(UNNORMALIZED_EVENT_SHAPE);
   }
   
   /// Description
   virtual std::string description() const;
   /// For copying purposes
   virtual ConicalMeasure* create() const {return new ConicalMeasure(*this);}
   
   /// fast jet distance
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      PseudoJet lightAxis = lightFrom(axis);
      return particle.squared_distance(lightAxis);
   }
   
   /// fast beam distance
   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      return _RcutoffSq;
   }
   
   
   /// true jet distance
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      PseudoJet lightAxis = lightFrom(axis);
      double jet_dist = particle.squared_distance(lightAxis)/_RcutoffSq;
      double jet_perp = particle.perp();
      
      if (_beta == 2.0) {
         return jet_perp * jet_dist;
      } else {
         return jet_perp * pow(jet_dist,_beta/2.0);
      }
   }
   
   /// true beam distance
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      return particle.perp();
   }
   
   /// no denominator used for this measure
   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
   
protected:
   double _beta;     ///< angular exponent
   double _Rcutoff;  ///< effective jet radius
   double _RcutoffSq;///< effective jet radius squared
};
   


///------------------------------------------------------------------------
/// \class OriginalGeometricMeasure
/// \brief Dimensionful event-shape measure, with dot-product distances
///
/// This class is the original (and hopefully now correctly coded) geometric measure.
/// This measure is defined by the Lorentz dot product between
/// the particle and the axis.  This class does not include normalization of tau_N.
/// New in Nsubjettiness v2.2
/// NOTE: This is defined differently from the DeprecatedGeometricMeasure which are now commented out.
///------------------------------------------------------------------------
class OriginalGeometricMeasure : public MeasureDefinition {

public:
   /// Constructor
   OriginalGeometricMeasure(double Rcutoff)
   : MeasureDefinition(), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
      if (Rcutoff <= 0) throw Error("OriginalGeometricMeasure:  You must choose Rcutoff > 0.");
      setTauMode(UNNORMALIZED_EVENT_SHAPE);
      setAxisScaling(false);  // No need to rescale axes (for speed up in one-pass minimization)
   }

   /// Description
   virtual std::string description() const;
   /// For copying purposes
   virtual OriginalGeometricMeasure* create() const {return new OriginalGeometricMeasure(*this);}
   
   // This class uses the default jet_distance_squared and beam_distance_squared

   /// true jet measure
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return dot_product(lightFrom(axis), particle)/_RcutoffSq;
   }

   /// true beam measure
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      fastjet::PseudoJet beam_a(0,0,1,1);
      fastjet::PseudoJet beam_b(0,0,-1,1);
      double min_perp = std::min(dot_product(beam_a, particle),dot_product(beam_b, particle));
      return min_perp;
   }

   /// no denominator needed for this measure.
   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
   
protected:
   double _Rcutoff;   ///< Effective jet radius (rho = R^2)
   double _RcutoffSq; ///< Effective jet radius squared

};


///------------------------------------------------------------------------
/// \class ModifiedGeometricMeasure
/// \brief Dimensionful event-shape measure, with dot-product distances, modified beam measure
///
/// This class is the Modified geometric measure.  This jet measure is defined by the Lorentz dot product between
/// the particle and the axis, as in the Original Geometric Measure. The beam measure is defined differently from
/// the above OriginalGeometric to allow for more conical jets. New in Nsubjettiness v2.2
///------------------------------------------------------------------------
class ModifiedGeometricMeasure : public MeasureDefinition {

public:
   /// Constructor
   ModifiedGeometricMeasure(double Rcutoff)
   :  MeasureDefinition(), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)) {
      if (Rcutoff <= 0) throw Error("ModifiedGeometricMeasure:  You must choose Rcutoff > 0.");
      setTauMode(UNNORMALIZED_EVENT_SHAPE);
      setAxisScaling(false);  // No need to rescale axes (for speed up in one-pass minimization)
   }

   /// Description
   virtual std::string description() const;
   /// For copying purposes
   virtual ModifiedGeometricMeasure* create() const {return new ModifiedGeometricMeasure(*this);}
   
   // This class uses the default jet_distance_squared and beam_distance_squared

   /// True jet measure
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return dot_product(lightFrom(axis), particle)/_RcutoffSq;
   }

   /// True beam measure
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      fastjet::PseudoJet lightParticle = lightFrom(particle);
      return 0.5*particle.mperp()*lightParticle.pt();
   }

   /// This measure does not require a denominator
   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
   
protected:
   double _Rcutoff;   ///< Effective jet radius (rho = R^2)
   double _RcutoffSq; ///< Effective jet radius squared


};

///------------------------------------------------------------------------
/// \class ConicalGeometricMeasure
/// \brief Dimensionful event-shape measure, basis for XCone jet algorithm
///
/// This class is the Conical Geometric measure.  This measure is defined by the Lorentz dot product between
/// the particle and the axis normalized by the axis and particle pT, as well as a factor of cosh(y) to vary
/// the rapidity depepdence of the beam. New in Nsubjettiness v2.2, and the basis for the XCone jet algorithm
///------------------------------------------------------------------------
class ConicalGeometricMeasure : public MeasureDefinition {

public:
   
   /// Constructor
   ConicalGeometricMeasure(double jet_beta, double beam_gamma, double Rcutoff)
   :  MeasureDefinition(),  _jet_beta(jet_beta), _beam_gamma(beam_gamma), _Rcutoff(Rcutoff), _RcutoffSq(sq(Rcutoff)){
      if (jet_beta <= 0)   throw Error("ConicalGeometricMeasure:  You must choose beta > 0.");
      if (beam_gamma <= 0) throw Error("ConicalGeometricMeasure:  You must choose gamma > 0.");
      if (Rcutoff <= 0)   throw Error("ConicalGeometricMeasure:  You must choose Rcutoff > 0.");
      setTauMode(UNNORMALIZED_EVENT_SHAPE);
   }

   /// Description
   virtual std::string description() const;
   /// For copying purposes
   virtual ConicalGeometricMeasure* create() const {return new ConicalGeometricMeasure(*this);}
   
   /// fast jet measure
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
      return pseudoRsquared;
   }

   /// fast beam measure
   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      return _RcutoffSq;
   }

   /// true jet measure
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      double jet_dist = jet_distance_squared(particle,axis)/_RcutoffSq;
      if (jet_dist > 0.0) {
         fastjet::PseudoJet lightParticle = lightFrom(particle);
         double weight = (_beam_gamma == 1.0) ? 1.0 : std::pow(0.5*lightParticle.pt(),_beam_gamma - 1.0);
         return particle.pt() * weight * std::pow(jet_dist,_jet_beta/2.0);
      } else {
         return 0.0;
      }
   }

   /// true beam measure
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      fastjet::PseudoJet lightParticle = lightFrom(particle);
      double weight = (_beam_gamma == 1.0) ? 1.0 : std::pow(0.5*lightParticle.pt(),_beam_gamma - 1.0);
      return particle.pt() * weight;
   }

   /// no denominator needed
   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
      
protected:
   double _jet_beta;   ///< jet angular exponent
   double _beam_gamma; ///< beam angular exponent (gamma = 1.0 is recommended)
   double _Rcutoff;    ///< effective jet radius
   double _RcutoffSq;  ///< effective jet radius squared
   
};

///------------------------------------------------------------------------
/// \class XConeMeasure
/// \brief Dimensionful event-shape measure used in XCone jet algorithm
///
/// This class is the XCone Measure.  This is the default measure for use with the
/// XCone algorithm. It is identical to the conical geometric measure but with gamma = 1.0.
///------------------------------------------------------------------------
class XConeMeasure : public ConicalGeometricMeasure {

public:
   /// Constructor
   XConeMeasure(double jet_beta, double R)
   :   ConicalGeometricMeasure(jet_beta,
                               1.0, // beam_gamma, hard coded at gamma = 1.0 default
                               R    // Rcutoff scale
                               ) { }

   /// Description
   virtual std::string description() const;
   /// For copying purposes
   virtual XConeMeasure* create() const {return new XConeMeasure(*this);}

};

///------------------------------------------------------------------------
/// \class LightLikeAxis
/// \brief Helper class to define light-like axes directions
///
/// This is a helper class for the minimization routines.
/// It creates a convenient way of defining axes in order to better facilitate calculations.
///------------------------------------------------------------------------
class LightLikeAxis {

public:
   /// Bare constructor
   LightLikeAxis() : _rap(0.0), _phi(0.0), _weight(0.0), _mom(0.0) {}
   /// Constructor
   LightLikeAxis(double my_rap, double my_phi, double my_weight, double my_mom) :
   _rap(my_rap), _phi(my_phi), _weight(my_weight), _mom(my_mom) {}
   
   /// Rapidity
   double rap() const {return _rap;}
   /// Azimuth
   double phi() const {return _phi;}
   /// weight factor
   double weight() const {return _weight;}
   /// pt momentum
   double mom() const {return _mom;}
   
   /// set rapidity
   void set_rap(double my_set_rap) {_rap = my_set_rap;}
   /// set azimuth
   void set_phi(double my_set_phi) {_phi = my_set_phi;}
   /// set weight factor
   void set_weight(double my_set_weight) {_weight = my_set_weight;}
   /// set pt momentum
   void set_mom(double my_set_mom) {_mom = my_set_mom;}
   /// set all kinematics
   void reset(double my_rap, double my_phi, double my_weight, double my_mom) {_rap=my_rap; _phi=my_phi; _weight=my_weight; _mom=my_mom;}

   /// Return PseudoJet version
   fastjet::PseudoJet ConvertToPseudoJet();
   
   /// Squared distance to PseudoJet
   double DistanceSq(const fastjet::PseudoJet& input) const {
      return DistanceSq(input.rap(),input.phi());
   }

   /// Distance to PseudoJet
   double Distance(const fastjet::PseudoJet& input) const {
      return std::sqrt(DistanceSq(input));
   }

   /// Squared distance to Lightlikeaxis
   double DistanceSq(const LightLikeAxis& input) const {
      return DistanceSq(input.rap(),input.phi());
   }

   /// Distance to Lightlikeaxis
   double Distance(const LightLikeAxis& input) const {
      return std::sqrt(DistanceSq(input));
   }

private:
   double _rap;    ///< rapidity
   double _phi;    ///< azimuth
   double _weight; ///< weight factor
   double _mom;    ///< pt momentum
   
   /// Internal squared distance calculation
   double DistanceSq(double rap2, double phi2) const {
      double rap1 = _rap;
      double phi1 = _phi;
      
      double distRap = rap1-rap2;
      double distPhi = std::fabs(phi1-phi2);
      if (distPhi > M_PI) {distPhi = 2.0*M_PI - distPhi;}
      return distRap*distRap + distPhi*distPhi;
   }
   
   /// Internal distance calculation
   double Distance(double rap2, double phi2) const {
      return std::sqrt(DistanceSq(rap2,phi2));
   }
      
};
   
   
////------------------------------------------------------------------------
///// \class DeprecatedGeometricCutoffMeasure
//// This class is the old, incorrectly coded geometric measure.
//// It is kept in case anyone wants to check old code, but should not be used for production purposes.
//class DeprecatedGeometricCutoffMeasure : public MeasureDefinition {
//
//public:
//
//   // Please, please don't use this.
//   DeprecatedGeometricCutoffMeasure(double jet_beta, double Rcutoff)
//   :   MeasureDefinition(),
//      _jet_beta(jet_beta),
//      _beam_beta(1.0), // This is hard coded, since alternative beta_beam values were never checked.
//      _Rcutoff(Rcutoff),
//      _RcutoffSq(sq(Rcutoff)) {
//         setTauMode(UNNORMALIZED_EVENT_SHAPE);
//         setAxisScaling(false);
//         if (jet_beta != 2.0) {
//         throw Error("Geometric minimization is currently only defined for beta = 2.0.");
//      }
//   }
//
//   virtual std::string description() const;
//
//   virtual DeprecatedGeometricCutoffMeasure* create() const {return new DeprecatedGeometricCutoffMeasure(*this);}
//
//   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
//      fastjet::PseudoJet lightAxis = lightFrom(axis);
//      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
//      return pseudoRsquared;
//   }
//
//   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
//      return _RcutoffSq;
//   }
//
//   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
//      fastjet::PseudoJet lightAxis = lightFrom(axis);
//      double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(lightAxis.pt(),_beam_beta - 1.0);
//      return particle.pt() * weight * std::pow(jet_distance_squared(particle,axis),_jet_beta/2.0);
//   }
//
//   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
//      double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(particle.pt()/particle.e(),_beam_beta - 1.0);
//      return particle.pt() * weight * std::pow(_Rcutoff,_jet_beta);
//   }
//
//   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
//      return std::numeric_limits<double>::quiet_NaN();
//   }
//
//   virtual std::vector<fastjet::PseudoJet> get_one_pass_axes(int n_jets,
//                                                             const std::vector<fastjet::PseudoJet>& inputs,
//                                                             const std::vector<fastjet::PseudoJet>& seedAxes,
//                                                             int nAttempts,    // cap number of iterations
//                                                             double accuracy   // cap distance of closest approach
//                                                            ) const;
//
//protected:
//   double _jet_beta;
//   double _beam_beta;
//   double _Rcutoff;
//   double _RcutoffSq;
//
//};
//
//// ------------------------------------------------------------------------
//// / \class DeprecatedGeometricMeasure
//// Same as DeprecatedGeometricMeasureCutoffMeasure, but with Rcutoff taken to infinity.
//// NOTE:  This class should not be used for production purposes.
//class DeprecatedGeometricMeasure : public DeprecatedGeometricCutoffMeasure {
//
//public:
//   DeprecatedGeometricMeasure(double beta)
//   : DeprecatedGeometricCutoffMeasure(beta,std::numeric_limits<double>::max()) {
//      _RcutoffSq = std::numeric_limits<double>::max();
//      setTauMode(UNNORMALIZED_JET_SHAPE);
//   }
//
//   virtual std::string description() const;
//
//   virtual DeprecatedGeometricMeasure* create() const {return new DeprecatedGeometricMeasure(*this);}
//};


} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_MEASUREDEFINITION_HH__
