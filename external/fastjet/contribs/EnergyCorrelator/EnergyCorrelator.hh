#ifndef __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__
#define __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__

//  EnergyCorrelator Package
//  Questions/Comments?  Email the authors.
//    larkoski@mit.edu, lnecib@mit.edu,
//    gavin.salam@cern.ch jthaler@jthaler.net
//
//  Copyright (c) 2013-2016
//  Andrew Larkoski, Lina Necib Gavin Salam, and Jesse Thaler
//
//  $Id: EnergyCorrelator.hh 1098 2018-01-07 20:17:52Z linoush $
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

#include <fastjet/internal/base.hh>
#include "fastjet/FunctionOfPseudoJet.hh"

#include <string>
#include <climits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

/// \mainpage EnergyCorrelator contrib
///
/// The EnergyCorrelator contrib provides an implementation of energy
/// correlators and their ratios as described in arXiv:1305.0007 by
/// Larkoski, Salam and Thaler.  Additionally, the ratio observable
/// D2 described in arXiv:1409.6298 by Larkoski, Moult and Neill
/// is also included in this contrib. Finally, a generalized version of
/// the energy correlation functions is added, defined in
/// arXiv:1609.07483 by Moult, Necib and Thaler, which allow the
/// definition of the M series, N series, and U series observables.
/// There is also a generalized version of D2.
///
///
/// <p>There are 4 main classes:
///
/// - EnergyCorrelator
/// - EnergyCorrelatorRatio
/// - EnergyCorrelatorDoubleRatio
/// - EnergyCorrelatorGeneralized
///
/// <p>There are five classes that define useful combinations of the ECFs.
///
/// - EnergyCorrelatorNseries
/// - EnergyCorrelatorMseries
/// - EnergyCorrelatorUseries
/// - EnergyCorrelatorD2
/// - EnergyCorrelatorGeneralizedD2
///
/// <p> There are also aliases for easier access:
/// - EnergyCorrelatorCseries (same as EnergyCorrelatorDoubleRatio)
/// - EnergyCorrelatorC1      (EnergyCorrelatorCseries with i=1)
/// - EnergyCorrelatorC2      (EnergyCorrelatorCseries with i=2)
/// - EnergyCorrelatorN2      (EnergyCorrelatorNseries with i=2)
/// - EnergyCorrelatorN3      (EnergyCorrelatorNseries with i=3)
/// - EnergyCorrelatorM2      (EnergyCorrelatorMseries with i=2)
/// - EnergyCorrelatorU1      (EnergyCorrelatorUseries with i=1)
/// - EnergyCorrelatorU2      (EnergyCorrelatorUseries with i=2)
/// - EnergyCorrelatorU3      (EnergyCorrelatorUseries with i=3)
///
/// Each of these classes is a FastJet FunctionOfPseudoJet.
/// EnergyCorrelatorDoubleRatio (which is equivalent to EnergyCorrelatorCseries)
/// is in particular is useful for quark/gluon discrimination and boosted
/// object tagging.
///
/// Using the original 2- and 3-point correlators, EnergyCorrelationD2 has
/// been shown to be the optimal combination for boosted 2-prong tagging.
///
/// The EnergyCorrelatorNseries and EnergyCorrelatorMseries use
/// generalized correlation functions with different angular scaling,
/// and are intended for use on 2-prong and 3-prong jets.
/// The EnergyCorrelatorUseries is useful for quark/gluon discrimimation.
///
/// See the file example.cc for an illustration of usage and
/// example_basic_usage.cc for the most commonly used functions.

//------------------------------------------------------------------------
/// \class EnergyCorrelator
/// ECF(N,beta) is the N-point energy correlation function, with an angular exponent beta.
///
/// It is defined as follows
///
///  - \f$ \mathrm{ECF}(1,\beta)  = \sum_i E_i \f$
///  - \f$ \mathrm{ECF}(2,\beta)  = \sum_{i<j} E_i E_j \theta_{ij}^\beta \f$
///  - \f$ \mathrm{ECF}(3,\beta)  = \sum_{i<j<k} E_i E_j E_k (\theta_{ij} \theta_{ik} \theta_{jk})^\beta \f$
///  - \f$ \mathrm{ECF}(4,\beta)  = \sum_{i<j<k<l} E_i E_j E_k E_l (\theta_{ij}  \theta_{ik} \theta_{il} \theta_{jk} \theta_{jl} \theta_{kl})^\beta \f$
///  - ...
///
/// The correlation can be determined with energies and angles (as
/// given above) or with transverse momenta and boost invariant angles
/// (the code's default). The choice is controlled by
/// EnergyCorrelator::Measure provided in the constructor.
///
/// The current implementation handles values of N up to and including 5.
/// Run times scale as n^N/N!, where n is the number of particles in a jet.

    class EnergyCorrelator : public FunctionOfPseudoJet<double> {
        friend class EnergyCorrelatorGeneralized;  ///< This allow ECFG to access the energy and angle definitions
                                                  ///< of this class, which are otherwise private.
    public:

        enum Measure {
            pt_R,     ///< use transverse momenta and boost-invariant angles,
            ///< eg \f$\mathrm{ECF}(2,\beta) = \sum_{i<j} p_{ti} p_{tj} \Delta R_{ij}^{\beta} \f$
            E_theta,   ///  use energies and angles,
            ///  eg \f$\mathrm{ECF}(2,\beta) = \sum_{i<j} E_{i} E_{j}   \theta_{ij}^{\beta} \f$
            E_inv     ///  use energies and invariant mass,
            ///  eg \f$\mathrm{ECF}(2,\beta) = \sum_{i<j} E_{i} E_{j}   (\frac{2 p_{i} \cdot p_{j}}{E_{i} E_{j}})^{\beta/2} \f$
        };

        enum Strategy {
            slow,          ///< interparticle angles are not cached.
            ///< For N>=3 this leads to many expensive recomputations,
            ///< but has only O(n) memory usage for n particles

            storage_array  /// the interparticle angles are cached. This gives a significant speed
            /// improvement for N>=3, but has a memory requirement of (4n^2) bytes.
        };

    public:

        /// constructs an N-point correlator with angular exponent beta,
        /// using the specified choice of energy and angular measure as well
        /// one of two possible underlying computational Strategy
        EnergyCorrelator(unsigned int N,
                         double beta,
                         Measure measure = pt_R,
                         Strategy strategy = storage_array) :
                _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};

        /// destructor
        virtual ~EnergyCorrelator(){}

        /// returns the value of the energy correlator for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

        /// returns the the part of the description related to the parameters
        std::string description_parameters() const;
        std::string description_no_N() const;

    private:

        unsigned int _N;
        double _beta;
        Measure _measure;
        Strategy _strategy;

        double energy(const PseudoJet& jet) const;
        double angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const;
        double multiply_angles(double angles[], int n_angles, unsigned int N_total) const;
        void precompute_energies_and_angles(std::vector<fastjet::PseudoJet> const &particles, double* energyStore, double** angleStore) const;
        double evaluate_n3(unsigned int nC, unsigned int n_angles, double* energyStore, double** angleStore) const;
        double evaluate_n4(unsigned int nC, unsigned int n_angles, double* energyStore, double** angleStore) const;
        double evaluate_n5(unsigned int nC, unsigned int n_angles, double* energyStore, double** angleStore) const;
    };

// core EnergyCorrelator::result code in .cc file.


//------------------------------------------------------------------------
/// \class EnergyCorrelatorRatio
/// A class to calculate the ratio of (N+1)-point to N-point energy correlators,
///     ECF(N+1,beta)/ECF(N,beta),
/// called \f$ r_N^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorRatio : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an (N+1)-point to N-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorRatio(unsigned int N,
                              double  beta,
                              EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                              EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorRatio() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        unsigned int _N;
        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };

    inline double EnergyCorrelatorRatio::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelator(_N + 1, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelator(_N, _beta, _measure, _strategy).result(jet);

        return numerator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorDoubleRatio
/// Calculates the double ratio of energy correlators, ECF(N-1,beta)*ECF(N+1)/ECF(N,beta)^2.
///
/// A class to calculate a double ratio of energy correlators,
///     ECF(N-1,beta)*ECF(N+1,beta)/ECF(N,beta)^2,
/// called \f$C_N^{(\beta)}\f$ in the publication, and equal to
/// \f$ r_N^{(\beta)}/r_{N-1}^{(\beta)} \f$.
///

    class EnergyCorrelatorDoubleRatio : public FunctionOfPseudoJet<double> {

    public:

        EnergyCorrelatorDoubleRatio(unsigned int N,
                                    double beta,
                                    EnergyCorrelator::Measure measure = EnergyCorrelator::pt_R,
                                    EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _N(N), _beta(beta), _measure(measure), _strategy(strategy) {

                    if (_N < 1) throw Error("EnergyCorrelatorDoubleRatio:  N must be 1 or greater.");

                };

        virtual ~EnergyCorrelatorDoubleRatio() {}


        /// returns the value of the energy correlator double-ratio for a
        /// jet's constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        unsigned int _N;
        double _beta;
        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorDoubleRatio::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelator(_N - 1, _beta, _measure, _strategy).result(jet) * EnergyCorrelator(_N + 1, _beta, _measure, _strategy).result(jet);
        double denominator = pow(EnergyCorrelator(_N, _beta, _measure, _strategy).result(jet), 2.0);

        return numerator/denominator;

    }

//------------------------------------------------------------------------
/// \class EnergyCorrelatorC1
/// A class to calculate the normalized 2-point energy correlators,
///     ECF(2,beta)/ECF(1,beta)^2,
/// called \f$ C_1^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorC1 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorC1(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorC1() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorC1::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelator(2, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelator(1, _beta, _measure, _strategy).result(jet);

        return numerator/denominator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorC2
/// A class to calculate the double ratio of 3-point to 2-point
/// energy correlators,
///     ECF(3,beta)*ECF(1,beta)/ECF(2,beta)^2,
/// called \f$ C_2^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorC2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a 3-point to 2-point correlator double ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorC2(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorC2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorC2::result(const PseudoJet& jet) const {

        double numerator3 = EnergyCorrelator(3, _beta, _measure, _strategy).result(jet);
        double numerator1 = EnergyCorrelator(1, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelator(2, _beta, _measure, _strategy).result(jet);

        return numerator3*numerator1/denominator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorD2
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     ECF(3,beta)*ECF(1,beta)^3/ECF(2,beta)^3,
/// called \f$ D_2^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorD2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorD2(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorD2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorD2::result(const PseudoJet& jet) const {

        double numerator3 = EnergyCorrelator(3, _beta, _measure, _strategy).result(jet);
        double numerator1 = EnergyCorrelator(1, _beta, _measure, _strategy).result(jet);
        double denominator2 = EnergyCorrelator(2, _beta, _measure, _strategy).result(jet);

        return numerator3*numerator1*numerator1*numerator1/denominator2/denominator2/denominator2;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorGeneralized
/// A generalized and normalized version of the N-point energy correlators, with
/// angular exponent beta and v number of pairwise angles.  When \f$v = {N \choose 2}\f$
/// (or, for convenience, \f$v = -1\f$), EnergyCorrelatorGeneralized just gives normalized
/// versions of EnergyCorrelator:
///  - \f$ \mathrm{ECFG}(-1,1,\beta) = \mathrm{ECFN}(N,\beta) = \mathrm{ECF}(N,\beta)/\mathrm{ECF}(1,\beta)\f$
///
/// Note that there is no separate class that implements ECFN, though it is a
/// notation that we will use in this documentation.  Examples of the low-point normalized
/// correlators are:
///  - \f$\mathrm{ECFN}(1,\beta)  = 1\f$
///  - \f$\mathrm{ECFN}(2,\beta)  = \sum_{i<j} z_i z_j \theta_{ij}^\beta \f$
///  - \f$\mathrm{ECFN}(3,\beta)  = \sum_{i<j<k} z_i z_j z_k (\theta_{ij} \theta_{ik} \theta_{jk})^\beta \f$
///  - \f$\mathrm{ECFN}(4,\beta)  = \sum_{i<j<k<l} z_i z_j z_k z_l (\theta_{ij}  \theta_{ik} \theta_{il} \theta_{jk} \theta_{jl} \theta_{kl})^\beta \f$
///  - ...
/// where the \f$z_i\f$'s are the energy fractions.
///
/// When a new value of v is given, the generalized energy correlators are defined as
///  - \f$\mathrm{ECFG}(0,1,\beta)  = 1\f$
///  - \f$\mathrm{ECFG}(1,2,\beta)  = \sum_{i<j} z_i z_j \theta_{ij}^\beta \f$
///  - \f$\mathrm{ECFG}(v,3,\beta)  = \sum_{i<j<k} z_i z_j z_k \prod_{m = 1}^{v} \min^{(m)} \{ \theta_{ij}, \theta_{jk}, \theta_{ki} \}^\beta \f$
///  - \f$\mathrm{ECFG}(v,4,\beta)  = \sum_{i<j<k<l} z_i z_j z_k z_l \prod_{m = 1}^{v} \min^{(m)} \{ \theta_{ij}, \theta_{ik}, \theta_{il}, \theta_{jk}, \theta_{jl}, \theta_{kl} \}^\beta \f$
///  - \f$\mathrm{ECFG}(v,n,\beta)  = \sum_{i_1 < i_2 < \dots < i_n} z_{i_1} z_{i_2} \dots z_{i_n} \prod_{m = 1}^{v} \min^{(m)}_{s < t \in \{i_1, i_2 , \dots, i_n \}} \{ \theta_{st}^{\beta} \}\f$,
///
/// where \f$\min^{(m)}\f$ means the m-th smallest element of the list.
///
/// The correlation can be determined with energies and angles (as
/// given above) or with transverse momenta and boost invariant angles
/// (the code's default). The choice is controlled by
/// EnergyCorrelator::Measure provided in the constructor.
///
/// The current implementation handles values of N up to and including 5.
///
    class EnergyCorrelatorGeneralized : public FunctionOfPseudoJet<double> {
    public:

        /// constructs an N-point correlator with v_angles pairwise angles
        /// and angular exponent beta,
        /// using the specified choice of energy and angular measure as well
        /// one of two possible underlying computational Strategy
        EnergyCorrelatorGeneralized(int v_angles,
                                    unsigned int N,
                                    double beta,
                                    EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                                    EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                :  _angles(v_angles), _N(N), _beta(beta), _measure(measure), _strategy(strategy),  _helper_correlator(1,_beta, _measure, _strategy) {};

        /// destructor
        virtual ~EnergyCorrelatorGeneralized(){}

        /// returns the value of the normalized energy correlator for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).

        double result(const PseudoJet& jet) const;
        std::vector<double> result_all_angles(const PseudoJet& jet) const;

    private:

        int _angles;
        unsigned int _N;
        double _beta;
        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;
        EnergyCorrelator _helper_correlator;

        double energy(const PseudoJet& jet) const;
        double angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const;
        double multiply_angles(double angles[], int n_angles, unsigned int N_total) const;
        void precompute_energies_and_angles(std::vector<fastjet::PseudoJet> const &particles, double* energyStore, double** angleStore) const;
        double evaluate_n3(unsigned int nC, unsigned int n_angles, double* energyStore, double** angleStore) const;
        double evaluate_n4(unsigned int nC, unsigned int n_angles, double* energyStore, double** angleStore) const;
        double evaluate_n5(unsigned int nC, unsigned int n_angles, double* energyStore, double** angleStore) const;
    };



//------------------------------------------------------------------------
/// \class EnergyCorrelatorGeneralizedD2
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     ECFN(3,alpha)/ECFN(2,beta)^3 alpha/beta,
/// called \f$ D_2^{(\alpha, \beta)} \f$ in the publication.
    class EnergyCorrelatorGeneralizedD2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorGeneralizedD2(
                double alpha,
                double  beta,
                EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _alpha(alpha), _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorGeneralizedD2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _alpha;
        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorGeneralizedD2::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelatorGeneralized(-1, 3, _alpha, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorGeneralized(-1, 2, _beta, _measure, _strategy).result(jet);

        return numerator/pow(denominator, 3.0*_alpha/_beta);

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorNseries
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     N_n = ECFG(2,n+1,beta)/ECFG(1,n,beta)^2,
/// called \f$ N_i^{(\alpha, \beta)} \f$ in the publication.
/// By definition, N_1^{beta} = ECFG(1, 2, 2*beta), where the angular exponent
/// is twice as big since the N series should involve two pairwise angles.
    class EnergyCorrelatorNseries : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a n 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorNseries(
                unsigned int n,
                double  beta,
                EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _n(n), _beta(beta), _measure(measure), _strategy(strategy) {

                    if (_n < 1) throw Error("EnergyCorrelatorNseries:  n must be 1 or greater.");

                };

        virtual ~EnergyCorrelatorNseries() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        unsigned int _n;
        double _beta;
        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;

    };


    inline double EnergyCorrelatorNseries::result(const PseudoJet& jet) const {

        if (_n == 1) return EnergyCorrelatorGeneralized(1, 2, 2*_beta, _measure, _strategy).result(jet);
        // By definition, N1 = ECFN(2, 2 beta)
        double numerator = EnergyCorrelatorGeneralized(2, _n + 1, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorGeneralized(1, _n, _beta, _measure, _strategy).result(jet);

        return numerator/denominator/denominator;

    }



//------------------------------------------------------------------------
/// \class EnergyCorrelatorN2
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     ECFG(2,3,beta)/ECFG(1,2,beta)^2,
/// called \f$ N_2^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorN2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorN2(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorN2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorN2::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelatorGeneralized(2, 3, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorGeneralized(1, 2, _beta, _measure, _strategy).result(jet);

        return numerator/denominator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorN3
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     ECFG(2,4,beta)/ECFG(1,3,beta)^2,
/// called \f$ N_3^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorN3 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorN3(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorN3() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorN3::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelatorGeneralized(2, 4, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorGeneralized(1, 3, _beta, _measure, _strategy).result(jet);

        return numerator/denominator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorMseries
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     M_n = ECFG(1,n+1,beta)/ECFG(1,n,beta),
/// called \f$ M_i^{(\alpha, \beta)} \f$ in the publication.
/// By definition, M_1^{beta} = ECFG(1,2,beta)
    class EnergyCorrelatorMseries : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a n 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorMseries(
                unsigned int n,
                double  beta,
                EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _n(n), _beta(beta), _measure(measure), _strategy(strategy) {

                    if (_n < 1) throw Error("EnergyCorrelatorMseries:  n must be 1 or greater.");

                };

        virtual ~EnergyCorrelatorMseries() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        unsigned int _n;
        double _beta;
        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;

    };


    inline double EnergyCorrelatorMseries::result(const PseudoJet& jet) const {

        if (_n == 1) return EnergyCorrelatorGeneralized(1, 2, _beta, _measure, _strategy).result(jet);

        double numerator = EnergyCorrelatorGeneralized(1, _n + 1, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorGeneralized(1, _n, _beta, _measure, _strategy).result(jet);

        return numerator/denominator;

    }

//------------------------------------------------------------------------
/// \class EnergyCorrelatorM2
/// A class to calculate the observable formed from the ratio of the
/// 3-point and 2-point energy correlators,
///     ECFG(1,3,beta)/ECFG(1,2,beta),
/// called \f$ M_2^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorM2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs an 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorM2(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorM2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorM2::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelatorGeneralized(1, 3, _beta, _measure, _strategy).result(jet);
        double denominator = EnergyCorrelatorGeneralized(1, 2, _beta, _measure, _strategy).result(jet);

        return numerator/denominator;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorCseries
/// Calculates the C series energy correlators, ECFN(N-1,beta)*ECFN(N+1,beta)/ECFN(N,beta)^2.
/// This is equivalent to EnergyCorrelatorDoubleRatio
///
/// A class to calculate a double ratio of energy correlators,
///     ECFN(N-1,beta)*ECFN(N+1,beta)/ECFN(N,beta)^2,
/// called \f$C_N^{(\beta)}\f$ in the publication, and equal to
/// \f$ r_N^{(\beta)}/r_{N-1}^{(\beta)} \f$.
///

    class EnergyCorrelatorCseries : public FunctionOfPseudoJet<double> {

    public:

        EnergyCorrelatorCseries(unsigned int N,
                                    double beta,
                                    EnergyCorrelator::Measure measure = EnergyCorrelator::pt_R,
                                    EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _N(N), _beta(beta), _measure(measure), _strategy(strategy) {

                    if (_N < 1) throw Error("EnergyCorrelatorCseries:  N must be 1 or greater.");

                };

        virtual ~EnergyCorrelatorCseries() {}


        /// returns the value of the energy correlator double-ratio for a
        /// jet's constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        unsigned int _N;
        double _beta;
        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorCseries::result(const PseudoJet& jet) const {

        double numerator = EnergyCorrelatorGeneralized(-1, _N - 1, _beta, _measure, _strategy).result(jet) * EnergyCorrelatorGeneralized(-1, _N + 1, _beta, _measure, _strategy).result(jet);
        double denominator = pow(EnergyCorrelatorGeneralized(-1, _N, _beta, _measure, _strategy).result(jet), 2.0);

        return numerator/denominator;

    }

//------------------------------------------------------------------------
/// \class EnergyCorrelatorUseries
/// A class to calculate the observable used for quark versus gluon discrimination
///     U_n = ECFG(1,n+1,beta),
/// called \f$ U_i^{(\beta)} \f$ in the publication.

    class EnergyCorrelatorUseries : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a n 3-point to 2-point correlator ratio with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorUseries(
                unsigned int n,
                double  beta,
                EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _n(n), _beta(beta), _measure(measure), _strategy(strategy) {

                    if (_n < 1) throw Error("EnergyCorrelatorUseries:  n must be 1 or greater.");

                };

        virtual ~EnergyCorrelatorUseries() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        unsigned int _n;
        double _beta;
        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;

    };


    inline double EnergyCorrelatorUseries::result(const PseudoJet& jet) const {

        double answer = EnergyCorrelatorGeneralized(1, _n + 1, _beta, _measure, _strategy).result(jet);
        return answer;

    }


//------------------------------------------------------------------------
/// \class EnergyCorrelatorU1
/// A class to calculate the observable formed from
///     ECFG(1,2,beta),
/// called \f$ U_1^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorU1 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a 2-point correlator with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorU1(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorU1() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorU1::result(const PseudoJet& jet) const {

        double answer = EnergyCorrelatorGeneralized(1, 2, _beta, _measure, _strategy).result(jet);

        return answer;

    }


    //------------------------------------------------------------------------
    /// \class EnergyCorrelatorU2
    /// A class to calculate the observable formed from
    ///     ECFG(1,3,beta),
    /// called \f$ U_2^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorU2 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a 3-point correlator with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorU2(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorU2() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorU2::result(const PseudoJet& jet) const {

        double answer = EnergyCorrelatorGeneralized(1, 3, _beta, _measure, _strategy).result(jet);

        return answer;

    }


    //------------------------------------------------------------------------
    /// \class EnergyCorrelatorU3
    /// A class to calculate the observable formed from
    ///     ECFG(1,4,beta),
    /// called \f$ U_3^{(\beta)} \f$ in the publication.
    class EnergyCorrelatorU3 : public FunctionOfPseudoJet<double> {

    public:

        /// constructs a 4-point correlator with
        /// angular exponent beta, using the specified choice of energy and
        /// angular measure as well one of two possible underlying
        /// computational strategies
        EnergyCorrelatorU3(double  beta,
                           EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                           EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
                : _beta(beta), _measure(measure), _strategy(strategy) {};

        virtual ~EnergyCorrelatorU3() {}

        /// returns the value of the energy correlator ratio for a jet's
        /// constituents. (Normally accessed by the parent class's
        /// operator()).
        double result(const PseudoJet& jet) const;

        std::string description() const;

    private:

        double _beta;

        EnergyCorrelator::Measure _measure;
        EnergyCorrelator::Strategy _strategy;


    };


    inline double EnergyCorrelatorU3::result(const PseudoJet& jet) const {

        double answer = EnergyCorrelatorGeneralized(1, 4, _beta, _measure, _strategy).result(jet);

        return answer;

    }



} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__
