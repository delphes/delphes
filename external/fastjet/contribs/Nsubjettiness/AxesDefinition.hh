//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: AxesDefinition.hh 833 2015-07-23 14:35:23Z jthaler $
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

#ifndef __FASTJET_CONTRIB_AXES_DEFINITION_HH__
#define __FASTJET_CONTRIB_AXES_DEFINITION_HH__


#include "MeasureDefinition.hh"
#include "ExtraRecombiners.hh"

#include "fastjet/PseudoJet.hh"
#include <fastjet/LimitedWarning.hh>

#include <iomanip>
#include <cmath>
#include <vector>
#include <list>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
   
// The following AxesDefinitions are currently available (and the relevant arguments, if needed)
class KT_Axes;
class CA_Axes;
class AntiKT_Axes;         // (R0)
class WTA_KT_Axes;
class WTA_CA_Axes;
class GenKT_Axes;          // (p, R0 = infinity)
class WTA_GenKT_Axes;      // (p, R0 = infinity)
class GenET_GenKT_Axes;    // (delta, p, R0 = infinity)
class Manual_Axes;
   
class OnePass_KT_Axes;
class OnePass_CA_Axes;
class OnePass_AntiKT_Axes;       // (R0)
class OnePass_WTA_KT_Axes;
class OnePass_WTA_CA_Axes;
class OnePass_GenKT_Axes;        // (p, R0 = infinity)
class OnePass_WTA_GenKT_Axes;    // (p, R0 = infinity)
class OnePass_GenET_GenKT_Axes;  // (delta, p, R0 = infinity)
class OnePass_Manual_Axes;
   
class MultiPass_Axes;            // (NPass) (currently only defined for KT_Axes)
class MultiPass_Manual_Axes;     // (NPass)

class Comb_GenKT_Axes;           // (nExtra, p, R0 = infinity)
class Comb_WTA_GenKT_Axes;       // (nExtra, p, R0 = infinity)
class Comb_GenET_GenKT_Axes;     // (nExtra, delta, p, R0 = infinity)

///////
//
// AxesDefinition
//
///////

///------------------------------------------------------------------------
/// \class AxesDefinition
/// \brief Base class for axes definitions
///
/// A generic AxesDefinition first finds a set of seed axes.
/// Then, if desired, uses measure information
/// (from MeasureDefinition) to refine those axes starting from those seed axes.
/// The AxesDefinitions are typically based on sequential jet algorithms.
///------------------------------------------------------------------------
class AxesDefinition {
   
public:
   
   /// This function should be overloaded in all derived classes, and defines how to find the seed axes.
   /// If desired, the measure information (which might be NULL) can be used to test multiple axes choices, but should
   /// not be used for iterative refining (since that is the job of MeasureDefinition).
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets,
                                                             const std::vector<fastjet::PseudoJet>& inputs,
                                                             const MeasureDefinition * measure) const = 0;
   
   /// Short description of AxesDefinitions (and any parameters)
   virtual std::string short_description() const = 0;

   /// Long description of AxesDefinitions (and any parameters)
   virtual std::string description() const = 0;
   
   /// This has to be defined in all derived classes, and allows these to be copied around.
   virtual AxesDefinition* create() const = 0;
   
public:
   
   /// Starting from seeds, refine axes using one or more passes.
   /// Note that in order to do >0 passes, we need information from the MeasureDefinition about how to do the appropriate minimization.
   std::vector<fastjet::PseudoJet> get_refined_axes(int n_jets,
                                                    const std::vector<fastjet::PseudoJet>& inputs,
                                                    const std::vector<fastjet::PseudoJet>& seedAxes,
                                                    const MeasureDefinition * measure = NULL) const {

      assert(n_jets == (int)seedAxes.size()); //added int casting to get rid of compiler warning
            
      if (_Npass == 0) {
         // no refining, just use seeds
         return seedAxes;
      } else if (_Npass == 1) {
         if (measure == NULL) throw Error("AxesDefinition:  One-pass minimization requires specifying a MeasureDefinition.");
         
         // do one pass minimum using measure definition
         return measure->get_one_pass_axes(n_jets, inputs, seedAxes,_nAttempts,_accuracy);
      } else {
         if (measure == NULL) throw Error("AxesDefinition:  Multi-pass minimization requires specifying a MeasureDefinition.");
         return get_multi_pass_axes(n_jets, inputs, seedAxes, measure);
      }
   }
   
   /// Combines get_starting_axes with get_refined_axes.
   /// In the Njettiness class, these two steps are done separately in order to store seed axes information.
   std::vector<fastjet::PseudoJet> get_axes(int n_jets,
                                            const std::vector<fastjet::PseudoJet>& inputs,
                                            const MeasureDefinition * measure = NULL) const {
      std::vector<fastjet::PseudoJet> seedAxes = get_starting_axes(n_jets, inputs, measure);
      return get_refined_axes(n_jets,inputs,seedAxes,measure);
   }

   
   /// Short-hand for the get_axes function.  Useful when trying to write terse code.
   inline std::vector<fastjet::PseudoJet> operator() (int n_jets,
                                               const std::vector<fastjet::PseudoJet>& inputs,
                                               const MeasureDefinition * measure = NULL) const {
      return get_axes(n_jets,inputs,measure);
   }
   
   /// \enum AxesRefiningEnum
   /// Defines the cases of zero pass and one pass for convenience
   enum AxesRefiningEnum {
      UNDEFINED_REFINE = -1, // added to create a default value
      NO_REFINING = 0,
      ONE_PASS = 1,
      MULTI_PASS = 100,
   };
   
   /// A integer that is used externally to decide how to do multi-pass minimization
   int nPass() const { return _Npass; }

   /// A flag that indicates whether results are deterministics.
   bool givesRandomizedResults() const {
      return (_Npass > 1);
   }
   
   /// A flag that indicates whether manual axes are being used.
   bool needsManualAxes() const {
      return _needsManualAxes; // if there is no starting axes finder
   }
   
   /// Allows user to change number of passes.  Also used internally to set nPass.
   /// Can also specify details of one/multi pass minimziation
   void setNPass(int nPass,
                 int nAttempts = 1000,
                 double accuracy  = 0.0001,
                 double noise_range = 1.0 // only needed for MultiPass minimization
                 )
   {
      _Npass = nPass;
      _nAttempts = nAttempts;
      _accuracy = accuracy;
      _noise_range = noise_range;
      if (nPass < 0) throw Error("AxesDefinition requires a nPass >= 0");
   }
   
   /// Destructor
   virtual ~AxesDefinition() {};
   
protected:
   
   /// Default constructor contains no information.  Number of passes has to be set
   /// manually by derived classes using setNPass function.
   AxesDefinition() : _Npass(UNDEFINED_REFINE),
                     _nAttempts(0),
                     _accuracy(0.0),
                     _noise_range(0.0),
                     _needsManualAxes(false) {}

   /// Does multi-pass minimization by randomly jiggling the axes within _noise_range
   std::vector<fastjet::PseudoJet> get_multi_pass_axes(int n_jets,
                                                       const std::vector<fastjet::PseudoJet>& inputs,
                                                       const std::vector<fastjet::PseudoJet>& seedAxes,
                                                       const MeasureDefinition* measure) const;
   
   /// Function to jiggle axes within _noise_range
   PseudoJet jiggle(const PseudoJet& axis) const;
   
   int _Npass;            ///< Number of passes (0 = no refining, 1 = one-pass, >1 multi-pass)
   int _nAttempts;        ///< Number of attempts per pass
   double _accuracy;      ///< Accuracy goal per pass
   double _noise_range;   ///< Noise in rapidity/phi (for multi-pass minimization only)
   bool _needsManualAxes; ///< Flag to indicate special case of manual axes
};
  
///------------------------------------------------------------------------
/// \class ExclusiveJetAxes
/// \brief Base class for axes defined from exclusive jet algorithm
///
/// This class finds axes by clustering particles with an exclusive jet definition.
/// This can be implemented with different jet algorithms.  The user can call this directly
/// using their favorite fastjet::JetDefinition
///------------------------------------------------------------------------
class ExclusiveJetAxes : public AxesDefinition {
   
public:
   /// Constructor takes JetDefinition as an argument
   ExclusiveJetAxes(fastjet::JetDefinition def)
   : AxesDefinition(), _def(def) {
      setNPass(NO_REFINING);    // default to no minimization
   }
   
   /// Starting axes obtained by creating a cluster sequenence and running exclusive_jets.
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets,
                                                             const std::vector <fastjet::PseudoJet> & inputs,
                                                             const MeasureDefinition * ) const {
      fastjet::ClusterSequence jet_clust_seq(inputs, _def);
      
      std::vector<fastjet::PseudoJet> axes = jet_clust_seq.exclusive_jets_up_to(n_jets);
      
      if ((int)axes.size() < n_jets) {
         _too_few_axes_warning.warn("ExclusiveJetAxes::get_starting_axes:  Fewer than N axes found; results are unpredictable.");
         axes.resize(n_jets);  // resize to make sure there are enough axes to not yield an error elsewhere
      }
      
      return axes;
   }
   
   /// Short description
   virtual std::string short_description() const { return "ExclAxes";}
   /// Long description
   virtual std::string description() const { return "ExclAxes: " + _def.description();}
   
   /// To make it possible to copy around.
   virtual ExclusiveJetAxes* create() const {return new ExclusiveJetAxes(*this);}

private:
   fastjet::JetDefinition _def; ///< Jet definition to use.
   static LimitedWarning _too_few_axes_warning;
};

///------------------------------------------------------------------------
/// \class ExclusiveCombinatorialJetAxes
/// \brief Base class for axes defined from exclusive jet algorithm, checking combinatorial options
///
/// This class finds axes by clustering particles with an exclusive jet definition.
/// It takes an extra number of jets (specificed by the user via nExtra), and then finds the set of N that minimizes N-jettiness.
/// WARNING: If one wants to be guarenteed that results improve by increasing nExtra, then one should use
/// winner-take-all-style recombination schemes
///------------------------------------------------------------------------
class ExclusiveCombinatorialJetAxes : public AxesDefinition {
   
public:
   /// Constructor takes JetDefinition and nExtra as options (nExtra=0 acts the same as ExclusiveJetAxes)
   ExclusiveCombinatorialJetAxes(fastjet::JetDefinition def, int nExtra = 0)
   : AxesDefinition(), _def(def), _nExtra(nExtra) {
      if (nExtra < 0) throw Error("Need nExtra >= 0");
      setNPass(NO_REFINING);   // default to no minimization
   }
   
    /// Find n_jets + _nExtra axes, and then choose the n_jets subset with the smallest N-(sub)jettiness value.
    virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets, 
                                                           const std::vector<fastjet::PseudoJet> & inputs,
                                                           const MeasureDefinition *measure) const {
      int starting_number = n_jets + _nExtra;
      fastjet::ClusterSequence jet_clust_seq(inputs, _def);
      std::vector<fastjet::PseudoJet> starting_axes = jet_clust_seq.exclusive_jets_up_to(starting_number);
       
       if ((int)starting_axes.size() < n_jets) {
          _too_few_axes_warning.warn("ExclusiveCombinatorialJetAxes::get_starting_axes:  Fewer than N + nExtra axes found; results are unpredictable.");
          starting_axes.resize(n_jets);  // resize to make sure there are enough axes to not yield an error elsewhere
       }
       
      std::vector<fastjet::PseudoJet> final_axes;

      // check so that no computation time is wasted if there are no extra axes
      if (_nExtra == 0) final_axes = starting_axes;

      else {
 
        // define string of 1's based on number of desired jets
        std::string bitmask(n_jets, 1);
        // expand the array size to the total number of jets with extra 0's at the end, makes string easy to permute
        bitmask.resize(starting_number, 0); 
 
        double min_tau = std::numeric_limits<double>::max();
        std::vector<fastjet::PseudoJet> temp_axes;
         
        do {

          temp_axes.clear();
 
          // only take an axis if it is listed as true (1) in the string
          for (int i = 0; i < (int)starting_axes.size(); ++i) {
            if (bitmask[i]) temp_axes.push_back(starting_axes[i]);
          }
 
          double temp_tau = measure->result(inputs, temp_axes);
          if (temp_tau < min_tau) {
            min_tau = temp_tau;
            final_axes = temp_axes;
          }
 
          // permutes string of 1's and 0's according to next lexicographic ordering and returns true
          // continues to loop through all possible lexicographic orderings
          // returns false and breaks the loop when there are no more possible orderings
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      }

      return final_axes;
    }

   /// Short description
   virtual std::string short_description() const { return "ExclCombAxes";}
   /// Long description
   virtual std::string description() const { return "ExclCombAxes: " + _def.description();}
   /// To make it possible to copy around.
   virtual ExclusiveCombinatorialJetAxes* create() const {return new ExclusiveCombinatorialJetAxes(*this);}

private:
   fastjet::JetDefinition _def;   ///< Jet definition to use
   int _nExtra;                   ///< Extra axes to find
   static LimitedWarning _too_few_axes_warning;
};
   
///------------------------------------------------------------------------
/// \class HardestJetAxes
/// \brief Base class for axes defined from an inclusive jet algorithm
///
/// This class finds axes by running an inclusive algorithm and then finding the n hardest jets.
/// This can be implemented with different jet algorithms, and can be called by the user.
///------------------------------------------------------------------------
class HardestJetAxes : public AxesDefinition {
public:
   /// Constructor takes JetDefinition
   HardestJetAxes(fastjet::JetDefinition def)
   : AxesDefinition(), _def(def) {
      setNPass(NO_REFINING);    // default to no minimization
   }
   
   /// Finds seed axes by running a ClusterSequence, running inclusive_jets, and finding the N hardest
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int n_jets,
                                                             const std::vector <fastjet::PseudoJet> & inputs,
                                                             const MeasureDefinition * ) const {
      fastjet::ClusterSequence jet_clust_seq(inputs, _def);
      std::vector<fastjet::PseudoJet> axes = sorted_by_pt(jet_clust_seq.inclusive_jets());
      
      if ((int)axes.size() < n_jets) {
         _too_few_axes_warning.warn("HardestJetAxes::get_starting_axes:  Fewer than N axes found; results are unpredictable.");
      }
      
      axes.resize(n_jets);  // only keep n hardest
      return axes;
   }
   
   /// Short description
   virtual std::string short_description() const { return "HardAxes";}
   /// Long description
   virtual std::string description() const { return "HardAxes: " + _def.description();}
   /// To make it possible to copy around.
   virtual HardestJetAxes* create() const {return new HardestJetAxes(*this);}
   
private:
   fastjet::JetDefinition _def;  ///< Jet Definition to use.
   
   static LimitedWarning _too_few_axes_warning;

};
   
///------------------------------------------------------------------------
/// \class KT_Axes
/// \brief Axes from exclusive kT
///
/// Axes from kT algorithm with E_scheme recombination.
///------------------------------------------------------------------------
class KT_Axes : public ExclusiveJetAxes {
public:
   /// Constructor
   KT_Axes()
   : ExclusiveJetAxes(fastjet::JetDefinition(fastjet::kt_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                             fastjet::E_scheme,
                                             fastjet::Best)
                      ) {
      setNPass(NO_REFINING);
   }

   /// Short description
   virtual std::string short_description() const {
      return "KT";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "KT Axes";
      return stream.str();
   };
   
   /// For copying purposes
   virtual KT_Axes* create() const {return new KT_Axes(*this);}

};

///------------------------------------------------------------------------
/// \class CA_Axes
/// \brief Axes from exclusive CA
///
/// Axes from CA algorithm with E_scheme recombination.
///------------------------------------------------------------------------
class CA_Axes : public ExclusiveJetAxes {
public:
   /// Constructor
   CA_Axes()
   : ExclusiveJetAxes(fastjet::JetDefinition(fastjet::cambridge_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                             fastjet::E_scheme,
                                             fastjet::Best)
                      ) {
      setNPass(NO_REFINING);
   }

   /// Short description
   virtual std::string short_description() const {
      return "CA";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "CA Axes";
      return stream.str();
   };
   
   /// For copying purposes
   virtual CA_Axes* create() const {return new CA_Axes(*this);}
   
};

   
///------------------------------------------------------------------------
/// \class AntiKT_Axes
/// \brief Axes from inclusive anti-kT
///
/// Axes from anti-kT algorithm and E_scheme.
/// The one parameter R0 is subjet radius
///------------------------------------------------------------------------
class AntiKT_Axes : public HardestJetAxes {

public:
   /// Constructor.  Takes jet radius as argument
   AntiKT_Axes(double R0)
   : HardestJetAxes(fastjet::JetDefinition(fastjet::antikt_algorithm,
                                           R0,
                                           fastjet::E_scheme,
                                           fastjet::Best)
                    ), _R0(R0) {
      setNPass(NO_REFINING);
   }

   /// Short description
   virtual std::string short_description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "AKT" << _R0;
      return stream.str();
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Anti-KT Axes (R0 = " << _R0 << ")";
      return stream.str();
   };
   
   /// For copying purposes
   virtual AntiKT_Axes* create() const {return new AntiKT_Axes(*this);}
   
protected:
   double _R0;  ///<  AKT jet radius

};

///------------------------------------------------------------------------
/// \class JetDefinitionWrapper
/// \brief Wrapper for jet definitions (for memory management)
///
/// This class was introduced to avoid issue of a FastJet bug when using genKT clustering
/// Now using this for all AxesDefinition with a manual recombiner to use the delete_recombiner_when_unused function
///------------------------------------------------------------------------
class JetDefinitionWrapper {

public: 
   
   /// Default Constructor
   JetDefinitionWrapper(JetAlgorithm jet_algorithm_in, double R_in, double xtra_param_in, const JetDefinition::Recombiner *recombiner) {
      jet_def = fastjet::JetDefinition(jet_algorithm_in, R_in, xtra_param_in);
      jet_def.set_recombiner(recombiner);
      jet_def.delete_recombiner_when_unused();        // added to prevent memory leaks
   }

   /// Additional constructor so that build-in FastJet algorithms can also be called
   JetDefinitionWrapper(JetAlgorithm jet_algorithm_in, double R_in, const JetDefinition::Recombiner *recombiner, fastjet::Strategy strategy_in) {
      jet_def = fastjet::JetDefinition(jet_algorithm_in, R_in, recombiner, strategy_in);
      jet_def.delete_recombiner_when_unused();
   }

   /// Return jet definition
   JetDefinition getJetDef() {
      return jet_def;
   }

private:
   JetDefinition jet_def;  ///< my jet definition
};

///------------------------------------------------------------------------
/// \class WTA_KT_Axes
/// \brief Axes from exclusive kT, winner-take-all recombination
///
/// Axes from kT algorithm and winner-take-all recombination
///------------------------------------------------------------------------
class WTA_KT_Axes : public ExclusiveJetAxes {
public:
   /// Constructor
   WTA_KT_Axes()
   : ExclusiveJetAxes(JetDefinitionWrapper(fastjet::kt_algorithm,
                                          fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                          _recomb = new WinnerTakeAllRecombiner(), // Needs to be explicitly declared (this will be deleted by JetDefinitionWrapper)
                                          fastjet::Best).getJetDef()
                      ) {
      setNPass(NO_REFINING);
    }

   /// Short description
   virtual std::string short_description() const {
      return "WTA KT";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Winner-Take-All KT Axes";
      return stream.str();
   };
   
   /// For copying purposes
   virtual WTA_KT_Axes* create() const {return new WTA_KT_Axes(*this);}

private:
   const WinnerTakeAllRecombiner *_recomb;  ///< Internal recombiner

};
   
///------------------------------------------------------------------------
/// \class WTA_CA_Axes
/// \brief Axes from exclusive CA, winner-take-all recombination
///
/// Axes from CA algorithm and winner-take-all recombination
///------------------------------------------------------------------------
class WTA_CA_Axes : public ExclusiveJetAxes {
public:
   /// Constructor
   WTA_CA_Axes()
   : ExclusiveJetAxes(JetDefinitionWrapper(fastjet::cambridge_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
                                             _recomb = new WinnerTakeAllRecombiner(), // Needs to be explicitly declared (this will be deleted by JetDefinitionWrapper)
                                             fastjet::Best).getJetDef()) {
    setNPass(NO_REFINING);
  }

   /// Short description
   virtual std::string short_description() const {
      return "WTA CA";
   };
   
   /// Long descriptions
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Winner-Take-All CA Axes";
      return stream.str();
   };
   
   /// For copying purposes
   virtual WTA_CA_Axes* create() const {return new WTA_CA_Axes(*this);}
   
private:
   const WinnerTakeAllRecombiner *_recomb;  ///< Internal recombiner

};


///------------------------------------------------------------------------
/// \class GenKT_Axes
/// \brief Axes from exclusive generalized kT
///
/// Axes from a general KT algorithm (standard E-scheme recombination)
/// Requires the power of the KT algorithm to be used and the radius parameter
///------------------------------------------------------------------------
class GenKT_Axes : public ExclusiveJetAxes {
   
public:
   /// Constructor
   GenKT_Axes(double p, double R0 = fastjet::JetDefinition::max_allowable_R)
   : ExclusiveJetAxes(fastjet::JetDefinition(fastjet::genkt_algorithm,
                                             R0,
                                             p)), _p(p), _R0(R0) {
      if (p < 0) throw Error("GenKT_Axes:  Currently only p >=0 is supported.");
      setNPass(NO_REFINING);
   }
   
   /// Short description
   virtual std::string short_description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "GenKT Axes";
      return stream.str();
   };
   
   /// Long descriptions
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "General KT (p = " << _p << "), R0 = " << _R0;
      return stream.str();
   };
   
   /// For copying purposes
   virtual GenKT_Axes* create() const {return new GenKT_Axes(*this);}
   
protected:
   double _p;   ///< genkT power
   double _R0;  ///< jet radius
};
   
   
///------------------------------------------------------------------------
/// \class WTA_GenKT_Axes
/// \brief Axes from exclusive generalized kT, winner-take-all recombination
///
/// Axes from a general KT algorithm with a Winner Take All Recombiner
/// Requires the power of the KT algorithm to be used and the radius parameter
///------------------------------------------------------------------------
class WTA_GenKT_Axes : public ExclusiveJetAxes {

public:
   /// Constructor
   WTA_GenKT_Axes(double p, double R0 = fastjet::JetDefinition::max_allowable_R)
   : ExclusiveJetAxes(JetDefinitionWrapper(fastjet::genkt_algorithm,
                                            R0,
                                            p,
                                            _recomb = new WinnerTakeAllRecombiner()
                                            ).getJetDef()), _p(p), _R0(R0) {
      if (p < 0) throw Error("WTA_GenKT_Axes:  Currently only p >=0 is supported.");
      setNPass(NO_REFINING);
   }

   /// Short description
   virtual std::string short_description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "WTA, GenKT Axes";
      return stream.str();
   };
   
   /// Long descriptions
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Winner-Take-All General KT (p = " << _p << "), R0 = " << _R0;
      return stream.str();
   };
   
   /// For copying purposes
   virtual WTA_GenKT_Axes* create() const {return new WTA_GenKT_Axes(*this);}
   
protected:
   double _p;   ///< genkT power
   double _R0;  ///< jet radius
   const WinnerTakeAllRecombiner *_recomb; ///< Internal recombiner
};
   
///------------------------------------------------------------------------
/// \class GenET_GenKT_Axes
/// \brief Axes from exclusive kT, generalized Et-scheme recombination
///
/// Class using general KT algorithm with a more general recombination scheme
/// Requires power of KT algorithm, power of recombination weights, and radius parameter
///------------------------------------------------------------------------
class GenET_GenKT_Axes : public ExclusiveJetAxes {

public:
   /// Constructor
   GenET_GenKT_Axes(double delta, double p, double R0 = fastjet::JetDefinition::max_allowable_R)
   : ExclusiveJetAxes((JetDefinitionWrapper(fastjet::genkt_algorithm, R0, p, _recomb = new GeneralEtSchemeRecombiner(delta))).getJetDef() ),
    _delta(delta), _p(p), _R0(R0) {
       if (p < 0) throw Error("GenET_GenKT_Axes:  Currently only p >=0 is supported.");
       if (delta <= 0) throw Error("GenET_GenKT_Axes:  Currently only delta >0 is supported.");
       setNPass(NO_REFINING);
   }

   /// Short description
   virtual std::string short_description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "GenET, GenKT Axes";
      return stream.str();
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2);
      // TODO: if _delta is huge, change to "WTA"
      if (_delta < std::numeric_limits<int>::max()) stream << "General Recombiner (delta = " << _delta << "), " << "General KT (p = " << _p << ") Axes, R0 = " << _R0;
      else stream << "Winner-Take-All General KT (p = " << _p << "), R0 = " << _R0;

      return stream.str();
   };
   
   /// For copying purposes
   virtual GenET_GenKT_Axes* create() const {return new GenET_GenKT_Axes(*this);}
   
protected:
   double _delta; ///< Recombination pT weighting
   double _p;     ///< GenkT power
   double _R0;    ///< jet radius
   const GeneralEtSchemeRecombiner *_recomb;   ///< Internal recombiner
};

///------------------------------------------------------------------------
/// \class OnePass_KT_Axes
/// \brief Axes from exclusive kT, with one-pass minimization
///
/// Onepass minimization from kt axes
///------------------------------------------------------------------------
class OnePass_KT_Axes : public KT_Axes {
public:
   /// Constructor
   OnePass_KT_Axes() : KT_Axes() {
      setNPass(ONE_PASS);
   }
   
   /// Short description
   virtual std::string short_description() const {
      return "OnePass KT";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from KT Axes";
      return stream.str();
   };
   
   /// For copying purposes
   virtual OnePass_KT_Axes* create() const {return new OnePass_KT_Axes(*this);}
   

};

///------------------------------------------------------------------------
/// \class OnePass_CA_Axes
/// \brief Axes from exclusive CA, with one-pass minimization
///
/// Onepass minimization from CA axes
///------------------------------------------------------------------------
class OnePass_CA_Axes : public CA_Axes {
public:
   /// Constructor
   OnePass_CA_Axes() : CA_Axes() {
      setNPass(ONE_PASS);
   }

   /// Short description
   virtual std::string short_description() const {
      return "OnePass CA";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from CA Axes";
      return stream.str();
   };
   
   /// For copying purposes
   virtual OnePass_CA_Axes* create() const {return new OnePass_CA_Axes(*this);}


};
   
///------------------------------------------------------------------------
/// \class OnePass_AntiKT_Axes
/// \brief Axes from inclusive anti-kT, with one-pass minimization
///
/// Onepass minimization from AntiKT axes, one parameter R0
///------------------------------------------------------------------------
class OnePass_AntiKT_Axes : public AntiKT_Axes {

public:
   /// Constructor
   OnePass_AntiKT_Axes(double R0) : AntiKT_Axes(R0) {
      setNPass(ONE_PASS);
   }
   
   /// Short Description
   virtual std::string short_description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "OnePassAKT" << _R0;
      return stream.str();
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Anti-KT Axes (R0 = " << _R0 << ")";
      return stream.str();
   };
   
   /// For copying purposes
   virtual OnePass_AntiKT_Axes* create() const {return new OnePass_AntiKT_Axes(*this);}

};

///------------------------------------------------------------------------
/// \class OnePass_WTA_KT_Axes
/// \brief Axes from exclusive kT, winner-take-all recombination, with one-pass minimization
///
/// Onepass minimization from winner-take-all kt axes
///------------------------------------------------------------------------
class OnePass_WTA_KT_Axes : public WTA_KT_Axes {
public:
   /// Constructor
   OnePass_WTA_KT_Axes() : WTA_KT_Axes() {
      setNPass(ONE_PASS);
   }
   
   /// Short description
   virtual std::string short_description() const {
      return "OnePass WTA KT";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Winner-Take-All KT Axes";
      return stream.str();
   };
   
   /// For copying purposes
   virtual OnePass_WTA_KT_Axes* create() const {return new OnePass_WTA_KT_Axes(*this);}
   

};

///------------------------------------------------------------------------
/// \class OnePass_WTA_CA_Axes
/// \brief Axes from exclusive CA, winner-take-all recombination, with one-pass minimization
///
/// Onepass minimization from winner-take-all CA axes
///------------------------------------------------------------------------
class OnePass_WTA_CA_Axes : public WTA_CA_Axes {
   
public:
   /// Constructor
   OnePass_WTA_CA_Axes() : WTA_CA_Axes() {
      setNPass(ONE_PASS);
   }

   /// Short description
   virtual std::string short_description() const {
      return "OnePass WTA CA";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Winner-Take-All CA Axes";
      return stream.str();
   };
   
   /// For copying purposes
   virtual OnePass_WTA_CA_Axes* create() const {return new OnePass_WTA_CA_Axes(*this);}
   
};

///------------------------------------------------------------------------
/// \class OnePass_GenKT_Axes
/// \brief Axes from exclusive generalized kT with one-pass minimization
///
/// Onepass minimization, General KT Axes (standard E-scheme recombination)
///------------------------------------------------------------------------
class OnePass_GenKT_Axes : public GenKT_Axes {
   
public:
   /// Constructor
   OnePass_GenKT_Axes(double p, double R0 = fastjet::JetDefinition::max_allowable_R) : GenKT_Axes(p, R0) {
      setNPass(ONE_PASS);
   }
   
   /// Short description
   virtual std::string short_description() const {
      return "OnePass GenKT";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from General KT (p = " << _p << "), R0 = " << _R0;
      return stream.str();
   };
   
   /// For copying purposes
   virtual OnePass_GenKT_Axes* create() const {return new OnePass_GenKT_Axes(*this);}
};
   
///------------------------------------------------------------------------
/// \class OnePass_WTA_GenKT_Axes
/// \brief Axes from exclusive generalized kT, winner-take-all recombination, with one-pass minimization
///
/// Onepass minimization from winner-take-all, General KT Axes
///------------------------------------------------------------------------
class OnePass_WTA_GenKT_Axes : public WTA_GenKT_Axes {
   
public:
   /// Constructor
   OnePass_WTA_GenKT_Axes(double p, double R0 = fastjet::JetDefinition::max_allowable_R) : WTA_GenKT_Axes(p, R0) {
      setNPass(ONE_PASS);
   }

   /// Short description
   virtual std::string short_description() const {
      return "OnePass WTA GenKT";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Winner-Take-All General KT (p = " << _p << "), R0 = " << _R0;
      return stream.str();
   };
   
   /// For copying purposes
   virtual OnePass_WTA_GenKT_Axes* create() const {return new OnePass_WTA_GenKT_Axes(*this);}
};

///------------------------------------------------------------------------
/// \class OnePass_GenET_GenKT_Axes
/// \brief Axes from exclusive generalized kT, generalized Et-scheme recombination, with one-pass minimization
///
/// Onepass minimization from General Recomb, General KT axes
///------------------------------------------------------------------------
class OnePass_GenET_GenKT_Axes : public GenET_GenKT_Axes {
   
public:
   /// Constructor
   OnePass_GenET_GenKT_Axes(double delta, double p, double R0 = fastjet::JetDefinition::max_allowable_R) : GenET_GenKT_Axes(delta, p, R0) {
      setNPass(ONE_PASS);
   }

   /// Short description
   virtual std::string short_description() const {
      return "OnePass GenET, GenKT";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2);
      if (_delta < std::numeric_limits<int>::max()) stream << "One-Pass Minimization from General Recombiner (delta = " 
        << _delta << "), " << "General KT (p = " << _p << ") Axes, R0 = " << _R0;
      else stream << "One-Pass Minimization from Winner-Take-All General KT (p = " << _p << "), R0 = " << _R0;
      return stream.str();
   };
   
   /// For copying purposes
   virtual OnePass_GenET_GenKT_Axes* create() const {return new OnePass_GenET_GenKT_Axes(*this);}
};


///------------------------------------------------------------------------
/// \class Manual_Axes
/// \brief Manual axes finding
///
/// Allows the user to set the axes manually
///------------------------------------------------------------------------
class Manual_Axes : public AxesDefinition {
public:
   /// Constructor.  Note that _needsManualAxes is set to true.
   Manual_Axes() : AxesDefinition() {
      setNPass(NO_REFINING);
      _needsManualAxes = true;
   }
   
   /// This is now a dummy function since this is manual mode
   virtual std::vector<fastjet::PseudoJet> get_starting_axes(int,
                                                             const std::vector<fastjet::PseudoJet>&,
                                                             const MeasureDefinition *) const;

   
   /// Short description
   virtual std::string short_description() const {
      return "Manual";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Manual Axes";
      return stream.str();
   };
   
   /// For copying purposes
   virtual Manual_Axes* create() const {return new Manual_Axes(*this);}


};

///------------------------------------------------------------------------
/// \class OnePass_Manual_Axes
/// \brief Manual axes finding, with one-pass minimization
///
/// One pass minimization from manual starting point
///------------------------------------------------------------------------
class OnePass_Manual_Axes : public Manual_Axes {
public:
   /// Constructor.  Note that _needsManualAxes is set to true.
   OnePass_Manual_Axes() : Manual_Axes() {
      setNPass(ONE_PASS);
   }
   
   /// Short description
   virtual std::string short_description() const {
      return "OnePass Manual";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "One-Pass Minimization from Manual Axes";
      return stream.str();
   };
   
   // For copying purposes
   virtual OnePass_Manual_Axes* create() const {return new OnePass_Manual_Axes(*this);}

};
   
///------------------------------------------------------------------------
/// \class MultiPass_Axes
/// \brief Manual axes finding, with multi-pass (randomized) minimization
///
/// Multi-pass minimization from kT starting point
///------------------------------------------------------------------------
class MultiPass_Axes : public KT_Axes {

public:
   
   /// Constructor
   MultiPass_Axes(unsigned int Npass) : KT_Axes() {
      setNPass(Npass);
   }

   /// Short description
   virtual std::string short_description() const {
      return "MultiPass";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Multi-Pass Axes (Npass = " << _Npass << ")";
      return stream.str();
   };
   
   /// For copying purposs
   virtual MultiPass_Axes* create() const {return new MultiPass_Axes(*this);}
   
};

///------------------------------------------------------------------------
/// \class MultiPass_Manual_Axes
/// \brief Axes finding from exclusive kT, with multi-pass (randomized) minimization
///
/// multi-pass minimization from kT starting point
///------------------------------------------------------------------------
class MultiPass_Manual_Axes : public Manual_Axes {

public:
   /// Constructor
   MultiPass_Manual_Axes(unsigned int Npass) : Manual_Axes() {
      setNPass(Npass);
   }

   /// Short Description
   virtual std::string short_description() const {
      return "MultiPass Manual";
   };
   
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Multi-Pass Manual Axes (Npass = " << _Npass << ")";
      return stream.str();
   };
   
   /// For copying purposes
   virtual MultiPass_Manual_Axes* create() const {return new MultiPass_Manual_Axes(*this);}
   
};

///------------------------------------------------------------------------
/// \class Comb_GenKT_Axes
/// \brief Axes from exclusive generalized kT with combinatorial testing
///
/// Axes from kT algorithm (standard E-scheme recombination)
/// Requires nExtra parameter and returns set of N that minimizes N-jettiness
/// Note that this method is not guaranteed to find a deeper minimum than GenKT_Axes
///------------------------------------------------------------------------
class Comb_GenKT_Axes : public ExclusiveCombinatorialJetAxes {
public:
   /// Constructor
   Comb_GenKT_Axes(int nExtra, double p, double R0 = fastjet::JetDefinition::max_allowable_R)
   : ExclusiveCombinatorialJetAxes(fastjet::JetDefinition(fastjet::genkt_algorithm, R0, p), nExtra),
      _p(p), _R0(R0) {
      if (p < 0) throw Error("Comb_GenKT_Axes:  Currently only p >=0 is supported.");
      setNPass(NO_REFINING);
   }
   
   /// Short description
   virtual std::string short_description() const {
      return "N Choose M GenKT";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "N Choose M Minimization (nExtra = " << _nExtra << ") from General KT (p = " << _p << "), R0 = " << _R0;
      return stream.str();
   };
   
   /// For copying purposes
   virtual Comb_GenKT_Axes* create() const {return new Comb_GenKT_Axes(*this);}
   
private:
   double _nExtra;   ///< Number of extra axes
   double _p;        ///< GenkT power
   double _R0;       ///< jet radius
};
   
   

///------------------------------------------------------------------------
/// \class Comb_WTA_GenKT_Axes
/// \brief Axes from exclusive generalized kT, winner-take-all recombination, with combinatorial testing
///
/// Axes from kT algorithm and winner-take-all recombination
/// Requires nExtra parameter and returns set of N that minimizes N-jettiness
///------------------------------------------------------------------------
class Comb_WTA_GenKT_Axes : public ExclusiveCombinatorialJetAxes {
public:
   /// Constructor
   Comb_WTA_GenKT_Axes(int nExtra, double p, double R0 = fastjet::JetDefinition::max_allowable_R)
   : ExclusiveCombinatorialJetAxes((JetDefinitionWrapper(fastjet::genkt_algorithm, R0, p, _recomb = new WinnerTakeAllRecombiner())).getJetDef(), nExtra),
    _p(p), _R0(R0) {
       if (p < 0) throw Error("Comb_WTA_GenKT_Axes:  Currently only p >=0 is supported.");
       setNPass(NO_REFINING);
    }

   /// Short description
   virtual std::string short_description() const {
      return "N Choose M WTA GenKT";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "N Choose M Minimization (nExtra = " << _nExtra << ") from Winner-Take-All General KT (p = " << _p << "), R0 = " << _R0;
      return stream.str();
   };
   
   /// For copying purposes
   virtual Comb_WTA_GenKT_Axes* create() const {return new Comb_WTA_GenKT_Axes(*this);}

private:
   double _nExtra;   ///< Number of extra axes
   double _p;        ///< GenkT power
   double _R0;       ///< jet radius
   const WinnerTakeAllRecombiner *_recomb;   ///< Internal recombiner
};
   
///------------------------------------------------------------------------
/// \class Comb_GenET_GenKT_Axes
/// \brief Axes from exclusive generalized kT, generalized Et-scheme recombination, with combinatorial testing
///
/// Axes from kT algorithm and General Et scheme recombination
/// Requires nExtra parameter and returns set of N that minimizes N-jettiness
///------------------------------------------------------------------------
class Comb_GenET_GenKT_Axes : public ExclusiveCombinatorialJetAxes {
public:
   /// Constructor
   Comb_GenET_GenKT_Axes(int nExtra, double delta, double p, double R0 = fastjet::JetDefinition::max_allowable_R)
   : ExclusiveCombinatorialJetAxes((JetDefinitionWrapper(fastjet::genkt_algorithm, R0, p, _recomb = new GeneralEtSchemeRecombiner(delta))).getJetDef(), nExtra),
    _delta(delta), _p(p), _R0(R0) {
       if (p < 0) throw Error("Comb_GenET_GenKT_Axes:  Currently only p >=0 is supported.");
       if (delta <= 0) throw Error("Comb_GenET_GenKT_Axes:  Currently only delta >=0 is supported.");
        setNPass(NO_REFINING);
    }

   /// Short description
   virtual std::string short_description() const {
      return "N Choose M GenET GenKT";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2);
      if (_delta < std::numeric_limits<int>::max()) stream << "N choose M Minimization (nExtra = " << _nExtra 
        << ") from General Recombiner (delta = " << _delta << "), " << "General KT (p = " << _p << ") Axes, R0 = " << _R0;
      else stream << "N choose M Minimization (nExtra = " << _nExtra << ") from Winner-Take-All General KT (p = " << _p << "), R0 = " << _R0;
      return stream.str();
   };
   
   /// For copying purposes
   virtual Comb_GenET_GenKT_Axes* create() const {return new Comb_GenET_GenKT_Axes(*this);}

private:
   double _nExtra;   ///< Number of extra axes
   double _delta;    ///< Recombination pT weighting exponent
   double _p;        ///< GenkT power
   double _R0;       ///< jet radius
   const GeneralEtSchemeRecombiner *_recomb;  ///<  Internal recombiner
};
   

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESS_HH__

