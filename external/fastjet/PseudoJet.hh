//STARTHEADER
// $Id: PseudoJet.hh 2728 2011-11-20 14:18:59Z salam $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER


#ifndef __FASTJET_PSEUDOJET_HH__
#define __FASTJET_PSEUDOJET_HH__

#include<valarray>
#include<vector>
#include<cassert>
#include<cmath>
#include<iostream>
#include "fastjet/internal/numconsts.hh"
#include "fastjet/internal/IsBase.hh"
#include "fastjet/SharedPtr.hh"
#include "fastjet/Error.hh"
#include "fastjet/PseudoJetStructureBase.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//using namespace std;

/// Used to protect against parton-level events where pt can be zero
/// for some partons, giving rapidity=infinity. KtJet fails in those cases.
const double MaxRap = 1e5;

/// default value for phi, meaning it (and rapidity) have yet to be calculated) 
const double pseudojet_invalid_phi = -100.0;

// forward definition
class ClusterSequenceAreaBase;

/// @ingroup basic_classes
/// \class PseudoJet
/// Class to contain pseudojets, including minimal information of use to
/// jet-clustering routines.
class PseudoJet {

 public:
  //----------------------------------------------------------------------
  /// @name Constructors and destructor
  //\{
  /// default constructor, which as of FJ3.0 provides an object for
  /// which all operations are now valid and which has zero momentum
  ///
  // (cf. this is actually OK from a timing point of view and in some
  // cases better than just having the default constructor for the
  // internal shared pointer: see PJtiming.cc and the notes therein)
  PseudoJet() : _px(0), _py(0), _pz(0), _E(0) {_finish_init(); _reset_indices();}
  /// construct a pseudojet from explicit components
  PseudoJet(const double px, const double py, const double pz, const double E);

  /// constructor from any object that has px,py,pz,E = some_four_vector[0--3],
  template <class L> PseudoJet(const L & some_four_vector);

  // Constructor that performs minimal initialisation (only that of
  // the shared pointers), of use in certain speed-critical contexts
  //
  // NB: "dummy" is commented to avoid unused-variable compiler warnings
  PseudoJet(bool /* dummy */) {}

  /// default (virtual) destructor
  virtual ~PseudoJet(){};
  //\} ---- end of constructors and destructors --------------------------

  //----------------------------------------------------------------------
  /// @name Kinematic access functions
  //\{
  //----------------------------------------------------------------------
  inline double E()   const {return _E;}
  inline double e()   const {return _E;} // like CLHEP
  inline double px()  const {return _px;}
  inline double py()  const {return _py;}
  inline double pz()  const {return _pz;}

  /// returns phi (in the range 0..2pi)
  inline double phi() const {return phi_02pi();}

  /// returns phi in the range -pi..pi
  inline double phi_std()  const {
    _ensure_valid_rap_phi();
    return _phi > pi ? _phi-twopi : _phi;}

  /// returns phi in the range 0..2pi
  inline double phi_02pi() const {
    _ensure_valid_rap_phi();
    return _phi;
  }

  /// returns the rapidity or some large value when the rapidity
  /// is infinite
  inline double rap() const {
    _ensure_valid_rap_phi();
    return _rap;
  }

  /// the same as rap()
  inline double rapidity() const {return rap();} // like CLHEP

  /// returns the pseudo-rapidity or some large value when the
  /// rapidity is infinite
  double pseudorapidity() const;
  double eta() const {return pseudorapidity();}

  /// returns the squared transverse momentum
  inline double pt2() const {return _kt2;}
  /// returns the scalar transverse momentum
  inline double  pt() const {return sqrt(_kt2);} 
  /// returns the squared transverse momentum
  inline double perp2() const {return _kt2;}  // like CLHEP
  /// returns the scalar transverse momentum
  inline double  perp() const {return sqrt(_kt2);}    // like CLHEP
  /// returns the squared transverse momentum
  inline double kt2() const {return _kt2;} // for bkwds compatibility

  /// returns the squared invariant mass // like CLHEP
  inline double  m2() const {return (_E+_pz)*(_E-_pz)-_kt2;}    
  /// returns the invariant mass 
  /// (If m2() is negative then -sqrt(-m2()) is returned, as in CLHEP)
  inline double  m() const;    

  /// returns the squared transverse mass = kt^2+m^2
  inline double mperp2() const {return (_E+_pz)*(_E-_pz);}
  /// returns the transverse mass = sqrt(kt^2+m^2)
  inline double mperp() const {return sqrt(std::abs(mperp2()));}
  /// returns the squared transverse mass = kt^2+m^2
  inline double mt2() const {return (_E+_pz)*(_E-_pz);}
  /// returns the transverse mass = sqrt(kt^2+m^2)
  inline double mt() const {return sqrt(std::abs(mperp2()));}

  /// return the squared 3-vector modulus = px^2+py^2+pz^2
  inline double modp2() const {return _kt2+_pz*_pz;}
  /// return the 3-vector modulus = sqrt(px^2+py^2+pz^2)
  inline double modp() const {return sqrt(_kt2+_pz*_pz);}

  /// return the transverse energy
  inline double Et() const {return (_kt2==0) ? 0.0 : _E/sqrt(1.0+_pz*_pz/_kt2);}
  /// return the transverse energy squared
  inline double Et2() const {return (_kt2==0) ? 0.0 : _E*_E/(1.0+_pz*_pz/_kt2);}

  /// returns component i, where X==0, Y==1, Z==2, E==3
  double operator () (int i) const ; 
  /// returns component i, where X==0, Y==1, Z==2, E==3
  inline double operator [] (int i) const { return (*this)(i); }; // this too



  /// returns kt distance (R=1) between this jet and another
  double kt_distance(const PseudoJet & other) const;

  /// returns squared cylinder (rap-phi) distance between this jet and another
  double plain_distance(const PseudoJet & other) const;
  /// returns squared cylinder (rap-phi) distance between this jet and
  /// another
  inline double squared_distance(const PseudoJet & other) const {
    return plain_distance(other);}

  /// return the cylinder (rap-phi) distance between this jet and another,
  /// \f$\Delta_R = \sqrt{\Delta y^2 + \Delta \phi^2}\f$.
  inline double delta_R(const PseudoJet & other) const {
    return sqrt(squared_distance(other));
  }

  /// returns other.phi() - this.phi(), constrained to be in 
  /// range -pi .. pi
  double delta_phi_to(const PseudoJet & other) const;

  //// this seemed to compile except if it was used
  //friend inline double 
  //  kt_distance(const PseudoJet & jet1, const PseudoJet & jet2) { 
  //                                      return jet1.kt_distance(jet2);}

  /// returns distance between this jet and the beam
  inline double beam_distance() const {return _kt2;}

  /// return a valarray containing the four-momentum (components 0-2
  /// are 3-mom, component 3 is energy).
  std::valarray<double> four_mom() const;

  //\}  ------- end of kinematic access functions

  // taken from CLHEP
  enum { X=0, Y=1, Z=2, T=3, NUM_COORDINATES=4, SIZE=NUM_COORDINATES };


  //----------------------------------------------------------------------
  /// @name Kinematic modification functions
  //\{
  //----------------------------------------------------------------------
  /// transform this jet (given in the rest frame of prest) into a jet
  /// in the lab frame [NOT FULLY TESTED]
  PseudoJet & boost(const PseudoJet & prest);
  /// transform this jet (given in lab) into a jet in the rest
  /// frame of prest  [NOT FULLY TESTED]
  PseudoJet & unboost(const PseudoJet & prest);

  void operator*=(double);
  void operator/=(double);
  void operator+=(const PseudoJet &);
  void operator-=(const PseudoJet &);

  /// reset the 4-momentum according to the supplied components and
  /// put the user and history indices back to their default values
  inline void reset(double px, double py, double pz, double E);

  /// reset the PseudoJet to be equal to psjet (including its
  /// indices); NB if the argument is derived from a PseudoJet then
  /// the "reset" used will be the templated version
  ///
  /// Note: this is included on top of the templated version because
  /// PseudoJet is not "derived" from PseudoJet, so the templated
  /// reset would not handle this case properly.
  inline void reset(const PseudoJet & psjet) {
    (*this) = psjet;
  }

  /// reset the 4-momentum according to the supplied generic 4-vector
  /// (accessible via indexing, [0]==px,...[3]==E) and put the user
  /// and history indices back to their default values.
  template <class L> inline void reset(const L & some_four_vector) {
    // check if some_four_vector can be cast to a PseudoJet
    //
    // Note that a regular dynamic_cast would not work here because
    // there is no guarantee that L is polymorphic. We use a more
    // complex construct here that works also in such a case. As for
    // dynamic_cast, NULL is returned if L is not derived from
    // PseudoJet
    const PseudoJet * pj = cast_if_derived<const PseudoJet>(&some_four_vector);

    if (pj){
      (*this) = *pj;
    } else {
      reset(some_four_vector[0], some_four_vector[1],
	    some_four_vector[2], some_four_vector[3]);
    }
  }

  /// reset the PseudoJet according to the specified pt, rapidity,
  /// azimuth and mass (also resetting indices, etc.)
  /// (phi should satisfy -2pi<phi<4pi)
  inline void reset_PtYPhiM(double pt_in, double y_in, double phi_in, double m_in=0.0) {
    reset_momentum_PtYPhiM(pt_in, y_in, phi_in, m_in);
    _reset_indices();
  }

  /// reset the 4-momentum according to the supplied components 
  /// but leave all other information (indices, user info, etc.)
  /// untouched
  inline void reset_momentum(double px, double py, double pz, double E);

  /// reset the 4-momentum according to the components of the supplied
  /// PseudoJet, including cached components; note that the template
  /// version (below) will be called for classes derived from PJ.
  inline void reset_momentum(const PseudoJet & pj);

  /// reset the 4-momentum according to the specified pt, rapidity,
  /// azimuth and mass (phi should satisfy -2pi<phi<4pi)
  void reset_momentum_PtYPhiM(double pt, double y, double phi, double m=0.0);

  /// reset the 4-momentum according to the supplied generic 4-vector
  /// (accessible via indexing, [0]==px,...[3]==E), but leave all
  /// other information (indices, user info, etc.)  untouched
  template <class L> inline void reset_momentum(const L & some_four_vector) {
    reset_momentum(some_four_vector[0], some_four_vector[1],
		   some_four_vector[2], some_four_vector[3]);
  }

  /// in some cases when setting a 4-momentum, the user/program knows
  /// what rapidity and azimuth are associated with that 4-momentum;
  /// by calling this routine the user can provide the information
  /// directly to the PseudoJet and avoid expensive rap-phi
  /// recalculations.
  ///
  /// - \param rap  rapidity
  /// - \param phi  (in range -twopi...4*pi)
  ///
  /// USE WITH CAUTION: there are no checks that the rapidity and
  /// azimuth supplied are sensible, nor does this reset the
  /// 4-momentum components if things don't match.
  void set_cached_rap_phi(double rap, double phi);


  //\} --- end of kin mod functions ------------------------------------

  //----------------------------------------------------------------------
  /// @name User index functions
  ///
  /// To allow the user to set and access an integer index which can
  /// be exploited by the user to associate extra information with a
  /// particle/jet (for example pdg id, or an indication of a
  /// particle's origin within the user's analysis)
  //
  //\{

  /// return the user_index, 
  inline int user_index() const {return _user_index;}
  /// set the user_index, intended to allow the user to add simple
  /// identifying information to a particle/jet
  inline void set_user_index(const int index) {_user_index = index;}

  //\} ----- end of use index functions ---------------------------------

  //----------------------------------------------------------------------
  /// @name User information types and functions
  ///
  /// Allows PseudoJet to carry extra user info (as an object derived from
  /// UserInfoBase).
  //\{

  /// @ingroup user_info
  /// \class UserInfoBase
  /// a base class to hold extra user information in a PseudoJet
  ///
  /// This is a base class to help associate extra user information
  /// with a jet. The user should store their information in a class
  /// derived from this. This allows information of arbitrary
  /// complexity to be easily associated with a PseudoJet (in contrast
  /// to the user index). For example, in a Monte Carlo simulation,
  /// the user information might include the PDG ID, and the position
  /// of the production vertex for the particle.
  ///
  /// The PseudoJet is able to store a shared pointer to any object
  /// derived from UserInfo. The use of a shared pointer frees the
  /// user of the need to handle the memory management associated with
  /// the information.
  ///
  /// Having the user information derive from a common base class also
  /// facilitates dynamic casting, etc.
  ///
  class UserInfoBase{
  public:
    // dummy ctor
    UserInfoBase(){};

    // dummy virtual dtor
    // makes it polymorphic to allow for dynamic_cast
    virtual ~UserInfoBase(){}; 
  };

  /// error class to be thrown if accessing user info when it doesn't
  /// exist
  class InexistentUserInfo : public Error {
  public:
    InexistentUserInfo();
  };

  /// sets the internal shared pointer to the user information.
  ///
  /// Note that the PseudoJet will now _own_ the pointer, and delete
  /// the corresponding object when it (the jet, and any copies of the jet)
  /// goes out of scope. 
  void set_user_info(UserInfoBase * user_info_in) {
    _user_info.reset(user_info_in);
  }

  /// returns a reference to the dynamic cast conversion of user_info
  /// to type L.
  ///
  /// Usage: suppose you have previously set the user info with a pointer
  /// to an object of type MyInfo, 
  ///
  ///   class MyInfo: public PseudoJet::UserInfoBase {
  ///      MyInfo(int id) : _pdg_id(id);
  ///      int pdg_id() const {return _pdg_id;}
  ///      int _pdg_id;
  ///   };
  ///
  ///   PseudoJet particle(...);
  ///   particle.set_user_info(new MyInfo(its_pdg_id));
  ///
  /// Then you would access that pdg_id() as
  ///
  ///   particle.user_info<MyInfo>().pdg_id();
  ///
  /// It's overkill for just a single integer, but scales easily to
  /// more extensive information.
  ///
  /// Note that user_info() throws an InexistentUserInfo() error if
  /// there is no user info; throws a std::bad_cast if the conversion
  /// doesn't work
  ///
  /// If this behaviour does not fit your needs, use instead the the
  /// user_info_ptr() or user_info_shared_ptr() member functions.
  template<class L>
  const L & user_info() const{
    if (_user_info.get() == 0) throw InexistentUserInfo();
    return dynamic_cast<const L &>(* _user_info.get());
  }

  /// returns true if the PseudoJet has user information
  bool has_user_info() const{
    return _user_info.get();
  }

  /// returns true if the PseudoJet has user information than can be
  /// cast to the template argument type.
  template<class L>
  bool has_user_info() const{
    return _user_info.get() && dynamic_cast<const L *>(_user_info.get());
  }

  /// retrieve a pointer to the (const) user information
  const UserInfoBase * user_info_ptr() const{
    if (!_user_info()) return NULL;
    return _user_info.get();
  }


  /// retrieve a (const) shared pointer to the user information
  const SharedPtr<UserInfoBase> & user_info_shared_ptr() const{
    return _user_info;
  }

  /// retrieve a (non-const) shared pointer to the user information;
  /// you can use this, for example, to set the shared pointer, eg
  ///
  /// \code
  ///   p2.user_info_shared_ptr() = p1.user_info_shared_ptr();
  /// \endcode
  ///
  /// or 
  ///
  /// \code
  ///   SharedPtr<PseudoJet::UserInfoBase> info_shared(new MyInfo(...));
  ///   p2.user_info_shared_ptr() = info_shared;
  /// \endcode
  SharedPtr<UserInfoBase> & user_info_shared_ptr(){
    return _user_info;
  }

  // \} --- end of extra info functions ---------------------------------

  //----------------------------------------------------------------------
  /// @name Description
  ///
  /// Since a PseudoJet can have a structure that contains a variety
  /// of information, we provide a description that allows one to check
  /// exactly what kind of PseudoJet we are dealing with
  //
  //\{

  /// return a string describing what kind of PseudoJet we are dealing with 
  std::string description() const;

  //\} ----- end of description functions ---------------------------------

  //-------------------------------------------------------------
  /// @name Access to the associated ClusterSequence object.
  ///
  /// In addition to having kinematic information, jets may contain a
  /// reference to an associated ClusterSequence (this is the case,
  /// for example, if the jet has been returned by a ClusterSequence
  /// member function).
  //\{
  //-------------------------------------------------------------
  /// returns true if this PseudoJet has an associated ClusterSequence.
  bool has_associated_cluster_sequence() const;
  /// shorthand for has_associated_cluster_sequence()
  bool has_associated_cs() const {return has_associated_cluster_sequence();}

  /// returns true if this PseudoJet has an associated and still
  /// valid(ated) ClusterSequence.
  bool has_valid_cluster_sequence() const;
  /// shorthand for has_valid_cluster_sequence()
  bool has_valid_cs() const {return has_valid_cluster_sequence();}

  /// get a (const) pointer to the parent ClusterSequence (NULL if
  /// inexistent)
  const ClusterSequence* associated_cluster_sequence() const;
  // shorthand for associated_cluster_sequence()
  const ClusterSequence* associated_cs() const {return associated_cluster_sequence();}

  /// if the jet has a valid associated cluster sequence then return a
  /// pointer to it; otherwise throw an error
  inline const ClusterSequence * validated_cluster_sequence() const {
    return validated_cs();
  }
  /// shorthand for validated_cluster_sequence()
  const ClusterSequence * validated_cs() const;

  /// if the jet has valid area information then return a pointer to
  /// the associated ClusterSequenceAreaBase object; otherwise throw an error
  inline const ClusterSequenceAreaBase * validated_cluster_sequence_area_base() const {
    return validated_csab();
  }

  /// shorthand for validated_cluster_sequence_area_base()
  const ClusterSequenceAreaBase * validated_csab() const;
  //\}

  //-------------------------------------------------------------
  /// @name Access to the associated PseudoJetStructureBase object.
  ///
  /// In addition to having kinematic information, jets may contain a
  /// reference to an associated ClusterSequence (this is the case,
  /// for example, if the jet has been returned by a ClusterSequence
  /// member function).
  //\{
  //-------------------------------------------------------------

  /// set the associated structure
  void set_structure_shared_ptr(const SharedPtr<PseudoJetStructureBase> &structure);

  /// return true if there is some structure associated with this PseudoJet
  bool has_structure() const;

  /// return a pointer to the structure (of type
  /// PseudoJetStructureBase*) associated with this PseudoJet.
  ///
  /// return NULL if there is no associated structure
  const PseudoJetStructureBase* structure_ptr() const;
  
  /// return a non-const pointer to the structure (of type
  /// PseudoJetStructureBase*) associated with this PseudoJet.
  ///
  /// return NULL if there is no associated structure
  ///
  /// Only use this if you know what you are doing. In any case,
  /// prefer the 'structure_ptr()' (the const version) to this method,
  /// unless you really need a write access to the PseudoJet's
  /// underlying structure.
  PseudoJetStructureBase* structure_non_const_ptr();
  
  /// return a pointer to the structure (of type
  /// PseudoJetStructureBase*) associated with this PseudoJet.
  ///
  /// throw an error if there is no associated structure
  const PseudoJetStructureBase* validated_structure_ptr() const;
  
  /// return a reference to the shared pointer to the
  /// PseudoJetStructureBase associated with this PseudoJet
  const SharedPtr<PseudoJetStructureBase> & structure_shared_ptr() const;

  /// returns a reference to the structure casted to the requested
  /// structure type
  ///
  /// If there is no sructure associated, an Error is thrown.
  /// If the type is not met, a std::bad_cast error is thrown.
  template<typename StructureType>
  const StructureType & structure() const;

  /// check if the PseudoJet has the structure resulting from a Transformer 
  /// (that is, its structure is compatible with a Transformer::StructureType).
  /// If there is no structure, false is returned.
  template<typename TransformerType>
  bool has_structure_of() const;

  /// this is a helper to access any structure created by a Transformer 
  /// (that is, of type Transformer::StructureType).
  ///
  /// If there is no structure, or if the structure is not compatible
  /// with TransformerType, an error is thrown.
  template<typename TransformerType>
  const typename TransformerType::StructureType & structure_of() const;

  //\}

  //-------------------------------------------------------------
  /// @name Methods for access to information about jet structure
  ///
  /// These allow access to jet constituents, and other jet
  /// subtructure information. They only work if the jet is associated
  /// with a ClusterSequence.
  //-------------------------------------------------------------
  //\{

  /// check if it has been recombined with another PseudoJet in which
  /// case, return its partner through the argument. Otherwise,
  /// 'partner' is set to 0.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual bool has_partner(PseudoJet &partner) const;

  /// check if it has been recombined with another PseudoJet in which
  /// case, return its child through the argument. Otherwise, 'child'
  /// is set to 0.
  /// 
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual bool has_child(PseudoJet &child) const;

  /// check if it is the product of a recombination, in which case
  /// return the 2 parents through the 'parent1' and 'parent2'
  /// arguments. Otherwise, set these to 0.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual bool has_parents(PseudoJet &parent1, PseudoJet &parent2) const;

  /// check if the current PseudoJet contains the one passed as
  /// argument.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual bool contains(const PseudoJet &constituent) const;

  /// check if the current PseudoJet is contained the one passed as
  /// argument.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual bool is_inside(const PseudoJet &jet) const;


  /// returns true if the PseudoJet has constituents
  virtual bool has_constituents() const;

  /// retrieve the constituents. 
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence or other substructure information
  virtual std::vector<PseudoJet> constituents() const;


  /// returns true if the PseudoJet has support for exclusive subjets
  virtual bool has_exclusive_subjets() const;

  /// return a vector of all subjets of the current jet (in the sense
  /// of the exclusive algorithm) that would be obtained when running
  /// the algorithm with the given dcut. 
  ///
  /// Time taken is O(m ln m), where m is the number of subjets that
  /// are found. If m gets to be of order of the total number of
  /// constituents in the jet, this could be substantially slower than
  /// just getting that list of constituents.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  std::vector<PseudoJet> exclusive_subjets (const double & dcut) const;

  /// return the size of exclusive_subjets(...); still n ln n with same
  /// coefficient, but marginally more efficient than manually taking
  /// exclusive_subjets.size()
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  int n_exclusive_subjets(const double & dcut) const;

  /// return the list of subjets obtained by unclustering the supplied
  /// jet down to nsub subjets. Throws an error if there are fewer than
  /// nsub particles in the jet.
  ///
  /// For ClusterSequence type jets, requires nsub ln nsub time
  ///
  /// An Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  std::vector<PseudoJet> exclusive_subjets (int nsub) const;

  /// return the list of subjets obtained by unclustering the supplied
  /// jet down to nsub subjets (or all constituents if there are fewer
  /// than nsub).
  ///
  /// For ClusterSequence type jets, requires nsub ln nsub time
  ///
  /// An Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  std::vector<PseudoJet> exclusive_subjets_up_to (int nsub) const;

  /// return the dij that was present in the merging nsub+1 -> nsub 
  /// subjets inside this jet.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  double exclusive_subdmerge(int nsub) const;

  /// return the maximum dij that occurred in the whole event at the
  /// stage that the nsub+1 -> nsub merge of subjets occurred inside 
  /// this jet.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  double exclusive_subdmerge_max(int nsub) const;


  /// returns true if a jet has pieces
  ///
  /// By default a single particle or a jet coming from a
  /// ClusterSequence have no pieces and this methos will return false.
  ///
  /// In practice, this is equivalent to have an structure of type
  /// CompositeJetStructure.
  virtual bool has_pieces() const;


  /// retrieve the pieces that make up the jet. 
  ///
  /// If the jet does not support pieces, an error is throw
  virtual std::vector<PseudoJet> pieces() const;


  // the following ones require a computation of the area in the
  // parent ClusterSequence (See ClusterSequenceAreaBase for details)
  //------------------------------------------------------------------

  /// check if it has a defined area
  virtual bool has_area() const;

  /// return the jet (scalar) area.
  /// throws an Error if there is no support for area in the parent CS
  virtual double area() const;

  /// return the error (uncertainty) associated with the determination
  /// of the area of this jet.
  /// throws an Error if there is no support for area in the parent CS
  virtual double area_error() const;

  /// return the jet 4-vector area.
  /// throws an Error if there is no support for area in the parent CS
  virtual PseudoJet area_4vector() const;

  /// true if this jet is made exclusively of ghosts.
  /// throws an Error if there is no support for area in the parent CS
  virtual bool is_pure_ghost() const;

  //\} --- end of jet structure -------------------------------------



  //----------------------------------------------------------------------
  /// @name Members mainly intended for internal use
  //----------------------------------------------------------------------
  //\{
  /// return the cluster_hist_index, intended to be used by clustering
  /// routines.
  inline int cluster_hist_index() const {return _cluster_hist_index;}
  /// set the cluster_hist_index, intended to be used by clustering routines.
  inline void set_cluster_hist_index(const int index) {_cluster_hist_index = index;}

  /// alternative name for cluster_hist_index() [perhaps more meaningful]
  inline int cluster_sequence_history_index() const {
    return cluster_hist_index();}
  /// alternative name for set_cluster_hist_index(...) [perhaps more
  /// meaningful]
  inline void set_cluster_sequence_history_index(const int index) {
    set_cluster_hist_index(index);}

  //\} ---- end of internal use functions ---------------------------

 protected:  

  SharedPtr<PseudoJetStructureBase> _structure;
  SharedPtr<UserInfoBase> _user_info;


 private: 
  // NB: following order must be kept for things to behave sensibly...
  double _px,_py,_pz,_E;
  mutable double _phi, _rap;
  double _kt2; 
  int    _cluster_hist_index, _user_index;

  /// calculate phi, rap, kt2 based on the 4-momentum components
  void _finish_init();
  /// set the indices to default values
  void _reset_indices();

  /// ensure that the internal values for rapidity and phi 
  /// correspond to 4-momentum structure
  inline void _ensure_valid_rap_phi() const {
    if (_phi == pseudojet_invalid_phi) _set_rap_phi();
  }

  /// set cached rapidity and phi values
  void _set_rap_phi() const;
};


//----------------------------------------------------------------------
// routines for basic binary operations

PseudoJet operator+(const PseudoJet &, const PseudoJet &);
PseudoJet operator-(const PseudoJet &, const PseudoJet &);
PseudoJet operator*(double, const PseudoJet &);
PseudoJet operator*(const PseudoJet &, double);
PseudoJet operator/(const PseudoJet &, double);

/// returns true if the 4 momentum components of the two PseudoJets
/// are identical and all the internal indices (user, cluster_history)
/// + structure and user-info shared pointers are too
bool operator==(const PseudoJet &, const PseudoJet &);

/// inequality test which is exact opposite of operator==
inline bool operator!=(const PseudoJet & a, const PseudoJet & b) {return !(a==b);}

/// Can only be used with val=0 and tests whether all four
/// momentum components are equal to val (=0.0)
bool operator==(const PseudoJet & jet, const double val);

/// Can only be used with val=0 and tests whether at least one of the
/// four momentum components is different from val (=0.0)
inline bool operator!=(const PseudoJet & a, const double & val) {return !(a==val);}

inline double dot_product(const PseudoJet & a, const PseudoJet & b) {
  return a.E()*b.E() - a.px()*b.px() - a.py()*b.py() - a.pz()*b.pz();
}

/// returns true if the momenta of the two input jets are identical
bool have_same_momentum(const PseudoJet &, const PseudoJet &);

/// return a pseudojet with the given pt, y, phi and mass
/// (phi should satisfy -2pi<phi<4pi)
PseudoJet PtYPhiM(double pt, double y, double phi, double m = 0.0);

//----------------------------------------------------------------------
// Routines to do with providing sorted arrays of vectors.

/// return a vector of jets sorted into decreasing transverse momentum
std::vector<PseudoJet> sorted_by_pt(const std::vector<PseudoJet> & jets);

/// return a vector of jets sorted into increasing rapidity
std::vector<PseudoJet> sorted_by_rapidity(const std::vector<PseudoJet> & jets);

/// return a vector of jets sorted into decreasing energy
std::vector<PseudoJet> sorted_by_E(const std::vector<PseudoJet> & jets);

/// return a vector of jets sorted into increasing pz
std::vector<PseudoJet> sorted_by_pz(const std::vector<PseudoJet> & jets);

//----------------------------------------------------------------------
// some code to help sorting

/// sort the indices so that values[indices[0->n-1]] is sorted
/// into increasing order 
void sort_indices(std::vector<int> & indices, 
		  const std::vector<double> & values);

/// given a vector of values with a one-to-one correspondence with the
/// vector of objects, sort objects into an order such that the
/// associated values would be in increasing order (but don't actually
/// touch the values vector in the process).
template<class T> std::vector<T> objects_sorted_by_values(const std::vector<T> & objects, 
					      const std::vector<double> & values);

/// \if internal_doc
/// @ingroup internal
/// \class IndexedSortHelper
/// a class that helps us carry out indexed sorting.
/// \endif
class IndexedSortHelper {
public:
  inline IndexedSortHelper (const std::vector<double> * reference_values) {
    _ref_values = reference_values;
  };
  inline int operator() (const int & i1, const int & i2) const {
    return  (*_ref_values)[i1] < (*_ref_values)[i2];
  };
private:
  const std::vector<double> * _ref_values;
};


//----------------------------------------------------------------------
/// constructor from any object that has px,py,pz,E = some_four_vector[0--3],
// NB: do not know if it really needs to be inline, but when it wasn't
//     linking failed with g++ (who knows what was wrong...)
template <class L> inline  PseudoJet::PseudoJet(const L & some_four_vector) {
  reset(some_four_vector);
}

//----------------------------------------------------------------------
inline void PseudoJet::_reset_indices() { 
  set_cluster_hist_index(-1);
  set_user_index(-1);
  _structure.reset();
  _user_info.reset();
}


// taken literally from CLHEP
inline double PseudoJet::m() const {
  double mm = m2();
  return mm < 0.0 ? -std::sqrt(-mm) : std::sqrt(mm);
}


inline void PseudoJet::reset(double px_in, double py_in, double pz_in, double E_in) {
  _px = px_in;
  _py = py_in;
  _pz = pz_in;
  _E  = E_in;
  _finish_init();
  _reset_indices();
}

inline void PseudoJet::reset_momentum(double px_in, double py_in, double pz_in, double E_in) {
  _px = px_in;
  _py = py_in;
  _pz = pz_in;
  _E  = E_in;
  _finish_init();
}

inline void PseudoJet::reset_momentum(const PseudoJet & pj) {
  _px  = pj._px ;
  _py  = pj._py ;
  _pz  = pj._pz ;
  _E   = pj._E  ;
  _phi = pj._phi;
  _rap = pj._rap;
  _kt2 = pj._kt2;
}

//-------------------------------------------------------------------------------
// implementation of the templated accesses to the underlying structyre
//-------------------------------------------------------------------------------

// returns a reference to the structure casted to the requested
// structure type
//
// If there is no sructure associated, an Error is thrown.
// If the type is not met, a std::bad_cast error is thrown.
template<typename StructureType>
const StructureType & PseudoJet::structure() const{
  return dynamic_cast<const StructureType &>(* validated_structure_ptr());
  
}

// check if the PseudoJet has the structure resulting from a Transformer 
// (that is, its structure is compatible with a Transformer::StructureType)
template<typename TransformerType>
bool PseudoJet::has_structure_of() const{
  if (!_structure()) return false;

  return dynamic_cast<const typename TransformerType::StructureType *>(_structure.get()) != 0;
}

// this is a helper to access a structure created by a Transformer 
// (that is, of type Transformer::StructureType)
// NULL is returned if the corresponding type is not met
template<typename TransformerType>
const typename TransformerType::StructureType & PseudoJet::structure_of() const{
  if (!_structure()) 
    throw Error("Trying to access the structure of a PseudoJet without an associated structure");

  return dynamic_cast<const typename TransformerType::StructureType &>(*_structure);
}



//-------------------------------------------------------------------------------
// helper functions to build a jet made of pieces
//
// Note that there are more complete versions of these functions, with
// an additional argument for a recombination scheme, in
// JetDefinition.hh
// -------------------------------------------------------------------------------

/// build a "CompositeJet" from the vector of its pieces
///
/// In this case, E-scheme recombination is assumed to compute the
/// total momentum
PseudoJet join(const std::vector<PseudoJet> & pieces);

/// build a MergedJet from a single PseudoJet
PseudoJet join(const PseudoJet & j1);

/// build a MergedJet from 2 PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2);

/// build a MergedJet from 3 PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3);

/// build a MergedJet from 4 PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, const PseudoJet & j4);



FASTJET_END_NAMESPACE

#endif // __FASTJET_PSEUDOJET_HH__
