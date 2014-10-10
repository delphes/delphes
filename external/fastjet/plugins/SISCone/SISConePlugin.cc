
// fastjet stuff
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"

// sisocne stuff
#include "momentum.h"
#include "siscone.h"

// other stuff
#include<sstream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;
using namespace siscone;

/// shortcut for converting siscone Cmomentum into PseudoJet
template<> PseudoJet::PseudoJet(const siscone::Cmomentum & four_vector) {
  (*this) = PseudoJet(four_vector.px,four_vector.py,four_vector.pz,
                      four_vector.E);
}

//======================================================================
// wrap-up around siscone's user-defined scales
namespace siscone_plugin_internal{
  /// @ingroup internal
  /// \class SISConeUserScale
  /// class that makes the transition between the internal SISCone
  /// user-defined scale choice (using SISCone's Cjet) and
  /// user-defined scale choices in the plugn above (using FastJet's
  /// PseudoJets)
  class SISConeUserScale  : public siscone::Csplit_merge::Cuser_scale_base{
  public:
    /// ctor takes the "fastjet-style" user-defined scale as well as a
    /// reference to the current cluster sequence (to access the
    /// particles if needed)
    SISConeUserScale(const SISConePlugin::UserScaleBase *user_scale,
		     const ClusterSequence &cs)
      : _user_scale(user_scale), _cs(cs){}

    /// returns the scale associated to a given jet
    virtual double operator()(const siscone::Cjet &jet) const{
      return _user_scale->result(_build_PJ_from_Cjet(jet));
    }

    /// returns true id the scasle associated to jet a is larger than
    /// the scale associated to jet b
    virtual bool is_larger(const siscone::Cjet &a, const siscone::Cjet &b) const{
      return _user_scale->is_larger(_build_PJ_from_Cjet(a), _build_PJ_from_Cjet(b));
    }

  private:
    /// constructs a PseudoJet from a siscone::Cjet
    ///
    /// Note that it is tempting to overload the PseudoJet ctor. This
    /// would not work because down the line we need to access the
    /// original PseudoJet through the ClusterSequence and therefore
    /// the PseudoJet structure needs to be aware of the
    /// ClusterSequence.
    PseudoJet _build_PJ_from_Cjet(const siscone::Cjet &jet) const{
      PseudoJet j(jet.v.px, jet.v.py, jet.v.pz, jet.v.E);
      j.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(
                                   new SISConePlugin::UserScaleBaseStructureType<siscone::Cjet>(jet,_cs)));      
      return j;
    }

    const SISConePlugin::UserScaleBase *_user_scale;
    const ClusterSequence & _cs;
  };
}

// end of the internal material
//======================================================================


/////////////////////////////////////////////
// static members declaration              //
/////////////////////////////////////////////
std::auto_ptr<SISConePlugin>           SISConePlugin::stored_plugin;
std::auto_ptr<std::vector<PseudoJet> > SISConePlugin::stored_particles;
std::auto_ptr<Csiscone>                SISConePlugin::stored_siscone;


/////////////////////////////////////////////
// now comes the implementation itself     //
/////////////////////////////////////////////
string SISConePlugin::description () const {
  ostringstream desc;
  
  const string on = "on";
  const string off = "off";

  string sm_scale_string = "split-merge uses " + 
    split_merge_scale_name(Esplit_merge_scale(split_merge_scale()));

  desc << "SISCone jet algorithm with " ;
  desc << "cone_radius = "       << cone_radius        () << ", ";
  if (_progressive_removal)
    desc << "progressive-removal mode, ";
  else  
    desc << "overlap_threshold = " << overlap_threshold  () << ", ";
  desc << "n_pass_max = "        << n_pass_max         () << ", ";
  desc << "protojet_ptmin = "    << protojet_ptmin()      << ", ";
  if (_progressive_removal && _user_scale) {
    desc << "using a user-defined scale for ordering of stable cones";
    string user_scale_desc = _user_scale->description();
    if (user_scale_desc != "") desc << " (" << user_scale_desc << ")";
  } else {
    desc <<  sm_scale_string;
  }
  if (!_progressive_removal){
    desc << ", caching turned "      << (caching() ? on : off);
    desc << ", SM stop scale = "     << _split_merge_stopping_scale;
  }

  // add a note to the description if we use the pt-weighted splitting
  if (_use_pt_weighted_splitting){
    desc << ", using pt-weighted splitting";
  }

  if (_use_jet_def_recombiner){
    desc << ", using jet-definition's own recombiner";
  }

  // create a fake siscone object so that we can find out more about it
  Csiscone siscone;
  if (siscone.merge_identical_protocones) {
    desc << ", and (IR unsafe) merge_indentical_protocones=true" ;
  }

  desc << ", SISCone code v" << siscone_version();

  return desc.str();
}


// overloading the base class implementation
void SISConePlugin::run_clustering(ClusterSequence & clust_seq) const {

  Csiscone::set_banner_stream(clust_seq.fastjet_banner_stream());

  Csiscone   local_siscone;
  Csiscone * siscone;
  
  unsigned n = clust_seq.jets().size();

  bool new_siscone = true; // by default we'll be running it

  if (caching() && !_progressive_removal) {

    // Establish if we have a cached run with the same R, npass and
    // particles. If not then do any tidying up / reallocation that's
    // necessary for the next round of caching, otherwise just set
    // relevant pointers so that we can reuse and old run.
    if (stored_siscone.get() != 0) {
      new_siscone = !(stored_plugin->cone_radius()   == cone_radius()
                      && stored_plugin->n_pass_max() == n_pass_max()  
                      && stored_particles->size()    == n);
      if (!new_siscone) {
        for(unsigned i = 0; i < n; i++) {
          // only check momentum because indices will be correctly dealt
          // with anyway when extracting the clustering order.
          new_siscone |= !have_same_momentum(clust_seq.jets()[i], 
                                             (*stored_particles)[i]);
        }
      }
    } 
      
    // allocate the new siscone, etc., if need be
    if (new_siscone) {
      stored_siscone  .reset( new Csiscone );
      stored_particles.reset( new std::vector<PseudoJet>(clust_seq.jets()));
      reset_stored_plugin();
    }

    siscone = stored_siscone.get();
  } else {
    siscone = &local_siscone;
  }

  // make sure stopping scale is set in siscone
  siscone->SM_var2_hardest_cut_off = _split_merge_stopping_scale*_split_merge_stopping_scale;

  // set the specific parameters
  // when running with ghosts for passive areas, do not put the
  // ghosts into the stable-cone search (not relevant)
  siscone->stable_cone_soft_pt2_cutoff = ghost_separation_scale()
                                         * ghost_separation_scale();
  // set the type of splitting we want (default=std one, true->pt-weighted split)
  siscone->set_pt_weighted_splitting(_use_pt_weighted_splitting);

  if (new_siscone) {
    // transfer fastjet initial particles into the siscone type
    std::vector<Cmomentum> siscone_momenta(n);
    for(unsigned i = 0; i < n; i++) {
      const PseudoJet & p = clust_seq.jets()[i]; // shorthand
      siscone_momenta[i] = Cmomentum(p.px(), p.py(), p.pz(), p.E());
    }
    
    // run the jet finding
    //cout << "plg sms: " << split_merge_scale() << endl;
    if (_progressive_removal){
      // handle the optional user-defined scale choice
      SharedPtr<siscone_plugin_internal::SISConeUserScale> internal_scale;
      if (_user_scale){
	internal_scale.reset(new siscone_plugin_internal::SISConeUserScale(_user_scale, clust_seq));
	siscone->set_user_scale(internal_scale.get());
      }
      siscone->compute_jets_progressive_removal(siscone_momenta, cone_radius(),
						n_pass_max(), protojet_or_ghost_ptmin(), 
						Esplit_merge_scale(split_merge_scale()));
    } else {
      siscone->compute_jets(siscone_momenta, cone_radius(), overlap_threshold(),
			    n_pass_max(), protojet_or_ghost_ptmin(), 
			    Esplit_merge_scale(split_merge_scale()));
    }
  } else {
    // rerun the jet finding
    // just run the overlap part of the jets.
    //cout << "plg rcmp sms: " << split_merge_scale() << endl;
    siscone->recompute_jets(overlap_threshold(), protojet_or_ghost_ptmin(), 
			    Esplit_merge_scale(split_merge_scale()));
  }

  // extract the jets [in reverse order -- to get nice ordering in pt at end]
  int njet = siscone->jets.size();

  // allocate space for the extras object
  SISConeExtras * extras = new SISConeExtras(n);

  // the ordering in which the inclusive jets are transfered here is
  // deliberate and ensures that when a user asks for
  // inclusive_jets(), they are provided in the order in which SISCone
  // created them.
  for (int ijet = njet-1; ijet >= 0; ijet--) {
    const Cjet & jet = siscone->jets[ijet]; // shorthand
    
    // Successively merge the particles that make up the cone jet
    // until we have all particles in it.  Start off with the zeroth
    // particle.
    int jet_k = jet.contents[0];
    for (unsigned ipart = 1; ipart < jet.contents.size(); ipart++) {
      // take the last result of the merge
      int jet_i = jet_k;
      // and the next element of the jet
      int jet_j = jet.contents[ipart];
      // and merge them (with a fake dij)
      double dij = 0.0;

      if (_use_jet_def_recombiner) {
	clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, jet_k);
      } else {
	// create the new jet by hand so that we can adjust its user index
	PseudoJet newjet = clust_seq.jets()[jet_i] + clust_seq.jets()[jet_j];
	// set the user index to be the pass in which the jet was discovered
       // *** The following line was commented for 3.0.1 ***
	//newjet.set_user_index(jet.pass);
	clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, newjet, jet_k);
      }

    }

    // we have merged all the jet's particles into a single object, so now
    // "declare" it to be a beam (inclusive) jet.
    // [NB: put a sensible looking d_iB just to be nice...]
    double d_iB = clust_seq.jets()[jet_k].perp2();
    clust_seq.plugin_record_iB_recombination(jet_k, d_iB);

    // now record the pass of the jet in the extras object
    extras->_pass[clust_seq.jets()[jet_k].cluster_hist_index()] = jet.pass;
  }

  // now copy the list of protocones into an "extras" objects
  for (unsigned ipass = 0; ipass < siscone->protocones_list.size(); ipass++) {
    for (unsigned ipc = 0; ipc < siscone->protocones_list[ipass].size(); ipc++) {
      PseudoJet protocone(siscone->protocones_list[ipass][ipc]);
      protocone.set_user_index(ipass);
      extras->_protocones.push_back(protocone);
    }
  }
  extras->_most_ambiguous_split = siscone->most_ambiguous_split;

  // tell it what the jet definition was
  extras->_jet_def_plugin = this;

  // give the extras object to the cluster sequence.
  // 
  // As of v3.1 of FastJet, extras are automatically owned (as
  // SharedPtr) by the ClusterSequence and auto_ptr is deprecated. So
  // we can use a simple pointer here
  //clust_seq.plugin_associate_extras(std::auto_ptr<ClusterSequence::Extras>(extras));
  clust_seq.plugin_associate_extras(extras);
}


void SISConePlugin::reset_stored_plugin() const{
  stored_plugin.reset( new SISConePlugin(*this));
}

// //======================================================================
// // material to handle user-defined scales
// 
// //--------------------------------------------------
// // SISCone structure type
// 
// // the textual descripotion
// std::string SISConePlugin::UserScaleBase::StructureType::description() const{
//   return "PseudoJet wrapping a siscone::Cjet stable cone"; 
// }
// 
// // retrieve the constituents
// //
// // if you simply need to iterate over the constituents, it will be
// // faster to access them via constituent(i)
// vector<PseudoJet> SISConePlugin::UserScaleBase::StructureType::constituents(const PseudoJet &) const{ 
//   vector<PseudoJet> constits;
//   constits.reserve(_jet.n);
//   for (unsigned int i=0; i<(unsigned int)_jet.n;i++)
//     constits.push_back(constituent(i));
//   return constits;
// }
// 
// // returns the number of constituents
// unsigned int SISConePlugin::UserScaleBase::StructureType::size() const{
//   return _jet.n;
// }
// 
// // returns the index (in the original particle list) of the ith
// // constituent
// int SISConePlugin::UserScaleBase::StructureType::constituent_index(unsigned int i) const{ 
//   return _jet.contents[i];
// }
// 
// // returns the ith constituent (as a PseusoJet)
// const PseudoJet & SISConePlugin::UserScaleBase::StructureType::constituent(unsigned int i) const{
//   return _cs.jets()[_jet.contents[i]];
// }
// 
// // returns the scalar pt of this stable cone
// double SISConePlugin::UserScaleBase::StructureType::pt_tilde() const{
//   return _jet.pt_tilde;
// }
// 
// // returns the sm_var2 (signed ordering variable squared) for this stable cone
// double SISConePlugin::UserScaleBase::StructureType::ordering_var2() const{
//   return _jet.sm_var2;
// }


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
