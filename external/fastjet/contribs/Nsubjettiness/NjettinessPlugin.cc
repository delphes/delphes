//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: NjettinessPlugin.cc 663 2014-06-03 21:26:41Z jthaler $
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

#include "NjettinessPlugin.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{



std::string NjettinessPlugin::description() const {return "N-jettiness jet finder";}


// Clusters the particles according to the Njettiness jet algorithm
// Apologies for the complication with this code, but we need to make
// a fake jet clustering tree.  The partitioning is done by getPartitionList
void NjettinessPlugin::run_clustering(ClusterSequence& cs) const
{
   std::vector<fastjet::PseudoJet> particles = cs.jets();

   // HACK: remove area information from particles (in case this is called by
   // a ClusterSequenceArea.  Will be fixed in a future FastJet release)
   for (unsigned i = 0; i < particles.size(); i++) {
      particles[i].set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>());
   }
   
   
   _njettinessFinder.getTau(_N, particles);

   std::vector<std::list<int> > partition = _njettinessFinder.getPartitionList(particles);

   std::vector<fastjet::PseudoJet> jet_indices_for_extras;

   // output clusterings for each jet
   for (size_t i0 = 0; i0 < partition.size(); ++i0) {
      size_t i = partition.size() - 1 - i0; // reversed order of reading to match axes order
      std::list<int>& indices = partition[i];
      if (indices.size() == 0) continue;
      while (indices.size() > 1) {
         int merge_i = indices.back(); indices.pop_back();
         int merge_j = indices.back(); indices.pop_back();
         int newIndex;
         double fakeDij = -1.0;
      
         cs.plugin_record_ij_recombination(merge_i, merge_j, fakeDij, newIndex);

         indices.push_back(newIndex);
      }
      double fakeDib = -1.0;
      
      int finalJet = indices.back();
      cs.plugin_record_iB_recombination(finalJet, fakeDib);
      jet_indices_for_extras.push_back(cs.jets()[finalJet]);  // Get the four vector for the final jets to compare later.
   }

   //HACK:  Re-reverse order of reading to match CS order
   reverse(jet_indices_for_extras.begin(),jet_indices_for_extras.end());

   NjettinessExtras * extras = new NjettinessExtras(_njettinessFinder.currentTauComponents(),jet_indices_for_extras,_njettinessFinder.currentAxes());
   cs.plugin_associate_extras(std::auto_ptr<ClusterSequence::Extras>(extras));
   
}


} // namespace contrib

FASTJET_END_NAMESPACE
