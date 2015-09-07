/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Delphes.h"

#include "modules/AngularSmearing.h"
#include "modules/PhotonConversions.h"
#include "modules/ParticlePropagator.h"
#include "modules/Efficiency.h"
#include "modules/IdentificationMap.h"
#include "modules/EnergySmearing.h"
#include "modules/MomentumSmearing.h"
#include "modules/ImpactParameterSmearing.h"
#include "modules/TimeSmearing.h"
#include "modules/SimpleCalorimeter.h"
#include "modules/Calorimeter.h"
#include "modules/Isolation.h"
#include "modules/EnergyScale.h"
#include "modules/UniqueObjectFinder.h"
#include "modules/TrackCountingBTagging.h"
#include "modules/BTagging.h"
#include "modules/TauTagging.h"
#include "modules/TreeWriter.h"
#include "modules/Merger.h"
#include "modules/LeptonDressing.h"
#include "modules/PileUpMerger.h"
#include "modules/JetPileUpSubtractor.h"
#include "modules/TrackPileUpSubtractor.h"
#include "modules/TaggingParticlesSkimmer.h"
#include "modules/PileUpJetID.h"
#include "modules/ConstituentFilter.h"
#include "modules/StatusPidFilter.h"
#include "modules/PdgCodeFilter.h"
#include "modules/Cloner.h"
#include "modules/Weighter.h"
#include "modules/Hector.h"
#include "modules/JetFlavorAssociation.h"
#include "modules/JetFakeParticle.h"
#include "modules/ExampleModule.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Delphes+;

#pragma link C++ class AngularSmearing+;
#pragma link C++ class PhotonConversions+;
#pragma link C++ class ParticlePropagator+;
#pragma link C++ class Efficiency+;
#pragma link C++ class IdentificationMap+;
#pragma link C++ class EnergySmearing+;
#pragma link C++ class MomentumSmearing+;
#pragma link C++ class ImpactParameterSmearing+;
#pragma link C++ class TimeSmearing+;
#pragma link C++ class SimpleCalorimeter+;
#pragma link C++ class Calorimeter+;
#pragma link C++ class Isolation+;
#pragma link C++ class EnergyScale+;
#pragma link C++ class UniqueObjectFinder+;
#pragma link C++ class TrackCountingBTagging+;
#pragma link C++ class BTagging+;
#pragma link C++ class TauTagging+;
#pragma link C++ class TreeWriter+;
#pragma link C++ class Merger+;
#pragma link C++ class LeptonDressing+;
#pragma link C++ class PileUpMerger+;
#pragma link C++ class JetPileUpSubtractor+;
#pragma link C++ class TrackPileUpSubtractor+;
#pragma link C++ class TaggingParticlesSkimmer+;
#pragma link C++ class PileUpJetID+;
#pragma link C++ class ConstituentFilter+;
#pragma link C++ class StatusPidFilter+;
#pragma link C++ class PdgCodeFilter+;
#pragma link C++ class Cloner+;
#pragma link C++ class Weighter+;
#pragma link C++ class Hector+;
#pragma link C++ class JetFlavorAssociation+;
#pragma link C++ class JetFakeParticle+;
#pragma link C++ class ExampleModule+;

#endif
