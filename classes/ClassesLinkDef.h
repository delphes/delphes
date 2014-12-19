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


/** ExRootAnalysisLinkDef
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"
#include "classes/DelphesFactory.h"

#include "classes/SortableObject.h"
#include "classes/DelphesClasses.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class DelphesModule+;
#pragma link C++ class DelphesFactory+;

#pragma link C++ class SortableObject+;

#pragma link C++ class Event+;
#pragma link C++ class LHCOEvent+;
#pragma link C++ class LHEFEvent+;
#pragma link C++ class HepMCEvent+;
#pragma link C++ class GenParticle+;
#pragma link C++ class Vertex+;
#pragma link C++ class MissingET+;
#pragma link C++ class ScalarHT+;
#pragma link C++ class Rho+;
#pragma link C++ class Weight+;
#pragma link C++ class Photon+;
#pragma link C++ class Electron+;
#pragma link C++ class Muon+;
#pragma link C++ class Jet+;
#pragma link C++ class Track+;
#pragma link C++ class Tower+;
#pragma link C++ class HectorHit+;

#pragma link C++ class Candidate+;

#endif

