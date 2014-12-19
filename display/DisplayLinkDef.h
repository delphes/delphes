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

#include "display/DelphesDisplay.h"
#include "display/DelphesCaloData.h"
#include "display/DelphesBranchElement.h"
#include "display/Delphes3DGeometry.h"
#include "display/DelphesEventDisplay.h"
#include "display/DelphesHtmlSummary.h"
#include "display/DelphesPlotSummary.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class DelphesDisplay+;
#pragma link C++ class DelphesCaloData+;
#pragma link C++ class DelphesBranchElement<DelphesCaloData>-!;
#pragma link C++ class DelphesBranchElement<TEveElementList>-!;
#pragma link C++ class DelphesBranchElement<TEveTrackList>-!;
#pragma link C++ class Delphes3DGeometry;
#pragma link C++ class DelphesEventDisplay;
#pragma link C++ class DelphesHtmlObjTable;
#pragma link C++ class DelphesHtmlSummary;
#pragma link C++ class DelphesPlotSummary;

#endif

