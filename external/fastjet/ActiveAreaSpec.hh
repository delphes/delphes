//STARTHEADER
// $Id: ActiveAreaSpec.hh 2577 2011-09-13 15:11:38Z salam $
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


#ifndef __FASTJET_ACTIVEAREASPEC_HH__
#define __FASTJET_ACTIVEAREASPEC_HH__

// for backwards compatibility fastjet/ActiveAreaSpec.hh provides the
// ActiveAreSpec class that is equivalent to GhostedAreaSpec. The
// latter should now be used (defined in the
// fastjet/GhostedAreaSpec.hh header). ActiveAreaSpec is not
// guaranteed to work in future release of FastJet
#warning This file includes fastjet/ActiveAreaSpec.hh, \
a deprecated FastJet header provided only for backward compatibility. \
This is not guaranteed to work in future releases of FastJet. \
Please consider including fastjet/GhostedAreaSpec.hh directly. \
Similarily, if you use the (deprecated) class ActiveAreaSpec, \
please switch to the equivalent GhostedAreaSpec instead.


#include "fastjet/GhostedAreaSpec.hh"

// NB: ActiveAreaSpec is left over from FJ 2.0 and 2.1;
//     In FJ 2.3 and 2.4 it was just typedefed to GhostedAreaSpec.
//     That's still the case in FJ 3.0, but the typedef has moved
//     from this header into the GhostedAreaSpec.hh header. This way,
//     anyone including ActiveAreaSpec.hh will get a deprecated warning, but
//     people who relied on getting ActiveAreaSpec indirectly will not 
//     have their code broken.


#endif // __FASTJET_ACTIVEAREASPEC_HH__
