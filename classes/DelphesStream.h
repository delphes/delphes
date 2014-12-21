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

#ifndef DelphesStream_h
#define DelphesStream_h

/** \class DelphesStream
 *
 *  Provides an interface to manipulate c strings as if they were input streams
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

class DelphesStream
{
public:

  DelphesStream(char *buffer);

  bool ReadDbl(double &value);
  bool ReadInt(int &value);

private:

  char *fBuffer;
  
  static bool fFirstLongMin;
  static bool fFirstLongMax;
  static bool fFirstHugePos;
  static bool fFirstHugeNeg;
  static bool fFirstZero;
};

#endif // DelphesStream_h


