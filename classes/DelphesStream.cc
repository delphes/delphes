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


/** \class DelphesStream
 *
 *  Provides an interface to manipulate c strings as if they were input streams
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesStream.h"

#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>

#include <iostream>

using namespace std;

//------------------------------------------------------------------------------

bool DelphesStream::fFirstLongMin = true;
bool DelphesStream::fFirstLongMax = true;
bool DelphesStream::fFirstHugePos = true;
bool DelphesStream::fFirstHugeNeg = true;
bool DelphesStream::fFirstZero = true;

//------------------------------------------------------------------------------

DelphesStream::DelphesStream(char *buffer) :
  fBuffer(buffer)
{
}

//------------------------------------------------------------------------------

bool DelphesStream::ReadDbl(double &value)
{
  char *start = fBuffer;
  errno = 0;
  value = strtod(start, &fBuffer);
  if(errno == ERANGE)
  {
    if(fFirstHugePos && value == HUGE_VAL)
    {
      fFirstHugePos = false; 
      cout << "** WARNING: too large positive value, return " << value << endl;
    }
    else if(fFirstHugeNeg && value == -HUGE_VAL)
    {
      fFirstHugeNeg = false; 
      cout << "** WARNING: too large negative value, return " << value << endl;
    }
    else if(fFirstZero)
    {
      fFirstZero = false; 
      value = 0.0;
      cout << "** WARNING: too small value, return " << value << endl;
    }
  }
  return start != fBuffer;
}

//------------------------------------------------------------------------------

bool DelphesStream::ReadInt(int &value)
{
  char *start = fBuffer;
  errno = 0;
  value = strtol(start, &fBuffer, 10);
  if(errno == ERANGE)
  {
    if(fFirstLongMin && value == LONG_MIN)
    {
      fFirstLongMin = false; 
      cout << "** WARNING: too large positive value, return " << value << endl;
    }
    else if(fFirstLongMax && value == LONG_MAX)
    {
      fFirstLongMax = false; 
      cout << "** WARNING: too large negative value, return " << value << endl;
    }
  }
  return start != fBuffer;
}

//------------------------------------------------------------------------------
