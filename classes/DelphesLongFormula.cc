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


#include "classes/DelphesLongFormula.h"

#include "TString.h"

#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

DelphesLongFormula::DelphesLongFormula() :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesLongFormula::DelphesLongFormula(const char *name, const char *expression) :
  TFormula()
{
}

//------------------------------------------------------------------------------

DelphesLongFormula::~DelphesLongFormula()
{
}

//------------------------------------------------------------------------------

Int_t DelphesLongFormula::Compile(const char *expression)
{
  TString buffer;
  const char *it;
  for(it = expression; *it; ++it)
  {
    if(*it == ' ' || *it == '\t' || *it == '\r' || *it == '\n' || *it == '\\' ) continue;
    buffer.Append(*it);
  }

  buffer.ReplaceAll("pt", "[pt]");
  buffer.ReplaceAll("eta", "[eta]");
  buffer.ReplaceAll("phi", "[phi]");
  buffer.ReplaceAll("energy", "[energy]");
  buffer.ReplaceAll("d0", "[d0]");
  buffer.ReplaceAll("dz", "[dz]");
  buffer.ReplaceAll("ctgTheta", "[ctgTheta]");

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    TFormula::SetMaxima(100000,1000,1000000);
  #endif
  
  if(TFormula::Compile(buffer) != 0)
  {
    throw runtime_error("Invalid Long Formula.");
  }

  return 0;
}

//------------------------------------------------------------------------------

Double_t DelphesLongFormula::Eval(Double_t pt, 
                                  Double_t eta, 
                                  Double_t phi, 
                                  Double_t energy,
                                  Double_t d0,
                                  Double_t dz,
                                  Double_t ctgTheta
                                 )
{

  TVarNameMap fVarNameMap;
  TVarValMap fVarValMap;

  fVarNameMap[this->GetParNumber("pt")]= "pt";
  fVarNameMap[this->GetParNumber("eta")]= "eta";
  fVarNameMap[this->GetParNumber("phi")]= "phi";
  fVarNameMap[this->GetParNumber("energy")]= "energy";
  fVarNameMap[this->GetParNumber("d0")]= "d0";
  fVarNameMap[this->GetParNumber("dz")]= "dz";
  fVarNameMap[this->GetParNumber("ctgTheta")]= "ctgTheta";

  fVarValMap["pt"]= pt;      
  fVarValMap["eta"]= eta;
  fVarValMap["phi"]= phi;
  fVarValMap["energy"]= energy;
  fVarValMap["d0"]= d0;
  fVarValMap["dz"]= dz;
  fVarValMap["ctgTheta"]= ctgTheta;

  Double_t vals[7];

  Int_t j = 0;
  for (Int_t i=0; i != 7; i++)
  {
     if ( fVarNameMap.find(i) != fVarNameMap.end() ) 
     {
        TString var_name = fVarNameMap[i];
        vals[i] = fVarValMap[var_name];
     }
     else
        vals[i] = 0.;
  }   
  return EvalPar(nullptr, vals);
}

//------------------------------------------------------------------------------
