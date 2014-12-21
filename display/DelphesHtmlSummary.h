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

// Delphes HTML table for the event display. 
// Based on the ROOT example "alice_esd_html_summary.C"

#ifndef DelphesHtmlSummary_h
#define DelphesHtmlSummary_h 

#include "TArrayF.h"
#include "TOrdCollection.h"


class DelphesHtmlObjTable : public TObject
{
public:                     // make them public for shorter code

   TString   fName;
   Int_t     fNValues;      // number of values
   Int_t     fNFields;      // number of fields
   TArrayF  *fValues;
   TString  *fLabels;
   Bool_t    fExpand;

   TString   fHtml;         // HTML output code

   void Build();
   void BuildTitle();
   void BuildLabels();
   void BuildTable();

public:
   DelphesHtmlObjTable(const char *name, Int_t nfields, Int_t nvals, Bool_t exp=kTRUE);
   virtual ~DelphesHtmlObjTable();

   void     SetLabel(Int_t col, const char *label) { fLabels[col] = label; }
   void     SetValue(Int_t col, Int_t row, Float_t val) { fValues[col].SetAt(val, row); }
   TString  Html() const { return fHtml; }

   ClassDef(DelphesHtmlObjTable, 0);
};

//==============================================================================

class DelphesHtmlSummary
{
public:                           // make them public for shorter code
   Int_t           fNTables;
   TOrdCollection *fObjTables;    // ->array of object tables
   TString         fHtml;         // output HTML string
   TString         fTitle;        // page title
   TString         fHeader;       // HTML header
   TString         fFooter;       // HTML footer

   void     MakeHeader();
   void     MakeFooter();

public:
   DelphesHtmlSummary(const char *title);
   virtual ~DelphesHtmlSummary();

   DelphesHtmlObjTable  *AddTable(const char *name, Int_t nfields, Int_t nvals, 
                           Bool_t exp=kTRUE, Option_t *opt="");
   DelphesHtmlObjTable  *GetTable(Int_t at) const { return (DelphesHtmlObjTable *)fObjTables->At(at); }
   void           Build();
   void           Clear(Option_t *option="");
   void           Reset(Option_t *option="");
   TString        Html() const { return fHtml; }

   ClassDef(DelphesHtmlSummary, 0);
};

#endif // DelphesHtmlSummary_h
