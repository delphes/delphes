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

#include "display/DelphesHtmlSummary.h"
#include "TGHtml.h"
#include "TEveElement.h"
#include "TEveManager.h"
#include "TEveEventManager.h"
#include "TEvePointSet.h"
#include "TEveTrack.h"

//==============================================================================

DelphesHtmlSummary *fgDelphesHtmlSummary = 0;
TGHtml      *fgHtml        = 0;

//==============================================================================

//______________________________________________________________________________
DelphesHtmlObjTable::DelphesHtmlObjTable(const char *name, Int_t nfields, Int_t nvals, Bool_t exp) : 
   fName(name), fNValues(nvals), fNFields(nfields), fExpand(exp)
{
   // Constructor.

   fValues = new TArrayF[fNFields];
   for (int i=0;i<fNFields;i++)
      fValues[i].Set(nvals);
   fLabels = new TString[fNFields];
}

//______________________________________________________________________________
DelphesHtmlObjTable::~DelphesHtmlObjTable()
{
   // Destructor.

   delete [] fValues;
   delete [] fLabels;
}

//______________________________________________________________________________
void DelphesHtmlObjTable::Build()
{
   // Build HTML code.

   fHtml = "<table width=100% border=1 cellspacing=0 cellpadding=0 bgcolor=f0f0f0> ",

   BuildTitle();
   if (fExpand && (fNFields > 0) && (fNValues > 0)) {
      BuildLabels();
      BuildTable();
   }

   fHtml += "</table>";
}

//______________________________________________________________________________
void DelphesHtmlObjTable::BuildTitle()
{
   // Build table title.
   
   fHtml += "<tr><td colspan=";
   fHtml += Form("%d>", fNFields+1);
   fHtml += "<table width=100% border=0 cellspacing=2 cellpadding=0 bgcolor=6e6ea0>";
   fHtml += "<tr><td align=left>";
   fHtml += "<font face=Verdana size=3 color=ffffff><b><i>";
   fHtml += fName;
   fHtml += "</i></b></font></td>";
   fHtml += "<td>";
   fHtml += "<td align=right> ";
   fHtml += "<font face=Verdana size=3 color=ffffff><b><i>";
   fHtml += Form("Size = %d", fNValues);
   fHtml += "</i></b></font></td></tr>";
   fHtml += "</table>";
   fHtml += "</td></tr>";
}

//______________________________________________________________________________
void DelphesHtmlObjTable::BuildLabels()
{
   // Build table labels.

   Int_t i;
   fHtml += "<tr bgcolor=c0c0ff>";
   fHtml += "<th> </th>"; // for the check boxes
   for (i=0;i<fNFields;i++) {
      fHtml += "<th> ";
      fHtml += fLabels[i];
      fHtml += " </th>"; // for the check boxes
   }
   fHtml += "</tr>";
}

//______________________________________________________________________________
void DelphesHtmlObjTable::BuildTable()
{
   // Build part of table with values.

   for (int i = 0; i < fNValues; i++) {
      if (i%2)
         fHtml += "<tr bgcolor=e0e0ff>";
      else
         fHtml += "<tr bgcolor=ffffff>";
      
      TString name = fName;
      name.ReplaceAll(" ", "_");
      // checkboxes
      fHtml += "<td bgcolor=d0d0ff align=\"center\">";
      fHtml += "<input type=\"checkbox\" name=\"";
      fHtml += name;
      fHtml += Form("[%d]\">",i);
      fHtml += "</td>";

      for (int j = 0; j < fNFields; j++) {
         fHtml += "<td width=";
         fHtml += Form("%d%%", 100/fNFields);
         fHtml += " align=\"center\"";
         fHtml += ">";
         fHtml += Form("%1.4f", fValues[j][i]);
         fHtml += "</td>";
      }
      fHtml += "</tr> ";
   }
}

//______________________________________________________________________________
DelphesHtmlSummary::DelphesHtmlSummary(const char *title) : fNTables(0), fTitle(title)
{
   // Constructor.

   fObjTables = new TOrdCollection();
}

//______________________________________________________________________________
DelphesHtmlSummary::~DelphesHtmlSummary()
{
   // Destructor.

   Reset();
}

//______________________________________________________________________________
DelphesHtmlObjTable *DelphesHtmlSummary::AddTable(const char *name, Int_t nfields, Int_t nvals,
                                    Bool_t exp, Option_t *option)
{
   // Add a new table in our list of tables.

   TString opt = option;
   opt.ToLower();
   DelphesHtmlObjTable *table = new DelphesHtmlObjTable(name, nfields, nvals, exp);
   fNTables++;
   if (opt.Contains("first"))
      fObjTables->AddFirst(table);
   else
      fObjTables->Add(table);
   return table;
}

//______________________________________________________________________________
void DelphesHtmlSummary::Clear(Option_t *option)
{
   // Clear the table list.

   if (option && option[0] == 'D')
      fObjTables->Delete(option);
   else
      fObjTables->Clear(option);
   fNTables = 0;
}

//______________________________________________________________________________
void DelphesHtmlSummary::Reset(Option_t *)
{
   // Reset (delete) the table list;

   delete fObjTables; fObjTables = 0;
   fNTables = 0;
}

//______________________________________________________________________________
void DelphesHtmlSummary::Build()
{
   // Build the summary.

   MakeHeader();
   for (int i=0;i<fNTables;i++) {
      GetTable(i)->Build();
      fHtml += GetTable(i)->Html();
   }
   MakeFooter();
}

//______________________________________________________________________________
void DelphesHtmlSummary::MakeHeader()
{
   // Make HTML header.

   fHeader  = "<html><head><title>";
   fHeader += fTitle;
   fHeader += "</title></head><body>";
   fHeader += "<center><h2><font color=#2222ee><i>";
   fHeader += fTitle;
   fHeader += "</i></font></h2></center>";
   fHtml    = fHeader;
}

//______________________________________________________________________________
void DelphesHtmlSummary::MakeFooter()
{
   // Make HTML footer.

   fFooter  = "<br><p><br><center><strong><font size=2 color=#2222ee>";
   fFooter += "Example of using Html widget to display tabular data";
   fFooter += "<br>";
   fFooter += "(c) 2007-2010 Bertrand Bellenot";
   fFooter += "</font></strong></center></body></html>";  
   fHtml   += fFooter;
}

