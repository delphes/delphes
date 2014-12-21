
/** \class ExRootProgressBar
 *
 *  Class showing progress bar
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootProgressBar.h"

#include "TSystem.h"

#include <iostream>

#include <string.h>
#include <stdio.h>

using namespace std;

ExRootProgressBar::ExRootProgressBar(Long64_t entries, Int_t width) :
  fEntries(entries), fEventCounter(0), fWidth(width), fTime(0), fHashes(-1), fBar(0)
{
  fBar = new char[width + 1];
  memset(fBar, '-', width);
  fBar[width] = 0;
}

//------------------------------------------------------------------------------

ExRootProgressBar::~ExRootProgressBar()
{
  if(fBar) delete[] fBar;
}

//------------------------------------------------------------------------------

void ExRootProgressBar::Update(Long64_t entry, Long64_t eventCounter, Bool_t last)
{
  ULong64_t time = gSystem->Now();

  if(time < fTime + 500 && entry < fEntries && !last) return;

  fTime = time;

  if(fEntries > 0)
  {
    Int_t hashes = Int_t(Double_t(entry)/fEntries*fWidth);
    if(hashes > fWidth) hashes = fWidth;
    if(hashes != fHashes)
    {
      memset(fBar, '#', hashes);
      memset(fBar + hashes, '-', fWidth - hashes);
      fHashes = hashes;
      fprintf(stderr, "** [%s] (%.2f%%)\r", fBar, Float_t(entry)/fEntries*100.0);
    }
  }
  else
  {
    if(eventCounter > fEventCounter)
    {
      fEventCounter = eventCounter;
      fprintf(stderr, "** %lld events processed\r", eventCounter);
    }
  }

  fflush(stderr);
}

//------------------------------------------------------------------------------

void ExRootProgressBar::Finish()
{
  fprintf(stderr, "\n");
  fflush(stderr);
}

//------------------------------------------------------------------------------

