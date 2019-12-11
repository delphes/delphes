/** \class HighMassVertexRecover
 *
 *  Try to assign a vertex also to tracks which have
 *  not been clusterized changing the mass hypotesis.
 *
 *  \author Olmo Cerri, Caltech
 *
 */

#include "modules/HighMassVertexRecover.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TObjArray.h"
#include "TLorentzVector.h"

#include <iostream>
#include <sstream>


using namespace std;
namespace vtx_DAZT
{
  static const Double_t c_light = 2.99792458e+8; // [m/s]
}
using namespace vtx_DAZT;

//------------------------------------------------------------------------------

HighMassVertexRecover::HighMassVertexRecover()
{
}

//------------------------------------------------------------------------------

HighMassVertexRecover::~HighMassVertexRecover()
{
}

//------------------------------------------------------------------------------

void HighMassVertexRecover::Init()
{
  // Get parameters from the card
  fVerbose         = GetInt("Verbose", 0);

  fSigmaCompatibility = GetDouble("SigmaCompatibility", 2.0);

  ExRootConfParam param = GetParam("MassHypot");
  UInt_t size = param.GetSize();
  if(size > 0)
  {
    for(UInt_t i = 0; i < size; ++i)
    {
      Double_t mass = GetDouble(param[i].GetString(), 1.);
      fMassList.push_back(mass);
    }
  }
  else
  {
    // Ordered by production probability in the SM
    fMassList.push_back(0.13957);   //pion (+/-)
    fMassList.push_back(0.49367);   // k (+/-)
    fMassList.push_back(0.93827);   // proton
    fMassList.push_back(0.00051);   //electron
    fMassList.push_back(0.10558);   //muon
  }

  if(fVerbose)
  {
    for(unsigned int i = 0; i < fMassList.size(); i++)
    {
      cout << "mass: " << fMassList[i] << endl;
    }
  }


  // import input array
  fTrackInputArray = ImportArray(GetString("TrackInputArray", "VertexFinderDAClusterizerZT/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  fVertexInputArray = ImportArray(GetString("VertexInputArray", "VertexFinderDAClusterizerZT/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();

  // create output array
  fTrackOutputArray = ExportArray(GetString("TrackOutputArray", "tracks"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void HighMassVertexRecover::Finish()
{
  delete fItTrackInputArray;
  delete fItVertexInputArray;
  delete fItVertexOutputArray;
}

//------------------------------------------------------------------------------

void HighMassVertexRecover::Process()
{
  if(fVerbose > 1)
  {
    cout << endl << endl << endl << "------------------ EVENT ----------------------------------------------" << endl;
  }
  // Make a copy of input VerticesfItVertexInputArray->Reset();
  Candidate * vertex;
  fItVertexInputArray->Reset();
  while((vertex = static_cast<Candidate*>(fItVertexInputArray->Next())))
  {
    fVertexOutputArray->Add(static_cast<Candidate*>(vertex->Clone()));
  }
  fItVertexOutputArray = fVertexOutputArray->MakeIterator();

  fItTrackInputArray->Reset();
  Candidate * track;
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {
    Candidate* mother = track;
    track = static_cast<Candidate*>(track->Clone());
    track->AddCandidate(mother);

    if(track->ClusterIndex == -1)
    {
      if(fVerbose > 1)
      {
        cout << endl << "PID: " << track->PID << endl;
        cout << "Pt: " << track->Momentum.Pt() << endl;
        cout << Form("Zd: %.3f +/- %.3f", track->Zd, track->ErrorDZ)  << endl;
      }
      vector<Candidate*> vtx;                       //vector on vertices compatible in z
      // Fill it with onyl the vertexes compatible
      fItVertexOutputArray->Reset();
      Candidate * vertex;
      while((vertex = static_cast<Candidate*>(fItVertexOutputArray->Next())))
      {
        double dz = vertex->Position.Z() - track->Zd;
        dz /= TMath::Hypot(track->ErrorDZ, vertex->PositionError.Z());

        Double_t sigd0 = fabs(track->D0)/track->ErrorD0;
        if(fabs(dz) < fSigmaCompatibility && sigd0 < 3) //Should be done better...is assuming v_xy = 0 and the stat threatment is poor/wrong
        {
          UInt_t i = 0;
          UInt_t nv = vtx.size();
          while(i<nv)
          {
            if(vtx[i]->SumPT2 > vertex->SumPT2) i++;
            else break;
          }
          vtx.insert(vtx.begin()+i, vertex);
        }
      }

      if(fVerbose > 2)
      {
        cout << "Compatible vertexes: " << vtx.size() << endl;
        for(unsigned int i = 0; i<vtx.size(); i++)
        {
          cout << Form("SumPT2: %.0f GeV^2      Z: %.3f +/- %.3f mm      T: %.3f +/- %.3f ps", vtx[i]->SumPT2, vtx[i]->Position.Z(), vtx[i]->PositionError.Z(), vtx[i]->Position.T(), vtx[i]->PositionError.T()) << endl;
        }
      }

      if(vtx.size() > 0)
      {
        // Try at first to see if the mass hypothesys can work
        // The vertexes are order by SumPT2 and the mass hypotesys for relevance
        UInt_t i = 0, match = 0;
        while(i<vtx.size() && match == 0)
        {
          UInt_t j = 0;
          while(j<fMassList.size() && match == 0)
          {
            auto Tpair = ComputeCATime(track, fMassList[j]);
            double dt = vtx[i]->Position.T()*1.E9/c_light - Tpair.first;  // [ps]
            // dt /= Tpair.second;
            dt /= TMath::Hypot(Tpair.second, vtx[i]->PositionError.T()*1.E9/c_light);

            if(fabs(dt)<fSigmaCompatibility)
            {
              match = 1;
              track->Mass = fMassList[j];
              track->ClusterIndex = vtx[i]->ClusterIndex;
              track->InitialPosition.SetT(vtx[i]->Position.T());
              track->InitialPosition.SetZ(vtx[i]->Position.Z());
              track->Td = Tpair.first * 1E-9 * c_light;

              vtx[i]->SumPT2 += track->Momentum.Pt()*track->Momentum.Pt();
              vtx[i]->SumPt += track->Momentum.Pt();
              if(fVerbose>5)
              {
                cout << "SM VTX found at " << i << ". dt = " << dt << "sigma" << endl;
                cout << "vtx=> SumPT2 " << vtx[i]->SumPT2 << " - DOF: " << vtx[i]->ClusterNDF << endl;
                cout << "Mass: " << fMassList[j] << endl;
              }
            }

            j++;
          }

          i++;
        }

        if(match == 0)
        {
          Double_t p = track->Momentum.Pt() * sqrt(1 + track->CtgTheta*track->CtgTheta);
          // Only Z
          // Double_t beta_z = vtx[0]->Position.Z() - track->Position.Z();
          // beta_z /= vtx[0]->Position.T() - track->Position.T();
          // Double_t e = track->Momentum.Pt() * track->CtgTheta / beta_z;

          //Full path length
          Double_t beta = track->L / (vtx[0]->Position.T() - track->Position.T());
          Double_t e = p / beta;

          track->Mass = sqrt(e*e - p*p);
          track->ClusterIndex = vtx[0]->ClusterIndex;
          track->InitialPosition.SetT(vtx[0]->Position.T());
          track->InitialPosition.SetZ(vtx[0]->Position.Z());
          track->Td = vtx[0]->Position.T();

          vtx[0]->SumPT2 += track->Momentum.Pt()*track->Momentum.Pt();
          vtx[0]->SumPt += track->Momentum.Pt();
          if(fVerbose>5)
          {
            cout << "BSM VTX fitted:" << endl;
            cout << "vtx=> SumPT2 " << vtx[0]->SumPT2 << " - DOF: " << vtx[0]->ClusterNDF << endl;
            cout << "Mass: " << track->Mass << endl;
          }
        }

      }
    }

    fTrackOutputArray->Add(track);
  }
}

//------------------------------------------------------------------------------
// Auxiliary function to compute the CA time and erro given the mass
pair<Double_t, Double_t> HighMassVertexRecover::ComputeCATime(Candidate * tk, Double_t m)
{
  Double_t p = tk->Momentum.Pt() * sqrt(1 + tk->CtgTheta*tk->CtgTheta);
  Double_t e = sqrt(p*p + m*m);

  Double_t t = tk->Position.T()*1.E9/c_light; // from [mm] to [ps]
  //Full path length
  t -= tk->L*1E9/(c_light*p/e);

  // Only Z
  // Double_t bz = tk->Momentum.Pt() * tk->CtgTheta/e;
  // t += (tk->Zd - tk->Position.Z())*1E9/(c_light*bz);

  pair<Double_t, Double_t> out;
  out.first = t;
  out.second = tk->ErrorT*1.E9/c_light; // [ps]

  return out;
}
