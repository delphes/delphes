/** \class PileUpJetID
 *
 *  CMS PileUp Jet ID Variables
 *
 *  \author S. Zenz
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

using namespace std;

class PileUpJetID: public DelphesModule
{
public:
  PileUpJetID() = default;

  void Init() override;
  void Process() override;

private:
  Double_t fJetPTMin;
  Double_t fParameterR;

  Double_t fMeanSqDeltaRMaxBarrel; // |eta| < 1.5
  Double_t fBetaMinBarrel; // |eta| < 2.5
  Double_t fMeanSqDeltaRMaxEndcap; // 1.5 < |eta| < 4.0
  Double_t fBetaMinEndcap; // 1.5 < |eta| < 4.0
  Double_t fMeanSqDeltaRMaxForward; // |eta| > 4.0

  Double_t fNeutralPTMin;
  Double_t fJetPTMinForNeutrals;

  /*
JAY
---

|Eta|<1.5

meanSqDeltaR betaStar SigEff BgdEff
0.13 0.92 96% 8%
0.13 0.95 97% 16%
0.13 0.97 98% 27%

|Eta|>1.5

meanSqDeltaR betaStar SigEff BgdEff
0.14 0.91 95% 15%
0.14 0.94 97% 19%
0.14 0.97 98% 29%

BRYAN
-----

Barrel (MeanSqDR, Beta, sig eff, bg eff):
0.10, 0.08, 90%, 8%
0.11, 0.12, 90%, 6%
0.13, 0.16, 89%, 5%

Endcap (MeanSqDR, Beta, sig eff, bg eff):
0.07, 0.06, 89%, 4%
0.08, 0.08, 92%, 6%
0.09, 0.08, 95%, 10%
0.10, 0.08, 97%, 13%

SETH GUESSES FOR |eta| > 4.0
----------------------------

MeanSqDeltaR
0.07
0.10
0.14
0.2
  */

  // If set to true, may have weird results for PFCHS
  // If set to false, uses everything within dR < fParameterR even if in other jets &c.
  // Results should be very similar for PF
  Int_t fUseConstituents;

  Bool_t fAverageEachTower;

  CandidatesCollection fJetInputArray; //!
  CandidatesCollection fTrackInputArray; // SCZ
  CandidatesCollection fNeutralInputArray;

  CandidatesCollection fOutputArray; //!
  CandidatesCollection fNeutralsInPassingJets; // SCZ
};

//------------------------------------------------------------------------------

void PileUpJetID::Init()
{
  fJetPTMin = GetDouble("JetPTMin", 20.0);
  fParameterR = GetDouble("ParameterR", 0.5);
  fUseConstituents = GetInt("UseConstituents", 0);

  fMeanSqDeltaRMaxBarrel = GetDouble("MeanSqDeltaRMaxBarrel", 0.1);
  fBetaMinBarrel = GetDouble("BetaMinBarrel", 0.1);
  fMeanSqDeltaRMaxEndcap = GetDouble("MeanSqDeltaRMaxEndcap", 0.1);
  fBetaMinEndcap = GetDouble("BetaMinEndcap", 0.1);
  fMeanSqDeltaRMaxForward = GetDouble("MeanSqDeltaRMaxForward", 0.1);
  fJetPTMinForNeutrals = GetDouble("JetPTMinForNeutrals", 20.0);
  fNeutralPTMin = GetDouble("NeutralPTMin", 2.0);

  fAverageEachTower = false; // for timing

  // import input arrays
  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fTrackInputArray = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fNeutralInputArray = ImportArray(GetString("NeutralInputArray", "ParticlePropagator/tracks"));

  // create output arrays
  fOutputArray = ExportArray(GetString("OutputArray", "jets"));
  fNeutralsInPassingJets = ExportArray(GetString("NeutralsInPassingJets", "eflowtowers"));
}

//------------------------------------------------------------------------------

void PileUpJetID::Process()
{
  fOutputArray->clear();
  fNeutralsInPassingJets->clear();

  TLorentzVector momentum, area;

  // loop over all input candidates
  for(Candidate *const &candidate : *fJetInputArray)
  {
    momentum = candidate->Momentum;
    area = candidate->Area;

    candidate->NTimeHits = 0;

    float sumpt = 0.;
    float sumptch = 0.;
    float sumptchpv = 0.;
    float sumptchpu = 0.;
    float sumdrsqptsq = 0.;
    float sumptsq = 0.;
    int nc = 0;
    int nn = 0;
    float pt_ann[5];

    for(int i = 0; i < 5; i++)
    {
      pt_ann[i] = 0.;
    }

    if(fUseConstituents)
    {
      for(Candidate *const &constituent : candidate->GetCandidates())
      {
        float pt = constituent->Momentum.Pt();
        float dr = candidate->Momentum.DeltaR(constituent->Momentum);
        //	cout << " There exists a constituent with dr=" << dr << endl;
        sumpt += pt;
        sumdrsqptsq += dr * dr * pt * pt;
        sumptsq += pt * pt;
        if(constituent->Charge == 0)
        {
          nn++;
        }
        else
        {
          if(constituent->IsRecoPU)
          {
            sumptchpu += pt;
          }
          else
          {
            sumptchpv += pt;
          }
          sumptch += pt;
          nc++;
        }
        for(int i = 0; i < 5; i++)
        {
          if(dr > 0.1 * i && dr < 0.1 * (i + 1))
          {
            pt_ann[i] += pt;
          }
        }
        float tow_sumW = 0;
        for(size_t i = 0; i < constituent->ECalEnergyTimePairs.size(); i++)
        {
          float w = std::sqrt(constituent->ECalEnergyTimePairs[i].first);
          if(fAverageEachTower)
          {
            tow_sumW += w;
          }
          else
          {
            candidate->NTimeHits++;
          }
        }
        if(fAverageEachTower && tow_sumW > 0.)
        {
          candidate->NTimeHits++;
        }
      }
    }
    else
    {
      // Not using constituents, using dr
      for(Candidate *const &trk : *fTrackInputArray)
      {
        if(trk->Momentum.DeltaR(candidate->Momentum) < fParameterR)
        {
          float pt = trk->Momentum.Pt();
          sumpt += pt;
          sumptch += pt;
          if(trk->IsRecoPU)
          {
            sumptchpu += pt;
          }
          else
          {
            sumptchpv += pt;
          }
          float dr = candidate->Momentum.DeltaR(trk->Momentum);
          sumdrsqptsq += dr * dr * pt * pt;
          sumptsq += pt * pt;
          nc++;
          for(int i = 0; i < 5; i++)
          {
            if(dr > 0.1 * i && dr < 0.1 * (i + 1))
            {
              pt_ann[i] += pt;
            }
          }
        }
      }
      for(Candidate *const &constituent : *fNeutralInputArray)
      {
        if(constituent->Momentum.DeltaR(candidate->Momentum) < fParameterR)
        {
          float pt = constituent->Momentum.Pt();
          sumpt += pt;
          float dr = candidate->Momentum.DeltaR(constituent->Momentum);
          sumdrsqptsq += dr * dr * pt * pt;
          sumptsq += pt * pt;
          nn++;
          for(int i = 0; i < 5; i++)
          {
            if(dr > 0.1 * i && dr < 0.1 * (i + 1))
            {
              pt_ann[i] += pt;
            }
          }
        }
      }
    }

    if(sumptch > 0.)
    {
      candidate->Beta = sumptchpv / sumptch;
      candidate->BetaStar = sumptchpu / sumptch;
    }
    else
    {
      candidate->Beta = -999.;
      candidate->BetaStar = -999.;
    }
    if(sumptsq > 0.)
    {
      candidate->MeanSqDeltaR = sumdrsqptsq / sumptsq;
    }
    else
    {
      candidate->MeanSqDeltaR = -999.;
    }
    candidate->NCharged = nc;
    candidate->NNeutrals = nn;
    if(sumpt > 0.)
    {
      candidate->PTD = std::sqrt(sumptsq) / sumpt;
      for(int i = 0; i < 5; i++)
      {
        candidate->FracPt[i] = pt_ann[i] / sumpt;
      }
    }
    else
    {
      candidate->PTD = -999.;
      for(int i = 0; i < 5; i++)
      {
        candidate->FracPt[i] = -999.;
      }
    }

    fOutputArray->emplace_back(candidate);

    // New stuff
    /*
    fMeanSqDeltaRMaxBarrel = GetDouble("MeanSqDeltaRMaxBarrel",0.1);
    fBetaMinBarrel = GetDouble("BetaMinBarrel",0.1);
    fMeanSqDeltaRMaxEndcap = GetDouble("MeanSqDeltaRMaxEndcap",0.1);
    fBetaMinEndcap = GetDouble("BetaMinEndcap",0.1);
    fMeanSqDeltaRMaxForward = GetDouble("MeanSqDeltaRMaxForward",0.1);
    */

    bool passId = false;
    if(candidate->Momentum.Pt() > fJetPTMinForNeutrals && candidate->MeanSqDeltaR > -0.1)
    {
      if(fabs(candidate->Momentum.Eta()) < 1.5)
      {
        passId = ((candidate->Beta > fBetaMinBarrel) && (candidate->MeanSqDeltaR < fMeanSqDeltaRMaxBarrel));
      }
      else if(fabs(candidate->Momentum.Eta()) < 4.0)
      {
        passId = ((candidate->Beta > fBetaMinEndcap) && (candidate->MeanSqDeltaR < fMeanSqDeltaRMaxEndcap));
      }
      else
      {
        passId = (candidate->MeanSqDeltaR < fMeanSqDeltaRMaxForward);
      }
    }

    //    cout << " Pt Eta MeanSqDeltaR Beta PassId " << candidate->Momentum.Pt()
    //	 << " " << candidate->Momentum.Eta() << " " << candidate->MeanSqDeltaR << " " << candidate->Beta << " " << passId << endl;

    if(passId)
    {
      if(fUseConstituents)
      {
        for(Candidate *const &constituent : candidate->GetCandidates())
        {
          if(constituent->Charge == 0 && constituent->Momentum.Pt() > fNeutralPTMin)
          {
            fNeutralsInPassingJets->emplace_back(constituent);
            //	    cout << "    Constitutent added Pt Eta Charge " << constituent->Momentum.Pt() << " " << constituent->Momentum.Eta() << " " << constituent->Charge << endl;
          }
        }
      }
      else
      { // use DeltaR
        for(Candidate *const &constituent : *fNeutralInputArray)
        {
          if(constituent->Momentum.DeltaR(candidate->Momentum) < fParameterR && constituent->Momentum.Pt() > fNeutralPTMin)
          {
            fNeutralsInPassingJets->emplace_back(constituent);
            //            cout << "    Constitutent added Pt Eta Charge " << constituent->Momentum.Pt() << " " << constituent->Momentum.Eta() << " " << constituent->Charge << endl;
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("PileUpJetID", PileUpJetID);
