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
  explicit PileUpJetID(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fJetPTMin(Steer<double>("JetPTMin", 20.0)),
    fParameterR(Steer<double>("ParameterR", 0.5)),
    //
    fUseConstituents(Steer<int>("UseConstituents", 0)),
    fMeanSqDeltaRMaxBarrel(Steer<double>("MeanSqDeltaRMaxBarrel", 0.1)),
    fBetaMinBarrel(Steer<double>("BetaMinBarrel", 0.1)),
    fMeanSqDeltaRMaxEndcap(Steer<double>("MeanSqDeltaRMaxEndcap", 0.1)),
    fBetaMinEndcap(Steer<double>("BetaMinEndcap", 0.1)),
    fMeanSqDeltaRMaxForward(Steer<double>("MeanSqDeltaRMaxForward", 0.1)),
    //
    fJetPTMinForNeutrals(Steer<double>("JetPTMinForNeutrals", 20.0)),
    fNeutralPTMin(Steer<double>("NeutralPTMin", 2.0))
  {
  }

  void Init() override
  {
    fJetInputArray = ImportArray(Steer<std::string>("JetInputArray", "FastJetFinder/jets"));
    fTrackInputArray = ImportArray(Steer<std::string>("TrackInputArray", "ParticlePropagator/tracks"));
    fNeutralInputArray = ImportArray(Steer<std::string>("NeutralInputArray", "ParticlePropagator/tracks"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "jets"));
    fNeutralsInPassingJets = ExportArray(Steer<std::string>("NeutralsInPassingJets", "eflowtowers"));
  }
  void Process() override;

private:
  const double fJetPTMin;
  const double fParameterR;

  // If set to true, may have weird results for PFCHS
  // If set to false, uses everything within dR < fParameterR even if in other jets &c.
  // Results should be very similar for PF
  const int fUseConstituents;

  const double fMeanSqDeltaRMaxBarrel; // |eta| < 1.5
  const double fBetaMinBarrel; // |eta| < 2.5
  const double fMeanSqDeltaRMaxEndcap; // 1.5 < |eta| < 4.0
  const double fBetaMinEndcap; // 1.5 < |eta| < 4.0
  const double fMeanSqDeltaRMaxForward; // |eta| > 4.0

  const double fJetPTMinForNeutrals;
  const double fNeutralPTMin;

  CandidatesCollection fJetInputArray; //!
  CandidatesCollection fTrackInputArray; // SCZ
  CandidatesCollection fNeutralInputArray;

  CandidatesCollection fOutputArray; //!
  CandidatesCollection fNeutralsInPassingJets; // SCZ

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

  const bool fAverageEachTower{false}; // for timing
};

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
    fMeanSqDeltaRMaxBarrel = Steer<double>("MeanSqDeltaRMaxBarrel",0.1);
    fBetaMinBarrel = Steer<double>("BetaMinBarrel",0.1);
    fMeanSqDeltaRMaxEndcap = Steer<double>("MeanSqDeltaRMaxEndcap",0.1);
    fBetaMinEndcap = Steer<double>("BetaMinEndcap",0.1);
    fMeanSqDeltaRMaxForward = Steer<double>("MeanSqDeltaRMaxForward",0.1);
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
