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

/** \class JetFakeParticle
 *
 *  Converts jet into particle with some PID,
 *  according to parametrized probability.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

using namespace std;

class JetFakeParticle: public DelphesModule
{
public:
  explicit JetFakeParticle(const DelphesParameters &moduleParams) : DelphesModule(moduleParams)
  {
    for(const std::pair<int, std::string> &efficiencyFormula : Steer<std::vector<std::pair<int, std::string> > >("EfficiencyFormula"))
    {
      if(const int pdgCode = efficiencyFormula.first; std::abs(pdgCode) != 11 && std::abs(pdgCode) != 13 && std::abs(pdgCode) != 22)
        throw runtime_error("Jets can only fake into electrons, muons or photons. Other particles are not authorized.");
      else
      {
        std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
        formula->Compile(efficiencyFormula.second);
        fEfficiencyMap[pdgCode] = std::move(formula);
      }
    }
    // set default efficiency formula
    if(fEfficiencyMap.count(0) == 0)
    {
      std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
      formula->Compile("0.0");
      fEfficiencyMap[0] = std::move(formula);
    }
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "FastJetFinder/jets"));
    fElectronOutputArray = ExportArray(Steer<std::string>("ElectronOutputArray", "fakeElectrons"));
    fMuonOutputArray = ExportArray(Steer<std::string>("MuonOutputArray", "fakeMuons"));
    fPhotonOutputArray = ExportArray(Steer<std::string>("PhotonOutputArray", "fakePhotons"));
    fJetOutputArray = ExportArray(Steer<std::string>("JetOutputArray", "jets"));
  }
  void Process() override;

private:
#if !defined(__CINT__) && !defined(__CLING__)
  typedef std::map<int, std::unique_ptr<DelphesFormula> > TFakeMap; //!
  TFakeMap fEfficiencyMap;
#endif

  CandidatesCollection fInputArray; //!

  CandidatesCollection fElectronOutputArray; //!
  CandidatesCollection fMuonOutputArray; //!
  CandidatesCollection fPhotonOutputArray; //!
  CandidatesCollection fJetOutputArray; //!
};

//------------------------------------------------------------------------------

void JetFakeParticle::Process()
{
  fElectronOutputArray->clear();
  fMuonOutputArray->clear();
  fPhotonOutputArray->clear();
  fJetOutputArray->clear();

  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    const double eta = candidateMomentum.Eta();
    const double phi = candidateMomentum.Phi();
    const double pt = candidateMomentum.Pt();
    const double e = candidateMomentum.E();

    const double r = gRandom->Uniform();
    double total = 0.;
    Candidate *fake = nullptr;

    // loop over map for this jet
    for(const std::pair<const int, std::unique_ptr<DelphesFormula> > &efficiencyMap : fEfficiencyMap)
    {
      const int pdgCodeOut = efficiencyMap.first;
      const double p = efficiencyMap.second->Eval(pt, eta, phi, e);

      if(total <= r && r < total + p)
      {
        fake = static_cast<Candidate *>(candidate->Clone());

        // convert jet

        if(std::abs(pdgCodeOut) == 11 || std::abs(pdgCodeOut) == 13)
        {
          if(candidate->Charge != 0)
            fake->Charge = candidate->Charge / std::abs(candidate->Charge);
          else
          {
            const double rs = gRandom->Uniform();
            fake->Charge = (rs < 0.5) ? -1 : 1;
          }
        }

        if(std::abs(pdgCodeOut) == 22) fake->PID = 22;

        if(std::abs(pdgCodeOut) == 11) fElectronOutputArray->emplace_back(fake);
        if(std::abs(pdgCodeOut) == 13) fMuonOutputArray->emplace_back(fake);
        if(std::abs(pdgCodeOut) == 22) fPhotonOutputArray->emplace_back(fake);

        break;
      }
      total += p;
    }

    if(!fake) fJetOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("JetFakeParticle", JetFakeParticle);
