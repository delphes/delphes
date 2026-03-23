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

/** \class PhotonID
 *
 *  Applies complex photon Id. Reconstructed photon candidtes are first separated into matched and non-matched to gen particles.
 *  Non-matched pass the "fake" efficiency. Matched photons get further splitted into isolated and non-isolated (user can choose criterion for isolation)
 *  Isolated photons pass the "prompt" efficiency while the non-isolated pass the "non-prompt" efficiency
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

using namespace std;

class PhotonID: public DelphesModule
{
public:
  explicit PhotonID(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fPTMin(Steer<double>("PTMin", 10.0)), // min pt to be considered, make sure this threshold is higher than threshold in particle filter
    fRelIsoMax(Steer<double>("fRelIsoMax", 0.3)), // to be tuned, since FS and delphes have different isolation profiles
    fPromptFormula(std::make_unique<DelphesFormula>()),
    fNonPromptFormula(std::make_unique<DelphesFormula>()),
    fFakeFormula(std::make_unique<DelphesFormula>())
  {
    // read PhotonID formulae
    fPromptFormula->Compile(Steer<std::string>("PromptFormula", "1.0"));
    fNonPromptFormula->Compile(Steer<std::string>("NonPromptFormula", "1.0"));
    fFakeFormula->Compile(Steer<std::string>("FakeFormula", "1.0"));
  }

  void Init() override
  {
    fInputPhotonArray = ImportArray(Steer<std::string>("InputPhotonArray", "PhotonIsolation/photons"));
    fInputGenArray = ImportArray(Steer<std::string>("InputGenArray", "GenParticleFilter/filteredParticles")); // use filtered collection for speed
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "photons"));
  }
  void Process() override;

private:
  bool isFake(const Candidate *obj);

  const double fPTMin;
  const double fRelIsoMax;

  const std::unique_ptr<DelphesFormula> fPromptFormula;
  const std::unique_ptr<DelphesFormula> fNonPromptFormula;
  const std::unique_ptr<DelphesFormula> fFakeFormula;

  CandidatesCollection fInputPhotonArray;
  CandidatesCollection fInputGenArray; // use filtered collection for speed

  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void PhotonID::Process()
{
  fOutputArray->clear();

  //cout<< "----  new event ---------"<<endl;

  for(Candidate *const &candidate : *fInputPhotonArray)
  {
    Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
    new_candidate->AddCandidate(candidate);

    const TLorentzVector &candidatePosition = new_candidate->Position;
    const TLorentzVector &candidateMomentum = new_candidate->Momentum;
    const double eta = candidatePosition.Eta();
    const double phi = candidatePosition.Phi();
    const double pt = candidateMomentum.Pt();
    const double e = candidateMomentum.E();

    if(pt < fPTMin) continue;

    //cout<< "              ---- photon -----: "<<pt<<","<<eta<<","<<phi<<endl;

    // find out if photon matches does not match photon in gen collection and apply fae efficiency
    if(isFake(new_candidate))
    {
      //cout<<"                    Fake!"<<endl;

      if(gRandom->Uniform() > fFakeFormula->Eval(pt, eta, phi, e)) continue;
      //cout<<"                    passed"<<endl;
      new_candidate->Status = 3;
      fOutputArray->emplace_back(new_candidate);
    }

    // if matches photon in gen collection
    else
    {
      const double relIso = new_candidate->IsolationVar;
      const bool isolated = (relIso < 0.3);
      //cout<<"                    Prompt!:   "<<relIso<<endl;

      // if isolated apply prompt formula
      if(isolated)
      {
        //cout<<"                       isolated!:   "<<relIso<<endl;
        if(gRandom->Uniform() > fPromptFormula->Eval(pt, eta, phi, e)) continue;
        //cout<<"                       passed"<<endl;
        new_candidate->Status = 1;
        fOutputArray->emplace_back(new_candidate);
      }

      // if non-isolated apply non-prompt formula
      else
      {
        //cout<<"                       non-isolated!:   "<<relIso<<endl;
        if(gRandom->Uniform() > fNonPromptFormula->Eval(pt, eta, phi, e)) continue;
        //cout<<"                       passed"<<endl;
        new_candidate->Status = 2;
        fOutputArray->emplace_back(new_candidate);
      }
    }
  }
}

//------------------------------------------------------------------------------

bool PhotonID::isFake(const Candidate *obj)
{
  const TLorentzVector &mom_rec = obj->Momentum;

  bool matches = false;
  for(Candidate *const &gen : *fInputGenArray)
  {
    const TLorentzVector &mom_gen = gen->Momentum;
    int status = gen->Status;
    int pdgCode = std::abs(gen->PID);
    float dPtOverPt = std::fabs((mom_gen.Pt() - mom_rec.Pt()) / mom_rec.Pt());
    float deltaR = mom_gen.DeltaR(mom_rec);

    if(status != 1) continue;
    if(pdgCode != 22) continue;
    if(dPtOverPt > 0.5) continue;
    if(deltaR > 0.1) continue;

    matches = true;
    break;
  }

  return !matches;
}

//------------------------------------------------------------------------------

REGISTER_MODULE("PhotonID", PhotonID);
