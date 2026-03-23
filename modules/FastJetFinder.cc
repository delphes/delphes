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

/** \class FastJetFinder
 *
 *  Finds jets using FastJet library.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

#include <vector>

#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh>
#include <fastjet/contribs/Nsubjettiness/Nsubjettiness.hh>
#include <fastjet/contribs/RecursiveTools/SoftDrop.hh>
#include <fastjet/contribs/ValenciaPlugin/ValenciaPlugin.hh>
#include <fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh>
#include <fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh>
#include <fastjet/plugins/SISCone/fastjet/SISConePlugin.hh>
#include <fastjet/tools/Filter.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Pruner.hh>

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

class FastJetFinder: public DelphesModule
{
public:
  explicit FastJetFinder(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fJetAlgorithm(Steer<int>("JetAlgorithm", 6)),
    fParameterR(Steer<double>("ParameterR", 0.5)),
    fParameterP(Steer<double>("ParameterP", -1.0)),
    //
    fJetPTMin(Steer<double>("JetPTMin", 10.0)),
    fConeRadius(Steer<double>("ConeRadius", 0.5)),
    fSeedThreshold(Steer<double>("SeedThreshold", 1.0)),
    fConeAreaFraction(Steer<double>("ConeAreaFraction", 1.0)),
    fMaxIterations(Steer<int>("MaxIterations", 100)),
    fMaxPairSize(Steer<int>("MaxPairSize", 2)),
    fIratch(Steer<int>("Iratch", 1)),
    fAdjacencyCut(Steer<int>("AdjacencyCut", 2)),
    fOverlapThreshold(Steer<double>("OverlapThreshold", 0.75)),
    //-- Exclusive clustering for e+e- collisions --
    fNJets(Steer<int>("NJets", 2)),
    fExclusiveClustering(Steer<bool>("ExclusiveClustering", false)),
    fDCut(Steer<double>("DCut", 0.0)),
    //-- Valencia Linear Collider algorithm
    fGamma(Steer<double>("Gamma", 1.0)),
    //-- N(sub)jettiness parameters --
    fComputeNsubjettiness(Steer<bool>("ComputeNsubjettiness", false)),
    fBeta(Steer<double>("Beta", 1.0)),
    fAxisMode(Steer<int>("AxisMode", 1)),
    fRcutOff(Steer<double>("RcutOff", 0.8)), // used only if Njettiness is used as jet clustering algo (case 8)
    fN(Steer<int>("N", 2)), // used only if Njettiness is used as jet clustering algo (case 8)
    //-- Trimming parameters --
    fComputeTrimming(Steer<bool>("ComputeTrimming", false)),
    fRTrim(Steer<double>("RTrim", 0.2)),
    fPtFracTrim(Steer<double>("PtFracTrim", 0.05)),
    //-- Pruning parameters --
    fComputePruning(Steer<bool>("ComputePruning", false)),
    fZcutPrun(Steer<double>("ZcutPrun", 0.1)),
    fRcutPrun(Steer<double>("RcutPrun", 0.5)),
    fRPrun(Steer<double>("RPrun", 0.8)),
    //-- SoftDrop parameters --
    fComputeSoftDrop(Steer<bool>("ComputeSoftDrop", false)),
    fBetaSoftDrop(Steer<double>("BetaSoftDrop", 0.0)),
    fSymmetryCutSoftDrop(Steer<double>("SymmetryCutSoftDrop", 0.1)),
    fR0SoftDrop(Steer<double>("R0SoftDrop=", 0.8)),
    // ---  Jet Area Parameters ---
    fAreaAlgorithm(Steer<int>("AreaAlgorithm", 0)),
    fComputeRho(Steer<bool>("ComputeRho", false)),
    // - ghost based areas -
    fGhostEtaMax(Steer<double>("GhostEtaMax", 5.0)),
    fRepeat(Steer<int>("Repeat", 1)),
    fGhostArea(Steer<double>("GhostArea", 0.01)),
    fGridScatter(Steer<double>("GridScatter", 1.0)),
    fPtScatter(Steer<double>("PtScatter", 0.1)),
    fMeanGhostPt(Steer<double>("MeanGhostPt", 1.0E-100)),
    // - voronoi based areas -
    fEffectiveRfact(Steer<double>("EffectiveRfact", 1.0)),
    //
    fRhoEtaRange(Steer<std::vector<std::pair<double, double> > >("RhoEtaRange", {}))
  {
  }

  void Init() override;
  void Process() override;
  void Finish() override
  {
    fEstimators.clear();
  }

private:
  const int fJetAlgorithm;
  const double fParameterR;
  const double fParameterP;

  const double fJetPTMin;
  const double fConeRadius;
  const double fSeedThreshold;
  const double fConeAreaFraction;
  const int fMaxIterations;
  const int fMaxPairSize;
  const int fIratch;
  const int fAdjacencyCut;
  const double fOverlapThreshold;

  //-- Exclusive clustering for e+e- collisions --
  const int fNJets;
  const bool fExclusiveClustering;
  const double fDCut;

  //-- Valencia Linear Collider algorithm
  const double fGamma;

  //-- N (sub)jettiness parameters --
  const bool fComputeNsubjettiness;
  const double fBeta;
  const int fAxisMode;
  const double fRcutOff;
  const int fN;

  //-- Trimming parameters --
  const bool fComputeTrimming;
  const double fRTrim;
  const double fPtFracTrim;

  //-- Pruning parameters --
  const bool fComputePruning;
  const double fZcutPrun;
  const double fRcutPrun;
  const double fRPrun;

  //-- SoftDrop parameters --
  const bool fComputeSoftDrop;
  const double fBetaSoftDrop;
  const double fSymmetryCutSoftDrop;
  const double fR0SoftDrop;

  // --- FastJet Area method --------
  std::unique_ptr<fastjet::AreaDefinition> fAreaDefinition;
  const int fAreaAlgorithm;
  const bool fComputeRho;

  // -- ghost based areas --
  const double fGhostEtaMax;
  const int fRepeat;
  const double fGhostArea;
  const double fGridScatter;
  const double fPtScatter;
  const double fMeanGhostPt;

  // -- voronoi areas --
  const double fEffectiveRfact;

  const std::vector<std::pair<double, double> > fRhoEtaRange;

#if !defined(__CINT__) && !defined(__CLING__)
  struct TEstimatorStruct
  {
    std::unique_ptr<fastjet::JetMedianBackgroundEstimator> estimator;
    double etaMin, etaMax;
  };

  std::vector<TEstimatorStruct> fEstimators; //!
#endif

  std::unique_ptr<JetDefinition::Plugin> fPlugin; //!
  std::unique_ptr<JetDefinition::Recombiner> fRecomb; //!

  std::unique_ptr<fastjet::contrib::AxesDefinition> fAxesDef;
  std::unique_ptr<fastjet::contrib::MeasureDefinition> fMeasureDef;

  std::unique_ptr<fastjet::contrib::NjettinessPlugin> fNjettinessPlugin; //!
  std::unique_ptr<fastjet::contrib::ValenciaPlugin> fValenciaPlugin; //!
  std::unique_ptr<fastjet::JetDefinition> fDefinition; //!

  CandidatesCollection fInputArray; //!

  CandidatesCollection fOutputArray; //!
  CandidatesCollection fRhoOutputArray; //!
  CandidatesCollection fConstituentsOutputArray; //!
};

//------------------------------------------------------------------------------

void FastJetFinder::Init()
{
  TEstimatorStruct estimatorStruct;

  // define algorithm

  //fBeta parameter see above
  fMeasureDef = std::make_unique<NormalizedMeasure>(fBeta, fParameterR);

  switch(fAxisMode)
  {
  default:
  case 1:
    fAxesDef = std::make_unique<WTA_KT_Axes>();
    break;
  case 2:
    fAxesDef = std::make_unique<OnePass_WTA_KT_Axes>();
    break;
  case 3:
    fAxesDef = std::make_unique<KT_Axes>();
    break;
  case 4:
    fAxesDef = std::make_unique<OnePass_KT_Axes>();
  }

  switch(fAreaAlgorithm)
  {
  default:
  case 0:
    break;
  case 1:
    fAreaDefinition = std::make_unique<AreaDefinition>(active_area_explicit_ghosts, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
    break;
  case 2:
    fAreaDefinition = std::make_unique<AreaDefinition>(one_ghost_passive_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
    break;
  case 3:
    fAreaDefinition = std::make_unique<AreaDefinition>(passive_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
    break;
  case 4:
    fAreaDefinition = std::make_unique<AreaDefinition>(VoronoiAreaSpec(fEffectiveRfact));
    break;
  case 5:
    fAreaDefinition = std::make_unique<AreaDefinition>(active_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
    break;
  }

  switch(fJetAlgorithm)
  {
  case 1:
    fPlugin = std::make_unique<CDFJetCluPlugin>(fSeedThreshold, fConeRadius, fAdjacencyCut, fMaxIterations, fIratch, fOverlapThreshold);
    fDefinition = std::make_unique<JetDefinition>(fPlugin.get());
    break;
  case 2:
    fPlugin = std::make_unique<CDFMidPointPlugin>(fSeedThreshold, fConeRadius, fConeAreaFraction, fMaxPairSize, fMaxIterations, fOverlapThreshold);
    fDefinition = std::make_unique<JetDefinition>(fPlugin.get());
    break;
  case 3:
    fPlugin = std::make_unique<SISConePlugin>(fConeRadius, fOverlapThreshold, fMaxIterations, fJetPTMin);
    fDefinition = std::make_unique<JetDefinition>(fPlugin.get());
    break;
  case 4:
    fDefinition = std::make_unique<JetDefinition>(kt_algorithm, fParameterR);
    break;
  case 5:
    fDefinition = std::make_unique<JetDefinition>(cambridge_algorithm, fParameterR);
    break;
  default:
  case 6:
    fDefinition = std::make_unique<JetDefinition>(antikt_algorithm, fParameterR);
    break;
  case 7:
    fRecomb = std::make_unique<WinnerTakeAllRecombiner>();
    fDefinition = std::make_unique<JetDefinition>(antikt_algorithm, fParameterR, fRecomb.get(), Best);
    break;
  case 8:
    fNjettinessPlugin = std::make_unique<NjettinessPlugin>(fN, Njettiness::wta_kt_axes, Njettiness::unnormalized_cutoff_measure, fBeta, fRcutOff);
    fDefinition = std::make_unique<JetDefinition>(fNjettinessPlugin.get());
    break;
  case 9:
    fValenciaPlugin = std::make_unique<ValenciaPlugin>(fParameterR, fBeta, fGamma);
    fDefinition = std::make_unique<JetDefinition>(fValenciaPlugin.get());
    break;
  case 10:
    fDefinition = std::make_unique<JetDefinition>(ee_genkt_algorithm, fParameterR, fParameterP);
    break;

  // kT durham algorithm, 2 options:
  // 1. njets mode: stop when reach predetermined n jet (optionally apply sqrt(ExclYmerge(n-1,n))*Evis) > cut offline)
  // 2. dcut mode: stop when all dij above some threshold dcut. Is applied if fDCut > 0.
  case 11:
    fDefinition = std::make_unique<JetDefinition>(ee_kt_algorithm);
    break;
  }

  ClusterSequence::print_banner();

  if(fComputeRho && fAreaDefinition)
  {
    // read eta ranges
    fEstimators.clear();
    for(const std::pair<double, double> &etaMinMax : fRhoEtaRange)
    {
      TEstimatorStruct &estimatorStruct = fEstimators.emplace_back();
      estimatorStruct.estimator = std::make_unique<JetMedianBackgroundEstimator>(
        SelectorRapRange(etaMinMax.first, etaMinMax.second),
        *fDefinition, *fAreaDefinition);
      estimatorStruct.etaMin = etaMinMax.first;
      estimatorStruct.etaMax = etaMinMax.second;
    }
  }

  fInputArray = ImportArray(Steer<std::string>("InputArray", "Calorimeter/towers"));
  fOutputArray = ExportArray(Steer<std::string>("OutputArray", "jets"));
  fRhoOutputArray = ExportArray(Steer<std::string>("RhoOutputArray", "rho"));
  fConstituentsOutputArray = ExportArray(Steer<std::string>("ConstituentsOutputArray", "constituents"));
}

//------------------------------------------------------------------------------

void FastJetFinder::Process()
{
  fOutputArray->clear();
  fRhoOutputArray->clear();
  fConstituentsOutputArray->clear();

  Candidate *candidate = nullptr, *constituent = nullptr;
  TLorentzVector momentum;

  double deta, dphi, detaMax, dphiMax;
  double time, timeWeight;
  double neutralEnergyFraction, chargedEnergyFraction;

  int number, ncharged, nneutrals;
  int charge;
  double rho = 0.0;
  PseudoJet jet, area;
  std::unique_ptr<ClusterSequence> sequence;
  vector<PseudoJet> inputList, outputList, subjets;
  vector<PseudoJet>::iterator itInputList, itOutputList;
  vector<TEstimatorStruct>::iterator itEstimators;
  double excl_ymerge12 = 0.0;
  double excl_ymerge23 = 0.0;
  double excl_ymerge34 = 0.0;
  double excl_ymerge45 = 0.0;
  double excl_ymerge56 = 0.0;

  DelphesFactory *factory = GetFactory();

  inputList.clear();

  // loop over input objects
  number = 0;
  for(Candidate *const &candidate : *fInputArray)
  {
    momentum = candidate->Momentum;
    jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
    jet.set_user_index(number);
    inputList.push_back(jet);
    ++number;
  }

  // construct jets
  if(fAreaDefinition)
    sequence = std::make_unique<ClusterSequenceArea>(inputList, *fDefinition, *fAreaDefinition);
  else
    sequence = std::make_unique<ClusterSequence>(inputList, *fDefinition);

  // compute rho and store it
  if(fComputeRho && fAreaDefinition)
  {
    for(itEstimators = fEstimators.begin(); itEstimators != fEstimators.end(); ++itEstimators)
    {
      itEstimators->estimator->set_particles(inputList);
      rho = itEstimators->estimator->rho();

      candidate = factory->NewCandidate();
      candidate->Momentum.SetPtEtaPhiE(rho, 0.0, 0.0, rho);
      candidate->Edges[0] = itEstimators->etaMin;
      candidate->Edges[1] = itEstimators->etaMax;
      fRhoOutputArray->emplace_back(candidate);
    }
  }

  outputList.clear();

  if(fExclusiveClustering)
  {
    try
    {
      if(fDCut > 0.0) // exclusive dcut mode
        outputList = sorted_by_pt(sequence->exclusive_jets(fDCut * fDCut));
      else // exclusive njet mode
        outputList = sorted_by_pt(sequence->exclusive_jets(fNJets));
    }
    catch(const fastjet::Error &)
    {
      outputList.clear();
    }

    excl_ymerge12 = sequence->exclusive_ymerge(1);
    excl_ymerge23 = sequence->exclusive_ymerge(2);
    excl_ymerge34 = sequence->exclusive_ymerge(3);
    excl_ymerge45 = sequence->exclusive_ymerge(4);
    excl_ymerge56 = sequence->exclusive_ymerge(5);
  }
  else
  {
    outputList = sorted_by_pt(sequence->inclusive_jets(fJetPTMin));
  }

  // loop over all jets and export them
  detaMax = 0.0;
  dphiMax = 0.0;

  for(itOutputList = outputList.begin(); itOutputList != outputList.end(); ++itOutputList)
  {
    jet = *itOutputList;
    if(fJetAlgorithm == 7) jet = join(jet.constituents());

    momentum.SetPxPyPzE(jet.px(), jet.py(), jet.pz(), jet.E());

    area.reset(0.0, 0.0, 0.0, 0.0);
    if(fAreaDefinition) area = itOutputList->area_4vector();

    candidate = factory->NewCandidate();

    time = 0.0;
    timeWeight = 0.0;

    charge = 0;

    ncharged = 0;
    nneutrals = 0;

    neutralEnergyFraction = 0.;
    chargedEnergyFraction = 0.;

    inputList.clear();
    inputList = sequence->constituents(*itOutputList);

    for(itInputList = inputList.begin(); itInputList != inputList.end(); ++itInputList)
    {
      if(itInputList->user_index() < 0) continue;
      constituent = static_cast<Candidate *>(fInputArray->at(itInputList->user_index()));

      deta = std::fabs(momentum.Eta() - constituent->Momentum.Eta());
      dphi = std::fabs(momentum.DeltaPhi(constituent->Momentum));
      if(deta > detaMax) detaMax = deta;
      if(dphi > dphiMax) dphiMax = dphi;

      if(constituent->Charge == 0)
      {
        nneutrals++;
        neutralEnergyFraction += constituent->Momentum.E();
      }
      else
      {
        ncharged++;
        chargedEnergyFraction += constituent->Momentum.E();
      }

      time += std::sqrt(constituent->Momentum.E()) * (constituent->Position.T());
      timeWeight += std::sqrt(constituent->Momentum.E());

      charge += constituent->Charge;

      fConstituentsOutputArray->emplace_back(constituent);
      candidate->AddCandidate(constituent);
    }

    candidate->Momentum = momentum;
    candidate->Position.SetT(time / timeWeight);
    candidate->Area.SetPxPyPzE(area.px(), area.py(), area.pz(), area.E());

    candidate->DeltaEta = detaMax;
    candidate->DeltaPhi = dphiMax;
    candidate->Charge = charge;
    candidate->NNeutrals = nneutrals;
    candidate->NCharged = ncharged;

    candidate->NeutralEnergyFraction = (momentum.E() > 0) ? neutralEnergyFraction / momentum.E() : 0.0;
    candidate->ChargedEnergyFraction = (momentum.E() > 0) ? chargedEnergyFraction / momentum.E() : 0.0;

    //for exclusive clustering, access y_n,n+1 as exclusive_ymerge (fNJets);
    candidate->ExclYmerge12 = excl_ymerge12;
    candidate->ExclYmerge23 = excl_ymerge23;
    candidate->ExclYmerge34 = excl_ymerge34;
    candidate->ExclYmerge45 = excl_ymerge45;
    candidate->ExclYmerge56 = excl_ymerge56;

    //------------------------------------
    // Trimming
    //------------------------------------

    if(fComputeTrimming)
    {

      fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, fRTrim), fastjet::SelectorPtFractionMin(fPtFracTrim));
      fastjet::PseudoJet trimmed_jet = trimmer(*itOutputList);

      candidate->TrimmedP4[0].SetPtEtaPhiM(trimmed_jet.pt(), trimmed_jet.eta(), trimmed_jet.phi(), trimmed_jet.m());

      // four hardest subjets
      subjets.clear();
      subjets = trimmed_jet.pieces();
      subjets = sorted_by_pt(subjets);

      candidate->NSubJetsTrimmed = subjets.size();

      for(size_t i = 0; i < subjets.size() and i < 4; i++)
      {
        if(subjets.at(i).pt() < 0) continue;
        candidate->TrimmedP4[i + 1].SetPtEtaPhiM(subjets.at(i).pt(), subjets.at(i).eta(), subjets.at(i).phi(), subjets.at(i).m());
      }
    }

    //------------------------------------
    // Pruning
    //------------------------------------

    if(fComputePruning)
    {

      fastjet::Pruner pruner(fastjet::JetDefinition(fastjet::cambridge_algorithm, fRPrun), fZcutPrun, fRcutPrun);
      fastjet::PseudoJet pruned_jet = pruner(*itOutputList);

      candidate->PrunedP4[0].SetPtEtaPhiM(pruned_jet.pt(), pruned_jet.eta(), pruned_jet.phi(), pruned_jet.m());

      // four hardest subjet
      subjets.clear();
      subjets = pruned_jet.pieces();
      subjets = sorted_by_pt(subjets);

      candidate->NSubJetsPruned = subjets.size();

      for(size_t i = 0; i < subjets.size() and i < 4; i++)
      {
        if(subjets.at(i).pt() < 0) continue;
        candidate->PrunedP4[i + 1].SetPtEtaPhiM(subjets.at(i).pt(), subjets.at(i).eta(), subjets.at(i).phi(), subjets.at(i).m());
      }
    }

    //------------------------------------
    // SoftDrop
    //------------------------------------

    if(fComputeSoftDrop)
    {
      contrib::SoftDrop softDrop(fBetaSoftDrop, fSymmetryCutSoftDrop, fR0SoftDrop);
      fastjet::PseudoJet softdrop_jet = softDrop(*itOutputList);

      candidate->SoftDroppedP4[0].SetPtEtaPhiM(softdrop_jet.pt(), softdrop_jet.eta(), softdrop_jet.phi(), softdrop_jet.m());

      // four hardest subjet

      subjets.clear();
      subjets = softdrop_jet.pieces();
      subjets = sorted_by_pt(subjets);
      candidate->NSubJetsSoftDropped = softdrop_jet.pieces().size();

      candidate->SoftDroppedJet = candidate->SoftDroppedP4[0];

      for(size_t i = 0; i < subjets.size() and i < 4; i++)
      {
        if(subjets.at(i).pt() < 0) continue;
        candidate->SoftDroppedP4[i + 1].SetPtEtaPhiM(subjets.at(i).pt(), subjets.at(i).eta(), subjets.at(i).phi(), subjets.at(i).m());
        if(i == 0) candidate->SoftDroppedSubJet1 = candidate->SoftDroppedP4[i + 1];
        if(i == 1) candidate->SoftDroppedSubJet2 = candidate->SoftDroppedP4[i + 1];
      }
    }

    // --- compute N-subjettiness with N = 1,2,3,4,5 ----

    if(fComputeNsubjettiness)
    {

      Nsubjettiness nSub1(1, *fAxesDef, *fMeasureDef);
      Nsubjettiness nSub2(2, *fAxesDef, *fMeasureDef);
      Nsubjettiness nSub3(3, *fAxesDef, *fMeasureDef);
      Nsubjettiness nSub4(4, *fAxesDef, *fMeasureDef);
      Nsubjettiness nSub5(5, *fAxesDef, *fMeasureDef);

      candidate->Tau[0] = nSub1(*itOutputList);
      candidate->Tau[1] = nSub2(*itOutputList);
      candidate->Tau[2] = nSub3(*itOutputList);
      candidate->Tau[3] = nSub4(*itOutputList);
      candidate->Tau[4] = nSub5(*itOutputList);
    }

    fOutputArray->emplace_back(candidate);
  }
}

REGISTER_MODULE("FastJetFinder", FastJetFinder);
