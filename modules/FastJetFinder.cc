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

#include "modules/FastJetFinder.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh"
#include "fastjet/plugins/SISCone/fastjet/SISConePlugin.hh"

#include "fastjet/contribs/Nsubjettiness/ExtraRecombiners.hh"
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"

#include "fastjet/contribs/ValenciaPlugin/ValenciaPlugin.hh"

#include "fastjet/contribs/RecursiveTools/SoftDrop.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

//------------------------------------------------------------------------------

FastJetFinder::FastJetFinder() :
  fPlugin(0), fRecomb(0), fAxesDef(0), fMeasureDef(0), fNjettinessPlugin(0), fValenciaPlugin(0),
  fDefinition(0), fAreaDefinition(0), fItInputArray(0)
{
}

//------------------------------------------------------------------------------

FastJetFinder::~FastJetFinder()
{
}

//------------------------------------------------------------------------------

void FastJetFinder::Init()
{
  JetDefinition::Plugin *plugin = 0;
  JetDefinition::Recombiner *recomb = 0;
  ExRootConfParam param;
  Long_t i, size;
  Double_t etaMin, etaMax;
  TEstimatorStruct estimatorStruct;

  // define algorithm

  fJetAlgorithm = GetInt("JetAlgorithm", 6);
  fParameterR = GetDouble("ParameterR", 0.5);
  fParameterP = GetDouble("ParameterP", -1.0);

  fConeRadius = GetDouble("ConeRadius", 0.5);
  fSeedThreshold = GetDouble("SeedThreshold", 1.0);
  fConeAreaFraction = GetDouble("ConeAreaFraction", 1.0);
  fMaxIterations = GetInt("MaxIterations", 100);
  fMaxPairSize = GetInt("MaxPairSize", 2);
  fIratch = GetInt("Iratch", 1);
  fAdjacencyCut = GetInt("AdjacencyCut", 2);
  fOverlapThreshold = GetDouble("OverlapThreshold", 0.75);

  fJetPTMin = GetDouble("JetPTMin", 10.0);

  //-- N(sub)jettiness parameters --

  fComputeNsubjettiness = GetBool("ComputeNsubjettiness", false);
  fBeta = GetDouble("Beta", 1.0);
  fAxisMode = GetInt("AxisMode", 1);
  fRcutOff = GetDouble("RcutOff", 0.8); // used only if Njettiness is used as jet clustering algo (case 8)
  fN = GetInt("N", 2); // used only if Njettiness is used as jet clustering algo (case 8)

  //-- Exclusive clustering for e+e- collisions --

  fNJets = GetInt("NJets", 2);
  fExclusiveClustering = GetBool("ExclusiveClustering", false);
  fDCut = GetDouble("DCut", 0.0);

  //-- Valencia Linear Collider algorithm

  fGamma = GetDouble("Gamma", 1.0);
  //fBeta parameter see above

  fMeasureDef = new NormalizedMeasure(fBeta, fParameterR);

  switch(fAxisMode)
  {
  default:
  case 1:
    fAxesDef = new WTA_KT_Axes();
    break;
  case 2:
    fAxesDef = new OnePass_WTA_KT_Axes();
    break;
  case 3:
    fAxesDef = new KT_Axes();
    break;
  case 4:
    fAxesDef = new OnePass_KT_Axes();
  }

  //-- Trimming parameters --

  fComputeTrimming = GetBool("ComputeTrimming", false);
  fRTrim = GetDouble("RTrim", 0.2);
  fPtFracTrim = GetDouble("PtFracTrim", 0.05);

  //-- Pruning parameters --

  fComputePruning = GetBool("ComputePruning", false);
  fZcutPrun = GetDouble("ZcutPrun", 0.1);
  fRcutPrun = GetDouble("RcutPrun", 0.5);
  fRPrun = GetDouble("RPrun", 0.8);

  //-- SoftDrop parameters --

  fComputeSoftDrop = GetBool("ComputeSoftDrop", false);
  fBetaSoftDrop = GetDouble("BetaSoftDrop", 0.0);
  fSymmetryCutSoftDrop = GetDouble("SymmetryCutSoftDrop", 0.1);
  fR0SoftDrop = GetDouble("R0SoftDrop=", 0.8);

  // ---  Jet Area Parameters ---

  fAreaAlgorithm = GetInt("AreaAlgorithm", 0);
  fComputeRho = GetBool("ComputeRho", false);

  // - ghost based areas -
  fGhostEtaMax = GetDouble("GhostEtaMax", 5.0);
  fRepeat = GetInt("Repeat", 1);
  fGhostArea = GetDouble("GhostArea", 0.01);
  fGridScatter = GetDouble("GridScatter", 1.0);
  fPtScatter = GetDouble("PtScatter", 0.1);
  fMeanGhostPt = GetDouble("MeanGhostPt", 1.0E-100);

  // - voronoi based areas -
  fEffectiveRfact = GetDouble("EffectiveRfact", 1.0);

  switch(fAreaAlgorithm)
  {
  default:
  case 0:
    fAreaDefinition = 0;
    break;
  case 1:
    fAreaDefinition = new AreaDefinition(active_area_explicit_ghosts, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
    break;
  case 2:
    fAreaDefinition = new AreaDefinition(one_ghost_passive_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
    break;
  case 3:
    fAreaDefinition = new AreaDefinition(passive_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
    break;
  case 4:
    fAreaDefinition = new AreaDefinition(VoronoiAreaSpec(fEffectiveRfact));
    break;
  case 5:
    fAreaDefinition = new AreaDefinition(active_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
    break;
  }

  switch(fJetAlgorithm)
  {
  case 1:
    plugin = new CDFJetCluPlugin(fSeedThreshold, fConeRadius, fAdjacencyCut, fMaxIterations, fIratch, fOverlapThreshold);
    fDefinition = new JetDefinition(plugin);
    break;
  case 2:
    plugin = new CDFMidPointPlugin(fSeedThreshold, fConeRadius, fConeAreaFraction, fMaxPairSize, fMaxIterations, fOverlapThreshold);
    fDefinition = new JetDefinition(plugin);
    break;
  case 3:
    plugin = new SISConePlugin(fConeRadius, fOverlapThreshold, fMaxIterations, fJetPTMin);
    fDefinition = new JetDefinition(plugin);
    break;
  case 4:
    fDefinition = new JetDefinition(kt_algorithm, fParameterR);
    break;
  case 5:
    fDefinition = new JetDefinition(cambridge_algorithm, fParameterR);
    break;
  default:
  case 6:
    fDefinition = new JetDefinition(antikt_algorithm, fParameterR);
    break;
  case 7:
    recomb = new WinnerTakeAllRecombiner();
    fDefinition = new JetDefinition(antikt_algorithm, fParameterR, recomb, Best);
    break;
  case 8:
    fNjettinessPlugin = new NjettinessPlugin(fN, Njettiness::wta_kt_axes, Njettiness::unnormalized_cutoff_measure, fBeta, fRcutOff);
    fDefinition = new JetDefinition(fNjettinessPlugin);
    break;
  case 9:
    fValenciaPlugin = new ValenciaPlugin(fParameterR, fBeta, fGamma);
    fDefinition = new JetDefinition(fValenciaPlugin);
    break;
  case 10:
    fDefinition = new JetDefinition(ee_genkt_algorithm,fParameterR,fParameterP);
    break;

  // kT durham algorithm, 2 options:
  // 1. njets mode: stop when reach predetermined n jet (optionally apply sqrt(ExclYmerge(n-1,n))*Evis) > cut offline)
  // 2. dcut mode: stop when all dij above some threshold dcut. Is applied if fDCut > 0.
  case 11:
    fDefinition = new JetDefinition(ee_kt_algorithm);
    break;

  }

  fPlugin = plugin;
  fRecomb = recomb;

  ClusterSequence::print_banner();

  if(fComputeRho && fAreaDefinition)
  {
    // read eta ranges

    param = GetParam("RhoEtaRange");
    size = param.GetSize();

    fEstimators.clear();
    for(i = 0; i < size / 2; ++i)
    {
      etaMin = param[i * 2].GetDouble();
      etaMax = param[i * 2 + 1].GetDouble();
      estimatorStruct.estimator = new JetMedianBackgroundEstimator(SelectorRapRange(etaMin, etaMax), *fDefinition, *fAreaDefinition);
      estimatorStruct.etaMin = etaMin;
      estimatorStruct.etaMax = etaMax;
      fEstimators.push_back(estimatorStruct);
    }
  }

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "Calorimeter/towers"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));
  fRhoOutputArray = ExportArray(GetString("RhoOutputArray", "rho"));
  fConstituentsOutputArray = ExportArray(GetString("ConstituentsOutputArray", "constituents"));
}

//------------------------------------------------------------------------------

void FastJetFinder::Finish()
{
  vector<TEstimatorStruct>::iterator itEstimators;

  for(itEstimators = fEstimators.begin(); itEstimators != fEstimators.end(); ++itEstimators)
  {
    if(itEstimators->estimator) delete itEstimators->estimator;
  }

  if(fItInputArray) delete fItInputArray;
  if(fDefinition) delete fDefinition;
  if(fAreaDefinition) delete fAreaDefinition;
  if(fPlugin) delete static_cast<JetDefinition::Plugin *>(fPlugin);
  if(fRecomb) delete static_cast<JetDefinition::Recombiner *>(fRecomb);
  if(fNjettinessPlugin) delete static_cast<JetDefinition::Plugin *>(fNjettinessPlugin);
  if(fAxesDef) delete fAxesDef;
  if(fMeasureDef) delete fMeasureDef;
  if(fValenciaPlugin) delete static_cast<JetDefinition::Plugin *>(fValenciaPlugin);
}

//------------------------------------------------------------------------------

void FastJetFinder::Process()
{
  Candidate *candidate, *constituent;
  TLorentzVector momentum;

  Double_t deta, dphi, detaMax, dphiMax;
  Double_t time, timeWeight;
  Double_t neutralEnergyFraction, chargedEnergyFraction;

  Int_t number, ncharged, nneutrals;
  Int_t charge;
  Double_t rho = 0.0;
  PseudoJet jet, area;
  ClusterSequence *sequence;
  vector<PseudoJet> inputList, outputList, subjets;
  vector<PseudoJet>::iterator itInputList, itOutputList;
  vector<TEstimatorStruct>::iterator itEstimators;
  Double_t excl_ymerge12 = 0.0;
  Double_t excl_ymerge23 = 0.0;
  Double_t excl_ymerge34 = 0.0;
  Double_t excl_ymerge45 = 0.0;
  Double_t excl_ymerge56 = 0.0;

  DelphesFactory *factory = GetFactory();

  inputList.clear();

  // loop over input objects
  fItInputArray->Reset();
  number = 0;
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    momentum = candidate->Momentum;
    jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
    jet.set_user_index(number);
    inputList.push_back(jet);
    ++number;
  }

  // construct jets
  if(fAreaDefinition)
  {
    sequence = new ClusterSequenceArea(inputList, *fDefinition, *fAreaDefinition);
  }
  else
  {
    sequence = new ClusterSequence(inputList, *fDefinition);
  }

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
      fRhoOutputArray->Add(candidate);
    }
  }

  outputList.clear();

  if(fExclusiveClustering)
  {
    try
    {
      // exclusive dcut mode
      if (fDCut > 0.0)
      {
        outputList = sorted_by_pt(sequence->exclusive_jets(fDCut*fDCut));
      }
      else
      {
        // exclusive njet mode
        outputList = sorted_by_pt(sequence->exclusive_jets(fNJets));
      }
    }
    catch(fastjet::Error &)
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

    neutralEnergyFraction =0.;
    chargedEnergyFraction =0.;

    inputList.clear();
    inputList = sequence->constituents(*itOutputList);

    for(itInputList = inputList.begin(); itInputList != inputList.end(); ++itInputList)
    {
      if(itInputList->user_index() < 0) continue;
      constituent = static_cast<Candidate *>(fInputArray->At(itInputList->user_index()));

      deta = TMath::Abs(momentum.Eta() - constituent->Momentum.Eta());
      dphi = TMath::Abs(momentum.DeltaPhi(constituent->Momentum));
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

      time += TMath::Sqrt(constituent->Momentum.E()) * (constituent->Position.T());
      timeWeight += TMath::Sqrt(constituent->Momentum.E());

      charge += constituent->Charge;

      fConstituentsOutputArray->Add(constituent);
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

    candidate->NeutralEnergyFraction = (momentum.E() > 0 ) ? neutralEnergyFraction/momentum.E() : 0.0;
    candidate->ChargedEnergyFraction = (momentum.E() > 0 ) ? chargedEnergyFraction/momentum.E() : 0.0;

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

    fOutputArray->Add(candidate);
  }
  delete sequence;
}
