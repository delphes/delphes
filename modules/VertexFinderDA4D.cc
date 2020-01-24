/** \class VertexFinderDA4D
 *
 *  Cluster vertices from tracks using deterministic annealing and timing information
 *
 *  \authors O. Cerri
 *
 */


#include "modules/VertexFinderDA4D.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMatrixT.h"
#include "TLatex.h"
#include "TVector3.h"

#include "TAxis.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TFile.h"
#include "TColor.h"
#include "TLegend.h"

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <vector>

using namespace std;

namespace vtx_DAZT
{
  static const Double_t c_light = 2.99792458e+8; // [m/s]
}
using namespace vtx_DAZT;

//------------------------------------------------------------------------------

VertexFinderDA4D::VertexFinderDA4D()
{
  fVerbose = 0;
  fMaxIterations = 0;
  fBetaMax = 0;
  fBetaStop = 0;
  fBetaPurge = 0;
  fVertexZSize = 0;
  fVertexTSize = 0;
  fCoolingFactor = 0;
  fDzCutOff = 0;
  fD0CutOff = 0;
  fDtCutOff = 0;
  fPtMin = 0;
  fPtMax = 0;
  fD2Merge = 0;
  fMuOutlayer = 0;
  fMinTrackProb = 0;
}

//------------------------------------------------------------------------------

VertexFinderDA4D::~VertexFinderDA4D()
{
}

//------------------------------------------------------------------------------

void VertexFinderDA4D::Init()
{
  fVerbose         = GetInt("Verbose", 0);

  fMaxIterations   = GetInt("MaxIterations", 100);
  fMaxVertexNumber = GetInt("MaxVertexNumber", 500);

  fBetaMax         = GetDouble("BetaMax", 1.5);
  fBetaPurge       = GetDouble("BetaPurge", 1.);
  fBetaStop        = GetDouble("BetaStop", 0.2);

  fVertexZSize     = GetDouble("VertexZSize", 0.1); //in mm
  fVertexTSize     = 1E12*GetDouble("VertexTimeSize", 15E-12); //Convert from [s] to [ps]

  fCoolingFactor   = GetDouble("CoolingFactor", 0.8); // Multiply T so to cooldown must be <1

  fDzCutOff        = GetDouble("DzCutOff", 40);      // For the moment 3*DzCutOff is hard cut off for the considered tracks
  fD0CutOff        = GetDouble("D0CutOff", .5);       // d0/sigma_d0, used to compute the pi (weight) of the track
  fDtCutOff        = GetDouble("DtCutOff", 160);     // [ps], 3*DtCutOff is hard cut off for tracks
  fPtMin           = GetDouble("PtMin", 0.5);        // Minimum pt accepted for tracks
  fPtMax           = GetDouble("PtMax", 50);        // Maximum pt accepted for tracks


  fD2UpdateLim     = GetDouble("D2UpdateLim", .5);   // ((dz/ZSize)^2+(dt/TSize)^2)/nv limit for merging vertices
  fD2Merge         = GetDouble("D2Merge", 4.0);      // (dz/ZSize)^2+(dt/TSize)^2 limit for merging vertices
  fMuOutlayer      = GetDouble("MuOutlayer", 4);     // Outlayer rejection exponent
  fMinTrackProb    = GetDouble("MinTrackProb", 0.6); // Minimum probability to be assigned at a vertex
  fMinNTrack       = GetInt("MinNTrack", 10);        // Minimum number of tracks per vertex

  fFigFolderPath   = GetString("DebugFigPath", ".");

  fInputArray = ImportArray(GetString("InputArray", "TrackSmearing/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  fTrackOutputArray = ExportArray(GetString("TrackOutputArray", "tracks"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));

  fInputGenVtx = ImportArray(GetString("InputGenVtx", "PileUpMerger/vertices"));
  fItInputGenVtx = fInputGenVtx->MakeIterator();

  if (fBetaMax < fBetaPurge)
  {
    fBetaPurge = fBetaMax;
    if (fVerbose)
    {
      cout << "BetaPurge set to " << fBetaPurge << endl;
    }
  }

  if (fBetaPurge < fBetaStop)
  {
    fBetaStop = fBetaPurge;
    if (fVerbose)
    {
      cout << "BetaPurge set to " << fBetaPurge << endl;
    }
  }
}

//------------------------------------------------------------------------------

void VertexFinderDA4D::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void VertexFinderDA4D::Process()
{
  fInputArray->Sort();

  if (fVerbose)
  {
     cout<< endl << "      Start processing vertices with VertexFinderDA4D" << endl;
     cout<<" Found "<<fInputArray->GetEntriesFast()<<" input tracks"<<endl;
  }

  // clusterize tracks
  TObjArray *ClusterArray = new TObjArray;
  clusterize(*ClusterArray);

  if(fVerbose>10)
  {
    unsigned int N = fEnergy_rec.size();
    TGraph* gr1 = new TGraph(N, &fBeta_rec[0], &fNvtx_rec[0]);
    gr1->SetName("gr1");
    gr1->GetXaxis()->SetTitle("beta");
    gr1->GetYaxis()->SetTitle("# Vtx");
    TGraph* gr2 = new TGraph(N, &fBeta_rec[0], &fEnergy_rec[0]);
    gr2->SetName("gr2");
    gr2->GetXaxis()->SetTitle("beta");
    gr2->GetYaxis()->SetTitle("Total Energy");
    TGraph* gr3 = new TGraph(N, &fNvtx_rec[0], &fEnergy_rec[0]);
    gr3->SetName("gr3");
    gr3->GetXaxis()->SetTitle("# Vtx");
    gr3->GetYaxis()->SetTitle("Total Energy");

    auto f = new TFile("~/Desktop/debug/EnergyStat.root", "recreate");
    gr1->Write("gr1");
    gr2->Write("gr2");
    gr3->Write("gr3");

    f->Close();
  }

  if (fVerbose){std::cout <<  " clustering returned  "<< ClusterArray->GetEntriesFast() << " clusters  from " << fInputArray->GetEntriesFast() << " input tracks" <<std::endl;}

  // //loop over vertex candidates
  TIterator * ItClusterArray = ClusterArray->MakeIterator();
  ItClusterArray->Reset();
  Candidate *candidate;
  unsigned int k = 0;
  while((candidate = static_cast<Candidate*>(ItClusterArray->Next())))
  {
    if(fVerbose)
    {
     cout << Form("Cluster %d has %d tracks ", k, candidate->GetCandidates()->GetEntriesFast()) << endl;
    }
    if(candidate->ClusterNDF>0)
    {
      // Estimate the vertex resolution
      // loop over tracks belonging to this vertex
      TIter it1(candidate->GetCandidates());
      it1.Reset();

      Candidate *track;
      double sum_Dt_2 = 0;
      double sum_Dz_2 = 0;
      double sum_wt = 0;
      double sum_wz = 0;
      while((track = static_cast<Candidate*>(it1.Next())))
      {
        double dz = candidate->Position.Z() - track->Zd;
        double dt = candidate->Position.T() - track->Td;

        double wz = track->VertexingWeight/(track->ErrorDZ*track->ErrorDZ);
        double wt = track->VertexingWeight/(track->ErrorT*track->ErrorT);

        sum_Dt_2 += wt*dt*dt;
        sum_Dz_2 += wz*dz*dz;
        sum_wt += wt;
        sum_wz += wz;
      }

      double sigma_z = sqrt(sum_Dz_2/sum_wz);
      double sigma_t = sqrt(sum_Dt_2/sum_wt);
      candidate->PositionError.SetXYZT(0.0, 0.0, sigma_z , sigma_t);
      if(fVerbose > 3)
      {
        cout << "k: " << k << endl;
        cout << "Sigma z: " << sigma_z*1E3 << " um" << endl;
        cout << "Sigma t: " << sigma_t*1E9/c_light << " ps" << endl;
      }

      fVertexOutputArray->Add(candidate);
      k++;
    }
   }// end of cluster loop

  delete ClusterArray;
}

//------------------------------------------------------------------------------

void VertexFinderDA4D::clusterize(TObjArray &clusters)
{
  tracks_t tks;
  fill(tks);
  unsigned int nt=tks.getSize();
  if(fVerbose)
  {
    cout << "Tracks added: " << nt << endl;
  }
  if (nt == 0) return;



  vertex_t vtx; // the vertex prototypes
  vtx.ZSize = fVertexZSize;
  vtx.TSize = fVertexTSize;
  // initialize:single vertex at infinite temperature
  vtx.addItem(0, 0, 1);

  // Fit the vertex at T=inf and return the starting temperature
  double beta=beta0(tks, vtx);

  if( fVerbose > 1 )
  {
    cout << "Cluster position at T=inf: z = " << vtx.z[0] << " mm , t = " << vtx.t[0] << " ps" << "  pk = " << vtx.pk[0] << endl;
    cout << Form("Beta Start = %2.1e", beta) << endl;
  }

  if( fVerbose > 2){cout << "Cool down untill reaching the temperature to finish increasing the number of vertexes" << endl;}

  double rho0=0.0;  // start with no outlier rejection

  unsigned int last_round = 0;
  while(last_round < 2)
  {

    unsigned int niter=0;
    double delta2 = 0;
    do  {
      delta2 = update(beta, tks, vtx, rho0);

      if (fVerbose > 3)
      {
        cout << "Update " << niter << " : " << delta2 << endl;
      }
      niter++;
    }
    while (delta2 > fD2UpdateLim &&  niter < fMaxIterations);


    unsigned int n_it = 0;
    while(merge(vtx, fD2Merge) && n_it < fMaxIterations)
    {
      unsigned int niter=0;
      double delta2 = 0;
      do  {
        delta2 = update(beta, tks, vtx, rho0);
        niter++;
      }
      while (delta2 > fD2UpdateLim &&  niter < fMaxIterations);
      n_it++;

    }

    beta /= fCoolingFactor;

    if( beta < fBetaStop )
    {
      split(beta, vtx, tks);
    }
    else
    {
      beta = fBetaStop;
      last_round++;
    }

    if(fVerbose > 3)
    {
      cout << endl << endl << " ----- Beta = " << beta << " --------" << endl;
      cout << "Nv: " << vtx.getSize() << endl;
    }
  }

  if( fVerbose > 4)
  {
    for(unsigned int k = 0; k < vtx.getSize(); k++)
    {
      cout << Form("Vertex %d next beta_c = %.3f", k, vtx.beta_c[k]) << endl;
    }
  }

  if(fVerbose > 2)  {cout << "Adiabatic switch on of outlayr rejection" << endl;}
  rho0 = 1./nt;
  const double N_cycles = 10;
  for(unsigned int f = 1; f <= N_cycles; f++)
  {
    unsigned int niter=0;
    double delta2 = 0;
    do  {
      delta2 = update(beta, tks, vtx, rho0 * f/N_cycles);
      niter++;
    }
    while (delta2 > 0.3*fD2UpdateLim &&  niter < fMaxIterations);
  }

  do {
    beta /= fCoolingFactor;
    if(beta > fBetaPurge) beta = fBetaPurge;
    unsigned int i_pu = 0;
    for(int min_trk = 2; min_trk<=fMinNTrack; min_trk++)
    {
      while( purge(vtx, tks, rho0, beta, fMinTrackProb, min_trk) )
      {
        unsigned int niter=0;
        double delta2 = 0;
        do  {
          delta2 = update(beta, tks, vtx, rho0);
          niter++;
        }
        while (delta2 > fD2UpdateLim &&  niter < fMaxIterations);
        i_pu++;
      }
    }

    unsigned int n_it = 0;
    while(merge(vtx, fD2Merge) && n_it < fMaxIterations)
    {
      unsigned int niter=0;
      double delta2 = 0;
      do  {
        delta2 = update(beta, tks, vtx, rho0);
        niter++;
      }
      while (delta2 > fD2UpdateLim &&  niter < fMaxIterations);
      n_it++;

    }
  } while( beta < fBetaPurge );


  if(fVerbose > 2){cout << "Cooldown untill the limit before assigning track to vertices" << endl;}
  last_round = 0;
  while(last_round < 2)
  {
    unsigned int niter=0;
    double delta2 = 0;
    do  {
      delta2 = update(beta, tks, vtx, rho0);
      niter++;
    }
    while (delta2 > 0.3*fD2UpdateLim &&  niter < fMaxIterations);

    beta /= fCoolingFactor;
    if ( beta >= fBetaMax )
    {
      beta = fBetaMax;
      last_round++;
    }
  }


  // Build the cluster candidates
  for(unsigned int k = 0; k < vtx.getSize(); k++)
  {
    DelphesFactory *factory = GetFactory();
    Candidate * candidate = factory->NewCandidate();

    candidate->ClusterIndex = k;
    candidate->Position.SetXYZT(0.0, 0.0, vtx.z[k] , vtx.t[k]*1E-9*c_light);
    candidate->InitialPosition.SetXYZT(0.0, 0.0, vtx.z[k] , vtx.t[k]*1E-9*c_light);    
    candidate->PositionError.SetXYZT(0.0, 0.0, fVertexZSize , fVertexTSize*1E-9*c_light);
    candidate->SumPT2 = 0;
    candidate->SumPt = 0;
    candidate->ClusterNDF = 0;

    clusters.Add(candidate);
  }


  // Assign each track to the most probable vertex
  double Z_init = rho0 * exp(-beta * fMuOutlayer * fMuOutlayer); // Add fDtCutOff here toghether  with this
  vector<double> pk_exp_mBetaE = Compute_pk_exp_mBetaE(beta, vtx, tks, Z_init);
  for(unsigned int i = 0; i< tks.getSize(); i++)
  {
    if(tks.w[i] <= 0) continue;

    double p_max = 0;
    unsigned int k_max = 0;

    for(unsigned int k = 0; k < vtx.getSize(); k++)
    {
      unsigned int idx = k*nt + i;
      if(pk_exp_mBetaE[idx] == 0 || tks.Z[i] == 0 || vtx.pk[k] == 0)
      {
        continue;
      }

      double pv_max = vtx.pk[k] / (vtx.pk[k] + rho0 * exp(-beta * fMuOutlayer* fMuOutlayer));
      double p = pk_exp_mBetaE[idx] / tks.Z[i];

      p /= pv_max;

      if(p > p_max)
      {
        p_max = p;
        k_max = k;
      }
    }

    if(p_max > fMinTrackProb)
    {
      tks.tt[i]->ClusterIndex = k_max;
      tks.tt[i]->InitialPosition.SetT(1E-9*vtx.t[k_max]*c_light);
      tks.tt[i]->InitialPosition.SetZ(vtx.z[k_max]);

      ((Candidate *) clusters.At(k_max))->AddCandidate(tks.tt[i]);
      ((Candidate *) clusters.At(k_max))->SumPT2 += tks.tt[i]->Momentum.Pt()*tks.tt[i]->Momentum.Pt();
      ((Candidate *) clusters.At(k_max))->SumPt += tks.tt[i]->Momentum.Pt();
      ((Candidate *) clusters.At(k_max))->ClusterNDF += 1;
    }
    else
    {
      tks.tt[i]->ClusterIndex = -1;
      tks.tt[i]->InitialPosition.SetT(1E3*1000000*c_light);
      tks.tt[i]->InitialPosition.SetZ(1E8);
    }
    fTrackOutputArray->Add(tks.tt[i]);
  }

}

//------------------------------------------------------------------------------
// Definition of the distance metrci between track and vertex
double VertexFinderDA4D::Energy(double t_z, double v_z, double dz2_o, double t_t, double v_t, double dt2_o)
{
  return (t_z - v_z)*(t_z - v_z)* dz2_o + (t_t - v_t)*(t_t - v_t)*dt2_o;
}

//------------------------------------------------------------------------------
// Fill tks with the input candidates array
void VertexFinderDA4D::fill(tracks_t &tks)
{
  tks.sum_w_o_dt2 = 0;
  tks.sum_w_o_dz2 = 0;
  tks.sum_w = 0;

  Candidate *candidate;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    unsigned int discard = 0;

    double pt = candidate->Momentum.Pt();
    if(pt<fPtMin || pt>fPtMax) discard = 1;

    // ------------- Compute cloasest approach Z ----------------
    double z = candidate->DZ; // [mm]

    candidate->Zd = candidate->DZ; //Set the cloasest approach z
    if(fabs(z) > 3*fDzCutOff) discard = 1;

    // ------------- Compute cloasest approach T ----------------
    //Asumme pion mass which is the most common particle
    double M = 0.139570;
    candidate->Mass = M;
    double p = pt * sqrt(1 + candidate->CtgTheta*candidate->CtgTheta);
    double e = sqrt(p*p + M*M);

    double t = candidate->Position.T()*1.E9/c_light; // from [mm] to [ps]
    if(t <= -9999) discard = 1;                    // Means that the time information has not been added

    // DEBUG Here backpropagete for the whole length and not noly for z. Could improve resolution
    // double bz = pt * candidate->CtgTheta/e;
    // t += (z - candidate->Position.Z())*1E9/(c_light*bz);

    // Use full path Length
    t -= candidate->L*1E9/(c_light*p/e);

    candidate->Td = t*1E-9*c_light;
    if(fabs(t) > 3*fDtCutOff) discard = 1;

    // auto genp = (Candidate*) candidate->GetCandidates()->At(0);
    // cout << "Eta: " << candidate->Position.Eta() << endl;
    // cout << genp->Momentum.Pt() << " -- " << candidate->Momentum.Pt() << endl;
    // cout << genp->Momentum.Pz() << " -- " << candidate->Momentum.Pz() << endl;
    // cout << genp->Momentum.P() << " -- " << p << endl;
    // cout << genp->Momentum.E() << " -- " << e << endl;
    // cout << Form("bz_true: %.4f -- bz_gen: %.4f", genp->Momentum.Pz()/genp->Momentum.E(), bz) << endl;

    double dz2_o = candidate->ErrorDZ*candidate->ErrorDZ;
    dz2_o += fVertexZSize*fVertexZSize;
    // when needed add beam spot width (x-y)?? mha?
    dz2_o = 1/dz2_o; //Multipling is faster than dividing all the times

    double dt2_o = candidate->ErrorT*1.E9/c_light; // [ps]
    dt2_o *= dt2_o;
    dt2_o += fVertexTSize*fVertexTSize; // [ps^2]
    // Ideally we should also add the induced uncertantiy from dz, z_out, pt, ctgthetaand all the other thing used above (total around 5ps). For the moment we compensae using a high value for vertex time.
    dt2_o = 1/dt2_o;

    double w;
    if(fD0CutOff > 0 && candidate->ErrorD0 > 0)
    {
      double d0_sig = fabs(candidate->D0/candidate->ErrorD0);
      w = exp(d0_sig*d0_sig - fD0CutOff*fD0CutOff);
      w = 1./(1. + w);
      if (w < 1E-4) discard = 1;
    }
    else
    {
      w = 1;
    }
    candidate->VertexingWeight = w;


    if(discard)
    {
      candidate->ClusterIndex = -1;
      candidate->InitialPosition.SetT(1E3*1000000*c_light);
      candidate->InitialPosition.SetZ(1E8);
      fTrackOutputArray->Add(candidate);
    }
    else
    {
      tks.sum_w_o_dt2 += w * dt2_o;
      tks.sum_w_o_dz2 += w * dz2_o;
      tks.sum_w += w;
      tks.addItem(z, t, dz2_o, dt2_o, &(*candidate), w, candidate->PID); //PROVA: rimuovi &(*---)
    }

  }

  if(fVerbose > 1)
  {
    cout << "----->Filled tracks" << endl;
    cout << "M        z           dz        t            dt        w" << endl;
    for(unsigned int i = 0; i < tks.getSize(); i++)
    {
      cout << Form("%d\t%1.1e\t%1.1e\t%1.1e\t%1.1e\t%1.1e", tks.PID[i], tks.z[i], 1/sqrt(tks.dz2_o[i]), tks.t[i], 1/sqrt(tks.dt2_o[i]), tks.w[i]) << endl;
    }
  }

  return;
}

//------------------------------------------------------------------------------
// Compute higher phase transition temperature
double VertexFinderDA4D::beta0(tracks_t & tks, vertex_t &vtx)
{
  if(vtx.getSize() != 1)
  {
    throw std::invalid_argument( "Unexpected number of vertices" );
  }

  unsigned int nt = tks.getSize();

  //Set vertex position at T=inf as the weighted average of the tracks
  double sum_wz = 0, sum_wt = 0;
  for(unsigned int i = 0; i < nt; i++)
  {
    sum_wz += tks.w[i] * tks.z[i] * tks.dz2_o[i];
    sum_wt += tks.w[i] * tks.t[i] * tks.dt2_o[i];
  }
  vtx.t[0] = sum_wt / tks.sum_w_o_dt2;
  vtx.z[0] = sum_wz / tks.sum_w_o_dz2;

  // Compute the posterior distribution covariance matrix elements
  double s_zz = 0, s_tt = 0, s_tz = 0;
  for(unsigned int i = 0; i < nt; i++)
  {
    double dz = (tks.z[i] - vtx.z[0]) * tks.dz_o[i];
    double dt = (tks.t[i] - vtx.t[0]) * tks.dt_o[i];

    s_zz += tks.w[i] * dz * dz;
    s_tt += tks.w[i] * dt * dt;
    s_tz += tks.w[i] * dt * dz;
  }
  s_tt /= tks.sum_w;
  s_zz /= tks.sum_w;
  s_tz /= tks.sum_w;

  // Copute the max eighenvalue
  double beta_c = (s_tt - s_zz)*(s_tt - s_zz) + 4*s_tz*s_tz;
  beta_c = 1. / (s_tt + s_zz + sqrt(beta_c));

  double out;
  if (beta_c < fBetaMax)
  {
    // Cool down up to a step before the phase transition
    out = beta_c * sqrt(fCoolingFactor);
  }
  else
  {
    out = fBetaMax * fCoolingFactor;
  }

  return out;
}

//------------------------------------------------------------------------------
// Compute the new vertexes position and mass (probability) -- mass constrained annealing without noise
// Compute and store the posterior covariance matrix elements
// Returns the squared sum of changes of vertexex position normalized by the vertex size declared in the init
double VertexFinderDA4D::update(double beta, tracks_t &tks, vertex_t &vtx, double rho0)
{
  unsigned int nt = tks.getSize();
  unsigned int nv = vtx.getSize();

  //initialize sums
  double Z_init = rho0 * exp(-beta * fMuOutlayer * fMuOutlayer);

  // Compute all the energies (aka distances) and normalization partition function
  vector<double> pk_exp_mBetaE = Compute_pk_exp_mBetaE(beta, vtx, tks, Z_init);

  double sum_pk = 0;
  double delta2_max = 0;
  for (unsigned int k = 0; k < nv; k++)
  {
    // Compute the new vertex positions and masses
    double pk_new = 0;
    double sw_z = 0, sw_t = 0;
    // Compute the posterior covariance matrix Elements
    double szz = 0, stt = 0, stz = 0;
    double sum_wt = 0, sum_wz = 0;
    double sum_ptt = 0, sum_pzz = 0, sum_ptz = 0;


    for (unsigned int i = 0; i < nt; i++)
    {
      unsigned int idx = k*nt + i;

      if(pk_exp_mBetaE[idx] == 0 || tks.Z[i] == 0)
      {
        continue;
      }

      double p_ygx = pk_exp_mBetaE[idx] / tks.Z[i];      //p(y|x), Gibbs distribution
      if(std::isnan(p_ygx) || std::isinf(p_ygx) || p_ygx > 1)
      {
        cout << Form("%1.6e    %1.6e", pk_exp_mBetaE[idx], tks.Z[i]);
        throw std::invalid_argument(Form("p_ygx is %.8f", p_ygx));
      }
      pk_new += tks.w[i] * p_ygx;

      double wt = tks.w[i] * p_ygx * tks.dt2_o[i];
      sw_t += wt * tks.t[i];
      sum_wt += wt;

      double wz = tks.w[i] * p_ygx * tks.dz2_o[i];
      sw_z += wz * tks.z[i];
      sum_wz += wz;

      // Add the track contribution to the covariance matrix
      double p_xgy = p_ygx * tks.w[i] / vtx.pk[k];
      double dt = (tks.t[i] - vtx.t[k]) * tks.dt_o[i];
      double dz = (tks.z[i] - vtx.z[k]) * tks.dz_o[i];

      double wtt = p_xgy * tks.dt2_o[i];
      double wzz = p_xgy * tks.dz2_o[i];
      double wtz = p_xgy * tks.dt_o[i] * tks.dz_o[i];

      stt += wtt * dt * dt;
      szz += wzz * dz * dz;
      stz += wtz * dt * dz;

      sum_ptt += wtt;
      sum_pzz += wzz;
      sum_ptz += wtz;
    }
    if(pk_new == 0)
    {
      vtx.removeItem(k);
      k--;
      // throw std::invalid_argument(Form("pk_new is %.8f", pk_new));
    }
    else
    {
      pk_new /= tks.sum_w;
      sum_pk += pk_new;

      stt /= sum_ptt;
      szz /= sum_pzz;
      stz /= sum_ptz;

      double new_t = sw_t/sum_wt;
      double new_z = sw_z/sum_wz;
      if(std::isnan(new_z) || std::isnan(new_t))
      {
        cout << endl << endl;
        cout << Form("t: %.3e   /   %.3e", sw_t, sum_wt) << endl;
        cout << Form("z: %.3e   /   %.3e", sw_z, sum_wz) << endl;
        cout << "pk " << k << "  " << vtx.pk[k] << endl;
        throw std::invalid_argument("new_z is nan");
      }

      double z_displ = (new_z - vtx.z[k])/fVertexZSize;
      double t_displ = (new_t - vtx.t[k])/fVertexTSize;
      double delta2 = z_displ*z_displ + t_displ*t_displ;

      if (delta2 > delta2_max) delta2_max =  delta2;

      vtx.z[k] = new_z;
      vtx.t[k] = new_t;
      vtx.pk[k] = pk_new;
      vtx.szz[k] = szz;
      vtx.stt[k] = stt;
      vtx.stz[k] = stz;
    }
  }

  if(fabs((sum_pk - 1.) > 1E-4))
  {
    cout << "sum_pk " << sum_pk << endl;
    for (unsigned int k = 0; k < nv; k++)
    {
      cout << Form("%d: %1.4e", k, vtx.pk[k]) << endl;
    }
    throw std::invalid_argument("Sum of masses not unitary");
  }
  // if(fVerbose > 3)
  // {
  //   cout << "===Update over" << endl;
  //   for (unsigned int k = 0; k < nv; k++)
  //   {
  //     cout << k << endl;
  //     cout << "z: " << vtx.z[k] << " , t: " << vtx.t[k] << " , p: " << vtx.pk[k] << endl;
  //     cout << " | " << vtx.szz[k] << "   " << vtx.stz[k] << "|" << endl;
  //     cout << " | " << vtx.stz[k] << "   " << vtx.stt[k] << "|" << endl << endl;
  //   }
  //   cout << "=======" << endl;
  // }

  return delta2_max;
}

//------------------------------------------------------------------------------
// Split critical vertices (beta_c < beta)
// Returns true if at least one cluster was split
bool VertexFinderDA4D::split(double &beta, vertex_t &vtx, tracks_t & tks)
{
  bool split = false;

  auto pair_bc_k = vtx.ComputeAllBeta_c(fVerbose);

  // If minimum beta_c is higher than beta, no split is necessaire
  if( pair_bc_k.first > beta )
  {
    split = false;
  }
  else
  {
    const unsigned int nv = vtx.getSize();
    for(unsigned int k = 0; k < nv; k++)
    {
      if( fVerbose > 3 )
      {
        cout << "vtx " << k << "  beta_c = " << vtx.beta_c[k] << endl;
      }
      if(vtx.beta_c[k] <= beta)
      {
        double z_old = vtx.z[k];
        double t_old = vtx.t[k];
        double pk_old = vtx.pk[k];

        // Compute splitting direction: given by the max eighenvalue eighenvector
        double zn = (vtx.szz[k] - vtx.stt[k])*(vtx.szz[k] - vtx.stt[k]) + 4*vtx.stz[k]*vtx.stz[k];
        zn = vtx.szz[k] - vtx.stt[k] + sqrt(zn);
        double tn = 2*vtx.stz[k];
        double norm = hypot(zn, tn);
        tn /= norm;
        zn /= norm;

        // Estimate subcluster positions and weight
        double p1=0, z1=0, t1=0, wz1=0, wt1=0;
        double p2=0, z2=0, t2=0, wz2=0, wt2=0;
        const unsigned int nt = tks.getSize();
        for(unsigned int i=0; i<nt; ++i)
        {
          if (tks.Z[i] > 0)
          {
            double lr = (tks.t[i] - vtx.t[k]) * tn + (tks.z[i]-vtx.z[k]) * zn;
            // winner-takes-all, usually overestimates splitting
            double tl = lr < 0 ? 1.: 0.;
            double tr = 1. - tl;

            // soften it, especially at low T
            // double arg = lr * sqrt(beta * ( zn*zn*tks.dz2_o[i] + tn*tn*tks.dt2_o[i] ) );
            // if(abs(arg) < 20)
            // {
            //   double t = exp(-arg);
            //   tl = t/(t+1.);
            //   tr = 1/(t+1.);
            // }

            double p = vtx.pk[k] * tks.w[i];
            p *= exp(-beta * Energy(tks.z[i], vtx.z[k], tks.dz2_o[i], tks.t[i], vtx.t[k], tks.dt2_o[i])) / tks.Z[i];
            double wt = p*tks.dt2_o[i];
            double wz = p*tks.dz2_o[i];
            p1 += p*tl;  z1 += wz*tl*tks.z[i]; t1 += wt*tl*tks.t[i]; wz1 += wz*tl; wt1 += wt*tl;
            p2 += p*tr;  z2 += wz*tr*tks.z[i]; t2 += wt*tr*tks.t[i]; wz2 += wz*tr; wt2 += wt*tr;
          }
        }

        if(wz1 > 0  && wt1 > 0 && wz2 > 0 && wt2 > 0)
        {
          t1 /= wt1;
          z1 /= wz1;
          t2 /= wt2;
          z2 /= wz2;

          if( fVerbose > 3 )
          {
            double aux = (z1-z2)*(z1-z2)/(fVertexZSize*fVertexZSize) + (t1-t2)*(t1-t2)/(fVertexTSize*fVertexTSize);
            cout << "weighted split:  delta = " << sqrt(aux) << endl;
          }
        }
        else
        {
          continue;
          // throw std::invalid_argument( "0 division" );
        }

        while(vtx.NearestCluster(t1, z1) != k || vtx.NearestCluster(t2, z2) != k)
        {
          t1 = 0.5 * (t1 + t_old);
          z1 = 0.5 * (z1 + z_old);
          t2 = 0.5 * (t2 + t_old);
          z2 = 0.5 * (z2 + z_old);
        }

        // Compute final distance and split if the distance is enough
        double delta2 = (z1-z2)*(z1-z2)/(fVertexZSize*fVertexZSize) + (t1-t2)*(t1-t2)/(fVertexTSize*fVertexTSize);
        if(delta2 > fD2Merge)
        {
          split = true;
          vtx.t[k] = t1;
          vtx.z[k] = z1;
          vtx.pk[k] = p1 * pk_old/(p1+p2);

          double new_t = t2;
          double new_z = z2;
          double new_pk = p2 * pk_old/(p1+p2);

          vtx.addItem(new_z, new_t, new_pk);

          if( fVerbose > 3 )
          {
            cout << "===Split happened on vtx " << k << endl;
            cout << "OLD     z: " << z_old << " , t: " << t_old << " , pk: " << pk_old << endl;
            cout << "NEW+    z: " << vtx.z[k] << " , t: " << vtx.t[k] << " , pk: " << vtx.pk[k] << endl;
            cout << "NEW-    z: " << new_z << " , t: " << new_t << " , pk: " << new_pk <<  endl;
          }
        }
      }
    }
  }
  return split;
}


//------------------------------------------------------------------------------
// Merge vertexes closer than declared dimensions
bool VertexFinderDA4D::merge(vertex_t & vtx, double d2_merge = 2)
{
  bool merged = false;

  if(vtx.getSize() < 2) return merged;

  bool last_merge = false;
  do {
    double min_d2 = d2_merge;
    unsigned int k1_min, k2_min;
    for(unsigned int k1 = 0; k1 < vtx.getSize(); k1++)
    {
      for(unsigned int k2 = k1+1; k2 < vtx.getSize();k2++)
      {
        double d2_tmp = vtx.DistanceSquare(k1, k2);
        if(d2_tmp < min_d2)
        {
          min_d2 = d2_tmp;
          k1_min = k1;
          k2_min = k2;
        }
      }
    }

    if(min_d2 < d2_merge)
    {
      vtx.mergeItems(k1_min, k2_min);
      last_merge = true;
      merged = true;
    }
    else last_merge = false;
  } while(last_merge);

  return merged;
}

// -----------------------------------------------------------------------------
// Compute all the energies and set the partition function normalization for each track
vector<double> VertexFinderDA4D::Compute_pk_exp_mBetaE(double beta, vertex_t &vtx, tracks_t &tks, double Z_init)
{
  unsigned int nt = tks.getSize();
  unsigned int nv = vtx.getSize();

  vector<double> pk_exp_mBetaE(nt * nv);
  for (unsigned int k = 0; k < nv; k++)
  {
    for (unsigned int i = 0; i < nt; i++)
    {
      if(k == 0) tks.Z[i] = Z_init;

      double aux = Energy(tks.z[i], vtx.z[k], tks.dz2_o[i], tks.t[i], vtx.t[k], tks.dt2_o[i]);
      aux = vtx.pk[k] * exp(-beta * aux);
      // if(aux < 1E-10) continue;
      tks.Z[i] += aux;

      unsigned int idx = k*nt + i;
      pk_exp_mBetaE[idx] = aux;
    }
  }
  return pk_exp_mBetaE;
}

//------------------------------------------------------------------------------
// Eliminate clusters with only one significant/unique track
bool VertexFinderDA4D::purge(vertex_t & vtx, tracks_t & tks, double & rho0, const double beta, double min_prob, double min_trk)
{
  const unsigned int nv = vtx.getSize();
  const unsigned int nt = tks.getSize();

  if (nv < 2)
    return false;

  double sumpmin = nt;
  unsigned int k0 = nv;

  int nUnique = 0;
  double sump = 0;

  double Z_init = rho0 * exp(-beta * fMuOutlayer * fMuOutlayer); // Add fDtCutOff here toghether  with this
  vector<double> pk_exp_mBetaE = Compute_pk_exp_mBetaE(beta, vtx, tks, Z_init);

  for (unsigned int k = 0; k < nv; ++k) {

    nUnique = 0;
    sump = 0;

    double pmax = vtx.pk[k] / (vtx.pk[k] + rho0 * exp(-beta * fMuOutlayer* fMuOutlayer));
    double pcut = min_prob * pmax;

    for (unsigned int i = 0; i < nt; ++i) {
      unsigned int idx = k*nt + i;

      if(pk_exp_mBetaE[idx] == 0 || tks.Z[i] == 0)
      {
        continue;
      }

      double p = pk_exp_mBetaE[idx] / tks.Z[i];
      sump += p;
      if( ( p > pcut ) & ( tks.w[i] > 0 ) ) nUnique++;
    }

    if ((nUnique < min_trk) && (sump < sumpmin)) {
      sumpmin = sump;
      k0 = k;
    }

  }

  if (k0 != nv) {
    if (fVerbose > 5) {
      std::cout  << Form("eliminating prototype at z = %.3f mm, t = %.0f ps", vtx.z[k0], vtx.t[k0]) << " with sump=" << sumpmin
		 << "  rho*nt =" << vtx.pk[k0]*nt
		 << endl;
    }
    vtx.removeItem(k0);
    return true;
  } else {
    return false;
  }
}