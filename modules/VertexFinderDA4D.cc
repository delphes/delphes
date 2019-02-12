/** \class VertexFinderDA4D
 *
 *  Cluster vertices from tracks using deterministic annealing and timing information
 *
 *  \authors M. Selvaggi, L. Gray
 *
 */

#include "modules/VertexFinderDA4D.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMatrixT.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVector3.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace std;

static const Double_t mm = 1.;
static const Double_t m = 1000. * mm;
static const Double_t ns = 1.;
static const Double_t s = 1.e+9 * ns;
static const Double_t c_light = 2.99792458e+8 * m / s;

struct track_t
{
  double z; // z-coordinate at point of closest approach to the beamline
  double t; // t-coordinate at point of closest approach to the beamline
  double dz2; // square of the error of z(pca)
  double dtz; // covariance of z-t
  double dt2; // square of the error of t(pca)
  Candidate *tt; // a pointer to the Candidate Track
  double Z; // Z[i]   for DA clustering
  double pi; // track weight
  double pt;
  double eta;
  double phi;
};

struct vertex_t
{
  double z;
  double t;
  double pk; // vertex weight for "constrained" clustering
  // --- temporary numbers, used during update
  double ei;
  double sw;
  double swz;
  double swt;
  double se;
  // ---for Tc
  double swE;
  double Tc;
};

static bool split(double beta, std::vector<track_t> &tks, std::vector<vertex_t> &y);
static double update1(double beta, std::vector<track_t> &tks, std::vector<vertex_t> &y);
static double update2(double beta, std::vector<track_t> &tks, std::vector<vertex_t> &y, double &rho0, const double dzCutOff);
static void dump(const double beta, const std::vector<vertex_t> &y, const std::vector<track_t> &tks);
static bool merge(std::vector<vertex_t> &);
static bool merge(std::vector<vertex_t> &, double &);
static bool purge(std::vector<vertex_t> &, std::vector<track_t> &, double &, const double, const double);
static void splitAll(std::vector<vertex_t> &y);
static double beta0(const double betamax, std::vector<track_t> &tks, std::vector<vertex_t> &y, const double coolingFactor);
static double Eik(const track_t &t, const vertex_t &k);

static bool recTrackLessZ1(const track_t &tk1, const track_t &tk2)
{
  return tk1.z < tk2.z;
}

using namespace std;

//------------------------------------------------------------------------------

VertexFinderDA4D::VertexFinderDA4D() :
  fVerbose(0), fMinPT(0), fVertexSpaceSize(0), fVertexTimeSize(0),
  fUseTc(0), fBetaMax(0), fBetaStop(0), fCoolingFactor(0),
  fMaxIterations(0), fDzCutOff(0), fD0CutOff(0), fDtCutOff(0)
{
}

//------------------------------------------------------------------------------

VertexFinderDA4D::~VertexFinderDA4D()
{
}

//------------------------------------------------------------------------------

void VertexFinderDA4D::Init()
{

  fVerbose = GetBool("Verbose", 1);
  fMinPT = GetDouble("MinPT", 0.1);
  fVertexSpaceSize = GetDouble("VertexSpaceSize", 0.5); //in mm
  fVertexTimeSize = GetDouble("VertexTimeSize", 10E-12); //in s
  fUseTc = GetBool("UseTc", 1);
  fBetaMax = GetDouble("BetaMax ", 0.1);
  fBetaStop = GetDouble("BetaStop", 1.0);
  fCoolingFactor = GetDouble("CoolingFactor", 0.8);
  fMaxIterations = GetInt("MaxIterations", 100);
  fDzCutOff = GetDouble("DzCutOff", 40); // Adaptive Fitter uses 30 mm but that appears to be a bit tight here sometimes
  fD0CutOff = GetDouble("D0CutOff", 30);
  fDtCutOff = GetDouble("DtCutOff", 100E-12); // dummy

  // convert stuff in cm, ns
  fVertexSpaceSize /= 10.0;
  fVertexTimeSize *= 1E9;
  fDzCutOff /= 10.0; // Adaptive Fitter uses 3.0 but that appears to be a bit tight here sometimes
  fD0CutOff /= 10.0;

  fInputArray = ImportArray(GetString("InputArray", "TrackSmearing/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void VertexFinderDA4D::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void VertexFinderDA4D::Process()
{
  Candidate *candidate, *track;
  TObjArray *ClusterArray;
  ClusterArray = new TObjArray;
  TIterator *ItClusterArray;
  Int_t ivtx = 0;

  fInputArray->Sort();

  TLorentzVector pos, mom;
  if(fVerbose)
  {
    cout << " start processing vertices ..." << endl;
    cout << " Found " << fInputArray->GetEntriesFast() << " input tracks" << endl;
    //loop over input tracks
    fItInputArray->Reset();
    while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
    {
      pos = candidate->InitialPosition;
      mom = candidate->Momentum;

      cout << "pt: " << mom.Pt() << ", eta: " << mom.Eta() << ", phi: " << mom.Phi() << ", z: " << candidate->DZ / 10 << endl;
    }
  }

  // clusterize tracks in Z
  clusterize(*fInputArray, *ClusterArray);

  if(fVerbose)
  {
    std::cout << " clustering returned  " << ClusterArray->GetEntriesFast() << " clusters  from " << fInputArray->GetEntriesFast() << " selected tracks" << std::endl;
  }

  //loop over vertex candidates
  ItClusterArray = ClusterArray->MakeIterator();
  ItClusterArray->Reset();
  while((candidate = static_cast<Candidate *>(ItClusterArray->Next())))
  {

    double meantime = 0.;
    double expv_x2 = 0.;
    double normw = 0.;
    double errtime = 0;

    double meanpos = 0.;
    double meanerr2 = 0.;
    double normpos = 0.;
    double errpos = 0.;

    double sumpt2 = 0.;

    int itr = 0;

    if(fVerbose) cout << "this vertex has: " << candidate->GetCandidates()->GetEntriesFast() << " tracks" << endl;

    // loop over tracks belonging to this vertex
    TIter it1(candidate->GetCandidates());
    it1.Reset();

    while((track = static_cast<Candidate *>(it1.Next())))
    {

      itr++;
      // TBC: the time is in ns for now TBC
      double t = track->InitialPosition.T() / c_light;
      double dt = track->ErrorT / c_light;
      const double time = t;
      const double inverr = 1.0 / dt;
      meantime += time * inverr;
      expv_x2 += time * time * inverr;
      normw += inverr;

      // compute error position TBC
      const double pt = track->Momentum.Pt();
      const double z = track->DZ / 10.0;
      const double err_pt = track->ErrorPT;
      const double err_z = track->ErrorDZ;

      const double wi = (pt / (err_pt * err_z)) * (pt / (err_pt * err_z));
      meanpos += z * wi;

      meanerr2 += err_z * err_z * wi;
      normpos += wi;
      sumpt2 += pt * pt;

      // while we are here store cluster index in tracks
      track->ClusterIndex = ivtx;
    }

    meantime = meantime / normw;
    expv_x2 = expv_x2 / normw;
    errtime = TMath::Sqrt((expv_x2 - meantime * meantime) / itr);
    meanpos = meanpos / normpos;
    meanerr2 = meanerr2 / normpos;
    errpos = TMath::Sqrt(meanerr2 / itr);

    candidate->Position.SetXYZT(0.0, 0.0, meanpos * 10.0, meantime * c_light);
    candidate->PositionError.SetXYZT(0.0, 0.0, errpos * 10.0, errtime * c_light);
    candidate->SumPT2 = sumpt2;
    candidate->ClusterNDF = itr;
    candidate->ClusterIndex = ivtx;

    fVertexOutputArray->Add(candidate);

    ivtx++;

    if(fVerbose)
    {
      std::cout << "x,y,z";
      std::cout << ",t";
      std::cout << "=" << candidate->Position.X() / 10.0 << " " << candidate->Position.Y() / 10.0 << " " << candidate->Position.Z() / 10.0;
      std::cout << " " << candidate->Position.T() / c_light;

      std::cout << std::endl;
      std::cout << "sumpt2 " << candidate->SumPT2 << endl;

      std::cout << "ex,ey,ez";
      std::cout << ",et";
      std::cout << "=" << candidate->PositionError.X() / 10.0 << " " << candidate->PositionError.Y() / 10.0 << " " << candidate->PositionError.Z() / 10.0;
      std::cout << " " << candidate->PositionError.T() / c_light;
      std::cout << std::endl;
    }
  } // end of cluster loop

  if(fVerbose)
  {
    std::cout << "PrimaryVertexProducerAlgorithm::vertices candidates =" << ClusterArray->GetEntriesFast() << std::endl;
  }

  //TBC maybe this can be done later
  // sort vertices by pt**2  vertex (aka signal vertex tagging)
  /*if(pvs.size()>1){
      sort(pvs.begin(), pvs.end(), VertexHigherPtSquared());
    }
     */

  delete ClusterArray;
}

//------------------------------------------------------------------------------

void VertexFinderDA4D::clusterize(const TObjArray &tracks, TObjArray &clusters)
{
  if(fVerbose)
  {
    cout << "###################################################" << endl;
    cout << "# VertexFinderDA4D::clusterize   nt=" << tracks.GetEntriesFast() << endl;
    cout << "###################################################" << endl;
  }

  vector<Candidate *> pv = vertices();

  if(fVerbose)
  {
    cout << "# VertexFinderDA4D::clusterize   pv.size=" << pv.size() << endl;
  }
  if(pv.size() == 0)
  {
    return;
  }

  // convert into vector of candidates
  //TObjArray *ClusterArray = pv.begin()->GetCandidates();
  //Candidate *aCluster = static_cast<Candidate*>(&(pv.at(0)));
  Candidate *aCluster = pv.at(0);

  // fill into clusters and merge

  if(fVerbose)
  {
    std::cout << '\t' << 0;
    std::cout << ' ' << (*pv.begin())->Position.Z() / 10.0 << ' ' << (*pv.begin())->Position.T() / c_light << std::endl;
  }

  for(vector<Candidate *>::iterator k = pv.begin() + 1; k != pv.end(); k++)
  {
    if(fVerbose)
    {
      std::cout << '\t' << std::distance(pv.begin(), k);
      std::cout << ' ' << (*k)->Position.Z() << ' ' << (*k)->Position.T() << std::endl;
    }

    // TBC - check units here
    if(std::abs((*k)->Position.Z() - (*(k - 1))->Position.Z()) / 10.0 > (2 * fVertexSpaceSize) || std::abs((*k)->Position.T() - (*(k - 1))->Position.Z()) / c_light > 2 * 0.010)
    {
      // close a cluster
      clusters.Add(aCluster);
      //aCluster.clear();
    }
    //for(unsigned int i=0; i<k->GetCandidates().GetEntriesFast(); i++){
    aCluster = *k;
    //}
  }
  clusters.Add(aCluster);

  if(fVerbose)
  {
    std::cout << "# VertexFinderDA4D::clusterize clusters.size=" << clusters.GetEntriesFast() << std::endl;
  }
}

//------------------------------------------------------------------------------

vector<Candidate *> VertexFinderDA4D::vertices()
{
  Candidate *candidate;
  UInt_t clusterIndex = 0;
  vector<Candidate *> clusters;

  vector<track_t> tks;
  track_t tr;
  Double_t z, dz, t, l, dt, d0, d0error;

  // loop over input tracks
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    //TBC everything in cm
    z = candidate->DZ / 10;
    tr.z = z;
    dz = candidate->ErrorDZ / 10;
    tr.dz2 = dz * dz // track error
      //TBC: beamspot size induced error, take 0 for now.
      // + (std::pow(beamspot.BeamWidthX()*cos(phi),2.)+std::pow(beamspot.BeamWidthY()*sin(phi),2.))/std::pow(tantheta,2.) // beam-width induced
      + fVertexSpaceSize * fVertexSpaceSize; // intrinsic vertex size, safer for outliers and short lived decays

    // TBC: the time is in ns for now TBC
    //t = candidate->Position.T()/c_light;
    t = candidate->InitialPosition.T() / c_light;
    l = candidate->L / c_light;
    double pt = candidate->Momentum.Pt();
    double eta = candidate->Momentum.Eta();
    double phi = candidate->Momentum.Phi();

    tr.pt = pt;
    tr.eta = eta;
    tr.phi = phi;
    tr.t = t; //
    tr.dtz = 0.;
    dt = candidate->ErrorT / c_light;
    tr.dt2 = dt * dt + fVertexTimeSize * fVertexTimeSize; // the ~injected~ timing error plus a small minimum vertex size in time
    if(fD0CutOff > 0)
    {

      d0 = TMath::Abs(candidate->D0) / 10.0;
      d0error = candidate->ErrorD0 / 10.0;

      tr.pi = 1. / (1. + exp((d0 * d0) / (d0error * d0error) - fD0CutOff * fD0CutOff)); // reduce weight for high ip tracks
    }
    else
    {
      tr.pi = 1.;
    }
    tr.tt = &(*candidate);
    tr.Z = 1.;

    // TBC now putting track selection here (> fPTMin)
    if(tr.pi > 1e-3 && tr.pt > fMinPT)
    {
      tks.push_back(tr);
    }
  }

  //print out input tracks

  if(fVerbose)
  {
    std::cout << " start processing vertices ..." << std::endl;
    std::cout << " Found " << tks.size() << " input tracks" << std::endl;
    //loop over input tracks

    for(std::vector<track_t>::const_iterator it = tks.begin(); it != tks.end(); it++)
    {
      double z = it->z;
      double pt = it->pt;
      double eta = it->eta;
      double phi = it->phi;
      double t = it->t;

      std::cout << "pt: " << pt << ", eta: " << eta << ", phi: " << phi << ", z: " << z << ", t: " << t << std::endl;
    }
  }

  unsigned int nt = tks.size();
  double rho0 = 0.0; // start with no outlier rejection

  if(tks.empty()) return clusters;

  vector<vertex_t> y; // the vertex prototypes

  // initialize:single vertex at infinite temperature
  vertex_t vstart;
  vstart.z = 0.;
  vstart.t = 0.;
  vstart.pk = 1.;
  y.push_back(vstart);
  int niter = 0; // number of iterations

  // estimate first critical temperature
  double beta = beta0(fBetaMax, tks, y, fCoolingFactor);
  niter = 0;
  while((update1(beta, tks, y) > 1.e-6) && (niter++ < fMaxIterations))
  {
  }

  // annealing loop, stop when T<Tmin  (i.e. beta>1/Tmin)
  while(beta < fBetaMax)
  {

    if(fUseTc)
    {
      update1(beta, tks, y);
      while(merge(y, beta))
      {
        update1(beta, tks, y);
      }
      split(beta, tks, y);
      beta = beta / fCoolingFactor;
    }
    else
    {
      beta = beta / fCoolingFactor;
      splitAll(y);
    }

    // make sure we are not too far from equilibrium before cooling further
    niter = 0;
    while((update1(beta, tks, y) > 1.e-6) && (niter++ < fMaxIterations))
    {
    }
  }

  if(fUseTc)
  {
    // last round of splitting, make sure no critical clusters are left
    update1(beta, tks, y);
    while(merge(y, beta))
    {
      update1(beta, tks, y);
    }
    unsigned int ntry = 0;
    while(split(beta, tks, y) && (ntry++ < 10))
    {
      niter = 0;
      while((update1(beta, tks, y) > 1.e-6) && (niter++ < fMaxIterations))
      {
      }
      merge(y, beta);
      update1(beta, tks, y);
    }
  }
  else
  {
    // merge collapsed clusters
    while(merge(y, beta))
    {
      update1(beta, tks, y);
    }
    if(fVerbose)
    {
      cout << "dump after 1st merging " << endl;
      dump(beta, y, tks);
    }
  }

  // switch on outlier rejection
  rho0 = 1. / nt;
  for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
  {
    k->pk = 1.;
  } // democratic
  niter = 0;
  while((update2(beta, tks, y, rho0, fDzCutOff) > 1.e-8) && (niter++ < fMaxIterations))
  {
  }
  if(fVerbose)
  {
    cout << "rho0=" << rho0 << " niter=" << niter << endl;
    dump(beta, y, tks);
  }

  // merge again  (some cluster split by outliers collapse here)
  while(merge(y))
  {
  }
  if(fVerbose)
  {
    cout << "dump after 2nd merging " << endl;
    dump(beta, y, tks);
  }

  // continue from freeze-out to Tstop (=1) without splitting, eliminate insignificant vertices
  while(beta <= fBetaStop)
  {
    while(purge(y, tks, rho0, beta, fDzCutOff))
    {
      niter = 0;
      while((update2(beta, tks, y, rho0, fDzCutOff) > 1.e-6) && (niter++ < fMaxIterations))
      {
      }
    }
    beta /= fCoolingFactor;
    niter = 0;
    while((update2(beta, tks, y, rho0, fDzCutOff) > 1.e-6) && (niter++ < fMaxIterations))
    {
    }
  }

  //   // new, one last round of cleaning at T=Tstop
  //   while(purge(y,tks,rho0, beta)){
  //     niter=0; while((update2(beta, tks,y,rho0, fDzCutOff) > 1.e-6)  && (niter++ < fMaxIterations)){  }
  //   }

  if(fVerbose)
  {
    cout << "Final result, rho0=" << rho0 << endl;
    dump(beta, y, tks);
  }

  // select significant tracks and use a TransientVertex as a container
  //GlobalError dummyError;

  // ensure correct normalization of probabilities, should make double assginment reasonably impossible
  for(unsigned int i = 0; i < nt; i++)
  {
    tks[i].Z = rho0 * exp(-beta * (fDzCutOff * fDzCutOff));
    for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
    {
      tks[i].Z += k->pk * exp(-beta * Eik(tks[i], *k));
    }
  }

  for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
  {

    DelphesFactory *factory = GetFactory();
    candidate = factory->NewCandidate();

    //cout<<"new vertex"<<endl;
    //GlobalPoint pos(0, 0, k->z);
    double time = k->t;
    double z = k->z;
    //vector< reco::TransientTrack > vertexTracks;
    //double max_track_time_err2 = 0;
    double mean = 0.;
    double expv_x2 = 0.;
    double normw = 0.;
    for(unsigned int i = 0; i < nt; i++)
    {
      const double invdt = 1.0 / std::sqrt(tks[i].dt2);
      if(tks[i].Z > 0)
      {
        double p = k->pk * exp(-beta * Eik(tks[i], *k)) / tks[i].Z;
        if((tks[i].pi > 0) && (p > 0.5))
        {
          //std::cout << "pushing back " << i << ' ' << tks[i].tt << std::endl;
          //vertexTracks.push_back(*(tks[i].tt)); tks[i].Z=0;

          candidate->AddCandidate(tks[i].tt);
          tks[i].Z = 0;

          mean += tks[i].t * invdt * p;
          expv_x2 += tks[i].t * tks[i].t * invdt * p;
          normw += invdt * p;
        } // setting Z=0 excludes double assignment
      }
    }

    mean = mean / normw;
    expv_x2 = expv_x2 / normw;
    const double time_var = expv_x2 - mean * mean;
    const double crappy_error_guess = std::sqrt(time_var);
    /*GlobalError dummyErrorWithTime(0,
                                   0,0,
                                   0,0,0,
                                   0,0,0,crappy_error_guess);*/
    //TransientVertex v(pos, time, dummyErrorWithTime, vertexTracks, 5);

    candidate->ClusterIndex = clusterIndex++;
    ;
    candidate->Position.SetXYZT(0.0, 0.0, z * 10.0, time * c_light);

    // TBC - fill error later ...
    candidate->PositionError.SetXYZT(0.0, 0.0, 0.0, crappy_error_guess * c_light);

    clusterIndex++;
    clusters.push_back(candidate);
  }

  return clusters;
}

//------------------------------------------------------------------------------

static double Eik(const track_t &t, const vertex_t &k)
{
  return std::pow(t.z - k.z, 2.) / t.dz2 + std::pow(t.t - k.t, 2.) / t.dt2;
}

//------------------------------------------------------------------------------

static void dump(const double beta, const vector<vertex_t> &y, const vector<track_t> &tks0)
{
  // copy and sort for nicer printout
  vector<track_t> tks;
  for(vector<track_t>::const_iterator t = tks0.begin(); t != tks0.end(); t++)
  {
    tks.push_back(*t);
  }
  std::stable_sort(tks.begin(), tks.end(), recTrackLessZ1);

  cout << "-----DAClusterizerInZT::dump ----" << endl;
  cout << " beta=" << beta << endl;
  cout << "                                                               z= ";
  cout.precision(4);
  for(vector<vertex_t>::const_iterator k = y.begin(); k != y.end(); k++)
  {
    //cout  <<  setw(8) << fixed << k->z;
  }
  cout << endl
       << "                                                               t= ";
  for(vector<vertex_t>::const_iterator k = y.begin(); k != y.end(); k++)
  {
    //cout  <<  setw(8) << fixed << k->t;
  }
  //cout << endl << "T=" << setw(15) << 1./beta <<"                                             Tc= ";
  for(vector<vertex_t>::const_iterator k = y.begin(); k != y.end(); k++)
  {
    //cout  <<  setw(8) << fixed << k->Tc ;
  }

  cout << endl
       << "                                                              pk=";
  double sumpk = 0;
  for(vector<vertex_t>::const_iterator k = y.begin(); k != y.end(); k++)
  {
    //cout <<  setw(8) <<  setprecision(3) <<  fixed << k->pk;
    sumpk += k->pk;
  }
  cout << endl;

  double E = 0, F = 0;
  cout << endl;
  cout << "----       z +/- dz        t +/- dt        ip +/-dip       pt    phi  eta    weights  ----" << endl;
  cout.precision(4);
  for(unsigned int i = 0; i < tks.size(); i++)
  {
    if(tks[i].Z > 0)
    {
      F -= log(tks[i].Z) / beta;
    }
    double tz = tks[i].z;
    double tt = tks[i].t;
    //cout <<  setw (3)<< i << ")" <<  setw (8) << fixed << setprecision(4)<<  tz << " +/-" <<  setw (6)<< sqrt(tks[i].dz2)
    //     << setw(8) << fixed << setprecision(4) << tt << " +/-" << setw(6) << std::sqrt(tks[i].dt2)  ;

    double sump = 0.;
    for(vector<vertex_t>::const_iterator k = y.begin(); k != y.end(); k++)
    {
      if((tks[i].pi > 0) && (tks[i].Z > 0))
      {
        //double p=pik(beta,tks[i],*k);
        double p = k->pk * std::exp(-beta * Eik(tks[i], *k)) / tks[i].Z;
        if(p > 0.0001)
        {
          //cout <<  setw (8) <<  setprecision(3) << p;
        }
        else
        {
          cout << "    .   ";
        }
        E += p * Eik(tks[i], *k);
        sump += p;
      }
      else
      {
        cout << "        ";
      }
    }
    cout << endl;
  }
  cout << endl
       << "T=" << 1 / beta << " E=" << E << " n=" << y.size() << "  F= " << F << endl
       << "----------" << endl;
}

//------------------------------------------------------------------------------

static double update1(double beta, vector<track_t> &tks, vector<vertex_t> &y)
{
  //update weights and vertex positions
  // mass constrained annealing without noise
  // returns the squared sum of changes of vertex positions

  unsigned int nt = tks.size();

  //initialize sums
  double sumpi = 0;
  for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); ++k)
  {
    k->sw = 0.;
    k->swz = 0.;
    k->swt = 0.;
    k->se = 0.;
    k->swE = 0.;
    k->Tc = 0.;
  }

  // loop over tracks
  for(unsigned int i = 0; i < nt; i++)
  {

    // update pik and Zi
    double Zi = 0.;
    for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); ++k)
    {
      k->ei = std::exp(-beta * Eik(tks[i], *k)); // cache exponential for one track at a time
      Zi += k->pk * k->ei;
    }
    tks[i].Z = Zi;

    // normalization for pk
    if(tks[i].Z > 0)
    {
      sumpi += tks[i].pi;
      // accumulate weighted z and weights for vertex update
      for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); ++k)
      {
        k->se += tks[i].pi * k->ei / Zi;
        const double w = k->pk * tks[i].pi * k->ei / (Zi * (tks[i].dz2 * tks[i].dt2));
        k->sw += w;
        k->swz += w * tks[i].z;
        k->swt += w * tks[i].t;
        k->swE += w * Eik(tks[i], *k);
      }
    }
    else
    {
      sumpi += tks[i].pi;
    }

  } // end of track loop

  // now update z and pk
  double delta = 0;
  for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
  {
    if(k->sw > 0)
    {
      const double znew = k->swz / k->sw;
      const double tnew = k->swt / k->sw;
      delta += std::pow(k->z - znew, 2.) + std::pow(k->t - tnew, 2.);
      k->z = znew;
      k->t = tnew;
      k->Tc = 2. * k->swE / k->sw;
    }
    else
    {
      // cout << " a cluster melted away ?  pk=" << k->pk <<  " sumw=" << k->sw <<  endl
      k->Tc = -1;
    }

    k->pk = k->pk * k->se / sumpi;
  }

  // return how much the prototypes moved
  return delta;
}

//------------------------------------------------------------------------------

static double update2(double beta, vector<track_t> &tks, vector<vertex_t> &y, double &rho0, double dzCutOff)
{
  // MVF style, no more vertex weights, update tracks weights and vertex positions, with noise
  // returns the squared sum of changes of vertex positions

  unsigned int nt = tks.size();

  //initialize sums
  for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
  {
    k->sw = 0.;
    k->swz = 0.;
    k->swt = 0.;
    k->se = 0.;
    k->swE = 0.;
    k->Tc = 0.;
  }

  // loop over tracks
  for(unsigned int i = 0; i < nt; i++)
  {

    // update pik and Zi and Ti
    double Zi = rho0 * std::exp(-beta * (dzCutOff * dzCutOff)); // cut-off (eventually add finite size in time)
    //double Ti = 0.; // dt0*std::exp(-beta*fDtCutOff);
    for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
    {
      k->ei = std::exp(-beta * Eik(tks[i], *k)); // cache exponential for one track at a time
      Zi += k->pk * k->ei;
    }
    tks[i].Z = Zi;

    // normalization
    if(tks[i].Z > 0)
    {
      // accumulate weighted z and weights for vertex update
      for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
      {
        k->se += tks[i].pi * k->ei / Zi;
        double w = k->pk * tks[i].pi * k->ei / (Zi * (tks[i].dz2 * tks[i].dt2));
        k->sw += w;
        k->swz += w * tks[i].z;
        k->swt += w * tks[i].t;
        k->swE += w * Eik(tks[i], *k);
      }
    }

  } // end of track loop

  // now update z
  double delta = 0;
  for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
  {
    if(k->sw > 0)
    {
      const double znew = k->swz / k->sw;
      const double tnew = k->swt / k->sw;
      delta += std::pow(k->z - znew, 2.) + std::pow(k->t - tnew, 2.);
      k->z = znew;
      k->t = tnew;
      k->Tc = 2 * k->swE / k->sw;
    }
    else
    {
      // cout << " a cluster melted away ?  pk=" << k->pk <<  " sumw=" << k->sw <<  endl;
      k->Tc = 0;
    }
  }

  // return how much the prototypes moved
  return delta;
}

//------------------------------------------------------------------------------

static bool merge(vector<vertex_t> &y)
{
  // merge clusters that collapsed or never separated, return true if vertices were merged, false otherwise

  if(y.size() < 2) return false;

  for(vector<vertex_t>::iterator k = y.begin(); (k + 1) != y.end(); k++)
  {
    if(std::abs((k + 1)->z - k->z) < 1.e-3 && std::abs((k + 1)->t - k->t) < 1.e-3)
    { // with fabs if only called after freeze-out (splitAll() at highter T)
      double rho = k->pk + (k + 1)->pk;
      if(rho > 0)
      {
        k->z = (k->pk * k->z + (k + 1)->z * (k + 1)->pk) / rho;
        k->t = (k->pk * k->t + (k + 1)->t * (k + 1)->pk) / rho;
      }
      else
      {
        k->z = 0.5 * (k->z + (k + 1)->z);
        k->t = 0.5 * (k->t + (k + 1)->t);
      }
      k->pk = rho;

      y.erase(k + 1);
      return true;
    }
  }

  return false;
}

//------------------------------------------------------------------------------

static bool merge(vector<vertex_t> &y, double &beta)
{
  // merge clusters that collapsed or never separated,
  // only merge if the estimated critical temperature of the merged vertex is below the current temperature
  // return true if vertices were merged, false otherwise
  if(y.size() < 2) return false;

  for(vector<vertex_t>::iterator k = y.begin(); (k + 1) != y.end(); k++)
  {
    if(std::abs((k + 1)->z - k->z) < 2.e-3 && std::abs((k + 1)->t - k->t) < 2.e-3)
    {
      double rho = k->pk + (k + 1)->pk;
      double swE = k->swE + (k + 1)->swE - k->pk * (k + 1)->pk / rho * (std::pow((k + 1)->z - k->z, 2.) + std::pow((k + 1)->t - k->t, 2.));
      double Tc = 2 * swE / (k->sw + (k + 1)->sw);

      if(Tc * beta < 1)
      {
        if(rho > 0)
        {
          k->z = (k->pk * k->z + (k + 1)->z * (k + 1)->pk) / rho;
          k->t = (k->pk * k->t + (k + 1)->t * (k + 1)->pk) / rho;
        }
        else
        {
          k->z = 0.5 * (k->z + (k + 1)->z);
          k->t = 0.5 * (k->t + (k + 1)->t);
        }
        k->pk = rho;
        k->sw += (k + 1)->sw;
        k->swE = swE;
        k->Tc = Tc;
        y.erase(k + 1);
        return true;
      }
    }
  }

  return false;
}

//------------------------------------------------------------------------------

static bool purge(vector<vertex_t> &y, vector<track_t> &tks, double &rho0, const double beta, const double dzCutOff)
{
  // eliminate clusters with only one significant/unique track
  if(y.size() < 2) return false;

  unsigned int nt = tks.size();
  double sumpmin = nt;
  vector<vertex_t>::iterator k0 = y.end();
  for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
  {
    int nUnique = 0;
    double sump = 0;
    double pmax = k->pk / (k->pk + rho0 * exp(-beta * dzCutOff * dzCutOff));
    for(unsigned int i = 0; i < nt; i++)
    {
      if(tks[i].Z > 0)
      {
        double p = k->pk * std::exp(-beta * Eik(tks[i], *k)) / tks[i].Z;
        sump += p;
        if((p > 0.9 * pmax) && (tks[i].pi > 0))
        {
          nUnique++;
        }
      }
    }

    if((nUnique < 2) && (sump < sumpmin))
    {
      sumpmin = sump;
      k0 = k;
    }
  }

  if(k0 != y.end())
  {
    //cout << "eliminating prototype at " << k0->z << "," << k0->t << " with sump=" << sumpmin << endl;
    //rho0+=k0->pk;
    y.erase(k0);
    return true;
  }
  else
  {
    return false;
  }
}

//------------------------------------------------------------------------------

static double beta0(double betamax, vector<track_t> &tks, vector<vertex_t> &y, const double coolingFactor)
{

  double T0 = 0; // max Tc for beta=0
  // estimate critical temperature from beta=0 (T=inf)
  unsigned int nt = tks.size();

  for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
  {

    // vertex fit at T=inf
    double sumwz = 0.;
    double sumwt = 0.;
    double sumw = 0.;
    for(unsigned int i = 0; i < nt; i++)
    {
      double w = tks[i].pi / (tks[i].dz2 * tks[i].dt2);
      sumwz += w * tks[i].z;
      sumwt += w * tks[i].t;
      sumw += w;
    }
    k->z = sumwz / sumw;
    k->t = sumwt / sumw;

    // estimate Tcrit, eventually do this in the same loop
    double a = 0, b = 0;
    for(unsigned int i = 0; i < nt; i++)
    {
      double dx = tks[i].z - (k->z);
      double dt = tks[i].t - (k->t);
      double w = tks[i].pi / (tks[i].dz2 * tks[i].dt2);
      a += w * (std::pow(dx, 2.) / tks[i].dz2 + std::pow(dt, 2.) / tks[i].dt2);
      b += w;
    }
    double Tc = 2. * a / b; // the critical temperature of this vertex
    if(Tc > T0) T0 = Tc;
  } // vertex loop (normally there should be only one vertex at beta=0)

  if(T0 > 1. / betamax)
  {
    return betamax / pow(coolingFactor, int(std::log(T0 * betamax) / std::log(coolingFactor)) - 1);
  }
  else
  {
    // ensure at least one annealing step
    return betamax / coolingFactor;
  }
}

//------------------------------------------------------------------------------

static bool split(double beta, vector<track_t> &tks, vector<vertex_t> &y)
{
  // split only critical vertices (Tc >~ T=1/beta   <==>   beta*Tc>~1)
  // an update must have been made just before doing this (same beta, no merging)
  // returns true if at least one cluster was split

  const double epsilon = 1e-3; // split all single vertices by 10 um
  bool split = false;

  // avoid left-right biases by splitting highest Tc first

  std::vector<std::pair<double, unsigned int> > critical;
  for(unsigned int ik = 0; ik < y.size(); ik++)
  {
    if(beta * y[ik].Tc > 1.)
    {
      critical.push_back(make_pair(y[ik].Tc, ik));
    }
  }
  std::stable_sort(critical.begin(), critical.end(), std::greater<std::pair<double, unsigned int> >());

  for(unsigned int ic = 0; ic < critical.size(); ic++)
  {
    unsigned int ik = critical[ic].second;
    // estimate subcluster positions and weight
    double p1 = 0, z1 = 0, t1 = 0, w1 = 0;
    double p2 = 0, z2 = 0, t2 = 0, w2 = 0;
    //double sumpi=0;
    for(unsigned int i = 0; i < tks.size(); i++)
    {
      if(tks[i].Z > 0)
      {
        //sumpi+=tks[i].pi;
        double p = y[ik].pk * exp(-beta * Eik(tks[i], y[ik])) / tks[i].Z * tks[i].pi;
        double w = p / (tks[i].dz2 * tks[i].dt2);
        if(tks[i].z < y[ik].z)
        {
          p1 += p;
          z1 += w * tks[i].z;
          t1 += w * tks[i].t;
          w1 += w;
        }
        else
        {
          p2 += p;
          z2 += w * tks[i].z;
          t2 += w * tks[i].t;
          w2 += w;
        }
      }
    }
    if(w1 > 0)
    {
      z1 = z1 / w1;
      t1 = t1 / w1;
    }
    else
    {
      z1 = y[ik].z - epsilon;
      t1 = y[ik].t - epsilon;
    }
    if(w2 > 0)
    {
      z2 = z2 / w2;
      t2 = t2 / w2;
    }
    else
    {
      z2 = y[ik].z + epsilon;
      t2 = y[ik].t + epsilon;
    }

    // reduce split size if there is not enough room
    if((ik > 0) && (y[ik - 1].z >= z1))
    {
      z1 = 0.5 * (y[ik].z + y[ik - 1].z);
      t1 = 0.5 * (y[ik].t + y[ik - 1].t);
    }
    if((ik + 1 < y.size()) && (y[ik + 1].z <= z2))
    {
      z2 = 0.5 * (y[ik].z + y[ik + 1].z);
      t2 = 0.5 * (y[ik].t + y[ik + 1].t);
    }

    // split if the new subclusters are significantly separated
    if((z2 - z1) > epsilon || std::abs(t2 - t1) > epsilon)
    {
      split = true;
      vertex_t vnew;
      vnew.pk = p1 * y[ik].pk / (p1 + p2);
      y[ik].pk = p2 * y[ik].pk / (p1 + p2);
      vnew.z = z1;
      vnew.t = t1;
      y[ik].z = z2;
      y[ik].t = t2;
      y.insert(y.begin() + ik, vnew);

      // adjust remaining pointers
      for(unsigned int jc = ic; jc < critical.size(); jc++)
      {
        if(critical[jc].second > ik)
        {
          critical[jc].second++;
        }
      }
    }
  }

  //  stable_sort(y.begin(), y.end(), clusterLessZ);
  return split;
}

//------------------------------------------------------------------------------

void splitAll(vector<vertex_t> &y)
{

  const double epsilon = 1e-3; // split all single vertices by 10 um
  const double zsep = 2 * epsilon; // split vertices that are isolated by at least zsep (vertices that haven't collapsed)
  const double tsep = 2 * epsilon; // check t as well

  vector<vertex_t> y1;

  for(vector<vertex_t>::iterator k = y.begin(); k != y.end(); k++)
  {
    if(((k == y.begin()) || (k - 1)->z < k->z - zsep) && (((k + 1) == y.end()) || (k + 1)->z > k->z + zsep))
    {
      // isolated prototype, split
      vertex_t vnew;
      vnew.z = k->z - epsilon;
      vnew.t = k->t - epsilon;
      (*k).z = k->z + epsilon;
      (*k).t = k->t + epsilon;
      vnew.pk = 0.5 * (*k).pk;
      (*k).pk = 0.5 * (*k).pk;
      y1.push_back(vnew);
      y1.push_back(*k);
    }
    else if(y1.empty() || (y1.back().z < k->z - zsep) || (y1.back().t < k->t - tsep))
    {
      y1.push_back(*k);
    }
    else
    {
      y1.back().z -= epsilon;
      y1.back().t -= epsilon;
      k->z += epsilon;
      k->t += epsilon;
      y1.push_back(*k);
    }
  } // vertex loop

  y = y1;
}
