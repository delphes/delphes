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

/** \class ARCTracing
 *
 *  Ray-tracing fast simulation of the ARC (see the header). Lengths in mm, energies in eV, angles in rad.
 *
 *  \author ALFA ARC study
 */

#include "modules/ARCTracing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TRandom3.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

using namespace std;

static const double kMassPi = 0.13957039, kMassK = 0.493677, kMassP = 0.93827209;
static const double kFrankTamm = 36.98; // Cherenkov photons per eV per mm at sin^2(theta) = 1: alpha / (hbar c)
static const double kDeg = TMath::Pi() / 180.0;

static inline double Dot(const double *a, const double *b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
static inline double Norm(const double *a) { return sqrt(Dot(a, a)); }
static inline void Unit(double *a)
{
  double n = Norm(a);
  if(n > 0) { a[0] /= n; a[1] /= n; a[2] /= n; }
}
static inline void Cross(const double *a, const double *b, double *c)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

//------------------------------------------------------------------------------

ARCTracing::ARCTracing() :
  fItInputArray(0), fInputArray(0), fOutputArray(0)
{
}

ARCTracing::~ARCTracing()
{
}

//------------------------------------------------------------------------------

double ARCTracing::Table::Eval(double x) const
{
  if(e.empty()) return 0.0;
  if(x <= e.front()) return v.front();
  if(x >= e.back()) return v.back();
  size_t i = upper_bound(e.begin(), e.end(), x) - e.begin();
  double f = (x - e[i - 1]) / (e[i] - e[i - 1]);
  return v[i - 1] + f * (v[i] - v[i - 1]);
}

//------------------------------------------------------------------------------

void ARCTracing::ReadTables(const char *file)
{
  ifstream in(file);
  if(!in) throw runtime_error(string("ARCTracing: cannot open the optics file ") + file);
  string line;
  while(getline(in, line))
  {
    if(line.empty() || line[0] == '#') continue;
    istringstream ss(line);
    string key, name;
    int n;
    ss >> key >> name >> n;
    if(key != "TABLE") continue;
    Table t;
    for(int i = 0; i < n && getline(in, line); ++i)
    {
      istringstream row(line);
      double e, v;
      row >> e >> v;
      t.e.push_back(e);
      t.v.push_back(v);
    }
    fTables.push_back(t);
    fTableNames.push_back(name);
  }
}

const ARCTracing::Table *ARCTracing::FindTable(const string &name) const
{
  for(size_t i = 0; i < fTableNames.size(); ++i)
    if(fTableNames[i] == name) return &fTables[i];
  throw runtime_error("ARCTracing: optical table not found: " + name);
}

//------------------------------------------------------------------------------

void ARCTracing::ReadCells(const char *file)
{
  ifstream in(file);
  if(!in) throw runtime_error(string("ARCTracing: cannot open the cell file ") + file);
  string line;
  fBarrelRin = 1e9;
  fEndcapZin = 1e9;
  while(getline(in, line))
  {
    if(line.empty() || line[0] == '#') continue;
    istringstream ss(line);
    string region, name;
    ss >> region >> name;
    Cell c;
    c.region = (region == "barrel") ? kBarrel : kEndcap;
    double Rmax, aeroN[3], aeroRmax, aero[3];
    for(int i = 0; i < 3; ++i) ss >> c.O[i];
    for(int i = 0; i < 3; ++i) ss >> c.ex[i];
    for(int i = 0; i < 3; ++i) ss >> c.ey[i];
    for(int i = 0; i < 3; ++i) ss >> c.ez[i];
    for(int i = 0; i < 3; ++i) ss >> c.C[i];
    ss >> c.R >> Rmax;
    for(int i = 0; i < 3; ++i) ss >> c.S[i];
    for(int i = 0; i < 3; ++i) ss >> c.sex[i];
    for(int i = 0; i < 3; ++i) ss >> c.sey[i];
    for(int i = 0; i < 3; ++i) ss >> c.sez[i];
    for(int i = 0; i < 3; ++i) ss >> c.hs[i];
    for(int i = 0; i < 3; ++i) ss >> aero[i];
    for(int i = 0; i < 3; ++i) ss >> aeroN[i];
    ss >> c.aeroHalf >> aeroRmax >> c.phi0 >> c.zmin >> c.zmax >> c.apothem;
    if(!ss) throw runtime_error("ARCTracing: malformed cell line: " + line);
    c.aeroZ = aero[2];
    // the frames are written with a finite number of digits: make them exactly orthonormal (the Cherenkov angle of the gas,
    // 50 mrad, is sensitive to 1e-4 in a direction's norm through cot(theta_c))
    Unit(c.ez);
    double pz = Dot(c.ex, c.ez);
    for(int i = 0; i < 3; ++i) c.ex[i] -= pz * c.ez[i];
    Unit(c.ex);
    Cross(c.ez, c.ex, c.ey);
    Unit(c.sez);
    double qz = Dot(c.sex, c.sez);
    for(int i = 0; i < 3; ++i) c.sex[i] -= qz * c.sez[i];
    Unit(c.sex);
    double sy[3];
    Cross(c.sez, c.sex, sy);
    if(Dot(sy, c.sey) < 0)
      for(int i = 0; i < 3; ++i) sy[i] = -sy[i];
    for(int i = 0; i < 3; ++i) c.sey[i] = sy[i];
    c.side = (c.O[2] + c.ez[2] * c.zmin < 0.0) ? -1 : 1;
    // the sensor normal pointing to the mirror surface, and the detection plane on that face
    double apex[3] = {c.C[0], c.C[1], c.C[2] + c.R};
    double toApex[3] = {apex[0] - c.S[0], apex[1] - c.S[1], apex[2] - c.S[2]};
    if(Dot(c.sez, toApex) < 0.0)
      for(int i = 0; i < 3; ++i) c.sez[i] = -c.sez[i];
    for(int i = 0; i < 3; ++i) c.S[i] += c.hs[2] * c.sez[i];
    // design knobs: sensor moved along the cell axis, mirror curvature changed with the apex fixed
    c.S[2] += fSensorShift;
    if(fMirrorRadiusScale != 1.0)
    {
      double Rnew = c.R * fMirrorRadiusScale;
      c.C[2] = apex[2] - Rnew;
      c.R = Rnew;
    }
    c.phiAxis = atan2(c.ez[1], c.ez[0]);
    c.zAxis = c.O[2];
    if(c.region == kBarrel)
      fBarrelRin = min(fBarrelRin, c.zmin);
    else
      fEndcapZin = min(fEndcapZin, fabs(c.O[2] + c.ez[2] * c.zmin));
    fCells.push_back(c);
  }
  if(fCells.empty()) throw runtime_error("ARCTracing: no cells read");
}

//------------------------------------------------------------------------------

void ARCTracing::BuildGrid(Radiator &r, const Table &pde, const Table &refl)
{
  double emin = max(r.n.e.front(), pde.e.front()), emax = min(r.n.e.back(), pde.e.back());
  const int N = 160;
  r.gridE.resize(N);
  r.gridN.resize(N);
  r.gridW.resize(N);
  for(int i = 0; i < N; ++i)
  {
    double e = emin + (emax - emin) * i / (N - 1.0);
    r.gridE[i] = e;
    r.gridN[i] = 1.0 + r.indexScale * (r.n.Eval(e) - 1.0);
    r.gridW[i] = fPDEScale * pde.Eval(e) * (refl.Ok() ? refl.Eval(e) : 1.0);
  }
}

// detected photons per mm of path (before the geometry), and the cumulative energy distribution of the detected photons
double ARCTracing::Yield(const Radiator &r, double beta, vector<double> *cdf) const
{
  const size_t N = r.gridE.size();
  double sum = 0.0;
  if(cdf) cdf->assign(N, 0.0);
  for(size_t i = 1; i < N; ++i)
  {
    double s0 = 1.0 - 1.0 / pow(r.gridN[i - 1] * beta, 2), s1 = 1.0 - 1.0 / pow(r.gridN[i] * beta, 2);
    double f = 0.5 * (max(s0, 0.0) * r.gridW[i - 1] + max(s1, 0.0) * r.gridW[i]) * (r.gridE[i] - r.gridE[i - 1]);
    sum += f;
    if(cdf) (*cdf)[i] = sum;
  }
  if(cdf && sum > 0)
    for(size_t i = 0; i < N; ++i) (*cdf)[i] /= sum;
  return kFrankTamm * fYieldScale * sum;
}

double ARCTracing::SampleEnergy(const Radiator &r, const vector<double> &cdf) const
{
  double u = gRandom->Uniform();
  size_t i = upper_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
  if(i == 0) return r.gridE.front();
  if(i >= cdf.size()) return r.gridE.back();
  double f = (u - cdf[i - 1]) / max(cdf[i] - cdf[i - 1], 1e-12);
  return r.gridE[i - 1] + f * (r.gridE[i] - r.gridE[i - 1]);
}

//------------------------------------------------------------------------------

void ARCTracing::Init()
{
  fBz = GetDouble("Bz", 3.0);
  string cellFile = GetString("CellFile", "cells_ALFA_o1_v01.txt");
  string opticsFile = GetString("OpticsFile", "arc_optics_tables.txt");
  ReadTables(opticsFile.c_str());

  fPDE = *FindTable(GetString("PDETable", "ARC_SiPM_QuantumEfficiency"));
  fPDEScale = GetDouble("PDEScale", 1.0);
  string reflName = GetString("MirrorReflectivityTable", "REFLECTIVITY__ARC_Mirror");
  if(!reflName.empty()) fRefl = *FindTable(reflName);

  fRad[kGas].n = *FindTable(GetString("GasIndexTable", "RINDEX__Xenon_ARC"));
  fRad[kGas].indexScale = GetDouble("GasIndexScale", 1.0);
  string gasAbs = GetString("GasAbsorptionTable", "");
  if(!gasAbs.empty()) fRad[kGas].abs = *FindTable(gasAbs);
  fRad[kGas].thickness = 0.0;

  fRad[kAerogel].n = *FindTable(GetString("AerogelIndexTable", "RINDEX__Aerogel_PFRICH"));
  fRad[kAerogel].indexScale = GetDouble("AerogelIndexScale", 1.0);
  // defaults as the Geant4 material of the ALFA full simulation (k4geo Aerogel_PFRICH: RINDEX and Rayleigh scattering with the
  // Belle II table; absorption deliberately not applied, see the comment in materials_arc_o1_v01.xml); scattered photons continue
  string aeroAbs = GetString("AerogelAbsorptionTable", "");
  string aeroRay = GetString("AerogelRayleighTable", "RAYLEIGH__Aerogel_BelleII");
  if(!aeroAbs.empty()) fRad[kAerogel].abs = *FindTable(aeroAbs);
  if(!aeroRay.empty()) fRad[kAerogel].ray = *FindTable(aeroRay);
  fRad[kAerogel].thickness = GetDouble("AerogelThickness", 0.0); // 0: from the geometry
  fAerogelAbsorption = GetBool("AerogelAbsorption", true);
  fAerogelRayleigh = GetBool("AerogelRayleigh", true);
  fAerogelSlabInPath = GetBool("AerogelSlabInPath", true);

  fPixelPitch = GetDouble("PixelPitch", 0.8);
  fSensorSizeScale = GetDouble("SensorSizeScale", 1.0);
  fSensorShift = GetDouble("SensorShift", 0.0);
  fMirrorRadiusScale = GetDouble("MirrorRadiusScale", 1.0);
  fYieldScale = GetDouble("YieldScale", 1.0);
  fDarkCounts = GetDouble("DarkCountsPerSensor", 0.0);

  fRecoEnergy = GetDouble("RecoEnergy", 3.0);
  fRecoRadius[kGas] = GetDouble("RecoRadiusGas", 2000.0);
  fRecoRadius[kAerogel] = GetDouble("RecoRadiusAerogel", 1925.0);
  fRecoZ[kGas] = GetDouble("RecoZGas", 2110.0);
  fRecoZ[kAerogel] = GetDouble("RecoZAerogel", 2025.0);
  fWindow = GetDouble("Window", 0.025);
  fMinPe = GetDouble("MinPe", 3.0);
  fSigmaCorr[kGas] = GetDouble("GasSigmaCorr", 0.0);
  fSigmaCorr[kAerogel] = GetDouble("AerogelSigmaCorr", 0.0);
  fSigmaSingleFloor[kGas] = GetDouble("GasSigmaSingleFloor", 1.0E-3);
  fSigmaSingleFloor[kAerogel] = GetDouble("AerogelSigmaSingleFloor", 1.5E-3);
  fSigmaTrack = GetDouble("SigmaTrack", 0.0);
  fUseTrueTrack = GetBool("UseTrueTrack", true);
  fUseTrueMomentumForPID = GetBool("UseTrueMomentumForPID", false);
  fTrueHypothesisAngle = GetBool("TrueHypothesisAngle", false);
  fSharedCells = GetBool("SharedCells", true); // reconstruct every track against all photons of its cell (overlapping rings in jets)
  fDebugMonoEnergy = GetDouble("DebugMonoEnergy", 0.0); // > 0: all photons at this energy (no chromatic error), for studies
  fDebugMidEmission = GetBool("DebugMidEmission", false); // all photons from the middle of the radiator path, for studies

  for(int i = 0; i < 2; ++i)
  {
    BuildGrid(fRad[i], fPDE, fRefl);
    fRecoIndex[i] = 1.0 + fRad[i].indexScale * (fRad[i].n.Eval(fRecoEnergy) - 1.0);
  }
  ReadCells(cellFile.c_str());

  fNtracks = fNcovered = 0;
  fNgen[0] = fNgen[1] = fNdet[0] = fNdet[1] = 0;
  fCount = kFALSE;
  fDebugFile = GetString("DebugFile", "");
  fDebug = fDebugFile.empty() ? 0 : new ofstream(fDebugFile.c_str());
  if(fDebug) *fDebug << "# T itrack region cell pdg p xEntry yEntry thetaTrueGas thetaTrueAero nGenGas nGenAero | H u v thetaGas thetaAero [mm, rad]\n";
  for(int i = 0; i < 2; ++i)
    for(int j = 0; j < 20; ++j) fLoss[i][j] = 0;

  fInputArray = ImportArray(GetString("InputArray", "TrackSmearing/tracks"));
  fItInputArray = fInputArray->MakeIterator();
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));

  cout << "ARCTracing: " << fCells.size() << " cells (barrel r_in " << fBarrelRin << " mm, endcap |z_in| " << fEndcapZin
       << " mm), n_gas(" << fRecoEnergy << " eV) = " << fRecoIndex[kGas] << ", n_aerogel = " << fRecoIndex[kAerogel]
       << ", gas yield at beta = 1: " << Yield(fRad[kGas], 1.0, 0) * 10 << " /cm, aerogel: " << Yield(fRad[kAerogel], 1.0, 0) * 10
       << " /cm (before the optics)" << endl;
}

void ARCTracing::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if(fDebug)
  {
    fDebug->close();
    delete fDebug;
    fDebug = 0;
  }
  cout << "ARCTracing: " << fNtracks << " tracks, " << fNcovered << " in the ARC; photons generated / detected: gas "
       << fNgen[kGas] << " / " << fNdet[kGas] << ", aerogel " << fNgen[kAerogel] << " / " << fNdet[kAerogel] << endl;
  const char *why[18] = {"below threshold", "aerogel: absorbed (up)", "aerogel: stuck (up)", "aerogel: total reflection (up)",
                         "no mirror intersection", "mirror behind", "gas absorption", "mirror outside the cell", "mirror reflectivity",
                         "after mirror: not downwards", "slab top behind", "slab: absorbed (back)", "slab: scattered back up",
                         "slab: total reflection (out)", "sensor: not facing", "sensor plane behind", "outside the SiPM array", "gas below slab: emitted downwards"};
  for(int rad = 0; rad < 2; ++rad)
  {
    cout << "ARCTracing: losses " << (rad == kGas ? "gas" : "aerogel") << ":";
    for(int j = 0; j < 18; ++j)
      if(fLoss[rad][j] > 0) cout << " [" << why[j] << " " << fLoss[rad][j] << "]";
    cout << endl;
  }
}

//------------------------------------------------------------------------------

void ARCTracing::Helix::Point(double sT, double *P, double *T) const
{
  double c, s;
  if(straight)
  {
    c = cos(phi0);
    s = sin(phi0);
    P[0] = x0 + sT * c;
    P[1] = y0 + sT * s;
  }
  else
  {
    double psi = -q * sT / rc;
    c = cos(phi0 + psi);
    s = sin(phi0 + psi);
    P[0] = x0 + (rc / -q) * (s - sin(phi0));
    P[1] = y0 + (rc / -q) * (-c + cos(phi0));
  }
  P[2] = z0 + sT * pz / pt;
  T[0] = pt * c / p;
  T[1] = pt * s / p;
  T[2] = pz / p;
}

void ARCTracing::ToLocal(const Cell &c, const double *P, double *L) const
{
  double d[3] = {P[0] - c.O[0], P[1] - c.O[1], P[2] - c.O[2]};
  L[0] = Dot(d, c.ex);
  L[1] = Dot(d, c.ey);
  L[2] = Dot(d, c.ez);
}

void ARCTracing::VecToLocal(const Cell &c, const double *V, double *L) const
{
  L[0] = Dot(V, c.ex);
  L[1] = Dot(V, c.ey);
  L[2] = Dot(V, c.ez);
}

double ARCTracing::HexDist(double x, double y, double a, double phi0) const
{
  double d = 0.0;
  for(int k = 0; k < 3; ++k)
  {
    double ang = (phi0 + 30.0 + 60.0 * k) * kDeg;
    d = max(d, fabs(x * cos(ang) + y * sin(ang)) / a);
  }
  return d;
}

bool ARCTracing::InHex(double x, double y, double a, double phi0) const
{
  return HexDist(x, y, a, phi0) <= 1.0;
}

//------------------------------------------------------------------------------

// transverse path length at which the helix reaches the cylindrical radius r (Newton from the straight line)
double ARCTracing::SolveRadius(const Helix &h, double r, bool &ok) const
{
  ok = false;
  double P[3], T[3];
  // straight-line estimate: |X0_T + s u_T| = r
  double ux = cos(h.phi0), uy = sin(h.phi0);
  double b = h.x0 * ux + h.y0 * uy, c = h.x0 * h.x0 + h.y0 * h.y0 - r * r;
  double disc = b * b - c;
  if(disc < 0) return 0.0;
  double s = -b + sqrt(disc);
  if(s < 0) return 0.0;
  if(h.straight)
  {
    ok = true;
    return s;
  }
  for(int it = 0; it < 30; ++it)
  {
    h.Point(s, P, T);
    double rho = hypot(P[0], P[1]);
    double g = rho - r;
    double dg = (P[0] * T[0] + P[1] * T[1]) / max(rho, 1e-9) * (h.p / h.pt);
    if(fabs(dg) < 1e-6) return 0.0;
    double ds = -g / dg;
    s += ds;
    if(s < 0) return 0.0;
    if(fabs(ds) < 1e-4)
    {
      ok = (s < TMath::Pi() * h.rc); // first crossing, before the helix turns back
      return s;
    }
  }
  return 0.0;
}

// transverse path length at which the helix reaches the cell-local plane z = zloc
double ARCTracing::SolvePlane(const Helix &h, const Cell &c, double zloc, double sT0, bool &ok) const
{
  ok = false;
  double P[3], T[3], L[3], TL[3], s = sT0;
  for(int it = 0; it < 30; ++it)
  {
    h.Point(s, P, T);
    ToLocal(c, P, L);
    VecToLocal(c, T, TL);
    double f = L[2] - zloc, df = TL[2] * h.p / h.pt;
    if(fabs(df) < 1e-6) return 0.0;
    double ds = -f / df;
    s += ds;
    if(fabs(ds) < 1e-4)
    {
      ok = (s > 0);
      return s;
    }
  }
  return 0.0;
}

// transverse path length at which the helix meets the mirror sphere
double ARCTracing::SolveSphere(const Helix &h, const Cell &c, double sT0, bool &ok) const
{
  ok = false;
  double P[3], T[3], L[3], TL[3], s = sT0;
  for(int it = 0; it < 30; ++it)
  {
    h.Point(s, P, T);
    ToLocal(c, P, L);
    VecToLocal(c, T, TL);
    double d[3] = {L[0] - c.C[0], L[1] - c.C[1], L[2] - c.C[2]};
    double r = Norm(d);
    double f = r - c.R, df = Dot(d, TL) / max(r, 1e-9) * h.p / h.pt;
    if(fabs(df) < 1e-6) return 0.0;
    double ds = -f / df;
    s += ds;
    if(fabs(ds) < 1e-4)
    {
      ok = (s > 0);
      return s;
    }
  }
  return 0.0;
}

//------------------------------------------------------------------------------

const ARCTracing::Cell *ARCTracing::FindCell(const Helix &h, double &sTentry) const
{
  bool okB = false;
  double sB = SolveRadius(h, fBarrelRin, okB);
  double sE = -1.0;
  int side = 0;
  if(h.pz != 0.0)
  {
    side = h.pz > 0 ? 1 : -1;
    sE = (side * fEndcapZin - h.z0) * h.pt / h.pz;
    if(sE <= 0) sE = -1.0;
  }
  int region;
  if(okB && (sE < 0 || sB <= sE))
  {
    region = kBarrel;
    sTentry = sB;
  }
  else if(sE > 0)
  {
    region = kEndcap;
    sTentry = sE;
  }
  else
    return 0;

  double P[3], T[3], L[3];
  h.Point(sTentry, P, T);
  const Cell *best = 0;
  double bestD = 1e9;
  double phiP = atan2(P[1], P[0]);
  for(size_t i = 0; i < fCells.size(); ++i)
  {
    const Cell &c = fCells[i];
    if(c.region != region) continue;
    if(region == kBarrel)
    {
      double dphi = TVector2::Phi_mpi_pi(phiP - c.phiAxis);
      if(fabs(dphi) > 0.25 || fabs(P[2] - c.zAxis) > 320.0) continue;
      ToLocal(c, P, L);
      if(L[2] < 0.5 * c.zmin) continue;
      double d = HexDist(L[0], L[1], c.apothem * L[2] / c.zmin, c.phi0);
      if(d < bestD) { bestD = d; best = &c; }
    }
    else
    {
      if(c.side != side) continue;
      ToLocal(c, P, L);
      double d = HexDist(L[0], L[1], c.apothem, c.phi0);
      if(d < bestD) { bestD = d; best = &c; }
    }
  }
  return (bestD < 1.2) ? best : 0;
}

//------------------------------------------------------------------------------

// Snell at a plane perpendicular to the local z axis, from index n1 into n2; false on total internal reflection
bool ARCTracing::RefractZ(double *k, double n1, double n2) const
{
  double f = n1 / n2;
  double sin2sq = f * f * (k[0] * k[0] + k[1] * k[1]);
  if(sin2sq >= 1.0) return false;
  k[0] *= f;
  k[1] *= f;
  k[2] = (k[2] < 0 ? -1.0 : 1.0) * sqrt(1.0 - sin2sq);
  return true;
}

static bool RefractN(double *k, const double *nrm, double n1, double n2)
{
  double cos1 = Dot(k, nrm);
  double tang[3] = {k[0] - cos1 * nrm[0], k[1] - cos1 * nrm[1], k[2] - cos1 * nrm[2]};
  double sin1 = Norm(tang);
  double sin2 = n1 / n2 * sin1;
  if(sin2 >= 1.0) return false;
  double sgn = cos1 < 0 ? -1.0 : 1.0;
  for(int i = 0; i < 3; ++i)
    k[i] = (sin1 > 1e-12 ? tang[i] / sin1 * sin2 : 0.0) + sgn * sqrt(1.0 - sin2 * sin2) * nrm[i];
  return true;
}

// Rayleigh scattering: new direction with the (1 + cos^2 theta) distribution about the old one
void ARCTracing::RayleighScatter(double *k) const
{
  double ct;
  do
  {
    ct = gRandom->Uniform(-1.0, 1.0);
  } while(gRandom->Uniform() > 0.5 * (1.0 + ct * ct));
  double st = sqrt(1.0 - ct * ct), phi = gRandom->Uniform(0.0, TMath::TwoPi());
  double a[3] = {fabs(k[0]) < 0.9 ? 1.0 : 0.0, fabs(k[0]) < 0.9 ? 0.0 : 1.0, 0.0}, u[3], v[3];
  double ak = Dot(a, k);
  for(int i = 0; i < 3; ++i) u[i] = a[i] - ak * k[i];
  Unit(u);
  Cross(k, u, v);
  for(int i = 0; i < 3; ++i) k[i] = ct * k[i] + st * (cos(phi) * u[i] + sin(phi) * v[i]);
}

// photon inside the aerogel slab (zBot < Q_z < zTop) with direction k: absorption and Rayleigh scattering until it leaves
// returns 0 absorbed, 1 left through the top, 2 left through the bottom, 3 stuck (too many scatterings)
int ARCTracing::SlabTransport(double *Q, double *k, double zBot, double zTop, double lAbs, double lRay) const
{
  for(int it = 0; it < 6; ++it)
  {
    if(fabs(k[2]) < 1e-9) return 3;
    double zExit = k[2] > 0 ? zTop : zBot;
    double dExit = (zExit - Q[2]) / k[2];
    double dAbs = lAbs > 0 ? -lAbs * log(gRandom->Uniform()) : 1e30;
    double dRay = lRay > 0 ? -lRay * log(gRandom->Uniform()) : 1e30;
    if(dExit <= dAbs && dExit <= dRay)
    {
      for(int i = 0; i < 3; ++i) Q[i] += dExit * k[i];
      Q[2] = zExit;
      return k[2] > 0 ? 1 : 2;
    }
    if(dAbs < dRay) return 0;
    for(int i = 0; i < 3; ++i) Q[i] += dRay * k[i];
    RayleighScatter(k);
  }
  return 3;
}

// one Cherenkov photon from the local point P (track direction T) with energy E: traced to a pixel centre D on the sensor face
bool ARCTracing::TracePhoton(const Cell &c, int rad, double beta, double *P, double *T, double energy, double *D)
{
  const Radiator &r = fRad[rad];
  auto lost = [&](int why) { if(fCount) ++fLoss[rad][why]; return false; };
  double nE = 1.0 + r.indexScale * (r.n.Eval(energy) - 1.0);
  double nGas = 1.0 + fRad[kGas].indexScale * (fRad[kGas].n.Eval(energy) - 1.0);
  double nAero = 1.0 + fRad[kAerogel].indexScale * (fRad[kAerogel].n.Eval(energy) - 1.0);
  double cosT = 1.0 / (nE * beta);
  if(cosT >= 1.0) return lost(0);
  double sinT = sqrt(1.0 - cosT * cosT), phi = gRandom->Uniform(0.0, TMath::TwoPi());
  // orthonormal frame around the track direction
  double a[3] = {fabs(T[0]) < 0.9 ? 1.0 : 0.0, fabs(T[0]) < 0.9 ? 0.0 : 1.0, 0.0};
  double u[3], v[3], k[3];
  double at = Dot(a, T);
  Unit(T);
  for(int i = 0; i < 3; ++i) u[i] = a[i] - at * T[i];
  Unit(u);
  Cross(T, u, v);
  for(int i = 0; i < 3; ++i) k[i] = cosT * T[i] + sinT * (cos(phi) * u[i] + sin(phi) * v[i]);
  Unit(k);

  double half = fRad[kAerogel].thickness > 0 ? 0.5 * fRad[kAerogel].thickness : c.aeroHalf;
  double zTop = c.aeroZ + half, zBot = c.aeroZ - half;
  const Radiator &aero = fRad[kAerogel];
  double lAbs = (fAerogelAbsorption && aero.abs.Ok()) ? aero.abs.Eval(energy) : 0.0;
  double lRay = (fAerogelRayleigh && aero.ray.Ok()) ? aero.ray.Eval(energy) : 0.0;

  double Q[3] = {P[0], P[1], P[2]};
  bool toMirror = true;
  if(rad == kGas && Q[2] < zBot)
  {
    // emitted in the gas below the slab (between the cell's inner face and the aerogel): up through the slab first
    if(k[2] <= 1e-6) return lost(17);
    double d0 = (zBot - Q[2]) / k[2];
    for(int i = 0; i < 3; ++i) Q[i] += d0 * k[i];
    Q[2] = zBot;
    RefractZ(k, nGas, nAero);
    int res = SlabTransport(Q, k, zBot, zTop, lAbs, lRay);
    if(res == 0) return lost(1);
    if(res == 3) return lost(2);
    if(!RefractZ(k, nAero, nGas)) return lost(3);
    toMirror = (res == 1);
  }
  else if(rad == kAerogel)
  {
    int res = SlabTransport(Q, k, zBot, zTop, lAbs, lRay);
    if(res == 0) return lost(1);
    if(res == 3) return lost(2);
    if(!RefractZ(k, nAero, nGas)) return lost(3);
    toMirror = (res == 1); // scattered out through the bottom: straight to the sensor (a background hit)
  }
  if(toMirror)
  {
    double dC[3] = {Q[0] - c.C[0], Q[1] - c.C[1], Q[2] - c.C[2]};
    double b = Dot(k, dC), cc = Dot(dC, dC) - c.R * c.R;
    double disc = b * b - cc;
    if(disc < 0) return lost(4);
    double d = -b + sqrt(disc);
    if(d <= 0) return lost(5);
    if(fRad[kGas].abs.Ok() && gRandom->Uniform() > exp(-d / fRad[kGas].abs.Eval(energy))) return lost(6);
    double M[3] = {Q[0] + d * k[0], Q[1] + d * k[1], Q[2] + d * k[2]};
    double aM = c.region == kBarrel ? c.apothem * M[2] / c.zmin : c.apothem;
    if(M[2] < c.zmin || !InHex(M[0], M[1], aM, c.phi0)) return lost(7);
    if(fRefl.Ok() && gRandom->Uniform() > fRefl.Eval(energy)) return lost(8);
    double nrm[3] = {(M[0] - c.C[0]) / c.R, (M[1] - c.C[1]) / c.R, (M[2] - c.C[2]) / c.R};
    double kn = Dot(k, nrm);
    for(int i = 0; i < 3; ++i) k[i] -= 2.0 * kn * nrm[i];
    for(int i = 0; i < 3; ++i) Q[i] = M[i];
    // back through the aerogel slab in front of the sensor
    if(fAerogelSlabInPath)
    {
      if(k[2] >= -1e-6) return lost(9);
      double d1 = (zTop - Q[2]) / k[2];
      if(d1 < 0) return lost(10);
      for(int i = 0; i < 3; ++i) Q[i] += d1 * k[i];
      Q[2] = zTop;
      RefractZ(k, nGas, nAero);
      int res = SlabTransport(Q, k, zBot, zTop, lAbs, lRay);
      if(res == 0) return lost(11);
      if(res != 2) return lost(12);
      if(!RefractZ(k, nAero, nGas)) return lost(13);
    }
  }
  // sensor face
  double denom = Dot(k, c.sez);
  if(denom >= -1e-9) return lost(14);
  double dS[3] = {c.S[0] - Q[0], c.S[1] - Q[1], c.S[2] - Q[2]};
  double ds = Dot(dS, c.sez) / denom;
  if(ds <= 0) return lost(15);
  double H[3] = {Q[0] + ds * k[0] - c.S[0], Q[1] + ds * k[1] - c.S[1], Q[2] + ds * k[2] - c.S[2]};
  double uu = Dot(H, c.sex), vv = Dot(H, c.sey);
  int nx = int(2.0 * c.hs[0] * fSensorSizeScale / fPixelPitch + 1e-6), ny = int(2.0 * c.hs[1] * fSensorSizeScale / fPixelPitch + 1e-6);
  double hx = 0.5 * nx * fPixelPitch, hy = 0.5 * ny * fPixelPitch;
  if(fabs(uu) > hx || fabs(vv) > hy) return lost(16);
  int ix = min(max(int((uu + hx) / fPixelPitch), 0), nx - 1), iy = min(max(int((vv + hy) / fPixelPitch), 0), ny - 1);
  double uc = -hx + (ix + 0.5) * fPixelPitch, vc = -hy + (iy + 0.5) * fPixelPitch;
  for(int i = 0; i < 3; ++i) D[i] = c.S[i] + uc * c.sex[i] + vc * c.sey[i];
  return true;
}

//------------------------------------------------------------------------------

// inverse ray tracing of one pixel hit D: the mirror point M where the reflection law holds for the emission point E,
// then the photon direction at E; for the aerogel the direction is refracted back into the aerogel (nominal indices)
double ARCTracing::RecoAngle(const Cell &c, const double *E, const double *T, const double *D, int rad) const
{
  double e1[3] = {E[0] - c.C[0], E[1] - c.C[1], E[2] - c.C[2]};
  Unit(e1);
  double dc[3] = {D[0] - c.C[0], D[1] - c.C[1], D[2] - c.C[2]};
  double p1 = Dot(dc, e1);
  double w[3] = {dc[0] - p1 * e1[0], dc[1] - p1 * e1[1], dc[2] - p1 * e1[2]};
  if(Norm(w) < 1e-9) return -1.0;
  double e2[3] = {w[0], w[1], w[2]};
  Unit(e2);
  double betaD = atan2(Dot(dc, e2), p1);
  auto f = [&](double al, double *M) {
    double tv[3];
    for(int i = 0; i < 3; ++i)
    {
      M[i] = c.C[i] + c.R * (cos(al) * e1[i] + sin(al) * e2[i]);
      tv[i] = -sin(al) * e1[i] + cos(al) * e2[i];
    }
    double a[3] = {E[0] - M[0], E[1] - M[1], E[2] - M[2]}, b[3] = {D[0] - M[0], D[1] - M[1], D[2] - M[2]};
    Unit(a);
    Unit(b);
    double s[3] = {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    return Dot(s, tv);
  };
  double M[3], lo = 0.0, hi = betaD, flo = f(lo, M);
  for(int it = 0; it < 45; ++it)
  {
    double mid = 0.5 * (lo + hi), fm = f(mid, M);
    if((fm > 0) == (flo > 0)) { lo = mid; flo = fm; }
    else hi = mid;
  }
  f(0.5 * (lo + hi), M);
  double k[3] = {M[0] - E[0], M[1] - E[1], M[2] - E[2]};
  Unit(k);
  double Tn[3] = {T[0], T[1], T[2]};
  Unit(Tn);
  if(rad == kAerogel)
  {
    // the analysis refracts with the local radial normal at E in the barrel, the beam axis in the endcaps
    double nrm[3];
    if(c.region == kBarrel)
    {
      double Ew[3] = {c.O[0] + E[0] * c.ex[0] + E[1] * c.ey[0] + E[2] * c.ez[0], c.O[1] + E[0] * c.ex[1] + E[1] * c.ey[1] + E[2] * c.ez[1], 0.0};
      double nw[3] = {Ew[0], Ew[1], 0.0};
      Unit(nw);
      VecToLocal(c, nw, nrm);
    }
    else
    {
      double nw[3] = {0.0, 0.0, c.side > 0 ? 1.0 : -1.0};
      VecToLocal(c, nw, nrm);
    }
    if(!RefractN(k, nrm, fRecoIndex[kGas], fRecoIndex[kAerogel])) return -1.0;
  }
  return acos(min(max(Dot(k, Tn), -1.0), 1.0));
}

// mean of the angles within the window around center, after two passes of 3-sigma (MAD) clipping around the median
void ARCTracing::RobustMean(const vector<double> &th, double center, double window, double &mean, double &rms, int &n) const
{
  vector<double> sel;
  for(size_t i = 0; i < th.size(); ++i)
    if(th[i] >= 0 && fabs(th[i] - center) < window) sel.push_back(th[i]);
  n = 0;
  mean = rms = -1.0;
  if(sel.empty()) return;
  for(int pass = 0; pass < 2 && !sel.empty(); ++pass)
  {
    vector<double> s = sel;
    sort(s.begin(), s.end());
    double med = s[s.size() / 2];
    double sig = window;
    if(s.size() > 3)
    {
      vector<double> dev;
      for(size_t i = 0; i < s.size(); ++i) dev.push_back(fabs(s[i] - med));
      sort(dev.begin(), dev.end());
      sig = 1.4826 * dev[dev.size() / 2];
    }
    double cut = max(3.0 * sig, 0.002);
    vector<double> keep;
    for(size_t i = 0; i < sel.size(); ++i)
      if(fabs(sel[i] - med) < cut) keep.push_back(sel[i]);
    sel = keep;
  }
  if(sel.empty()) return;
  double s1 = 0, s2 = 0;
  for(size_t i = 0; i < sel.size(); ++i) { s1 += sel[i]; s2 += sel[i] * sel[i]; }
  n = sel.size();
  mean = s1 / n;
  rms = sqrt(max(s2 / n - mean * mean, 0.0));
}

//------------------------------------------------------------------------------

// per-track state kept between the two passes of Process(): photons are generated for every track first and pooled per
// cell, then each track is reconstructed against all hits of its cell (SharedCells): overlapping rings in jets confuse each other
struct ARCTracingWork
{
  Candidate *candidate, *mother, *particle;
  double pTrue, mass, charge, beta, recoP;
  ARCTracing::Helix hTrue, hReco;
  const ARCTracing::Cell *cell;
  double sEntry, sM, sLim[2][2], sLow[2], L[2];
  bool pathOk;
  int nGen[2], nDet[2];
  Long64_t index;
  std::vector<double> hits; // the track's own detected photons (3 numbers each)
};

void ARCTracing::Process()
{
  Candidate *candidate, *mother, *particle;
  const double masses[3] = {kMassPi, kMassK, kMassP};
  const int pdg[3] = {211, 321, 2212};

  vector<ARCTracingWork> works;
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));
    const TLorentzVector &trueMomentum = particle->Momentum;
    const TLorentzVector &recoMomentum = candidate->Momentum;
    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());
    candidate->CherenkovAngleGas = candidate->CherenkovAngleAerogel = -1.0;
    candidate->CherenkovRecoAngleGas = candidate->CherenkovRecoAngleAerogel = -1.0;
    candidate->SigmaThetaCGas = candidate->SigmaThetaCAerogel = -1.0;
    candidate->NpeGas = candidate->NpeAerogel = 0.0;
    candidate->ARCchi2Pi = candidate->ARCchi2K = candidate->ARCchi2P = -1.0;
    candidate->ARCpid = 0;
    candidate->ARCprobPi = candidate->ARCprobK = candidate->ARCprobP = -1.0;
    candidate->ARCmass = -1.0;
    ++fNtracks;

    double pTrue = trueMomentum.P(), mass = trueMomentum.M(), charge = candidate->Charge;
    double beta = pTrue > 0 ? pTrue / sqrt(pTrue * pTrue + mass * mass) : 0.0;
    double nref[2] = {fRecoIndex[kGas], fRecoIndex[kAerogel]};
    for(int i = 0; i < 2; ++i)
    {
      double cs = beta > 0 ? 1.0 / (nref[i] * beta) : 2.0;
      double ang = cs < 1.0 ? acos(cs) : -1.0;
      if(i == kGas) candidate->CherenkovAngleGas = ang;
      else candidate->CherenkovAngleAerogel = ang;
    }

    // the helix of the particle (photon generation) and of the track used by the reconstruction
    auto makeHelix = [&](const TLorentzVector &mom, Helix &h) {
      h.x0 = particle->Position.X();
      h.y0 = particle->Position.Y();
      h.z0 = particle->Position.Z();
      h.pt = mom.Pt();
      h.pz = mom.Pz();
      h.p = mom.P();
      h.phi0 = mom.Phi();
      h.q = charge;
      h.straight = (fabs(fBz) < 1e-6 || charge == 0.0);
      h.rc = h.straight ? 1e12 : h.pt / (0.3 * fabs(fBz)) * 1000.0;
    };
    Helix hTrue, hReco;
    makeHelix(trueMomentum, hTrue);
    makeHelix(fUseTrueTrack ? trueMomentum : recoMomentum, hReco);
    bool usable = (pTrue > 0 && hTrue.pt > 1e-3 && (hTrue.straight || 2.0 * hTrue.rc > 300.0));

    double sEntry = 0.0;
    const Cell *cell = usable ? FindCell(hTrue, sEntry) : 0;
    if(cell)
    {
      const Cell &c = *cell;
      ++fNcovered;
      double half = fRad[kAerogel].thickness > 0 ? 0.5 * fRad[kAerogel].thickness : c.aeroHalf;
      bool ok0, ok1, ok2, ok3;
      double sA0 = SolvePlane(hTrue, c, c.aeroZ - half, sEntry, ok0);
      double sA1 = SolvePlane(hTrue, c, c.aeroZ + half, sEntry, ok1);
      double sApex = SolvePlane(hTrue, c, c.C[2] + c.R, sEntry, ok2);
      double sM = ok2 ? SolveSphere(hTrue, c, sApex, ok3) : 0.0;
      double sLim[2][2] = {{sA1, sM}, {sA0, sA1}};
      double sLow[2] = {sEntry, sA0}; // gas below the aerogel slab (the sensor and cooling plates are ignored)
      bool pathOk = ok0 && ok1 && ok2 && ok3 && sA1 > sA0 && sM > sA1 && sA0 >= sEntry;

      vector<double> hitsLocal; // 3 numbers per detected photon
      double L[2] = {0.0, 0.0};
      int nGen[2] = {0, 0}, nDet[2] = {0, 0};
      if(pathOk)
      {
        vector<double> cdf;
        fCount = kTRUE;
        for(int rad = 0; rad < 2; ++rad)
        {
          double lowFrac = 0.0;
          L[rad] = (sLim[rad][1] - sLim[rad][0]) * hTrue.p / hTrue.pt;
          if(rad == kGas)
          {
            double Llow = (sLow[1] - sLow[0]) * hTrue.p / hTrue.pt;
            lowFrac = Llow / (L[rad] + Llow);
            L[rad] += Llow;
          }
          double lambda = Yield(fRad[rad], beta, &cdf) * L[rad];
          nGen[rad] = lambda > 0 ? gRandom->Poisson(lambda) : 0;
          for(int ip = 0; ip < nGen[rad]; ++ip)
          {
            double sT, Pw[3], Tw[3], P[3], T[3], D[3];
            if(fDebugMidEmission)
            {
              // study: all photons from the point the reconstruction assumes (barrel: world radius RecoRadius), else the path midpoint
              bool okm = false;
              sT = c.region == kBarrel ? SolveRadius(hTrue, fRecoRadius[rad], okm) : 0.0;
              if(!okm) sT = 0.5 * (sLim[rad][0] + sLim[rad][1]);
            }
            else if(rad == kGas && gRandom->Uniform() < lowFrac) sT = gRandom->Uniform(sLow[0], sLow[1]);
            else sT = gRandom->Uniform(sLim[rad][0], sLim[rad][1]);
            hTrue.Point(sT, Pw, Tw);
            ToLocal(c, Pw, P);
            VecToLocal(c, Tw, T);
            if(TracePhoton(c, rad, beta, P, T, fDebugMonoEnergy > 0 ? fDebugMonoEnergy : SampleEnergy(fRad[rad], cdf), D))
            {
              ++nDet[rad];
              hitsLocal.insert(hitsLocal.end(), D, D + 3);
            }
          }
          fNgen[rad] += nGen[rad];
          fNdet[rad] += nDet[rad];
        }
        fCount = kFALSE;
        if(fDarkCounts > 0)
        {
          int nx = int(2.0 * c.hs[0] * fSensorSizeScale / fPixelPitch + 1e-6), ny = int(2.0 * c.hs[1] * fSensorSizeScale / fPixelPitch + 1e-6);
          int nd = gRandom->Poisson(fDarkCounts);
          for(int id = 0; id < nd; ++id)
          {
            double uc = (gRandom->Integer(nx) + 0.5 - 0.5 * nx) * fPixelPitch, vc = (gRandom->Integer(ny) + 0.5 - 0.5 * ny) * fPixelPitch;
            for(int i = 0; i < 3; ++i) hitsLocal.push_back(c.S[i] + uc * c.sex[i] + vc * c.sey[i]);
          }
        }
      }


      ARCTracingWork w;
      w.candidate = candidate; w.mother = mother; w.particle = particle;
      w.pTrue = pTrue; w.mass = mass; w.charge = charge; w.beta = beta; w.recoP = recoMomentum.P();
      w.hTrue = hTrue; w.hReco = hReco; w.cell = cell; w.sEntry = sEntry; w.sM = sM;
      for(int i = 0; i < 2; ++i) { w.sLim[i][0] = sLim[i][0]; w.sLim[i][1] = sLim[i][1]; w.sLow[i] = sLow[i]; w.L[i] = L[i]; w.nGen[i] = nGen[i]; w.nDet[i] = nDet[i]; }
      w.pathOk = pathOk; w.index = fNtracks; w.hits = hitsLocal;
      works.push_back(w);
    }
    else
    {
      ARCTracingWork w;
      w.candidate = candidate; w.mother = mother; w.particle = particle; w.cell = 0;
      works.push_back(w);
    }
  }

  // photons of all tracks of the event pooled per cell
  map<const Cell *, vector<double> > pooled;
  if(fSharedCells)
    for(size_t i = 0; i < works.size(); ++i)
      if(works[i].cell) pooled[works[i].cell].insert(pooled[works[i].cell].end(), works[i].hits.begin(), works[i].hits.end());

  // pass 2: reconstruction and particle identification, track by track, against the pooled hits of its cell
  for(size_t iw = 0; iw < works.size(); ++iw)
  {
    ARCTracingWork &w = works[iw];
    candidate = w.candidate; mother = w.mother; particle = w.particle;
    if(w.cell)
    {
      const Cell &c = *w.cell;
      const Helix &hTrue = w.hTrue, &hReco = w.hReco;
      double pTrue = w.pTrue, charge = w.charge, beta = w.beta, sEntry = w.sEntry, sM = w.sM;
      double nref[2] = {fRecoIndex[kGas], fRecoIndex[kAerogel]};
      double (&sLim)[2][2] = w.sLim; double (&sLow)[2] = w.sLow; double (&L)[2] = w.L;
      int (&nGen)[2] = w.nGen; int (&nDet)[2] = w.nDet;
      bool pathOk = w.pathOk;
      const vector<double> &hitsLocal = fSharedCells ? pooled[w.cell] : w.hits;
      (void)beta; (void)sEntry; (void)sM;

      // reconstruction: emission point at mid-radiator on the (true or smeared) track
      double E[2][3], Tt[2][3];
      bool recoOk[2] = {false, false};
      for(int rad = 0; rad < 2; ++rad)
      {
        bool ok = false;
        double sT = 0.0;
        if(c.region == kBarrel)
          sT = SolveRadius(hReco, fRecoRadius[rad], ok);
        else if(hReco.pz != 0.0)
        {
          sT = (c.side * fRecoZ[rad] - hReco.z0) * hReco.pt / hReco.pz;
          ok = sT > 0;
        }
        if(!ok) continue;
        double Pw[3], Tw[3];
        hReco.Point(sT, Pw, Tw);
        ToLocal(c, Pw, E[rad]);
        VecToLocal(c, Tw, Tt[rad]);
        recoOk[rad] = true;
      }
      vector<double> theta[2];
      size_t nHits = hitsLocal.size() / 3;
      for(int rad = 0; rad < 2; ++rad)
      {
        if(!recoOk[rad]) continue;
        for(size_t ih = 0; ih < nHits; ++ih)
          theta[rad].push_back(RecoAngle(c, E[rad], Tt[rad], &hitsLocal[3 * ih], rad));
      }

      if(fDebug)
      {
        double Pw[3], Tw[3], Pl[3];
        hTrue.Point(sEntry, Pw, Tw);
        ToLocal(c, Pw, Pl);
        double Pm[3], Tm[3], Px[3], Tx[3];
        hTrue.Point(0.5 * (sEntry + sM), Pm, Tm);
        hTrue.Point(sM, Px, Tx);
        *fDebug << "T " << w.index << " " << c.region << " " << (&c - &fCells[0]) << " " << particle->PID << " " << pTrue << " " << Pl[0] << " " << Pl[1]
                << " " << candidate->CherenkovAngleGas << " " << candidate->CherenkovAngleAerogel << " " << nGen[kGas] << " " << nGen[kAerogel]
                << " " << hypot(Pm[0], Pm[1]) << " " << hypot(Px[0], Px[1]) << "\n";
        for(size_t ih = 0; ih < nHits; ++ih)
        {
          const double *D = &hitsLocal[3 * ih];
          double dS[3] = {D[0] - c.S[0], D[1] - c.S[1], D[2] - c.S[2]};
          *fDebug << "H " << Dot(dS, c.sex) << " " << Dot(dS, c.sey) << " " << (recoOk[kGas] ? theta[kGas][ih] : -1.0) << " "
                  << (recoOk[kAerogel] ? theta[kAerogel][ih] : -1.0) << "\n";
        }
      }

      // acceptance of the optics for the hypotheses: from this track's photons, or from test photons at the hypothesis beta
      auto acceptance = [&](int rad, double betaH) {
        if(nGen[rad] >= 10) return double(nDet[rad]) / nGen[rad];
        if(!pathOk) return 0.0;
        vector<double> cdf;
        if(Yield(fRad[rad], betaH, &cdf) <= 0) return 0.0;
        int nt = 0, nd = 0;
        double Llow = rad == kGas ? (sLow[1] - sLow[0]) : 0.0, lowFrac = Llow / (Llow + sLim[rad][1] - sLim[rad][0]);
        for(; nt < 20; ++nt)
        {
          double sT = gRandom->Uniform() < lowFrac ? gRandom->Uniform(sLow[0], sLow[1]) : gRandom->Uniform(sLim[rad][0], sLim[rad][1]), Pw[3], Tw[3], P[3], T[3], D[3];
          hTrue.Point(sT, Pw, Tw);
          ToLocal(c, Pw, P);
          VecToLocal(c, Tw, T);
          if(TracePhoton(c, rad, betaH, P, T, SampleEnergy(fRad[rad], cdf), D)) ++nd;
        }
        return double(nd) / nt;
      };

      // the measured rings: the window with most photons among the hypotheses (and the true one, for validation)
      double pHyp = fUseTrueMomentumForPID ? pTrue : w.recoP;
      if(pHyp <= 0) pHyp = pTrue;
      struct Ring { double mean, rms, sigma, thetaH, nExp; int n; bool expected, found; };
      Ring ring[3][2];
      Ring meas[2];
      for(int rad = 0; rad < 2; ++rad) { meas[rad].n = 0; meas[rad].found = false; meas[rad].mean = meas[rad].sigma = -1.0; }
      for(int h = 0; h < 3; ++h)
      {
        double betaH = pHyp / sqrt(pHyp * pHyp + masses[h] * masses[h]);
        for(int rad = 0; rad < 2; ++rad)
        {
          Ring &r = ring[h][rad];
          double cs = 1.0 / (nref[rad] * betaH);
          r.thetaH = cs < 1.0 ? acos(cs) : -1.0;
          r.nExp = r.thetaH > 0 && pathOk ? Yield(fRad[rad], betaH, 0) * L[rad] * acceptance(rad, betaH) : 0.0;
          r.expected = r.thetaH > 0 && r.nExp >= fMinPe;
          r.n = 0;
          r.mean = r.rms = r.sigma = -1.0;
          r.found = false;
          if(r.thetaH > 0 && recoOk[rad])
          {
            RobustMean(theta[rad], r.thetaH, fWindow, r.mean, r.rms, r.n);
            r.found = r.n >= fMinPe;
            if(r.found)
            {
              double s1sq = max(r.rms * r.rms, fSigmaSingleFloor[rad] * fSigmaSingleFloor[rad]);
              r.sigma = sqrt(s1sq / r.n + fSigmaCorr[rad] * fSigmaCorr[rad] + fSigmaTrack * fSigmaTrack);
              if(r.n > meas[rad].n) meas[rad] = r;
            }
          }
        }
      }
      bool anyRing = meas[kGas].found || meas[kAerogel].found;
      if(anyRing)
      {
        double chi2[3];
        for(int h = 0; h < 3; ++h)
        {
          double c2 = 0.0;
          for(int rad = 0; rad < 2; ++rad)
          {
            const Ring &r = ring[h][rad];
            if(r.expected && r.found)
              c2 += pow((r.mean - r.thetaH) / r.sigma, 2);
            else if(r.expected && !r.found)
              c2 += 2.0 * r.nExp; // no ring where nExp photons were expected: -2 ln P(0 | nExp)
            else if(!r.expected && meas[rad].found)
              c2 += pow(meas[rad].mean / meas[rad].sigma, 2); // a ring measured where this hypothesis predicts none
          }
          chi2[h] = c2;
        }
        int best = 0;
        for(int h = 1; h < 3; ++h)
          if(chi2[h] < chi2[best]) best = h;
        candidate->ARCchi2Pi = chi2[0];
        candidate->ARCchi2K = chi2[1];
        candidate->ARCchi2P = chi2[2];
        candidate->ARCpid = (charge < 0.0 ? -1 : 1) * pdg[best];
        candidate->ARCmass = masses[best];
        double w[3], wsum = 0.0;
        for(int h = 0; h < 3; ++h)
        {
          w[h] = exp(-0.5 * (chi2[h] - chi2[best]));
          wsum += w[h];
        }
        candidate->ARCprobPi = w[0] / wsum;
        candidate->ARCprobK = w[1] / wsum;
        candidate->ARCprobP = w[2] / wsum;

        for(int rad = 0; rad < 2; ++rad)
        {
          Ring out = ring[best][rad].found ? ring[best][rad] : meas[rad];
          if(fTrueHypothesisAngle)
          {
            // validation: the window around the true angle, as the full-simulation analysis
            double thetaTrue = rad == kGas ? candidate->CherenkovAngleGas : candidate->CherenkovAngleAerogel;
            out.found = false;
            if(thetaTrue > 0 && recoOk[rad])
            {
              RobustMean(theta[rad], thetaTrue, fWindow, out.mean, out.rms, out.n);
              out.found = out.n >= fMinPe;
              double s1sq = max(out.rms * out.rms, fSigmaSingleFloor[rad] * fSigmaSingleFloor[rad]);
              out.sigma = out.found ? sqrt(s1sq / out.n + fSigmaCorr[rad] * fSigmaCorr[rad] + fSigmaTrack * fSigmaTrack) : -1.0;
            }
          }
          if(!out.found) continue;
          if(rad == kGas)
          {
            candidate->CherenkovRecoAngleGas = out.mean;
            candidate->SigmaThetaCGas = out.sigma;
            candidate->NpeGas = out.n;
          }
          else
          {
            candidate->CherenkovRecoAngleAerogel = out.mean;
            candidate->SigmaThetaCAerogel = out.sigma;
            candidate->NpeAerogel = out.n;
          }
        }
      }
    }

    candidate->AddCandidate(mother);
    fOutputArray->Add(candidate);
  }
}
