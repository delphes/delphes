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

#ifndef ARCTracing_h
#define ARCTracing_h

/** \class ARCTracing
 *
 *  Ray-tracing fast simulation of the ARC (Array of RICH Cells) for the ALFA detector.
 *
 *  Unlike a parametrised response fitted to the full simulation of one geometry, this module models the
 *  detector: for every charged track it
 *   1. finds the ARC cell the helix enters (barrel or endcap) and works in that cell's optical frame, read from a table
 *      dumped from the k4geo geometry (mirror centre and radius, SiPM array position, tilt and size, aerogel position);
 *   2. generates the Cherenkov photons along the helix in the aerogel slab and in the gas with the Frank-Tamm yield,
 *      the SiPM PDE and the mirror reflectivity, with the photon energy sampled from the detected spectrum, so that the
 *      dispersion n(E) of each radiator (chromatic error), the thresholds and the yield vs beta come out by themselves;
 *   3. traces every photon: refraction out of the aerogel, absorption and Rayleigh scattering in it, reflection on the
 *      spherical mirror (lost if outside the cell), the way back through the aerogel, the SiPM array (lost if outside),
 *      and the pixel: a ring that falls off the array, a ring larger than the sensor, a track that bends, are all modelled;
 *   4. reconstructs the Cherenkov angle of every pixel hit exactly as the full-simulation analysis does (inverse ray
 *      tracing from the pixel centre on the sensor face, emission assumed at mid-radiator on the track), and
 *   5. computes the per-track angle per radiator, its resolution and the pi/K/p probabilities (chi2 over the two radiators).
 *
 *  Gas, aerogel (index, dispersion, thickness, scattering), magnetic field, sensor size, pixel pitch, mirror curvature and
 *  photon-detection efficiency are parameters; another ARC geometry is another cell table (dump_arc_cells.py).
 *
 *  \author ALFA ARC study
 */

#include "classes/DelphesModule.h"

#include <fstream>
#include <string>
#include <vector>

class TIterator;
class TObjArray;
class Candidate;

class ARCTracing: public DelphesModule
{
public:
  ARCTracing();
  ~ARCTracing();

  void Init();
  void Process();
  void Finish();

  // types (public: the ROOT dictionary of the module sees them)
  enum { kGas = 0, kAerogel = 1 };
  enum { kBarrel = 0, kEndcap = 1 };

  // tabulated optical property vs photon energy [eV], linear interpolation, no extrapolation
  struct Table
  {
    std::vector<double> e, v;
    double Eval(double x) const;
    bool Ok() const { return e.size() >= 2; }
  };

  // one placed cell, with its optics in the cell frame (mm): z along the cell axis, x along the mirror offset
  struct Cell
  {
    int region, side;
    double O[3], ex[3], ey[3], ez[3]; // frame in world coordinates
    double C[3], R; // mirror centre, reflecting radius
    double S[3], sex[3], sey[3], sez[3], hs[3]; // sensor face centre, axes (sez: normal pointing to the mirror), half sizes
    double aeroZ, aeroHalf; // aerogel slab
    double phi0, zmin, zmax, apothem; // hexagonal cell: vertex direction, axial extent, apothem (at zmin in the barrel)
    double phiAxis, zAxis; // for the barrel pre-selection
  };

  // helix from the production vertex, parametrised by the transverse path length sT (mm)
  struct Helix
  {
    double x0, y0, z0, phi0, rc, q, pt, pz, p;
    bool straight;
    void Point(double sT, double *P, double *T) const;
  };

  struct Radiator
  {
    Table n, abs, ray;
    double indexScale, thickness;
    std::vector<double> gridE, gridN, gridW; // energy grid for the yield integral and the energy sampling
  };

  struct Photon
  {
    double thetaGas, thetaAero; // reconstructed angle under each radiator hypothesis [rad]
  };

private:
  void ReadTables(const char *file);
  void ReadCells(const char *file);
  const Table *FindTable(const std::string &name) const;
  void BuildGrid(Radiator &r, const Table &pde, const Table &refl);
  double Yield(const Radiator &r, double beta, std::vector<double> *cdf) const; // detected photons per mm of path
  double SampleEnergy(const Radiator &r, const std::vector<double> &cdf) const;

  void ToLocal(const Cell &c, const double *P, double *L) const;
  void VecToLocal(const Cell &c, const double *V, double *L) const;
  bool InHex(double x, double y, double a, double phi0) const;
  double HexDist(double x, double y, double a, double phi0) const;
  const Cell *FindCell(const Helix &h, double &sTentry) const;
  double SolvePlane(const Helix &h, const Cell &c, double zloc, double sT0, bool &ok) const;
  double SolveSphere(const Helix &h, const Cell &c, double sT0, bool &ok) const;
  double SolveRadius(const Helix &h, double r, bool &ok) const;

  bool TracePhoton(const Cell &c, int rad, double beta, double *P, double *T, double energy, double *D);
  bool RefractZ(double *k, double n1, double n2) const;
  int SlabTransport(double *Q, double *k, double zBot, double zTop, double lAbs, double lRay) const;
  void RayleighScatter(double *k) const;
  double RecoAngle(const Cell &c, const double *E, const double *T, const double *D, int rad) const;
  void RobustMean(const std::vector<double> &th, double center, double window, double &mean, double &rms, int &n) const;

  Double_t fBz;
  Double_t fPixelPitch, fSensorSizeScale, fSensorShift, fMirrorRadiusScale;
  Double_t fYieldScale, fPDEScale, fDarkCounts;
  Double_t fRecoEnergy, fRecoRadius[2], fRecoZ[2], fWindow, fMinPe;
  Double_t fSigmaCorr[2], fSigmaSingleFloor[2], fSigmaTrack;
  Double_t fDebugMonoEnergy;
  Bool_t fDebugMidEmission;
  Bool_t fSharedCells;
  Bool_t fAerogelSlabInPath, fAerogelRayleigh, fAerogelAbsorption, fUseTrueTrack, fUseTrueMomentumForPID, fTrueHypothesisAngle;

  std::vector<Table> fTables; //!
  std::vector<std::string> fTableNames; //!
  Table fPDE, fRefl; //!
  Radiator fRad[2]; //!
  double fRecoIndex[2];
  std::vector<Cell> fCells; //!
  double fBarrelRin, fEndcapZin;

  // statistics
  Long64_t fNtracks, fNcovered, fNgen[2], fNdet[2];
  Long64_t fLoss[2][20]; //! photons lost per radiator and reason (see Finish)
  Bool_t fCount; //! count the losses (off for the hypothesis test photons)
  std::string fDebugFile;
  std::ofstream *fDebug; //! per-track / per-photon dump for validation studies

  TIterator *fItInputArray; //!
  const TObjArray *fInputArray; //!
  TObjArray *fOutputArray; //!

  ClassDef(ARCTracing, 1)
};

#endif
