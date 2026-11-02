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
#include <iomanip>
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
  while(getline(in, line))
  {
    if(line.empty() || line[0] == '#') continue;
    istringstream ss(line);
    string region, name;
    ss >> region >> name;
    double raw[44];
    for(int i = 0; i < 44; ++i) ss >> raw[i];
    if(!ss) throw runtime_error("ARCTracing: malformed cell line: " + line);
    AddCell(region == "barrel" ? kBarrel : kEndcap, raw);
  }
  if(fCells.empty()) throw runtime_error("ARCTracing: no cells read");
}

// one cell from the 44 numbers of a cell-table line (after region and name): O, ex, ey, ez (world), mirror centre (cell
// frame), Rmin, Rmax, sensor centre, its axes and half sizes, aerogel centre and normal, aerogel half thickness and radius,
// hexagon orientation [deg], cell z range, apothem
void ARCTracing::AddCell(int region, const double *raw)
{
  Cell c;
  c.region = region;
  fCellRaw.push_back(vector<double>(raw, raw + 44));
  int k = 0;
  for(int i = 0; i < 3; ++i) c.O[i] = raw[k++];
  for(int i = 0; i < 3; ++i) c.ex[i] = raw[k++];
  for(int i = 0; i < 3; ++i) c.ey[i] = raw[k++];
  for(int i = 0; i < 3; ++i) c.ez[i] = raw[k++];
  for(int i = 0; i < 3; ++i) c.C[i] = raw[k++];
  c.R = raw[k++];
  k++; // Rmax
  for(int i = 0; i < 3; ++i) c.S[i] = raw[k++];
  for(int i = 0; i < 3; ++i) c.sex[i] = raw[k++];
  for(int i = 0; i < 3; ++i) c.sey[i] = raw[k++];
  for(int i = 0; i < 3; ++i) c.sez[i] = raw[k++];
  for(int i = 0; i < 3; ++i) c.hs[i] = raw[k++];
  c.aeroZ = raw[k + 2];
  k += 6; // aerogel centre and normal
  c.aeroHalf = raw[k++];
  k++; // aerogel radius
  c.phi0 = raw[k++];
  c.zmin = raw[k++];
  c.zmax = raw[k++];
  c.apothem = raw[k++];
  if(fCells.empty())
  {
    fBarrelRin = 1e9;
    fEndcapZin = 1e9;
  }
  // frames given with a finite number of digits: make them exactly orthonormal (the Cherenkov angle of the gas, 50 mrad, is
  // sensitive to 1e-4 in a direction's norm through cot(theta_c))
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

//------------------------------------------------------------------------------

namespace
{
// rigid transform x -> R x + t (the DD4hep placements of the k4geo ARC driver)
struct ArcTf
{
  double R[3][3], t[3];
};
ArcTf ArcIdentity()
{
  ArcTf a;
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j) a.R[i][j] = (i == j);
    a.t[i] = 0.0;
  }
  return a;
}
ArcTf ArcRot(int axis, double ang)
{
  ArcTf a = ArcIdentity();
  double c = cos(ang), s = sin(ang);
  int i = (axis + 1) % 3, j = (axis + 2) % 3;
  a.R[i][i] = c;
  a.R[i][j] = -s;
  a.R[j][i] = s;
  a.R[j][j] = c;
  return a;
}
ArcTf ArcShift(double x, double y, double z)
{
  ArcTf a = ArcIdentity();
  a.t[0] = x;
  a.t[1] = y;
  a.t[2] = z;
  return a;
}
ArcTf operator*(const ArcTf &a, const ArcTf &b)
{
  ArcTf c;
  for(int i = 0; i < 3; ++i)
  {
    c.t[i] = a.t[i];
    for(int j = 0; j < 3; ++j)
    {
      c.R[i][j] = 0.0;
      for(int k = 0; k < 3; ++k) c.R[i][j] += a.R[i][k] * b.R[k][j];
      c.t[i] += a.R[i][j] * b.t[j];
    }
  }
  return c;
}
void ArcPoint(const ArcTf &a, const double *p, double *q)
{
  for(int i = 0; i < 3; ++i) q[i] = a.t[i] + a.R[i][0] * p[0] + a.R[i][1] * p[1] + a.R[i][2] * p[2];
}
void ArcVector(const ArcTf &a, const double *v, double *w)
{
  for(int i = 0; i < 3; ++i) w[i] = a.R[i][0] * v[0] + a.R[i][1] * v[1] + a.R[i][2] * v[2];
}
} // namespace

// The ARC cells from parameters, as the k4geo ARC_geo_o1_v01 driver builds them: hexagonal pyramid cells repeated in phi in the
// barrel, hexagonal prisms on a grid of 60-degree sectors (with mirrored copies) in the two endcaps; each cell type has its own
// mirror (radius, centre) and sensor (offset, tilt) from the ARC optics optimisation. Lengths in mm, angles in rad (phi step
// in deg). The cells equal those of the DD4hep geometry of the k4geo driver (checked to 1e-6 mm for all 1203 cells).
void ARCTracing::BuildCells()
{
  const double barrelRin = GetDouble("BarrelInnerRadius", 1900.0), barrelRout = GetDouble("BarrelOuterRadius", 2100.0);
  const double endcapLength = GetDouble("EndcapLength", 200.0), endcapZ = GetDouble("EndcapZ", 2100.0);
  const double wall = GetDouble("VesselWallThickness", 10.0), mirrorThick = GetDouble("MirrorThickness", 2.0);
  const double sensorThick = GetDouble("SensorThickness", 2.0), sensorOffset = GetDouble("SensorZOffset", 5.0);
  const double aerogelThick = GetDouble("CellAerogelThickness", 10.0);
  const double mirrorShrink = GetDouble("BarrelMirrorZShrink", 6.0), mirrorOrigin = GetDouble("MirrorZOrigin", 370.0);
  const double side = GetDouble("HexagonSide", 148.15), apo = side * cos(TMath::Pi() / 6.0);
  const double phiStep = GetDouble("BarrelPhiStep", 13.333) * TMath::DegToRad();
  const int nPhiBarrel = GetInt("BarrelSectors", 27), nPhiEndcap = GetInt("EndcapSectors", 6);
  const double endcapPhiStep = GetDouble("EndcapPhiStep", 60.0) * TMath::DegToRad();
  const int maxCell = GetInt("BarrelMaxCell", 17);
  const double sensorB = GetDouble("BarrelSensorSize", 80.0), sensorE = GetDouble("EndcapSensorSize", 80.0);
  const double halfPi = 0.5 * TMath::Pi();

  // the sensor depth and the aerogel offset (the formula of the optimisation program, as in the driver)
  auto sensorDepth = [&](double origin, double tilt, double sideX, double &aero) {
    double z = origin, a = fabs(tilt), h = 0.5 * sideX;
    if(a > asin(sensorOffset / h)) z += sin(a) * h - sensorOffset;
    aero = sin(a) * sideX + aerogelThick + sensorThick;
    return z;
  };
  // one cell in the cell-table convention: frame with z = cell axis, x = mirror-centre offset
  auto addRow = [&](int region, const ArcTf &cellTf, const double *Cloc, double R, const ArcTf &sensTf, double sideX, double sideY,
                    const double *Aloc, const double *nloc, const double *hexloc, double zmin, double zmax) {
    double O[3], ez[3], C[3], Sl[3], S[3], A[3], nA[3], hv[3], zero[3] = {0, 0, 0}, z1[3] = {0, 0, 1};
    ArcPoint(cellTf, zero, O);
    ArcVector(cellTf, z1, ez);
    ArcPoint(cellTf, Cloc, C);
    ArcPoint(sensTf, zero, Sl);
    ArcPoint(cellTf, Sl, S);
    double v[3] = {C[0] - O[0], C[1] - O[1], C[2] - O[2]}, pv = Dot(v, ez);
    for(int i = 0; i < 3; ++i) v[i] -= pv * ez[i];
    if(Norm(v) < 0.5)
    {
      for(int i = 0; i < 3; ++i) v[i] = S[i] - O[i];
      pv = Dot(v, ez);
      for(int i = 0; i < 3; ++i) v[i] -= pv * ez[i];
    }
    Unit(v);
    double ex[3] = {v[0], v[1], v[2]}, ey[3];
    Cross(ez, ex, ey);
    auto loc = [&](const double *P, double *out) {
      double d[3] = {P[0] - O[0], P[1] - O[1], P[2] - O[2]};
      out[0] = Dot(d, ex);
      out[1] = Dot(d, ey);
      out[2] = Dot(d, ez);
    };
    auto lvec = [&](const double *w, double *out) {
      out[0] = Dot(w, ex);
      out[1] = Dot(w, ey);
      out[2] = Dot(w, ez);
    };
    ArcTf sw = cellTf * sensTf;
    double raw[44], tmp[3], ax[3];
    int k = 0;
    for(int i = 0; i < 3; ++i) raw[k++] = O[i];
    for(int i = 0; i < 3; ++i) raw[k++] = ex[i];
    for(int i = 0; i < 3; ++i) raw[k++] = ey[i];
    for(int i = 0; i < 3; ++i) raw[k++] = ez[i];
    loc(C, tmp);
    for(int i = 0; i < 3; ++i) raw[k++] = tmp[i];
    raw[k++] = R - mirrorThick;
    raw[k++] = R;
    loc(S, tmp);
    for(int i = 0; i < 3; ++i) raw[k++] = tmp[i];
    for(int a = 0; a < 3; ++a)
    {
      double u[3] = {double(a == 0), double(a == 1), double(a == 2)};
      ArcVector(sw, u, ax);
      lvec(ax, tmp);
      for(int i = 0; i < 3; ++i) raw[k++] = tmp[i];
    }
    raw[k++] = 0.5 * sideX;
    raw[k++] = 0.5 * sideY;
    raw[k++] = 0.5 * sensorThick;
    ArcPoint(cellTf, Aloc, A);
    loc(A, tmp);
    for(int i = 0; i < 3; ++i) raw[k++] = tmp[i];
    ArcVector(cellTf, nloc, nA);
    lvec(nA, tmp);
    for(int i = 0; i < 3; ++i) raw[k++] = tmp[i];
    raw[k++] = 0.5 * aerogelThick;
    raw[k++] = 1.5 * side;
    ArcVector(cellTf, hexloc, hv);
    lvec(hv, tmp);
    raw[k++] = atan2(tmp[1], tmp[0]) * TMath::RadToDeg();
    raw[k++] = zmin;
    raw[k++] = zmax;
    raw[k++] = apo;
    AddCell(region, raw);
  };

  // ---- barrel: hexagonal pyramids, cell number n = 1..maxCell at z > 0 and 2..maxCell mirrored at z < 0, each repeated in phi
  ExRootConfParam types = GetParam("BarrelCellTypes"); // {column row curvature xposition zposition sensorposition sensortilt}
  int nTypes = types.GetSize();
  const double rin = barrelRin + wall, rout = barrelRout - wall, routShape = rout / sqrt(1.0 + pow(side / rin, 2));
  const double mirrorOriginB = barrelRout - wall - mirrorOrigin - mirrorThick - mirrorShrink;
  const double sensorOriginB = barrelRin + wall + sensorOffset, zstep = 2.0 * apo;
  const ArcTf pyramid = ArcRot(1, -halfPi);
  for(int sgn = -1; sgn <= 1 && nTypes > 0; sgn += 2) // no BarrelCellTypes: no barrel
    for(int n = (sgn < 0 ? 2 : 1); n <= maxCell; ++n)
    {
      int col = n / 2, row = (n % 2) ? 1 : 2, it = -1;
      for(int i = 0; i < nTypes; ++i)
        if(types[i][0].GetInt() == col && types[i][1].GetInt() == row) it = i;
      if(it < 0) throw runtime_error(Form("ARCTracing: BarrelCellTypes has no column %d row %d", col, row));
      double R = types[it][2].GetDouble(), cx = types[it][3].GetDouble(), cz = mirrorOriginB + types[it][4].GetDouble();
      double dpos = types[it][5].GetDouble(), tilt = types[it][6].GetDouble();
      if(sgn < 0)
      {
        cx = -cx;
        dpos = -dpos;
        tilt = -tilt;
      }
      double Cloc[3] = {cx, 0.0, cz - mirrorShrink}, aoff;
      double szp = sensorDepth(sensorOriginB, tilt, sensorB, aoff);
      // sensor: Transform3D(Ry(90 deg - tilt), (-depth, 0, offset)) placed with a further Ry(90 deg); aerogel likewise
      ArcTf sensInner = ArcRot(1, halfPi - tilt), aerInner = ArcRot(1, halfPi);
      sensInner.t[0] = -szp;
      sensInner.t[2] = dpos;
      aerInner.t[0] = -szp - aoff;
      const ArcTf sensTf = ArcRot(1, halfPi) * sensInner, aerTf = ArcRot(1, halfPi) * aerInner;
      double zero[3] = {0, 0, 0}, z1[3] = {0, 0, 1}, Aloc[3], nloc[3], hexloc[3] = {0, 1, 0};
      ArcPoint(aerTf, zero, Aloc);
      ArcVector(aerTf, z1, nloc);
      double zabs = col * zstep - 0.5 * zstep * (row == 2);
      if(sgn < 0) zabs = -zabs;
      double phiOff = 0.5 * phiStep * (row == 2);
      for(int phin = 0; phin < nPhiBarrel; ++phin)
      {
        ArcTf cellTf = ArcRot(2, phiStep * phin + phiOff) * ArcShift(0, 0, zabs) * pyramid;
        addRow(kBarrel, cellTf, Cloc, R, sensTf, sensorB, sensorB, Aloc, nloc, hexloc, rin, routShape);
      }
    }

  // ---- endcaps: hexagonal prisms at grid positions of one 60-degree sector, some with a mirrored copy, 6 sectors, two sides
  ExRootConfParam ec = GetParam("EndcapCells"); // {row column x/side y/apothem mirrored curvature xposition zposition sensorposition sensortilt}
  int nEc = ec.GetSize();
  const double mirrorOriginE = 0.5 * endcapLength - wall - mirrorOrigin, sensorOriginE = -0.5 * endcapLength + wall + sensorOffset;
  const double halfLen = 0.5 * endcapLength - wall;
  for(int i = 0; i < nEc; ++i)
  {
    double x = ec[i][2].GetDouble() * side, y = ec[i][3].GetDouble() * apo;
    bool mirrored = ec[i][4].GetInt() != 0;
    double R = ec[i][5].GetDouble(), cxr = ec[i][6].GetDouble(), cz = mirrorOriginE + ec[i][7].GetDouble();
    double dpos = ec[i][8].GetDouble(), tilt = ec[i][9].GetDouble();
    double alpha = (x != 0.0) ? atan(y / x) : halfPi;
    if(alpha < 0) alpha += TMath::Pi();
    double aoff, szp = sensorDepth(sensorOriginE, tilt, sensorE, aoff);
    double Aloc[3] = {0, 0, sensorOriginE + aoff}, nloc[3] = {0, 0, 1}, hexloc[3] = {1, 0, 0};
    for(int side_ = 0; side_ < 2; ++side_)
    {
      ArcTf env = side_ == 0 ? ArcShift(0, 0, endcapZ) : ArcRot(0, TMath::Pi()) * ArcIdentity();
      if(side_ == 1) env.t[2] = -endcapZ;
      for(int phin = 0; phin < nPhiEndcap; ++phin)
        for(int ref = 0; ref < (mirrored ? 2 : 1); ++ref)
        {
          double sgn = ref ? -1.0 : 1.0;
          double Cloc[3] = {sgn * cxr * cos(alpha), cxr * sin(alpha), cz};
          ArcTf sensTf = ArcRot(2, ref ? halfPi - alpha : alpha - halfPi) * ArcRot(0, tilt) * ArcShift(0, dpos, szp);
          ArcTf cellTf = env * ArcRot(2, endcapPhiStep * phin) * ArcShift(sgn * x, y, 0);
          addRow(kEndcap, cellTf, Cloc, R, sensTf, sensorE, sensorE, Aloc, nloc, hexloc, -halfLen, halfLen);
        }
    }
  }
}

// cells given one by one in the card, for any layout: add Cell {region numbers...} with region barrel (cell z axis outward
// from the beam line, hexagonal pyramid: the apothem grows with the distance along the axis) or endcap (z axis along the
// beam, hexagonal prism), followed by the 44 numbers of a cell-table line (DumpCellsFile writes them in this form)
void ARCTracing::ReadInlineCells()
{
  ExRootConfParam param = GetParam("Cell");
  int n = param.GetSize();
  for(int i = 0; i < n; ++i)
  {
    if(param[i].GetSize() != 45) throw runtime_error(Form("ARCTracing: Cell entry %d has %d values instead of 45 (region and 44 numbers)", i, param[i].GetSize()));
    string region = param[i][0].GetString();
    if(region != "barrel" && region != "endcap") throw runtime_error("ARCTracing: Cell region must be barrel or endcap, not " + region);
    double raw[44];
    for(int k = 0; k < 44; ++k) raw[k] = param[i][k + 1].GetDouble();
    AddCell(region == "barrel" ? kBarrel : kEndcap, raw);
  }
}

// optical tables from the card: add OpticsTable {NAME e1 v1 e2 v2 ...} (energy [eV], value)
void ARCTracing::ReadTablesFromCard()
{
  ExRootConfParam param = GetParam("OpticsTable");
  int n = param.GetSize();
  for(int i = 0; i < n; ++i)
  {
    Table t;
    int m = param[i].GetSize();
    for(int j = 1; j + 1 < m; j += 2)
    {
      t.e.push_back(param[i][j].GetDouble());
      t.v.push_back(param[i][j + 1].GetDouble());
    }
    fTables.push_back(t);
    fTableNames.push_back(param[i][0].GetString());
  }
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
  // geometry and optics come from the card (parameters and tables); a cell file and an optics file can replace them
  string cellFile = GetString("CellFile", "");
  string opticsFile = GetString("OpticsFile", "");
  if(opticsFile.empty())
    ReadTablesFromCard();
  else
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
  fRecoMidPathFallback = GetBool("RecoMidPathFallback", false); // study option: gas emission point beyond the main mirror -> mid-path (default off: as the full-simulation analysis; better reach, fewer gas rings at 45-65 deg)
  // the mirrors are spherical caps (k4geo ARC driver: Sphere(R - t, R, 0, 3.14 / 3.5) around the cell axis, cut to the hexagon):
  // the sphere reflects only within that polar angle of its axis
  fMirrorCapCos = cos(GetDouble("MirrorCapAngle", 3.14 / 3.5));
  fMultiCell = GetBool("MultiCell", true); // as Geant4: a photon reflects off the mirror of the cell it reaches; reconstruction in the main cell (most hits), as the full-simulation analysis
  fSharedCells = GetBool("SharedCells", true); // reconstruct every track against all photons of its cell (overlapping rings in jets)
  // hypothesis test: expected angle of each hypothesis from test photons reconstructed like the data (includes the bias of the
  // reconstruction method, the bending and the cell geometry), and a missing ring weighted by P(N < MinPe | expected photons)
  fCalibratedHypotheses = GetBool("CalibratedHypotheses", true);
  fTestPhotons = GetInt("TestPhotons", 40);
  fDebugMonoEnergy = GetDouble("DebugMonoEnergy", 0.0); // > 0: all photons at this energy (no chromatic error), for studies
  fDebugMidEmission = GetBool("DebugMidEmission", false); // all photons from the middle of the radiator path, for studies

  for(int i = 0; i < 2; ++i)
  {
    BuildGrid(fRad[i], fPDE, fRefl);
    fRecoIndex[i] = 1.0 + fRad[i].indexScale * (fRad[i].n.Eval(fRecoEnergy) - 1.0);
  }
  // cells: a cell file, the parametric barrel and endcaps (BarrelCellTypes, EndcapCells) and single cells (add Cell); each
  // source is optional and they add up
  if(!cellFile.empty()) ReadCells(cellFile.c_str());
  BuildCells();
  ReadInlineCells();
  if(fCells.empty()) throw runtime_error("ARCTracing: no cells (give BarrelCellTypes / EndcapCells, add Cell entries or a CellFile)");
  string dumpCells = GetString("DumpCellsFile", ""); // write every cell as an "add Cell" card line (to start a new layout)
  if(!dumpCells.empty())
  {
    ofstream out(dumpCells.c_str());
    out << "## " << fCells.size() << " ARC cells: region, origin(3), x y z axes (3 each), mirror centre(3), mirror radius inner/outer,\n"
        << "## sensor centre(3), sensor axes x y z (3 each), sensor half sizes(3), aerogel centre(3) and normal(3), aerogel half\n"
        << "## thickness and radius, hexagon vertex angle [deg], cell z range (2), hexagon apothem; mm, in the cell frame after the origin and axes\n"
        << setprecision(10);
    for(size_t i = 0; i < fCells.size(); ++i)
    {
      out << "  add Cell {" << (fCells[i].region == kBarrel ? "barrel" : "endcap");
      for(size_t k = 0; k < fCellRaw[i].size(); ++k) out << " " << fCellRaw[i][k];
      out << "}\n";
    }
  }
  // neighbours: cells of the same region and side whose axis points at mid-depth are within 2.6 apothems (first ring: 2)
  fNeighbours.assign(fCells.size(), vector<int>());
  for(size_t i = 0; i < fCells.size(); ++i)
  {
    const Cell &a = fCells[i];
    double zi = 0.5 * (a.zmin + a.zmax), Pi[3] = {a.O[0] + zi * a.ez[0], a.O[1] + zi * a.ez[1], a.O[2] + zi * a.ez[2]};
    double ai = a.region == kBarrel ? a.apothem * zi / a.zmin : a.apothem;
    for(size_t j = 0; j < fCells.size(); ++j)
    {
      const Cell &b = fCells[j];
      if(j == i || b.region != a.region || b.side != a.side) continue;
      double zj = 0.5 * (b.zmin + b.zmax), d[3] = {b.O[0] + zj * b.ez[0] - Pi[0], b.O[1] + zj * b.ez[1] - Pi[1], b.O[2] + zj * b.ez[2] - Pi[2]};
      if(Norm(d) < 2.6 * ai) fNeighbours[i].push_back(j);
    }
  }

  fNtracks = fNcovered = 0;
  fNcrossed = 0;
  fNsensorNeighbour = 0;
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
       << fNgen[kGas] << " / " << fNdet[kGas] << ", aerogel " << fNgen[kAerogel] << " / " << fNdet[kAerogel]
       << "; photons on the mirror of a neighbouring cell: " << fNcrossed << ", on the array of a neighbouring cell: " << fNsensorNeighbour << endl;
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
bool ARCTracing::TracePhoton(const Cell &c0, int rad, double beta, double *P, double *T, double energy, double *D, int *cellDet)
{
  const Cell *cp = &c0; // the cell whose mirror and sensor the photon reaches (c0 unless MultiCell sends it to a neighbour)
  const Cell &c = c0;
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
    // the mirror the photon meets: its own cell's, if the forward intersection lies inside that cell's hexagon; otherwise (MultiCell)
    // the nearest neighbour mirror met inside the neighbour's hexagon: the photon is above a neighbouring cell when it reaches the
    // mirrors, or was emitted there by an inclined track (beyond its own cell's sphere); it then continues in that cell
    double dC[3] = {Q[0] - c.C[0], Q[1] - c.C[1], Q[2] - c.C[2]};
    double b = Dot(k, dC), cc = Dot(dC, dC) - c.R * c.R;
    double disc = b * b - cc;
    double d = disc >= 0 ? -b + sqrt(disc) : -1.0;
    double M[3] = {Q[0] + d * k[0], Q[1] + d * k[1], Q[2] + d * k[2]};
    int ownWhy = 0;
    if(disc < 0) ownWhy = 4;
    else if(d <= 0) ownWhy = 5;
    else
    {
      double aM = c.region == kBarrel ? c.apothem * M[2] / c.zmin : c.apothem;
      if(M[2] < c.zmin || !InHex(M[0], M[1], aM, c.phi0) || (M[2] - c.C[2]) < fMirrorCapCos * c.R) ownWhy = 7;
    }
    if(ownWhy)
    {
      if(!fMultiCell) return lost(ownWhy);
      double Qw[3], kw[3];
      for(int i = 0; i < 3; ++i)
      {
        Qw[i] = c.O[i] + Q[0] * c.ex[i] + Q[1] * c.ey[i] + Q[2] * c.ez[i];
        kw[i] = k[0] * c.ex[i] + k[1] * c.ey[i] + k[2] * c.ez[i];
      }
      const Cell *nb = 0;
      double bestD = 1e30, Qn[3], kn[3], Qb[3], kb[3];
      for(int j : fNeighbours[&c - &fCells[0]])
      {
        const Cell &n = fCells[j];
        ToLocal(n, Qw, Qn);
        VecToLocal(n, kw, kn);
        double dCn[3] = {Qn[0] - n.C[0], Qn[1] - n.C[1], Qn[2] - n.C[2]};
        double bn = Dot(kn, dCn), ccn = Dot(dCn, dCn) - n.R * n.R, discn = bn * bn - ccn;
        if(discn < 0) continue;
        double dn = -bn + sqrt(discn);
        if(dn <= 0 || dn >= bestD) continue;
        double Mn[3] = {Qn[0] + dn * kn[0], Qn[1] + dn * kn[1], Qn[2] + dn * kn[2]};
        double aMn = n.region == kBarrel ? n.apothem * Mn[2] / n.zmin : n.apothem;
        if(Mn[2] < n.zmin || !InHex(Mn[0], Mn[1], aMn, n.phi0) || (Mn[2] - n.C[2]) < fMirrorCapCos * n.R) continue;
        nb = &n; bestD = dn;
        for(int i = 0; i < 3; ++i) { Qb[i] = Qn[i]; kb[i] = kn[i]; M[i] = Mn[i]; }
      }
      if(!nb) return lost(ownWhy);
      cp = nb; d = bestD;
      for(int i = 0; i < 3; ++i) { Q[i] = Qb[i]; k[i] = kb[i]; }
      if(fCount) ++fNcrossed;
    }
    if(fRad[kGas].abs.Ok() && gRandom->Uniform() > exp(-d / fRad[kGas].abs.Eval(energy))) return lost(6);
    if(fRefl.Ok() && gRandom->Uniform() > fRefl.Eval(energy)) return lost(8);
    double nrm[3] = {(M[0] - cp->C[0]) / cp->R, (M[1] - cp->C[1]) / cp->R, (M[2] - cp->C[2]) / cp->R};
    double kn = Dot(k, nrm);
    for(int i = 0; i < 3; ++i) k[i] -= 2.0 * kn * nrm[i];
    for(int i = 0; i < 3; ++i) Q[i] = M[i];
    // back through the aerogel slab in front of the sensor (of the cell that reflected the photon)
    if(cp != &c0)
    {
      double halfP = fRad[kAerogel].thickness > 0 ? 0.5 * fRad[kAerogel].thickness : cp->aeroHalf;
      zTop = cp->aeroZ + halfP;
      zBot = cp->aeroZ - halfP;
    }
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
  // sensor face: the array of the cell that reflected the photon; with MultiCell, a photon that misses it can land on the array of a
  // neighbouring cell (an inclined track's light reflected off one mirror and focused next to its own sensor), as in Geant4
  auto onSensor = [&](const Cell &sc, const double *Qs, const double *ks, double &dist, double *Dout) {
    double denom = Dot(ks, sc.sez);
    if(denom >= -1e-9) return 14;
    double dS[3] = {sc.S[0] - Qs[0], sc.S[1] - Qs[1], sc.S[2] - Qs[2]};
    double ds = Dot(dS, sc.sez) / denom;
    if(ds <= 0) return 15;
    double H[3] = {Qs[0] + ds * ks[0] - sc.S[0], Qs[1] + ds * ks[1] - sc.S[1], Qs[2] + ds * ks[2] - sc.S[2]};
    double uu = Dot(H, sc.sex), vv = Dot(H, sc.sey);
    int nx = int(2.0 * sc.hs[0] * fSensorSizeScale / fPixelPitch + 1e-6), ny = int(2.0 * sc.hs[1] * fSensorSizeScale / fPixelPitch + 1e-6);
    double hx = 0.5 * nx * fPixelPitch, hy = 0.5 * ny * fPixelPitch;
    if(fabs(uu) > hx || fabs(vv) > hy) return 16;
    int ix = min(max(int((uu + hx) / fPixelPitch), 0), nx - 1), iy = min(max(int((vv + hy) / fPixelPitch), 0), ny - 1);
    double uc = -hx + (ix + 0.5) * fPixelPitch, vc = -hy + (iy + 0.5) * fPixelPitch;
    for(int i = 0; i < 3; ++i) Dout[i] = sc.S[i] + uc * sc.sex[i] + vc * sc.sey[i];
    dist = ds;
    return 0;
  };
  double dist = 0.0;
  int why = onSensor(*cp, Q, k, dist, D);
  if(why == 0)
  {
    if(cellDet) *cellDet = int(cp - &fCells[0]);
    return true;
  }
  if(!fMultiCell) return lost(why);
  double Qw[3], kw[3];
  for(int i = 0; i < 3; ++i)
  {
    Qw[i] = cp->O[i] + Q[0] * cp->ex[i] + Q[1] * cp->ey[i] + Q[2] * cp->ez[i];
    kw[i] = k[0] * cp->ex[i] + k[1] * cp->ey[i] + k[2] * cp->ez[i];
  }
  const Cell *best = 0;
  double bestDist = 1e30, Db[3];
  for(int j : fNeighbours[cp - &fCells[0]])
  {
    const Cell &n = fCells[j];
    double Qn[3], kn[3], Dn[3], dn = 0.0;
    ToLocal(n, Qw, Qn);
    VecToLocal(n, kw, kn);
    if(onSensor(n, Qn, kn, dn, Dn) == 0 && dn < bestDist) { best = &n; bestDist = dn; for(int i = 0; i < 3; ++i) Db[i] = Dn[i]; }
  }
  if(!best) return lost(why);
  for(int i = 0; i < 3; ++i) D[i] = Db[i];
  if(cellDet) *cellDet = int(best - &fCells[0]);
  if(fCount) ++fNsensorNeighbour;
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
  // the reflection point must lie between the two ends of the arc: without a sign change there is no specular solution in
  // [0, betaD] and the bisection would stop at an end point (a meaningless direction, angles up to pi)
  if((f(hi, M) > 0) == (flo > 0)) return -1.0;
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
  int nGen[2], nDet[2]; // nDet: detected in the main cell
  Long64_t index;
  int mainCell; // cell with most of the track's photons (the entry cell without MultiCell or without photons)
  std::map<int, std::vector<double> > hits; // the track's own detected photons per cell (3 numbers each, that cell's frame)
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
      if(fMultiCell && ok3)
      {
        // an inclined track can be above a neighbouring cell when it reaches the mirrors: the gas path ends on that cell's mirror
        // (the entry cell's sphere, continued beyond its hexagon, curves away below it and would cut the path short)
        double Pm[3], Tm[3], Lm[3];
        hTrue.Point(sM, Pm, Tm);
        ToLocal(c, Pm, Lm);
        double aM = c.region == kBarrel ? c.apothem * Lm[2] / c.zmin : c.apothem;
        if(!InHex(Lm[0], Lm[1], aM, c.phi0))
        {
          for(int j : fNeighbours[&c - &fCells[0]])
          {
            const Cell &n = fCells[j];
            bool okn = false;
            double sn = SolveSphere(hTrue, n, sM, okn);
            if(!okn || sn <= sA1) continue;
            hTrue.Point(sn, Pm, Tm);
            ToLocal(n, Pm, Lm);
            double an = n.region == kBarrel ? n.apothem * Lm[2] / n.zmin : n.apothem;
            if(Lm[2] >= n.zmin && InHex(Lm[0], Lm[1], an, n.phi0)) { sM = sn; break; }
          }
        }
      }
      double sLim[2][2] = {{sA1, sM}, {sA0, sA1}};
      double sLow[2] = {sEntry, sA0}; // gas below the aerogel slab (the sensor and cooling plates are ignored)
      bool pathOk = ok0 && ok1 && ok2 && ok3 && sA1 > sA0 && sM > sA1 && sA0 >= sEntry;

      map<int, vector<double> > hitsLocal; // per cell: 3 numbers per detected photon
      map<int, int> nDetCell[2];
      const int entryIdx = int(&c - &fCells[0]);
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
            int cd = entryIdx;
            if(TracePhoton(c, rad, beta, P, T, fDebugMonoEnergy > 0 ? fDebugMonoEnergy : SampleEnergy(fRad[rad], cdf), D, &cd))
            {
              ++nDet[rad];
              ++nDetCell[rad][cd];
              hitsLocal[cd].insert(hitsLocal[cd].end(), D, D + 3);
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
            for(int i = 0; i < 3; ++i) hitsLocal[entryIdx].push_back(c.S[i] + uc * c.sex[i] + vc * c.sey[i]);
          }
        }
      }


      // main cell: the one with most of the track's photons, as the full-simulation analysis; nDet counts that cell only
      int mainIdx = entryIdx;
      size_t mostHits = 0;
      for(map<int, vector<double> >::const_iterator it = hitsLocal.begin(); it != hitsLocal.end(); ++it)
        if(it->second.size() > mostHits) { mostHits = it->second.size(); mainIdx = it->first; }
      if(!fMultiCell) mainIdx = entryIdx;
      for(int rad = 0; rad < 2; ++rad) nDet[rad] = nDetCell[rad].count(mainIdx) ? nDetCell[rad][mainIdx] : 0;

      ARCTracingWork w;
      w.candidate = candidate; w.mother = mother; w.particle = particle; w.mainCell = mainIdx;
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
  map<int, vector<double> > pooled;
  if(fSharedCells)
    for(size_t i = 0; i < works.size(); ++i)
      if(works[i].cell)
        for(map<int, vector<double> >::const_iterator it = works[i].hits.begin(); it != works[i].hits.end(); ++it)
          pooled[it->first].insert(pooled[it->first].end(), it->second.begin(), it->second.end());

  // pass 2: reconstruction and particle identification, track by track, against the pooled hits of its cell
  for(size_t iw = 0; iw < works.size(); ++iw)
  {
    ARCTracingWork &w = works[iw];
    candidate = w.candidate; mother = w.mother; particle = w.particle;
    if(w.cell)
    {
      const Cell &cEntry = *w.cell;          // where the track enters: photon generation (test photons) starts here
      const Cell &c = fCells[w.mainCell];    // the main cell: mirror, sensor and frame of the reconstruction
      const int mainIdx = w.mainCell;
      const Helix &hTrue = w.hTrue, &hReco = w.hReco;
      double pTrue = w.pTrue, charge = w.charge, beta = w.beta, sEntry = w.sEntry, sM = w.sM;
      double nref[2] = {fRecoIndex[kGas], fRecoIndex[kAerogel]};
      double (&sLim)[2][2] = w.sLim; double (&sLow)[2] = w.sLow; double (&L)[2] = w.L;
      int (&nGen)[2] = w.nGen; int (&nDet)[2] = w.nDet;
      bool pathOk = w.pathOk;
      const vector<double> &hitsLocal = fSharedCells ? pooled[mainIdx] : w.hits[mainIdx];
      (void)beta; (void)sEntry; (void)sM;

      // reconstruction: emission point at mid-radiator on the (true or smeared) track
      double E[2][3], Tt[2][3];
      // largest physical Cherenkov angle (beta = 1) plus a margin of 4 windows: anything above is a failed reconstruction
      double thetaMaxReco[2];
      for(int rad = 0; rad < 2; ++rad) thetaMaxReco[rad] = acos(1.0 / fRecoIndex[rad]) + 4.0 * fWindow;
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
        if(rad == kGas && fRecoMidPathFallback)
        {
          // the fixed reconstruction radius can lie beyond the main cell's mirror for an inclined track (its mirror is lower near
          // the cell edge): then use the middle of the reconstructed track between the aerogel and that mirror
          double dE[3] = {E[rad][0] - c.C[0], E[rad][1] - c.C[1], E[rad][2] - c.C[2]};
          if(Norm(dE) >= c.R)
          {
            double halfA = fRad[kAerogel].thickness > 0 ? 0.5 * fRad[kAerogel].thickness : c.aeroHalf;
            bool okA, okP, okS;
            double sA = SolvePlane(hReco, c, c.aeroZ + halfA, sT, okA);
            double sP = SolvePlane(hReco, c, c.C[2] + c.R, sT, okP);
            double sMr = okP ? SolveSphere(hReco, c, sP, okS) : 0.0;
            if(!(okA && okP && okS && sMr > sA)) continue;
            sT = 0.5 * (sA + sMr);
            hReco.Point(sT, Pw, Tw);
            ToLocal(c, Pw, E[rad]);
          }
        }
        VecToLocal(c, Tw, Tt[rad]);
        recoOk[rad] = true;
      }
      vector<double> theta[2];
      size_t nHits = hitsLocal.size() / 3;
      for(int rad = 0; rad < 2; ++rad)
      {
        if(!recoOk[rad]) continue;
        for(size_t ih = 0; ih < nHits; ++ih)
          {
            double tr = RecoAngle(c, E[rad], Tt[rad], &hitsLocal[3 * ih], rad);
            if(tr >= 0 && tr < thetaMaxReco[rad]) theta[rad].push_back(tr);
          }
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
          ToLocal(cEntry, Pw, P);
          VecToLocal(cEntry, Tw, T);
          int cd = -1;
          if(TracePhoton(cEntry, rad, betaH, P, T, SampleEnergy(fRad[rad], cdf), D, &cd) && cd == mainIdx) ++nd;
        }
        return double(nd) / nt;
      };

      // expected reconstructed angle of a hypothesis: test photons at betaH emitted along the reconstructed helix, traced through the
      // cell and reconstructed exactly as the hits; median of the detected ones and their single-photon spread (1.4826 MAD) in sig1
      // (-1 if fewer than 5 are detected)
      auto expectedAngle = [&](int rad, double betaH, double &sig1) {
        sig1 = -1.0;
        if(!recoOk[rad]) return -1.0;
        vector<double> cdf;
        if(Yield(fRad[rad], betaH, &cdf) <= 0) return -1.0;
        double Llow = rad == kGas ? (sLow[1] - sLow[0]) : 0.0, lowFrac = Llow / (Llow + sLim[rad][1] - sLim[rad][0]);
        vector<double> th;
        for(int it = 0; it < fTestPhotons; ++it)
        {
          double sT = gRandom->Uniform() < lowFrac ? gRandom->Uniform(sLow[0], sLow[1]) : gRandom->Uniform(sLim[rad][0], sLim[rad][1]), Pw[3], Tw[3], P[3], T[3], Dh[3];
          hReco.Point(sT, Pw, Tw);
          ToLocal(cEntry, Pw, P);
          VecToLocal(cEntry, Tw, T);
          int cd = -1;
          if(TracePhoton(cEntry, rad, betaH, P, T, SampleEnergy(fRad[rad], cdf), Dh, &cd) && cd == mainIdx)
          {
            double tr = RecoAngle(c, E[rad], Tt[rad], Dh, rad);
            if(tr >= 0 && tr < thetaMaxReco[rad]) th.push_back(tr);
          }
        }
        if(th.size() < 5) return -1.0;
        std::nth_element(th.begin(), th.begin() + th.size() / 2, th.end());
        double med = th[th.size() / 2];
        vector<double> dev(th.size());
        for(size_t i = 0; i < th.size(); ++i) dev[i] = fabs(th[i] - med);
        std::nth_element(dev.begin(), dev.begin() + dev.size() / 2, dev.end());
        sig1 = 1.4826 * dev[dev.size() / 2];
        return med;
      };

      // the measured rings: the window with most photons among the hypotheses (and the true one, for validation)
      double pHyp = fUseTrueMomentumForPID ? pTrue : w.recoP;
      if(pHyp <= 0) pHyp = pTrue;
      struct Ring { double mean, rms, sigma, thetaH, nExp, sig1H; int n; bool expected, found; };
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
          r.sig1H = -1.0;
          if(fCalibratedHypotheses && r.thetaH > 0 && pathOk)
          {
            double te = expectedAngle(rad, betaH, r.sig1H);
            if(te > 0) r.thetaH = te;
          }
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
            {
              // calibrated: the per-track uncertainty from the single-photon spread of the test photons of this hypothesis, not from
              // the scatter of the few measured photons (which is often accidentally small)
              double sg = r.sigma;
              if(fCalibratedHypotheses && r.sig1H > 0)
              {
                double s1 = max(r.sig1H, fSigmaSingleFloor[rad]);
                sg = sqrt(s1 * s1 / r.n + fSigmaCorr[rad] * fSigmaCorr[rad] + fSigmaTrack * fSigmaTrack);
              }
              c2 += pow((r.mean - r.thetaH) / sg, 2);
            }
            else if(r.expected && !r.found)
            {
              if(fCalibratedHypotheses)
              {
                // no ring where nExp photons were expected: -2 ln P(N < MinPe | nExp), a ring needs MinPe photons
                double pk = exp(-r.nExp), sum = 0.0;
                for(int k = 0; k < int(ceil(fMinPe)); ++k) { sum += pk; pk *= r.nExp / (k + 1); }
                c2 += -2.0 * log(max(sum, 1e-300));
              }
              else
                c2 += 2.0 * r.nExp; // -2 ln P(0 | nExp)
            }
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
