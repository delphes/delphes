# ALFA cards with the ARC ray-tracing fast simulation (ARCTracing)

`ARCTracing` is a Delphes module that simulates the ARC (Array of RICH Cells) of the ALFA detector by ray tracing
instead of parametrised formulas: for every charged track it generates the Cherenkov photons along the helix in the
aerogel and in the gas from the optical tables of the Geant4 materials, traces them through the cell optics of the
k4geo geometry (aerogel slab, spherical mirror, SiPM array, pixels) and reconstructs the ring as the full-simulation
analysis does. Thresholds, photon yields, chromatic error, ring size on the array, bending and scattering come out of
the model; nothing is fitted to the full simulation.

Validated against the Geant4 full simulation (k4geo `ALFA_o1_v01`) for the three ARC configurations of the card, in the 14 cell-sized polar-angle bins
from 12 to 90 deg: of the 174 three-sigma pi/K and K/p reaches (gas and aerogel), 131 agree within 5 %, 35 within 10 %, and the
8 others are mostly in the barrel-endcap gap (39-44 deg). Validated range of the aerogel index: n <= 1.030.

## Files

| File | What |
|---|---|
| `modules/ARCTracing.h`, `modules/ARCTracing.cc` | the module |
| `classes/DelphesClasses.h/.cc`, `modules/TreeWriter.cc` | 16 ARC members on `Candidate` and `Track` (below) |
| `modules/ModulesLinkDef.h`, `Makefile` | registration and build |
| `cards/ARCTracing-ALFA/delphes_card_ALFA.tcl` | the ALFA card: 3 T, all-MAPS tracker, ARC, GRAiNITA ECAL. The ARC configuration (gas, SiPM arrays) is chosen with two switches at the top; the ARC geometry, the cell optics and the optical tables are parameters of the ARCTracing module in the card. No other file is read. |

## The ARC configurations of the card

The geometry is that of k4geo ALFA_o1_v01: the ALFA_o1_v00 detector of the ALFA group with the ARC cells of the CLD design
(CLD_o3_v01) as barrel and endcaps, 10 mm aerogel with n = 1.030 (raised from the CLD value 1.019), 0.8 mm pixels. Two
switches at the top of the card choose the configuration:

```tcl
set ARCGas    Xe2bar     ;# Xe2bar: xenon at 2 bar (n = 1.001313 at 3 eV) | C4F10: 1 bar (n = 1.001331)
set ARCArrays baseline   ;# baseline: 80 x 80 mm SiPM arrays (configuration B) | casestudy: 89.6 mm barrel / 108.8 mm endcap (configuration C)
```

They set the gas index table, the SiPM array size and the ARC material in the track fit (Xe 2 bar: 6.0 % X0 at 90 deg with
the baseline arrays, 6.2 % with the case-study arrays; C4F10: 4.4 %).

| label | aerogel | SiPM arrays | status |
|---|---|---|---|
| A | n = 1.019, 10 mm | 80 x 80 mm (CLD ARC cells) | validated in the central cell with an earlier version of the module (not in the card) |
| **B, baseline** (`ARCArrays baseline`, both gases) | **n = 1.030, 10 mm** | **80 x 80 mm** | **validated in 14 polar-angle bins (this version)** |
| **C, case study** (`ARCArrays casestudy`, Xe 2 bar) | **n = 1.030, 10 mm** | **89.6 mm barrel / 108.8 mm endcap, the largest that fit the cells** | **validated in 14 polar-angle bins (this version)** |

Configuration C is the ARC case study by Apranik. Its gain is in the number of rings, not in the separation per ring: for central
tracks an aerogel ring is found for 42 % of the pions at 2 GeV and 61 % of the kaons at 3 GeV (baseline: 30 % and 50 %).

## Another ARC, another gas

Everything the module needs is in the card:

- the vessel and cell layout: `BarrelInnerRadius`, `BarrelOuterRadius`, `EndcapZ`, `EndcapLength`, `VesselWallThickness`,
  `HexagonSide`, `BarrelPhiStep`, `BarrelSectors`, `BarrelMaxCell`, `EndcapSectors`;
- the optics of every cell type, one line each: `add BarrelCellTypes {column row mirror-radius mirror-x mirror-z sensor-offset
  sensor-tilt}` and `add EndcapCells {row column x/side y/apothem mirrored-copy mirror-radius mirror-x mirror-z sensor-offset
  sensor-tilt}` (mm, rad), placed as in the k4geo ARC_geo_o1_v01 construction;
- the mirrors, sensors and aerogel: `MirrorThickness`, `MirrorZOrigin`, `BarrelMirrorZShrink`, `MirrorCapAngle`,
  `SensorThickness`, `SensorZOffset`, `BarrelSensorSize`, `EndcapSensorSize`, `CellAerogelThickness`;
- the materials: `add OpticsTable {name energy value energy value ...}` (energy in eV): refractive indices, absorption and
  Rayleigh lengths, SiPM detection efficiency, mirror reflectivity. A new gas is one more table and `set GasIndexTable name`
  (or an existing table with `GasIndexScale` for another pressure).

The cells built from the parameters of this card reproduce the cells of the k4geo geometry to 1e-6 mm (checked against the
DD4hep geometry for all 1203 cells and all three configurations). A table of cells (`CellFile`) and of optics (`OpticsFile`)
can still replace the parameters, for a geometry that does not follow this construction.

## Build and run

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh   # ROOT, Pythia8, HepMC
./configure                                           # regenerates the Makefile with the new module
make -j8 HAS_PYTHIA8=true
./DelphesPythia8 cards/ARCTracing-ALFA/delphes_card_ALFA.tcl examples/Pythia8/ee_zh.cmnd zh.root
```

Single particles, as used for the validation: `examples/Pythia8/configParticleGun.cmnd` with `Main:spareMode1 = 211`
(or 321, 2212), `Main:spareParm1 = 60` (maximum momentum, GeV) and `Main:spareParm2 = 2.6` (maximum |eta|).

Speed: the ARC adds about 5 ms per e+e- -> ZH event at 240 GeV (about 210 us per charged track): 8.8 ms per event with the
baseline arrays (9.3 ms with the case-study arrays) against 3.8 ms without, Pythia and ROOT output included (10k events, AMD EPYC
7302). About half of it is the test photons of the calibrated hypothesis test (TestPhotons = 40); TestPhotons = 20 halves that part.

## Module parameters (the values of the cards)

| Parameter | Meaning | Card |
|---|---|---|
| geometry and optics parameters | see "Another ARC, another gas" above | the ALFA_o1_v01 values |
| `CellFile`, `OpticsFile` | optional: a cell table and an optics table instead of the parameters | "" (not used) |
| `Bz` | magnetic field [T] | `$::B` (3) |
| `GasIndexTable`, `GasIndexScale` | n(E) of the gas (an `OpticsTable`); n - 1 scaled (pressure) | from `ARCGas`: `RINDEX__Xenon_ARC` or `RINDEX__C4F10_PFRICH`; 1.0 |
| `AerogelIndexTable`, `AerogelIndexScale` | n(E) of the aerogel; n - 1 scaled | `RINDEX__Aerogel_PFRICH`, 1.0 (n = 1.030) |
| `AerogelRayleighTable`, `AerogelAbsorptionTable`, `GasAbsorptionTable` | attenuation tables; "" = none | `RAYLEIGH__Aerogel_BelleII`, "", "" (as the Geant4 material) |
| `AerogelThickness`, `AerogelSlabInPath` | slab thickness [mm] (0: from the geometry, 10 mm); gas photons cross the slab on the way back | 0, true |
| `PDETable`, `PDEScale`, `MirrorReflectivityTable`, `YieldScale` | photon detection | `ARC_SiPM_QuantumEfficiency`, 1, `REFLECTIVITY__ARC_Mirror`, 1 |
| `PixelPitch`, `SensorSizeScale`, `SensorShift`, `MirrorRadiusScale`, `DarkCountsPerSensor` | readout and design knobs | 0.8 mm, 1.0, 0, 1.0, 0 |
| `RecoEnergy`, `RecoRadiusGas/Aerogel`, `RecoZGas/Aerogel` | reconstruction: index at 3 eV, emission assumed at mid-radiator (barrel r, endcap z, mm) | 3.0, 2000/1925, 2110/2025 |
| `Window`, `MinPe` | window around the hypothesis angle [rad]; photons needed for a ring | 0.025, 3 |
| `GasSigmaSingleFloor`, `AerogelSigmaSingleFloor`, `GasSigmaCorr`, `AerogelSigmaCorr`, `SigmaTrack` | per-track resolution: sigma^2 = max(rms^2, floor^2)/N + corr^2 + track^2 [rad] | 1.0e-3, 1.5e-3, 0, 0, 0 |
| `SharedCells` | all tracks of a cell are reconstructed against the pooled photons of the cell | true |
| `MultiCell` | a photon reaches the mirror (and array) of a neighbouring cell when it is above that cell; each track is reconstructed in its main cell (most of its photons), as the full-simulation analysis | true |
| `MirrorCapAngle` | the mirrors are spherical caps of this polar half-angle around the cell axis, as in the k4geo driver [rad] | 3.14/3.5 |
| `RecoMidPathFallback` | study option: gas emission point at mid-path when the fixed radius is beyond the main cell's mirror | false |
| `CalibratedHypotheses`, `TestPhotons` | hypothesis test with the expected angle and single-photon spread of each hypothesis taken from test photons reconstructed like the data, and a missing ring weighted by P(N < MinPe); number of test photons per hypothesis and radiator | true, 40 |
| `UseTrueTrack`, `UseTrueMomentumForPID`, `TrueHypothesisAngle` | validation switches (true track and momentum; window around the true angle) | false, false, false |

Another ARC geometry is another cell table; another gas, pressure or aerogel is another table or scale.

## Track members

```
Track.CherenkovAngleGas / CherenkovAngleAerogel         true Cherenkov angle [rad], -1 if no ring
Track.CherenkovRecoAngleGas / CherenkovRecoAngleAerogel reconstructed angle [rad], -1 if no ring
Track.SigmaThetaCGas / SigmaThetaCAerogel               per-track angular resolution [rad]
Track.NpeGas / NpeAerogel                               photons in the reconstructed ring
Track.ARCchi2Pi / ARCchi2K / ARCchi2P                   per-hypothesis chi2 over both radiators, -1 if no ring
Track.ARCpid                                            +-211 / +-321 / +-2212 (smallest chi2), 0 if no ring
Track.ARCprobPi / ARCprobK / ARCprobP                   exp(-chi2/2) normalised over the three hypotheses, -1 if no ring
Track.ARCmass                                           mass of the best hypothesis [GeV], -1 if no ring
```

For hypothesis h the chi2 adds, per radiator, ((reco - theta_h)/sigma_h)^2 when a ring is expected and found,
-2 ln P(N < MinPe | N_exp,h) when a ring is expected and none is found (a ring needs MinPe photons), and (reco/sigma)^2
when a ring is found where h predicts none. theta_h, sigma_h and N_exp,h come from test photons emitted at beta_h along
the reconstructed track, traced through the cell and reconstructed exactly like the hits: theta_h is their median, so the
bias of the reconstruction method (about -2.5 mrad in the aerogel, as in the full simulation), the bending and the cell
geometry are part of the expectation; sigma_h = sigma_1/sqrt(N) with sigma_1 their single-photon spread. Without this
calibration (CalibratedHypotheses false) the test uses acos(1/(n beta)) and the scatter of the measured photons, which
ignores the bias of the method and underestimates the uncertainty of rings with few photons.

## The algorithm (per event, two passes)

1. For every track: the cell the helix enters (barrel or endcap); photons along the helix in aerogel and gas up to the
   mirror of the cell the track is in, N ~ Poisson(L dN/dL) with dN/dL = 36.98 /eV/mm * integral of sin^2 theta_c(E) PDE(E)
   R(E) dE, energies sampled from the integrand; each photon refracted out of the slab, Rayleigh-scattered (1 + cos^2) with
   the tabulated length, reflected by the first mirror cap it meets inside its cell's hexagon (its own cell's or a
   neighbour's: inclined tracks cross from one cell to the next), traced back through that cell's slab to its SiPM array
   (or a neighbouring array) and binned in pixels. Photons are kept per cell and pooled over the tracks of the event.
2. For every track: its main cell (most of its photons); inverse ray tracing of every pixel hit of that cell (pixel
   centre -> mirror point -> emission point at the fixed mid-radiator radius) to an angle, angles above the physical
   maximum dropped; window of +-25 mrad around the hypothesis angle; two passes of 3-sigma clipping with sigma =
   1.4826 MAD; mean, rms, N; resolution as above; chi2 and probabilities for pi, K, p.

3-sigma reach in the central barrel cell (82.6-90 deg), fast / full simulation, GeV (gas pi/K, K/p; aerogel pi/K, K/p):
case study C 24.6 / 25.0, 43.6 / 43.6, 7.9 / 7.6, 14.1 / 13.9; baseline B Xe 2 bar 24.9 / 24.7, 44.7 / 44.1, 7.7 / 7.4,
13.4 / 13.2; baseline B C4F10 25.7 / 24.4, 45.4 / 44.3, 7.6 / 7.4, 13.6 / 12.8. Known difference: at 45-60 deg the fast
simulation finds 80-90 % of the gas rings of the full simulation (light reaching a cell through more than one neighbour).
