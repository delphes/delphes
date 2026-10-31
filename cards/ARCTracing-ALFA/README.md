# ALFA cards with the ARC ray-tracing fast simulation (ARCTracing)

`ARCTracing` is a Delphes module that simulates the ARC (Array of RICH Cells) of the ALFA detector by ray tracing
instead of parametrised formulas: for every charged track it generates the Cherenkov photons along the helix in the
aerogel and in the gas from the optical tables of the Geant4 materials, traces them through the cell optics of the
k4geo geometry (aerogel slab, spherical mirror, SiPM array, pixels) and reconstructs the ring as the full-simulation
analysis does. Thresholds, photon yields, chromatic error, ring size on the array, bending and scattering come out of
the model; nothing is fitted to the full simulation.

Validated against the Geant4 full simulation (k4geo `ALFA_o1_v01`) on three configurations and both gases: 23 of 24
three-sigma reaches within 5 %, two configurations predicted blind. Validated range of the aerogel index: n <= 1.030.

## Files

| File | What |
|---|---|
| `modules/ARCTracing.h`, `modules/ARCTracing.cc` | the module |
| `classes/DelphesClasses.h/.cc`, `modules/TreeWriter.cc` | 16 ARC members on `Candidate` and `Track` (below) |
| `modules/ModulesLinkDef.h`, `Makefile` | registration and build |
| `cards/ARCTracing-ALFA/delphes_card_ALFA_Xe2bar.tcl` | the ALFA card: 3 T, all-MAPS tracker, ARC with xenon at 2 bar and aerogel n = 1.030, GRAiNITA ECAL |
| `cards/ARCTracing-ALFA/delphes_card_ALFA_C4F10.tcl` | the same with C4F10 at 1 bar in the ARC |
| `cards/ARCTracing-ALFA/cells_ALFA_o1_v01.txt` | the 1203 ARC cells of `ALFA_o1_v01` (mirror, SiPM array, aerogel per cell; 80 x 80 mm arrays) |
| `cards/ARCTracing-ALFA/cells_ALFA_o1_v01_bigSiPM.txt` | the same cells with the largest arrays that fit (89.6 mm barrel, 108.8 mm endcap) |
| `cards/ARCTracing-ALFA/arc_optics_tables.txt` | n(E) of Xe, C4F10 and aerogel, aerogel Rayleigh length, SiPM PDE, mirror reflectivity (from the Geant4 materials) |
| `cards/ARCTracing-ALFA/dump_arc_cells.py` | writes a cell table from a k4geo compact file (`python3 dump_arc_cells.py <compact.xml> <out.txt>`) |
| `cards/ARCTracing-ALFA/gun_{pip,kp,p}.cmnd` | Pythia8 single-particle guns (pi+, K+, p), log-uniform p in 0.1-60 GeV, abs(eta) < 2.6 |

## Which ARC configuration the cards describe

The validation and the design study use three labels; the cards in this repository are **B**:

| label | aerogel | SiPM arrays | status |
|---|---|---|---|
| A | n = 1.019, 10 mm | 80 x 80 mm (as built) | validated (blind) |
| **B** | **n = 1.030, 10 mm** | **80 x 80 mm (as built in k4geo ALFA_o1_v01)** | **validated; the two cards here** |
| C | n = 1.030, 10 mm | 89.6 mm barrel / 108.8 mm endcap, the largest that fit the cells | validated (blind); the optimised design |

Both cards state this in their first lines. To run the optimised design, change one line in either card:
`set CellFile cards/ARCTracing-ALFA/cells_ALFA_o1_v01_bigSiPM.txt` (the table is included). Everything else -- gas, pixels,
mirrors, reconstruction -- is identical. A dedicated `delphes_card_ALFA_optimised.tcl` will follow.

## Build and run

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh   # ROOT, Pythia8, HepMC
./configure                                           # regenerates the Makefile with the new module
make -j8 HAS_PYTHIA8=true
./DelphesPythia8 cards/ARCTracing-ALFA/delphes_card_ALFA_Xe2bar.tcl cards/ARCTracing-ALFA/gun_kp.cmnd kp.root
```

The cards refer to `cards/ARCTracing-ALFA/...` relative to the Delphes directory: run from there or make the two paths absolute.
About 500 tracks per second per core.

## Module parameters (the values of the cards)

| Parameter | Meaning | Card |
|---|---|---|
| `CellFile`, `OpticsFile` | cell table, optical tables | `cards/ARCTracing-ALFA/...` |
| `Bz` | magnetic field [T] | `$::B` (3) |
| `GasIndexTable`, `GasIndexScale` | n(E) of the gas; n - 1 scaled (pressure) | `RINDEX__Xenon_ARC`, 1.0 (Xe 2 bar); `RINDEX__C4F10_PFRICH` |
| `AerogelIndexTable`, `AerogelIndexScale` | n(E) of the aerogel; n - 1 scaled | `RINDEX__Aerogel_PFRICH`, 1.0 (n = 1.030) |
| `AerogelRayleighTable`, `AerogelAbsorptionTable`, `GasAbsorptionTable` | attenuation tables; "" = none | `RAYLEIGH__Aerogel_BelleII`, "", "" (as the Geant4 material) |
| `AerogelThickness`, `AerogelSlabInPath` | slab thickness [mm] (0: from the geometry, 10 mm); gas photons cross the slab on the way back | 0, true |
| `PDETable`, `PDEScale`, `MirrorReflectivityTable`, `YieldScale` | photon detection | `ARC_SiPM_QuantumEfficiency`, 1, `REFLECTIVITY__ARC_Mirror`, 1 |
| `PixelPitch`, `SensorSizeScale`, `SensorShift`, `MirrorRadiusScale`, `DarkCountsPerSensor` | readout and design knobs | 0.8 mm, 1.0, 0, 1.0, 0 |
| `RecoEnergy`, `RecoRadiusGas/Aerogel`, `RecoZGas/Aerogel` | reconstruction: index at 3 eV, emission assumed at mid-radiator (barrel r, endcap z, mm) | 3.0, 2000/1925, 2110/2025 |
| `Window`, `MinPe` | window around the hypothesis angle [rad]; photons needed for a ring | 0.025, 3 |
| `GasSigmaSingleFloor`, `AerogelSigmaSingleFloor`, `GasSigmaCorr`, `AerogelSigmaCorr`, `SigmaTrack` | per-track resolution: sigma^2 = max(rms^2, floor^2)/N + corr^2 + track^2 [rad] | 1.0e-3, 1.5e-3, 0, 0, 0 |
| `SharedCells` | all tracks of a cell are reconstructed against the pooled photons of the cell | true |
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

For hypothesis h the chi2 adds, per radiator, ((reco - theta_h)/sigma)^2 when a ring is expected and found,
2 N_exp,h (= -2 ln P(0 | N_exp,h)) when a ring is expected and none is found, and (reco/sigma)^2 when a ring is found
where h predicts none; N_exp,h is obtained by tracing test photons at beta_h through the cell, so the acceptance is
in the test.

## The algorithm (per event, two passes)

1. For every track: the cell the helix enters (barrel or endcap), in that cell's optical frame; photons along the helix
   in aerogel and gas, N ~ Poisson(L dN/dL) with dN/dL = 36.98 /eV/mm * integral of sin^2 theta_c(E) PDE(E) R(E) dE,
   energies sampled from the integrand; each photon refracted out of the slab, Rayleigh-scattered (1 + cos^2) with the
   tabulated length, reflected on the mirror sphere (lost outside the hexagonal cell), traced back through the slab to
   the SiPM array and binned in pixels. Photons are pooled per cell.
2. For every track: inverse ray tracing of every pixel hit of its cell (pixel centre -> mirror point -> emission point at
   mid-radiator) to an angle; window of +-25 mrad around the hypothesis angle; two passes of 3-sigma clipping with
   sigma = 1.4826 MAD; mean, rms, N; resolution as above; chi2 and probabilities for pi, K, p.

Performance of the two cards in the central barrel cell (3-sigma reach, fast / full simulation, GeV): Xe 2 bar
gas pi/K 23.7 / 24.7, K/p 40.8 / 44.1, aerogel pi/K 7.3 / 7.4, K/p 13.0 / 13.2; C4F10 gas 23.8 / 24.4, 41.3 / 44.3,
aerogel 7.3 / 7.4, 13.0 / 12.8.
