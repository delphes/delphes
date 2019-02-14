set MaxEvents 100

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  ParticlePropagator

  MuonTrackingEfficiency
  MuonSmearing

  TreeWriter
}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.00

  # magnetic field
  set Bz 3.8
}

##############################
# Muon tracking efficiency
##############################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons


  # tracking efficiency formula (efficiency formula as a function of pt, eta, phi, energy, d0 (mm), dz (mm), ctgTheta)
  set EfficiencyFormula {
                          (d0 < 10) *                                                     (pt <= 0.1)   * (0.00) +
                          (d0 < 10) *                (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)      * (0.73) +
                          (d0 < 10) *                (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e2)    * (0.95) +
                          (d0 < 10) *                 (abs(eta) <= 1.5) * (pt > 1.0e2)                  * (0.99) +
                          (d0 < 10) * (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.50) +
                          (d0 < 10) * (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e2) * (0.83) +
                          (d0 < 10) * (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e2)                * (0.90) +
                          (d0 < 10) * (abs(eta) > 2.5)                                                  * (0.00) +
                          (d0 > 10) *                                                                     (0.00)
                         }
}



################################                                                                    
# Muon track parameter smearing                                                                  
################################                                                                    

module TrackSmearing MuonSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # magnetic field
  set Bz 3.8


  set PResolutionFormula { (abs(eta) <= 0.5) * (pt > 0.1)      * sqrt(0.03^2 + pt^2*1.3e-3^2) +
                           (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 0.1) * sqrt(0.05^2 + pt^2*1.7e-3^2) +
                           (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1) * sqrt(0.15^2 + pt^2*3.1e-3^2)
                          }

  set CtgThetaResolutionFormula { 0.0 }

  set PhiResolutionFormula { 0.0 }

  set D0ResolutionFormula {
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 0.3543 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 0.2809 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 0.2304 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 0.1917 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 0.1737 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 0.1434 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.1060 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.0893 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.0753 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.0670 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.0577 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.0524 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.0452 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.0376 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.0350 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.0324 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.0283 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.0258 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.0237 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.0211 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.0191 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.0164 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.0150 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0143 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0130 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0130 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0116 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0116 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0110 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0110 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0110 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0104 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0109 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 134.7130 ) * 0.0110 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 0.4564 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 0.3580 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 0.3010 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 0.2353 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 0.2026 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 0.1595 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.1383 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.1119 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.0926 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.0816 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.0663 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.0553 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.0488 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.0431 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.0399 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.0357 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.0313 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.0277 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.0233 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.0221 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.0214 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.0180 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.0155 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0141 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0128 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0134 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0121 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0108 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0101 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0101 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0101 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0102 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0088 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 134.7130 ) * 0.0095 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 0.6970 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 0.6046 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 0.5315 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 0.4306 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 0.3398 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 0.2788 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.2387 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.1814 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.1557 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.1230 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.1009 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.0914 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.0767 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.0638 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.0544 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.0468 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.0425 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.0385 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.0331 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.0278 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.0256 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.0236 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.0217 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0196 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0176 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0165 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0157 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0150 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0144 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0144 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0137 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0130 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0137 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 134.7130 ) * 0.0137 
  }
  set DZResolutionFormula {
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 0.3693 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 0.3135 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 0.3125 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 0.2578 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 0.2221 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 0.1936 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.1686 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.1351 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.1113 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.0983 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.0882 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.0786 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.0684 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.0615 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.0551 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.0516 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.0484 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.0450 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.0416 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.0416 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.0382 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.0350 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.0317 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0316 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0316 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0316 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0348 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0316 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0316 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0316 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0284 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0283 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0315 +\
      ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 134.7130 ) * 0.0318 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 0.7062 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 0.6010 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 0.5992 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 0.4959 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 0.3877 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 0.3199 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.2649 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.2518 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.1982 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.1587 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.1399 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.1199 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.1031 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.0967 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.0805 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.0736 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.0707 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.0603 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.0609 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.0541 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.0511 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.0443 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.0409 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0408 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0409 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0377 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0375 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0377 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0342 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0342 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0343 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0343 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0309 +\
      ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 134.7130 ) * 0.0310 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 2.1717 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 2.0715 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 2.0679 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 1.7679 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 1.4393 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 1.1997 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.9800 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.8251 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.6695 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.5545 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.4366 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.3711 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.3319 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.2721 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.2443 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.2085 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.1816 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.1641 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.1451 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.1317 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.1117 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.1020 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.1017 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0983 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0882 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0847 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0814 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0784 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0817 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0750 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0816 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0820 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0814 +\
      ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 134.7130 ) * 0.0850 
  }
}


##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle
  add Branch MuonSmearing/muons Muon Muon
}
