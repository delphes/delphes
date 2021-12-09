module TauTagging TauTagging_R02N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
    set JetInputArray JetMomentumSmearing_VLCR02N2/JER_VLCjetsR02N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
    # default efficiency formula (misidentification rate)
    add EfficiencyFormula {0} {0.02}
    add EfficiencyFormula {11} {0.001}
    # efficiency formula for tau-jets

    add EfficiencyFormula {15} {
	(pt < 5) * (0.0) +
	(pt >=5) * (0.80)
    }
}

module TauTagging TauTagging_R02N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR02N3/JER_VLCjetsR02N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R02N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR02N4/JER_VLCjetsR02N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R02N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR02N5/JER_VLCjetsR02N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R02N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR02N6/JER_VLCjetsR02N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }


module TauTagging TauTagging_R05N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
    set JetInputArray JetMomentumSmearing_VLCR05N2/JER_VLCjetsR05N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
    # default efficiency formula (misidentification rate)
    add EfficiencyFormula {0} {0.02}
    add EfficiencyFormula {11} {0.001}
    # efficiency formula for tau-jets

    add EfficiencyFormula {15} {
	(pt < 5) * (0.0) +
	(pt >=5) * (0.80)
    }
}

module TauTagging TauTagging_R05N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR05N3/JER_VLCjetsR05N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R05N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR05N4/JER_VLCjetsR05N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R05N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR05N5/JER_VLCjetsR05N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R05N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR05N6/JER_VLCjetsR05N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R07N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR07N2/JER_VLCjetsR07N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R07N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR07N3/JER_VLCjetsR07N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R07N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR07N4/JER_VLCjetsR07N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R07N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR07N5/JER_VLCjetsR07N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R07N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR07N6/JER_VLCjetsR07N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R10N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR10N2/JER_VLCjetsR10N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R10N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR10N3/JER_VLCjetsR10N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R10N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR10N4/JER_VLCjetsR10N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R10N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR10N5/JER_VLCjetsR10N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R10N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR10N6/JER_VLCjetsR10N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R12N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR12N2/JER_VLCjetsR12N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R12N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR12N3/JER_VLCjetsR12N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R12N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR12N4/JER_VLCjetsR12N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R12N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR12N5/JER_VLCjetsR12N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R12N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR12N6/JER_VLCjetsR12N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R15N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR15N2/JER_VLCjetsR15N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R15N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR15N3/JER_VLCjetsR15N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R15N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR15N4/JER_VLCjetsR15N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R15N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR15N5/JER_VLCjetsR15N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }

module TauTagging TauTagging_R15N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR15N6/JER_VLCjetsR15N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }


 module TauTagging TauTagging_R02_inclusive {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetMomentumSmearing_VLCR02_inclusive/JER_VLCjetsR02_inclusive
  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 2.5
  add EfficiencyFormula {0} {0.02}
  add EfficiencyFormula {11} {0.001}
  add EfficiencyFormula {15} {
  (pt < 10) * (0.0) +
  (pt >=10) * (0.80)
  }
  }

  module TauTagging TauTagging_R05_inclusive {
   set ParticleInputArray Delphes/allParticles
   set PartonInputArray Delphes/partons
   set JetInputArray JetMomentumSmearing_VLCR05_inclusive/JER_VLCjetsR05_inclusive
   set DeltaR 0.5
   set TauPTMin 1.0
   set TauEtaMax 2.5
   add EfficiencyFormula {0} {0.02}
   add EfficiencyFormula {11} {0.001}
   add EfficiencyFormula {15} {
   (pt < 10) * (0.0) +
   (pt >=10) * (0.80)
   }
   }


module TauTagging TauTagging_R07_inclusive {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR07_inclusive/JER_VLCjetsR07_inclusive
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }


module TauTagging TauTagging_R10_inclusive {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR10_inclusive/JER_VLCjetsR10_inclusive
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }


module TauTagging TauTagging_R12_inclusive {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR12_inclusive/JER_VLCjetsR12_inclusive
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }


module TauTagging TauTagging_R15_inclusive {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray JetMomentumSmearing_VLCR15_inclusive/JER_VLCjetsR15_inclusive
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 2.5
 add EfficiencyFormula {0} {0.02}
 add EfficiencyFormula {11} {0.001}
 add EfficiencyFormula {15} {
 (pt < 10) * (0.0) +
 (pt >=10) * (0.80)
 }
 }
