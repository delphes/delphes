module TauTagging TauTagging_R05N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R05_N2/VLCjetsR05N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
    # default efficiency formula (misidentification rate)
    add EfficiencyFormula {0} {0.03}
    # efficiency formula for tau-jets
    
    add EfficiencyFormula {15} { 
	(pt < 5) * (0.0) +
	(pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
    } 
} 

module TauTagging TauTagging_R05N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R05_N3/VLCjetsR05N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R05N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R05_N4/VLCjetsR05N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R05N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R05_N5/VLCjetsR05N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R05N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R05_N6/VLCjetsR05N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R07N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R07_N2/VLCjetsR07N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R07N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R07_N3/VLCjetsR07N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R07N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R07_N4/VLCjetsR07N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R07N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R07_N5/VLCjetsR07N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R07N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R07_N6/VLCjetsR07N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R10N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R10_N2/VLCjetsR10N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R10N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R10_N3/VLCjetsR10N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R10N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R10_N4/VLCjetsR10N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R10N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R10_N5/VLCjetsR10N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R10N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R10_N6/VLCjetsR10N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R12N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R12_N2/VLCjetsR12N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R12N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R12_N3/VLCjetsR12N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R12N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R12_N4/VLCjetsR12N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R12N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R12_N5/VLCjetsR12N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R12N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R12_N6/VLCjetsR12N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R15N2 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R15_N2/VLCjetsR15N2
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R15N3 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R15_N3/VLCjetsR15N3
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R15N4 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R15_N4/VLCjetsR15N4
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R15N5 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R15_N5/VLCjetsR15N5
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

module TauTagging TauTagging_R15N6 {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R15_N6/VLCjetsR15N6
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 


module TauTagging TauTagging_R05_inclusive {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R05_inclusive/VLCjetsR05_inclusive
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 


module TauTagging TauTagging_R07_inclusive {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R07_inclusive/VLCjetsR07_inclusive
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 


module TauTagging TauTagging_R10_inclusive {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R10_inclusive/VLCjetsR10_inclusive
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 


module TauTagging TauTagging_R12_inclusive {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R12_inclusive/VLCjetsR12_inclusive
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 


module TauTagging TauTagging_R15_inclusive {
 set ParticleInputArray Delphes/allParticles
 set PartonInputArray Delphes/partons
 set JetInputArray FastJetFinderVLC_R15_inclusive/VLCjetsR15_inclusive
 set DeltaR 0.5
 set TauPTMin 1.0
 set TauEtaMax 4.0
 add EfficiencyFormula {0} {0.03} 
 add EfficiencyFormula {15} { 
 (pt < 5) * (0.0) +
 (pt >=5 && pt < 12.5) * (0.84) +
	(pt >=12.5 && pt < 25) * (0.79) +
	(pt >=25 && pt < 50) * (0.74) +
	(pt >=50 && pt < 75) * (0.66) +
	(pt >=75 && pt < 125) * (0.61) +
	(pt >=125 && pt < 250) * (0.51) +
	(pt >=250 ) * (0.36)
 } 
 } 

