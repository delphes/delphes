module BTagging BTagging_JER_WP50_R05N2 {
	set JetInputArray JetMomentumSmearing_VLCR05N2/JER_VLCjetsR05N2
	set BitNumber 0

	# 50% efficiency working point
    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
    # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
    # gluon's PDG code has the lowest priority

	# based on CLICdp-Note-2014-002

	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R05N2 {
	set JetInputArray JetMomentumSmearing_VLCR05N2/JER_VLCjetsR05N2
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R05N2 {
	set JetInputArray JetMomentumSmearing_VLCR05N2/JER_VLCjetsR05N2
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R05N3 {
	set JetInputArray JetMomentumSmearing_VLCR05N3/JER_VLCjetsR05N3
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R05N3 {
	set JetInputArray JetMomentumSmearing_VLCR05N3/JER_VLCjetsR05N3
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R05N3 {
	set JetInputArray JetMomentumSmearing_VLCR05N3/JER_VLCjetsR05N3
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R05N4 {
	set JetInputArray JetMomentumSmearing_VLCR05N4/JER_VLCjetsR05N4
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R05N4 {
	set JetInputArray JetMomentumSmearing_VLCR05N4/JER_VLCjetsR05N4
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R05N4 {
	set JetInputArray JetMomentumSmearing_VLCR05N4/JER_VLCjetsR05N4
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R05N5 {
	set JetInputArray JetMomentumSmearing_VLCR05N5/JER_VLCjetsR05N5
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R05N5 {
	set JetInputArray JetMomentumSmearing_VLCR05N5/JER_VLCjetsR05N5
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R05N5 {
	set JetInputArray JetMomentumSmearing_VLCR05N5/JER_VLCjetsR05N5
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R05N6 {
	set JetInputArray JetMomentumSmearing_VLCR05N6/JER_VLCjetsR05N6
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R05N6 {
	set JetInputArray JetMomentumSmearing_VLCR05N6/JER_VLCjetsR05N6
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R05N6 {
	set JetInputArray JetMomentumSmearing_VLCR05N6/JER_VLCjetsR05N6
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R07N2 {
	set JetInputArray JetMomentumSmearing_VLCR07N2/JER_VLCjetsR07N2
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R07N2 {
	set JetInputArray JetMomentumSmearing_VLCR07N2/JER_VLCjetsR07N2
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R07N2 {
	set JetInputArray JetMomentumSmearing_VLCR07N2/JER_VLCjetsR07N2
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R07N3 {
	set JetInputArray JetMomentumSmearing_VLCR07N3/JER_VLCjetsR07N3
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R07N3 {
	set JetInputArray JetMomentumSmearing_VLCR07N3/JER_VLCjetsR07N3
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R07N3 {
	set JetInputArray JetMomentumSmearing_VLCR07N3/JER_VLCjetsR07N3
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R07N4 {
	set JetInputArray JetMomentumSmearing_VLCR07N4/JER_VLCjetsR07N4
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R07N4 {
	set JetInputArray JetMomentumSmearing_VLCR07N4/JER_VLCjetsR07N4
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R07N4 {
	set JetInputArray JetMomentumSmearing_VLCR07N4/JER_VLCjetsR07N4
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R07N5 {
	set JetInputArray JetMomentumSmearing_VLCR07N5/JER_VLCjetsR07N5
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R07N5 {
	set JetInputArray JetMomentumSmearing_VLCR07N5/JER_VLCjetsR07N5
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R07N5 {
	set JetInputArray JetMomentumSmearing_VLCR07N5/JER_VLCjetsR07N5
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R07N6 {
	set JetInputArray JetMomentumSmearing_VLCR07N6/JER_VLCjetsR07N6
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R07N6 {
	set JetInputArray JetMomentumSmearing_VLCR07N6/JER_VLCjetsR07N6
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R07N6 {
	set JetInputArray JetMomentumSmearing_VLCR07N6/JER_VLCjetsR07N6
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R10N2 {
	set JetInputArray JetMomentumSmearing_VLCR10N2/JER_VLCjetsR10N2
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R10N2 {
	set JetInputArray JetMomentumSmearing_VLCR10N2/JER_VLCjetsR10N2
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R10N2 {
	set JetInputArray JetMomentumSmearing_VLCR10N2/JER_VLCjetsR10N2
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R10N3 {
	set JetInputArray JetMomentumSmearing_VLCR10N3/JER_VLCjetsR10N3
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R10N3 {
	set JetInputArray JetMomentumSmearing_VLCR10N3/JER_VLCjetsR10N3
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R10N3 {
	set JetInputArray JetMomentumSmearing_VLCR10N3/JER_VLCjetsR10N3
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R10N4 {
	set JetInputArray JetMomentumSmearing_VLCR10N4/JER_VLCjetsR10N4
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R10N4 {
	set JetInputArray JetMomentumSmearing_VLCR10N4/JER_VLCjetsR10N4
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R10N4 {
	set JetInputArray JetMomentumSmearing_VLCR10N4/JER_VLCjetsR10N4
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R10N5 {
	set JetInputArray JetMomentumSmearing_VLCR10N5/JER_VLCjetsR10N5
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R10N5 {
	set JetInputArray JetMomentumSmearing_VLCR10N5/JER_VLCjetsR10N5
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R10N5 {
	set JetInputArray JetMomentumSmearing_VLCR10N5/JER_VLCjetsR10N5
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R10N6 {
	set JetInputArray JetMomentumSmearing_VLCR10N6/JER_VLCjetsR10N6
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R10N6 {
	set JetInputArray JetMomentumSmearing_VLCR10N6/JER_VLCjetsR10N6
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R10N6 {
	set JetInputArray JetMomentumSmearing_VLCR10N6/JER_VLCjetsR10N6
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R12N2 {
	set JetInputArray JetMomentumSmearing_VLCR12N2/JER_VLCjetsR12N2
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R12N2 {
	set JetInputArray JetMomentumSmearing_VLCR12N2/JER_VLCjetsR12N2
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R12N2 {
	set JetInputArray JetMomentumSmearing_VLCR12N2/JER_VLCjetsR12N2
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R12N3 {
	set JetInputArray JetMomentumSmearing_VLCR12N3/JER_VLCjetsR12N3
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R12N3 {
	set JetInputArray JetMomentumSmearing_VLCR12N3/JER_VLCjetsR12N3
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R12N3 {
	set JetInputArray JetMomentumSmearing_VLCR12N3/JER_VLCjetsR12N3
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R12N4 {
	set JetInputArray JetMomentumSmearing_VLCR12N4/JER_VLCjetsR12N4
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R12N4 {
	set JetInputArray JetMomentumSmearing_VLCR12N4/JER_VLCjetsR12N4
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R12N4 {
	set JetInputArray JetMomentumSmearing_VLCR12N4/JER_VLCjetsR12N4
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R12N5 {
	set JetInputArray JetMomentumSmearing_VLCR12N5/JER_VLCjetsR12N5
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R12N5 {
	set JetInputArray JetMomentumSmearing_VLCR12N5/JER_VLCjetsR12N5
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R12N5 {
	set JetInputArray JetMomentumSmearing_VLCR12N5/JER_VLCjetsR12N5
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R12N6 {
	set JetInputArray JetMomentumSmearing_VLCR12N6/JER_VLCjetsR12N6
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R12N6 {
	set JetInputArray JetMomentumSmearing_VLCR12N6/JER_VLCjetsR12N6
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R12N6 {
	set JetInputArray JetMomentumSmearing_VLCR12N6/JER_VLCjetsR12N6
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R15N2 {
	set JetInputArray JetMomentumSmearing_VLCR15N2/JER_VLCjetsR15N2
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R15N2 {
	set JetInputArray JetMomentumSmearing_VLCR15N2/JER_VLCjetsR15N2
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R15N2 {
	set JetInputArray JetMomentumSmearing_VLCR15N2/JER_VLCjetsR15N2
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R15N3 {
	set JetInputArray JetMomentumSmearing_VLCR15N3/JER_VLCjetsR15N3
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R15N3 {
	set JetInputArray JetMomentumSmearing_VLCR15N3/JER_VLCjetsR15N3
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R15N3 {
	set JetInputArray JetMomentumSmearing_VLCR15N3/JER_VLCjetsR15N3
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R15N4 {
	set JetInputArray JetMomentumSmearing_VLCR15N4/JER_VLCjetsR15N4
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R15N4 {
	set JetInputArray JetMomentumSmearing_VLCR15N4/JER_VLCjetsR15N4
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R15N4 {
	set JetInputArray JetMomentumSmearing_VLCR15N4/JER_VLCjetsR15N4
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R15N5 {
	set JetInputArray JetMomentumSmearing_VLCR15N5/JER_VLCjetsR15N5
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R15N5 {
	set JetInputArray JetMomentumSmearing_VLCR15N5/JER_VLCjetsR15N5
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R15N5 {
	set JetInputArray JetMomentumSmearing_VLCR15N5/JER_VLCjetsR15N5
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R15N6 {
	set JetInputArray JetMomentumSmearing_VLCR15N6/JER_VLCjetsR15N6
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTagging_JER_WP70_R15N6 {
	set JetInputArray JetMomentumSmearing_VLCR15N6/JER_VLCjetsR15N6
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R15N6 {
	set JetInputArray JetMomentumSmearing_VLCR15N6/JER_VLCjetsR15N6
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}

########################
# inclusive clustering
########################

module BTagging BTagging_JER_WP50_R05_inclusive {                                                   
 set JetInputArray JetMomentumSmearing_VLCR05_inclusive/JER_VLCjetsR05_inclusive                         
 set BitNumber 0                                                                               
 source CLIC/CLICdet_BTag_50.tcl                                                                    
 }                                                                                             
module BTagging BTagging_JER_WP70_R05_inclusive {                                                   
 set JetInputArray JetMomentumSmearing_VLCR05_inclusive/JER_VLCjetsR05_inclusive                         
 set BitNumber 1                                                                               
 source CLIC/CLICdet_BTag_70.tcl                                                                    
}                                                                                              
module BTagging BTagging_JER_WP90_R05_inclusive {                                                   
 set JetInputArray JetMomentumSmearing_VLCR05_inclusive/JER_VLCjetsR05_inclusive                         
 set BitNumber 2                                                                               
 source CLIC/CLICdet_BTag_90.tcl                                                                    
}                                                                                              
module BTagging BTagging_JER_WP50_R07_inclusive {                                                   
 set JetInputArray JetMomentumSmearing_VLCR07_inclusive/JER_VLCjetsR07_inclusive                         
 set BitNumber 0                                                                               
 source CLIC/CLICdet_BTag_50.tcl                                                                    
 }                                                                                             
module BTagging BTagging_JER_WP70_R07_inclusive {                                                   
 set JetInputArray JetMomentumSmearing_VLCR07_inclusive/JER_VLCjetsR07_inclusive                         
 set BitNumber 1                                                                               
 source CLIC/CLICdet_BTag_70.tcl                                                                    
}                                                                                              
module BTagging BTagging_JER_WP90_R07_inclusive {                                                   
 set JetInputArray JetMomentumSmearing_VLCR07_inclusive/JER_VLCjetsR07_inclusive                         
 set BitNumber 2                                                                               
 source CLIC/CLICdet_BTag_90.tcl                                                                    
}                                                                                              
module BTagging BTagging_JER_WP50_R10_inclusive {                                                   
 set JetInputArray JetMomentumSmearing_VLCR10_inclusive/JER_VLCjetsR10_inclusive                         
 set BitNumber 0                                                                               
 source CLIC/CLICdet_BTag_50.tcl                                                                    
 }                                                                                             
module BTagging BTagging_JER_WP70_R10_inclusive {                                                   
 set JetInputArray JetMomentumSmearing_VLCR10_inclusive/JER_VLCjetsR10_inclusive                         
 set BitNumber 1                                                                               
 source CLIC/CLICdet_BTag_70.tcl                                                                    
}
module BTagging BTagging_JER_WP90_R10_inclusive {
 set JetInputArray JetMomentumSmearing_VLCR10_inclusive/JER_VLCjetsR10_inclusive
 set BitNumber 2
 source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R12_inclusive {
 set JetInputArray JetMomentumSmearing_VLCR12_inclusive/JER_VLCjetsR12_inclusive
 set BitNumber 0
 source CLIC/CLICdet_BTag_50.tcl
 }
module BTagging BTagging_JER_WP70_R12_inclusive {
 set JetInputArray JetMomentumSmearing_VLCR12_inclusive/JER_VLCjetsR12_inclusive
 set BitNumber 1
 source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R12_inclusive {
 set JetInputArray JetMomentumSmearing_VLCR12_inclusive/JER_VLCjetsR12_inclusive
 set BitNumber 2
 source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTagging_JER_WP50_R15_inclusive {
 set JetInputArray JetMomentumSmearing_VLCR15_inclusive/JER_VLCjetsR15_inclusive
 set BitNumber 0
 source CLIC/CLICdet_BTag_50.tcl
 }
module BTagging BTagging_JER_WP70_R15_inclusive {
 set JetInputArray JetMomentumSmearing_VLCR15_inclusive/JER_VLCjetsR15_inclusive
 set BitNumber 1
 source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTagging_JER_WP90_R15_inclusive {
 set JetInputArray JetMomentumSmearing_VLCR15_inclusive/JER_VLCjetsR15_inclusive
 set BitNumber 2
 source CLIC/CLICdet_BTag_90.tcl
}
