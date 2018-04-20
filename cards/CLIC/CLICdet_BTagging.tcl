module BTagging BTaggingWP50_R05N2 {
	set JetInputArray FastJetFinderVLC_R05_N2/VLCjetsR05N2
	set BitNumber 0

	# 50% efficiency working point
    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
    # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
    # gluon's PDG code has the lowest priority

	# based on CLICdp-Note-2014-002

	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R05N2 {
	set JetInputArray FastJetFinderVLC_R05_N2/VLCjetsR05N2
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R05N2 {
	set JetInputArray FastJetFinderVLC_R05_N2/VLCjetsR05N2
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R05N3 {
	set JetInputArray FastJetFinderVLC_R05_N3/VLCjetsR05N3
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R05N3 {
	set JetInputArray FastJetFinderVLC_R05_N3/VLCjetsR05N3
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R05N3 {
	set JetInputArray FastJetFinderVLC_R05_N3/VLCjetsR05N3
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R05N4 {
	set JetInputArray FastJetFinderVLC_R05_N4/VLCjetsR05N4
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R05N4 {
	set JetInputArray FastJetFinderVLC_R05_N4/VLCjetsR05N4
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R05N4 {
	set JetInputArray FastJetFinderVLC_R05_N4/VLCjetsR05N4
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R05N5 {
	set JetInputArray FastJetFinderVLC_R05_N5/VLCjetsR05N5
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R05N5 {
	set JetInputArray FastJetFinderVLC_R05_N5/VLCjetsR05N5
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R05N5 {
	set JetInputArray FastJetFinderVLC_R05_N5/VLCjetsR05N5
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R05N6 {
	set JetInputArray FastJetFinderVLC_R05_N6/VLCjetsR05N6
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R05N6 {
	set JetInputArray FastJetFinderVLC_R05_N6/VLCjetsR05N6
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R05N6 {
	set JetInputArray FastJetFinderVLC_R05_N6/VLCjetsR05N6
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R07N2 {
	set JetInputArray FastJetFinderVLC_R07_N2/VLCjetsR07N2
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R07N2 {
	set JetInputArray FastJetFinderVLC_R07_N2/VLCjetsR07N2
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R07N2 {
	set JetInputArray FastJetFinderVLC_R07_N2/VLCjetsR07N2
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R07N3 {
	set JetInputArray FastJetFinderVLC_R07_N3/VLCjetsR07N3
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R07N3 {
	set JetInputArray FastJetFinderVLC_R07_N3/VLCjetsR07N3
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R07N3 {
	set JetInputArray FastJetFinderVLC_R07_N3/VLCjetsR07N3
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R07N4 {
	set JetInputArray FastJetFinderVLC_R07_N4/VLCjetsR07N4
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R07N4 {
	set JetInputArray FastJetFinderVLC_R07_N4/VLCjetsR07N4
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R07N4 {
	set JetInputArray FastJetFinderVLC_R07_N4/VLCjetsR07N4
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R07N5 {
	set JetInputArray FastJetFinderVLC_R07_N5/VLCjetsR07N5
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R07N5 {
	set JetInputArray FastJetFinderVLC_R07_N5/VLCjetsR07N5
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R07N5 {
	set JetInputArray FastJetFinderVLC_R07_N5/VLCjetsR07N5
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R07N6 {
	set JetInputArray FastJetFinderVLC_R07_N6/VLCjetsR07N6
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R07N6 {
	set JetInputArray FastJetFinderVLC_R07_N6/VLCjetsR07N6
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R07N6 {
	set JetInputArray FastJetFinderVLC_R07_N6/VLCjetsR07N6
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R10N2 {
	set JetInputArray FastJetFinderVLC_R10_N2/VLCjetsR10N2
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R10N2 {
	set JetInputArray FastJetFinderVLC_R10_N2/VLCjetsR10N2
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R10N2 {
	set JetInputArray FastJetFinderVLC_R10_N2/VLCjetsR10N2
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R10N3 {
	set JetInputArray FastJetFinderVLC_R10_N3/VLCjetsR10N3
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R10N3 {
	set JetInputArray FastJetFinderVLC_R10_N3/VLCjetsR10N3
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R10N3 {
	set JetInputArray FastJetFinderVLC_R10_N3/VLCjetsR10N3
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R10N4 {
	set JetInputArray FastJetFinderVLC_R10_N4/VLCjetsR10N4
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R10N4 {
	set JetInputArray FastJetFinderVLC_R10_N4/VLCjetsR10N4
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R10N4 {
	set JetInputArray FastJetFinderVLC_R10_N4/VLCjetsR10N4
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R10N5 {
	set JetInputArray FastJetFinderVLC_R10_N5/VLCjetsR10N5
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R10N5 {
	set JetInputArray FastJetFinderVLC_R10_N5/VLCjetsR10N5
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R10N5 {
	set JetInputArray FastJetFinderVLC_R10_N5/VLCjetsR10N5
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R10N6 {
	set JetInputArray FastJetFinderVLC_R10_N6/VLCjetsR10N6
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R10N6 {
	set JetInputArray FastJetFinderVLC_R10_N6/VLCjetsR10N6
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R10N6 {
	set JetInputArray FastJetFinderVLC_R10_N6/VLCjetsR10N6
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R12N2 {
	set JetInputArray FastJetFinderVLC_R12_N2/VLCjetsR12N2
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R12N2 {
	set JetInputArray FastJetFinderVLC_R12_N2/VLCjetsR12N2
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R12N2 {
	set JetInputArray FastJetFinderVLC_R12_N2/VLCjetsR12N2
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R12N3 {
	set JetInputArray FastJetFinderVLC_R12_N3/VLCjetsR12N3
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R12N3 {
	set JetInputArray FastJetFinderVLC_R12_N3/VLCjetsR12N3
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R12N3 {
	set JetInputArray FastJetFinderVLC_R12_N3/VLCjetsR12N3
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R12N4 {
	set JetInputArray FastJetFinderVLC_R12_N4/VLCjetsR12N4
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R12N4 {
	set JetInputArray FastJetFinderVLC_R12_N4/VLCjetsR12N4
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R12N4 {
	set JetInputArray FastJetFinderVLC_R12_N4/VLCjetsR12N4
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R12N5 {
	set JetInputArray FastJetFinderVLC_R12_N5/VLCjetsR12N5
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R12N5 {
	set JetInputArray FastJetFinderVLC_R12_N5/VLCjetsR12N5
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R12N5 {
	set JetInputArray FastJetFinderVLC_R12_N5/VLCjetsR12N5
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R12N6 {
	set JetInputArray FastJetFinderVLC_R12_N6/VLCjetsR12N6
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R12N6 {
	set JetInputArray FastJetFinderVLC_R12_N6/VLCjetsR12N6
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R12N6 {
	set JetInputArray FastJetFinderVLC_R12_N6/VLCjetsR12N6
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R15N2 {
	set JetInputArray FastJetFinderVLC_R15_N2/VLCjetsR15N2
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R15N2 {
	set JetInputArray FastJetFinderVLC_R15_N2/VLCjetsR15N2
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R15N2 {
	set JetInputArray FastJetFinderVLC_R15_N2/VLCjetsR15N2
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R15N3 {
	set JetInputArray FastJetFinderVLC_R15_N3/VLCjetsR15N3
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R15N3 {
	set JetInputArray FastJetFinderVLC_R15_N3/VLCjetsR15N3
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R15N3 {
	set JetInputArray FastJetFinderVLC_R15_N3/VLCjetsR15N3
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R15N4 {
	set JetInputArray FastJetFinderVLC_R15_N4/VLCjetsR15N4
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R15N4 {
	set JetInputArray FastJetFinderVLC_R15_N4/VLCjetsR15N4
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R15N4 {
	set JetInputArray FastJetFinderVLC_R15_N4/VLCjetsR15N4
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R15N5 {
	set JetInputArray FastJetFinderVLC_R15_N5/VLCjetsR15N5
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R15N5 {
	set JetInputArray FastJetFinderVLC_R15_N5/VLCjetsR15N5
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R15N5 {
	set JetInputArray FastJetFinderVLC_R15_N5/VLCjetsR15N5
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R15N6 {
	set JetInputArray FastJetFinderVLC_R15_N6/VLCjetsR15N6
	set BitNumber 0
	source CLIC/CLICdet_BTag_50.tcl
}
module BTagging BTaggingWP70_R15N6 {
	set JetInputArray FastJetFinderVLC_R15_N6/VLCjetsR15N6
	set BitNumber 1
	source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R15N6 {
	set JetInputArray FastJetFinderVLC_R15_N6/VLCjetsR15N6
	set BitNumber 2
	source CLIC/CLICdet_BTag_90.tcl
}

########################
# inclusive clustering
########################

module BTagging BTaggingWP50_R05_inclusive {                                                   
 set JetInputArray FastJetFinderVLC_R05_inclusive/VLCjetsR05_inclusive                         
 set BitNumber 0                                                                               
 source CLIC/CLICdet_BTag_50.tcl                                                                    
 }                                                                                             
module BTagging BTaggingWP70_R05_inclusive {                                                   
 set JetInputArray FastJetFinderVLC_R05_inclusive/VLCjetsR05_inclusive                         
 set BitNumber 1                                                                               
 source CLIC/CLICdet_BTag_70.tcl                                                                    
}                                                                                              
module BTagging BTaggingWP90_R05_inclusive {                                                   
 set JetInputArray FastJetFinderVLC_R05_inclusive/VLCjetsR05_inclusive                         
 set BitNumber 2                                                                               
 source CLIC/CLICdet_BTag_90.tcl                                                                    
}                                                                                              
module BTagging BTaggingWP50_R07_inclusive {                                                   
 set JetInputArray FastJetFinderVLC_R07_inclusive/VLCjetsR07_inclusive                         
 set BitNumber 0                                                                               
 source CLIC/CLICdet_BTag_50.tcl                                                                    
 }                                                                                             
module BTagging BTaggingWP70_R07_inclusive {                                                   
 set JetInputArray FastJetFinderVLC_R07_inclusive/VLCjetsR07_inclusive                         
 set BitNumber 1                                                                               
 source CLIC/CLICdet_BTag_70.tcl                                                                    
}                                                                                              
module BTagging BTaggingWP90_R07_inclusive {                                                   
 set JetInputArray FastJetFinderVLC_R07_inclusive/VLCjetsR07_inclusive                         
 set BitNumber 2                                                                               
 source CLIC/CLICdet_BTag_90.tcl                                                                    
}                                                                                              
module BTagging BTaggingWP50_R10_inclusive {                                                   
 set JetInputArray FastJetFinderVLC_R10_inclusive/VLCjetsR10_inclusive                         
 set BitNumber 0                                                                               
 source CLIC/CLICdet_BTag_50.tcl                                                                    
 }                                                                                             
module BTagging BTaggingWP70_R10_inclusive {                                                   
 set JetInputArray FastJetFinderVLC_R10_inclusive/VLCjetsR10_inclusive                         
 set BitNumber 1                                                                               
 source CLIC/CLICdet_BTag_70.tcl                                                                    
}
module BTagging BTaggingWP90_R10_inclusive {
 set JetInputArray FastJetFinderVLC_R10_inclusive/VLCjetsR10_inclusive
 set BitNumber 2
 source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R12_inclusive {
 set JetInputArray FastJetFinderVLC_R12_inclusive/VLCjetsR12_inclusive
 set BitNumber 0
 source CLIC/CLICdet_BTag_50.tcl
 }
module BTagging BTaggingWP70_R12_inclusive {
 set JetInputArray FastJetFinderVLC_R12_inclusive/VLCjetsR12_inclusive
 set BitNumber 1
 source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R12_inclusive {
 set JetInputArray FastJetFinderVLC_R12_inclusive/VLCjetsR12_inclusive
 set BitNumber 2
 source CLIC/CLICdet_BTag_90.tcl
}
module BTagging BTaggingWP50_R15_inclusive {
 set JetInputArray FastJetFinderVLC_R15_inclusive/VLCjetsR15_inclusive
 set BitNumber 0
 source CLIC/CLICdet_BTag_50.tcl
 }
module BTagging BTaggingWP70_R15_inclusive {
 set JetInputArray FastJetFinderVLC_R15_inclusive/VLCjetsR15_inclusive
 set BitNumber 1
 source CLIC/CLICdet_BTag_70.tcl
}
module BTagging BTaggingWP90_R15_inclusive {
 set JetInputArray FastJetFinderVLC_R15_inclusive/VLCjetsR15_inclusive
 set BitNumber 2
 source CLIC/CLICdet_BTag_90.tcl
}
