############
# Jet finder VLC
############
#R05 N2
module FastJetFinder FastJetFinderVLC_R05_N2 {
    #  set InputArray Calorimeter/towers
    set InputArray EFlowFilter/eflow

    set OutputArray VLCjetsR05N2

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt, 7 anti-kt with winner-take-all axis (for N-subjettiness), 8 N-jettiness, 9 Valencia
    set NJets 2
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 0.5
    set Beta 1.0
    set Gamma 1.0

    set JetPTMin 20.0
}
#R05 N3
module FastJetFinder FastJetFinderVLC_R05_N3 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR05N3
    set NJets 3
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 0.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R05 N4
module FastJetFinder FastJetFinderVLC_R05_N4 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR05N4
    set NJets 4
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 0.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R05 N5
module FastJetFinder FastJetFinderVLC_R05_N5 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR05N5
    set NJets 5
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 0.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R05 N6
module FastJetFinder FastJetFinderVLC_R05_N6 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR05N6
    set NJets 6
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 0.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R07 N2
module FastJetFinder FastJetFinderVLC_R07_N2 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR07N2
    set NJets 2
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 0.7
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R07 N3
module FastJetFinder FastJetFinderVLC_R07_N3 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR07N3
    set NJets 3
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 0.7
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R07 N4
module FastJetFinder FastJetFinderVLC_R07_N4 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR07N4
    set NJets 4
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 0.7
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R07 N5
module FastJetFinder FastJetFinderVLC_R07_N5 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR07N5
    set NJets 5
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 0.7
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R07 N6
module FastJetFinder FastJetFinderVLC_R07_N6 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR07N6
    set NJets 6
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 0.7
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}

#R10N2
module FastJetFinder FastJetFinderVLC_R10_N2 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR10N2
    set NJets 2
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.0
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R10 N3
module FastJetFinder FastJetFinderVLC_R10_N3 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR10N3
    set NJets 3
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.0
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R10 N4
module FastJetFinder FastJetFinderVLC_R10_N4 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR10N4
    set NJets 4
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.0
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R10 N5
module FastJetFinder FastJetFinderVLC_R10_N5 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR10N5
    set NJets 5
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.0
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R10 N6
module FastJetFinder FastJetFinderVLC_R10_N6 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR10N6
    set NJets 6
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.0
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}

#R12 N2
module FastJetFinder FastJetFinderVLC_R12_N2 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR12N2
    set NJets 2
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.2
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R12 N3
module FastJetFinder FastJetFinderVLC_R12_N3 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR12N3
    set NJets 3
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.2
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R12 N4
module FastJetFinder FastJetFinderVLC_R12_N4 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR12N4
    set NJets 4
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.2
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R12 N5
module FastJetFinder FastJetFinderVLC_R12_N5 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR12N5
    set NJets 5
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.2
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R12 N6
module FastJetFinder FastJetFinderVLC_R12_N6 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR12N6
    set NJets 6
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.2
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}


#R15 N2
module FastJetFinder FastJetFinderVLC_R15_N2 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR15N2
    set NJets 2
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R15 N3
module FastJetFinder FastJetFinderVLC_R15_N3 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR15N3
    set NJets 3
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R15 N4
module FastJetFinder FastJetFinderVLC_R15_N4 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR15N4
    set NJets 4
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R15 N5
module FastJetFinder FastJetFinderVLC_R15_N5 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR15N5
    set NJets 5
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R15 N6
module FastJetFinder FastJetFinderVLC_R15_N6 {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR15N6
    set NJets 6
    set ExclusiveClustering true
    set JetAlgorithm 9
    set ParameterR 1.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}

########################
# inclusive clustering
# as means of comparison
########################
#R05
module FastJetFinder FastJetFinderVLC_R05_inclusive {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR05_inclusive
    set ExclusiveClustering false
    set JetAlgorithm 9
    set ParameterR 0.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R07
module FastJetFinder FastJetFinderVLC_R07_inclusive {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR07_inclusive
    set ExclusiveClustering false
    set JetAlgorithm 9
    set ParameterR 0.7
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R10
module FastJetFinder FastJetFinderVLC_R10_inclusive {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR10_inclusive
    set ExclusiveClustering false
    set JetAlgorithm 9
    set ParameterR 1.0
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R12
module FastJetFinder FastJetFinderVLC_R12_inclusive {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR12_inclusive
    set ExclusiveClustering false
    set JetAlgorithm 9
    set ParameterR 1.2
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
#R15
module FastJetFinder FastJetFinderVLC_R15_inclusive {
    set InputArray EFlowFilter/eflow
    set OutputArray VLCjetsR15_inclusive
    set ExclusiveClustering false
    set JetAlgorithm 9
    set ParameterR 1.5
    set Beta 1.0
    set Gamma 1.0
    set JetPTMin 20.0
}
