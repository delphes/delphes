##################################
# Jet finder in inclusive mode N=2
##################################
module FastJetFinder JetFinder_N2 {
    set InputArray EFlowFilter/eflow
    set OutputArray jets
    set JetAlgorithm 9

    set NJets 2
    set ExclusiveClustering true

    set ParameterR 2.0
    set Beta 1.0
    set Gamma 0.0
    set JetPTMin 0.0
}
##################################
# Jet finder in inclusive mode N=3
##################################
module FastJetFinder JetFinder_N3 {
    set InputArray EFlowFilter/eflow
    set OutputArray jets
    set JetAlgorithm 9

    set NJets 3
    set ExclusiveClustering true

    set ParameterR 2.0
    set Beta 1.0
    set Gamma 0.0
    set JetPTMin 0.0
}
##################################
# Jet finder in inclusive mode N=4
##################################
module FastJetFinder JetFinder_N4 {
    set InputArray EFlowFilter/eflow
    set OutputArray jets
    set JetAlgorithm 9

    set NJets 4
    set ExclusiveClustering true

    set ParameterR 2.0
    set Beta 1.0
    set Gamma 0.0
    set JetPTMin 0.0
}
##################################
# Jet finder in inclusive mode N=5
##################################
module FastJetFinder JetFinder_N5 {
    set InputArray EFlowFilter/eflow
    set OutputArray jets
    set JetAlgorithm 9

    set NJets 5
    set ExclusiveClustering true

    set ParameterR 2.0
    set Beta 1.0
    set Gamma 0.0
    set JetPTMin 0.0
}
##################################
# Jet finder in inclusive mode N=6
##################################
module FastJetFinder JetFinder_N6 {
    set InputArray EFlowFilter/eflow
    set OutputArray jets
    set JetAlgorithm 9

    set NJets 6
    set ExclusiveClustering true

    set ParameterR 2.0
    set Beta 1.0
    set Gamma 0.0
    set JetPTMin 0.0
}
