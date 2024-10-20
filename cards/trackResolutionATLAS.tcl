set PResolutionFormula { 0.0 }
# ctg theta and phi resolution from atlas ID TDR
# pT < 5 based on 1 GeV line
# pT > 15 based on 20 GeV line
# pT between is linear interpolation with w=0.3 from 1 GeV -> 20 GeV
# interpolation values chosen based on sqrt(s) = 13 TeV ttbar tracks, which are largely 
# low pT
set CtgThetaResolutionFormula { 
    ( abs(eta) > 0.0 && abs(eta) <= 0.25 ) * ( pt <= 5 ) * 0.001753554502369671 +\
    ( abs(eta) > 0.25 && abs(eta) <= 0.5 ) * ( pt <= 5 ) * 0.0018483412322274906 +\
    ( abs(eta) > 0.5 && abs(eta) <= 0.75 ) * ( pt <= 5 ) * 0.0020853080568720393 +\
    ( abs(eta) > 0.75 && abs(eta) <= 1.0 ) * ( pt <= 5 ) * 0.0024644549763033208 +\
    ( abs(eta) > 1.0 && abs(eta) <= 1.25 ) * ( pt <= 5 ) * 0.0030331753554502378 +\
    ( abs(eta) > 1.25 && abs(eta) <= 1.5 ) * ( pt <= 5 ) * 0.004360189573459721 +\
    ( abs(eta) > 1.5 && abs(eta) <= 1.75 ) * ( pt <= 5 ) * 0.006303317535545024 +\
    ( abs(eta) > 1.75 && abs(eta) <= 2.0 ) * ( pt <= 5 ) * 0.009099526066350713 +\
    ( abs(eta) > 2.0 && abs(eta) <= 2.25 ) * ( pt <= 5 ) * 0.013933649289099528 +\
    ( abs(eta) > 2.25 && abs(eta) <= 2.5 ) * ( pt <= 5 ) * 0.019478672985781993 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.25 ) * ( pt > 5 && pt <= 15 ) * 0.001447867298578201 +\
    ( abs(eta) > 0.25 && abs(eta) <= 0.5 ) * ( pt > 5 && pt <= 15 ) * 0.0014957345971564 +\
    ( abs(eta) > 0.5 && abs(eta) <= 0.75 ) * ( pt > 5 && pt <= 15 ) * 0.0016402843601895744 +\
    ( abs(eta) > 0.75 && abs(eta) <= 1.0 ) * ( pt > 5 && pt <= 15 ) * 0.001890047393364931 +\
    ( abs(eta) > 1.0 && abs(eta) <= 1.25 ) * ( pt > 5 && pt <= 15 ) * 0.0022838862559241713 +\
    ( abs(eta) > 1.25 && abs(eta) <= 1.5 ) * ( pt > 5 && pt <= 15 ) * 0.003246919431279625 +\
    ( abs(eta) > 1.5 && abs(eta) <= 1.75 ) * ( pt > 5 && pt <= 15 ) * 0.004690995260663507 +\
    ( abs(eta) > 1.75 && abs(eta) <= 2.0 ) * ( pt > 5 && pt <= 15 ) * 0.006718009478672987 +\
    ( abs(eta) > 2.0 && abs(eta) <= 2.25 ) * ( pt > 5 && pt <= 15 ) * 0.010293838862559242 +\
    ( abs(eta) > 2.25 && abs(eta) <= 2.5 ) * ( pt > 5 && pt <= 15 ) * 0.014371563981042657+\
    ( abs(eta) > 0.0 && abs(eta) <= 0.25 ) * ( pt > 15 ) * 0.0007345971563981041 +\
    ( abs(eta) > 0.25 && abs(eta) <= 0.5 ) * ( pt > 15 ) * 0.0006729857819905215 +\
    ( abs(eta) > 0.5 && abs(eta) <= 0.75 ) * ( pt > 15 ) * 0.0006018957345971564 +\
    ( abs(eta) > 0.75 && abs(eta) <= 1.0 ) * ( pt > 15 ) * 0.0005497630331753552 +\
    ( abs(eta) > 1.0 && abs(eta) <= 1.25 ) * ( pt > 15 ) * 0.0005355450236966823 +\
    ( abs(eta) > 1.25 && abs(eta) <= 1.5 ) * ( pt > 15 ) * 0.0006492890995260666 +\
    ( abs(eta) > 1.5 && abs(eta) <= 1.75 ) * ( pt > 15 ) * 0.000928909952606635 +\
    ( abs(eta) > 1.75 && abs(eta) <= 2.0 ) * ( pt > 15 ) * 0.0011611374407582938 +\
    ( abs(eta) > 2.0 && abs(eta) <= 2.25 ) * ( pt > 15 ) * 0.0018009478672985782 +\
    ( abs(eta) > 2.25 && abs(eta) <= 2.5 ) * ( pt > 15 ) * 0.002454976303317536 
}


set PhiResolutionFormula { 
    ( abs(eta) > 0.0 && abs(eta) <= 0.25 ) * ( pt <= 5 ) * 0.0015096525096525092 +\
    ( abs(eta) > 0.25 && abs(eta) <= 0.5 ) * ( pt <= 5 ) * 0.001625482625482625 +\
    ( abs(eta) > 0.5 && abs(eta) <= 0.75 ) * ( pt <= 5 ) * 0.0017181467181467181 +\
    ( abs(eta) > 0.75 && abs(eta) <= 1.0 ) * ( pt <= 5 ) * 0.0017876447876447872 +\
    ( abs(eta) > 1.0 && abs(eta) <= 1.25 ) * ( pt <= 5 ) * 0.0021003861003861 +\
    ( abs(eta) > 1.25 && abs(eta) <= 1.5 ) * ( pt <= 5 ) * 0.0022625482625482623 +\
    ( abs(eta) > 1.5 && abs(eta) <= 1.75 ) * ( pt <= 5 ) * 0.0025868725868725866 +\
    ( abs(eta) > 1.75 && abs(eta) <= 2.0 ) * ( pt <= 5 ) * 0.002934362934362934 +\
    ( abs(eta) > 2.0 && abs(eta) <= 2.25 ) * ( pt <= 5 ) * 0.0031428571428571426 +\
    ( abs(eta) > 2.25 && abs(eta) <= 2.5 ) * ( pt <= 5 ) * 0.003722007722007722+\
    
    ( abs(eta) > 0.0 && abs(eta) <= 0.25 ) * ( pt > 5 && pt <= 15 ) * 0.001100647431997914 +\
    ( abs(eta) > 0.25 && abs(eta) <= 0.5 ) * ( pt > 5 && pt <= 15 ) * 0.001184333014686712 +\
    ( abs(eta) > 0.5 && abs(eta) <= 0.75 ) * ( pt > 5 && pt <= 15 ) * 0.0012486191014165291 +\
    ( abs(eta) > 0.75 && abs(eta) <= 1.0 ) * ( pt > 5 && pt <= 15 ) * 0.001299004084470322 +\
    ( abs(eta) > 1.0 && abs(eta) <= 1.25 ) * ( pt > 5 && pt <= 15 ) * 0.0015234213956721995 +\
    ( abs(eta) > 1.25 && abs(eta) <= 1.5 ) * ( pt > 5 && pt <= 15 ) * 0.0016395394107934298 +\
    ( abs(eta) > 1.5 && abs(eta) <= 1.75 ) * ( pt > 5 && pt <= 15 ) * 0.0018697497175632221 +\
    ( abs(eta) > 1.75 && abs(eta) <= 2.0 ) * ( pt > 5 && pt <= 15 ) * 0.0021135717389415137 +\
    ( abs(eta) > 2.0 && abs(eta) <= 2.25 ) * ( pt > 5 && pt <= 15 ) * 0.002260385852090032 +\
    ( abs(eta) > 2.25 && abs(eta) <= 2.5 ) * ( pt > 5 && pt <= 15 ) * 0.002673604762318589+\

    ( abs(eta) > 0.0 && abs(eta) <= 0.25 ) * ( pt > 15 ) * 0.00014630225080385849 +\
    ( abs(eta) > 0.25 && abs(eta) <= 0.5 ) * ( pt > 15 ) * 0.000154983922829582 +\
    ( abs(eta) > 0.5 && abs(eta) <= 0.75 ) * ( pt > 15 ) * 0.00015305466237942118 +\
    ( abs(eta) > 0.75 && abs(eta) <= 1.0 ) * ( pt > 15 ) * 0.00015884244372990354 +\
    ( abs(eta) > 1.0 && abs(eta) <= 1.25 ) * ( pt > 15 ) * 0.0001771704180064309 +\
    ( abs(eta) > 1.25 && abs(eta) <= 1.5 ) * ( pt > 15 ) * 0.00018585209003215433 +\
    ( abs(eta) > 1.5 && abs(eta) <= 1.75 ) * ( pt > 15 ) * 0.00019646302250803858 +\
    ( abs(eta) > 1.75 && abs(eta) <= 2.0 ) * ( pt > 15 ) * 0.00019839228295819939 +\
    ( abs(eta) > 2.0 && abs(eta) <= 2.25 ) * ( pt > 15 ) * 0.00020128617363344054 +\
    ( abs(eta) > 2.25 && abs(eta) <= 2.5 ) * ( pt > 15 ) * 0.00022733118971061093

 }
# taken from arXiv:1405.6569 fig. 15
set D0ResolutionFormula {
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 0.5 && pt <= 1.0) * 0.08 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 1.0 && pt <= 2.0 ) * 0.045 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 2.0 && pt <= 4.0 ) * 0.028 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 4.0 && pt <= 6.0 ) * 0.020 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 6.0 && pt <= 10.0 ) * 0.016 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 10.0 && pt <= 30.0 ) * 0.014 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 30.0 ) * 0.0115
}
set DZResolutionFormula {
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 0.5 && pt <= 1.0) * 0.140 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 1.0 && pt <= 2.0 ) * 0.096 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 2.0 && pt <= 4.0 ) * 0.075 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 4.0 && pt <= 6.0 ) * 0.058 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 6.0 && pt <= 15.0 ) * 0.055 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 15.0 && pt <= 30.0 ) * 0.050 +\
    ( abs(eta) > 0.0 && abs(eta) <= 2.5 ) * ( pt > 30.0 ) * 0.047
}
