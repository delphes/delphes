  # misidentification rate (uds)
  # average efficiency: 0.00182
  add EfficiencyFormula {0} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.0038 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.00379 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.00361 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.00394 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.00424 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.00388 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.00269 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.00324 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.00291 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.00269 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.00182 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.00262 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.00213 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.00276 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.00227 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 3.89e-06 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.00135 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.000845 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.00191 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.00269 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.000414 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.000177 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.000379 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.000363 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.000169 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( -1.89e-06 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( -1.86e-05 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 3.06e-08 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( -2.29e-07 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 8.83e-07 )
  }
  
  # c-tagging efficiency
  # average efficiency: 0.198
  add EfficiencyFormula {4} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.306 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.4 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.451 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.485 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.491 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.198 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.301 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.377 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.421 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.449 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.0907 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.18 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.272 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.347 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.385 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.0387 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.0714 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.125 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.186 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.23 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.00354 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.0139 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.0165 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.025 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.0392 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.000778 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.00199 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.00779 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.00916 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.00666 )
  }
  
  # misidentification rate (b)
  # average efficiency: 0.0283
  add EfficiencyFormula {5} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.0861 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.0538 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.0404 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.0321 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.0258 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.0779 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.0514 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.0393 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.032 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.0277 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.0419 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.0486 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.0371 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.0331 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.0275 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.0247 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.0282 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.0319 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.0297 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.031 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.0065 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.00828 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.00793 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.00827 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.00596 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.0045 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.00185 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.00299 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.00243 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.0001 )
  }
