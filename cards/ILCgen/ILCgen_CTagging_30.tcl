  # misidentification rate (uds)
  # average efficiency: 0.0168
  add EfficiencyFormula {0} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.0294 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.0285 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.0286 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.0301 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.0311 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.0262 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.0239 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.0235 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.0244 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.0273 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.0197 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.0165 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.0183 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.0216 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.0221 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.0103 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.0125 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.0116 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.0138 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.0167 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.00848 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.00687 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.00718 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.00868 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.00706 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.00685 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.00459 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.00485 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.00657 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.0083 )
  }
  
  # c-tagging efficiency
  # average efficiency: 0.315
  add EfficiencyFormula {4} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.497 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.598 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.652 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.686 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.702 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.351 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.48 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.565 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.614 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.649 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.194 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.321 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.432 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.517 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.567 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.0942 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.161 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.239 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.313 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.37 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.0216 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.0524 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.0616 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.075 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.0988 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.00502 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.0135 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.0276 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.0428 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.0402 )
  }
  
  # misidentification rate (b)
  # average efficiency: 0.105
  add EfficiencyFormula {5} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.212 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.143 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.112 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.0862 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.0749 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.215 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.148 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.117 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.0992 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.0827 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.182 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.161 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.119 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.107 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.0884 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.137 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.136 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.133 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.108 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.101 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.0688 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.0724 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.0761 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.0751 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.064 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.0264 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.0337 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.0504 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.0501 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.056 )
  }
