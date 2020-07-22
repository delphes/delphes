  # misidentification rate (uds)
  # average efficiency: 0.00583
  add EfficiencyFormula {0} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.00885 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.00773 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.00734 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.00757 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.00884 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.00835 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.00764 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.00768 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.00768 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.00878 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.00703 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.00554 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.00637 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.00662 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.0079 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.00491 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.00465 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.00463 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.00421 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.0064 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.0057 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.00348 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.00316 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.00263 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.00401 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.00289 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.0054 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.00384 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.00187 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.00327 )
  }
  
  # misidentification rate (c)
  # average efficiency: 0.0544
  add EfficiencyFormula {4} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.0627 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.074 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.0869 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.096 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.119 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.054 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.0692 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.0835 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.0952 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.108 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.0369 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.0605 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.0745 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.0815 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.108 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.027 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.0387 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.0561 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.0653 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.0775 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.0125 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.0197 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.0195 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.026 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.0356 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.00479 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.00805 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.00797 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.0106 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.0144 )
  }
  
  # b-tagging efficiency
  # average efficiency: 0.529
  add EfficiencyFormula {5} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.682 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.803 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.852 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.887 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.899 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.568 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.743 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.816 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.854 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.881 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.432 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.632 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.748 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.801 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.851 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.28 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.443 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.563 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.643 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.717 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.143 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.206 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.226 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.281 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.325 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.0336 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.0865 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.152 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.166 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.162 )
  }
