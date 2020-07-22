  # misidentification rate (uds)
  # average efficiency: 0.165
  add EfficiencyFormula {0} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.193 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.148 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.139 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.14 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.147 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.216 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.159 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.146 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.147 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.156 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.222 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.159 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.153 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.154 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.165 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.226 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.163 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.148 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.155 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.162 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.266 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.162 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.138 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.13 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.14 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.223 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.19 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.139 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.124 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.142 )
  }
  
  # c-tagging efficiency
  # average efficiency: 0.547
  add EfficiencyFormula {4} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.704 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.766 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.811 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.842 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.858 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.604 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.69 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.761 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.802 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.826 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.505 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.579 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.673 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.742 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.786 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.436 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.46 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.525 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.606 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.66 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.328 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.304 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.307 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.33 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.395 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.251 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.209 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.21 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.203 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.233 )
  }
  
  # misidentification rate (b)
  # average efficiency: 0.343
  add EfficiencyFormula {5} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.358 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.258 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.217 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.184 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.175 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.41 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.298 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.248 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.218 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.184 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.465 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.36 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.3 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.261 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.222 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.49 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.454 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.437 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.393 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.342 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.441 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.44 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.484 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.497 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.528 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.279 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.337 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.305 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.354 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.337 )
  }
