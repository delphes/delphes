  # misidentification rate (uds)
  # average efficiency: 0.182
  add EfficiencyFormula {0} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.179 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.139 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.136 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.14 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.155 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.211 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.157 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.145 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.151 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.168 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.225 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.168 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.156 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.158 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.173 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.253 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.175 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.155 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.159 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.168 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.337 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.192 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.146 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.131 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.139 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.349 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.273 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.194 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.163 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.162 )
  }
  
  # misidentification rate (c)
  # average efficiency: 0.504
  add EfficiencyFormula {4} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.619 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.639 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.648 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.659 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.67 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.563 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.603 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.632 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.656 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.67 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.5 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.548 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.602 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.642 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.672 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.447 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.459 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.506 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.571 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.616 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.391 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.326 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.328 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.338 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.396 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.376 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.294 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.271 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.227 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.257 )
  }
  
  # b-tagging efficiency
  # average efficiency: 0.801
  add EfficiencyFormula {5} {
    ( abs(eta)<0.867 )*( energy<30 )*( 0.908 )+ 
    ( abs(eta)<0.867 )*( 30<=energy && energy<60 )*( 0.948 )+ 
    ( abs(eta)<0.867 )*( 60<=energy && energy<100 )*( 0.971 )+ 
    ( abs(eta)<0.867 )*( 100<=energy && energy<150 )*( 0.982 )+ 
    ( abs(eta)<0.867 )*( 150<=energy )*( 0.988 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( energy<30 )*( 0.85 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 30<=energy && energy<60 )*( 0.927 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 60<=energy && energy<100 )*( 0.96 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 100<=energy && energy<150 )*( 0.978 )+ 
    ( 0.867<=abs(eta) && abs(eta)<1.47 )*( 150<=energy )*( 0.982 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( energy<30 )*( 0.79 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 30<=energy && energy<60 )*( 0.88 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 60<=energy && energy<100 )*( 0.935 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 100<=energy && energy<150 )*( 0.963 )+ 
    ( 1.47<=abs(eta) && abs(eta)<1.83 )*( 150<=energy )*( 0.979 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( energy<30 )*( 0.708 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 30<=energy && energy<60 )*( 0.802 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 60<=energy && energy<100 )*( 0.882 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 100<=energy && energy<150 )*( 0.927 )+ 
    ( 1.83<=abs(eta) && abs(eta)<2.3 )*( 150<=energy )*( 0.944 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( energy<30 )*( 0.59 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 30<=energy && energy<60 )*( 0.639 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 60<=energy && energy<100 )*( 0.668 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 100<=energy && energy<150 )*( 0.721 )+ 
    ( 2.3<=abs(eta) && abs(eta)<2.99 )*( 150<=energy )*( 0.756 )+ 
    ( 2.99<=abs(eta) )*( energy<30 )*( 0.411 )+ 
    ( 2.99<=abs(eta) )*( 30<=energy && energy<60 )*( 0.477 )+ 
    ( 2.99<=abs(eta) )*( 60<=energy && energy<100 )*( 0.472 )+ 
    ( 2.99<=abs(eta) )*( 100<=energy && energy<150 )*( 0.5 )+ 
    ( 2.99<=abs(eta) )*( 150<=energy )*( 0.488 )
  }
