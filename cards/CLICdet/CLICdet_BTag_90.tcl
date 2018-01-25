# efficiency formula for b-jets
add EfficiencyFormula {5} {0.9}
# default efficiency formula (misidentification rate)
add EfficiencyFormula {0} {
	(energy >= 500 )*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-1 )+ \
		(energy >= 500 )*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 5e-2 )+ \
		(energy >= 500 )*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 4e-2 )+ \
		(energy >= 500 )*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 5e-2 )+ \
		(energy >= 500 )*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 6e-2 )+ \
		(energy >= 500 )*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 7e-2 )+ \
		(energy >= 500 )*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 7e-2 )+ \
		(energy >= 500 )*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 8e-2 )+ \
		(energy >= 500 )*( abs(eta) <=0.09                     ) * ( 8e-2 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-1 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-2 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-2 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 3e-2 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 2e-2 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 2e-2 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 2e-2 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 2e-2 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=0.09                     ) * ( 2e-2 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-1 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 5e-2 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-2 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 2e-2 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 2e-2 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 1e-2 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 1e-2 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 1e-2 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=0.09                     ) * ( 1e-2 )+ \
		(energy < 100 )*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-1 )+ \
		(energy < 100 )*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 8e-2 )+ \
		(energy < 100 )*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 4e-2 )+ \
		(energy < 100 )*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 3e-2 )+ \
		(energy < 100 )*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 2e-2 )+ \
		(energy < 100 )*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 1e-2 )+ \
		(energy < 100 )*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 1e-2 )+ \
		(energy < 100 )*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 1e-2 )+ \
		(energy < 100 )*( abs(eta) <=0.09                     ) * ( 1e-2 )
}

# efficiency formula for c-jets (misidentification rate)
add EfficiencyFormula {4} {
	(energy >= 500 )*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 8e-1 )+ \
		(energy >= 500 )*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-1 )+ \
		(energy >= 500 )*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-1 )+ \
		(energy >= 500 )*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 3e-1 )+ \
		(energy >= 500 )*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 3e-1 )+ \
		(energy >= 500 )*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 3e-1 )+ \
		(energy >= 500 )*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 4e-1 )+ \
		(energy >= 500 )*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 4e-1 )+ \
		(energy >= 500 )*( abs(eta) <=0.09                     ) * ( 4e-1 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 6e-1 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 5e-1 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 4e-1 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 3e-1 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 2e-1 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 2e-1 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 2e-1 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 2e-1 )+ \
		(energy < 500 && energy >= 250)*( abs(eta) <=0.09                     ) * ( 2e-1 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 5e-1 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-1 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-1 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 2e-1 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 2e-1 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 2e-1 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 2e-1 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 2e-1 )+ \
		(energy < 250 && energy >= 100)*( abs(eta) <=0.09                     ) * ( 2e-1 )+ \
		(energy < 100 )*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 5e-1 )+ \
		(energy < 100 )*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 5e-1 )+ \
		(energy < 100 )*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 4e-1 )+ \
		(energy < 100 )*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 3e-1 )+ \
		(energy < 100 )*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 3e-1 )+ \
		(energy < 100 )*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 2e-1 )+ \
		(energy < 100 )*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 2e-1 )+ \
		(energy < 100 )*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 2e-1 )+ \
		(energy < 100 )*( abs(eta) <=0.09                     ) * ( 2e-1 )
}
