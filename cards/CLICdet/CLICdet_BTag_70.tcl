# efficiency formula for b-jets
add EfficiencyFormula {5} {0.7}
# default efficiency formula (misidentification rate)
add EfficiencyFormula {0} {(energy >= 500 )*
    (energy >= 500 )*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-2 )+ \
			       (energy >= 500 )*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 9e-3 )+ \
			       (energy >= 500 )*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 8e-3 )+ \
			       (energy >= 500 )*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 8e-3 )+ \
			       (energy >= 500 )*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 1e-2 )+ \
			       (energy >= 500 )*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 9e-3 )+ \
			       (energy >= 500 )*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 1e-2 )+ \
			       (energy >= 500 )*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 1e-2 )+ \
			       (energy >= 500 )*( abs(eta) <=0.09                     ) * ( 1e-2 )+ \
			       (energy < 500 && energy >= 250)*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-2 )+ \
			       (energy < 500 && energy >= 250)*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 7e-3 )+ \
			       (energy < 500 && energy >= 250)*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 6e-3 )+ \
			       (energy < 500 && energy >= 250)*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 6e-3 )+ \
			       (energy < 500 && energy >= 250)*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 5e-3 )+ \
			       (energy < 500 && energy >= 250)*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 5e-3 )+ \
			       (energy < 500 && energy >= 250)*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 5e-3 )+ \
			       (energy < 500 && energy >= 250)*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 5e-3 )+ \
			       (energy < 500 && energy >= 250)*( abs(eta) <=0.09                     ) * ( 5e-3 )+ \
			       (energy < 250 && energy >= 100)*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-2 )+ \
			       (energy < 250 && energy >= 100)*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-3 )+ \
			       (energy < 250 && energy >= 100)*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-3 )+ \
			       (energy < 250 && energy >= 100)*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 3e-3 )+ \
			       (energy < 250 && energy >= 100)*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 2e-3 )+ \
			       (energy < 250 && energy >= 100)*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 2e-3 )+ \
			       (energy < 250 && energy >= 100)*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 2e-3 )+ \
			       (energy < 250 && energy >= 100)*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 2e-3 )+ \
			       (energy < 250 && energy >= 100)*( abs(eta) <=0.09                     ) * ( 2e-3 )+ \
			       (energy < 100 )*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-2 )+ \
			       (energy < 100 )*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 6e-3 )+ \
			       (energy < 100 )*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-3 )+ \
			       (energy < 100 )*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 2e-3 )+ \
			       (energy < 100 )*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 1e-3 )+ \
			       (energy < 100 )*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 1e-3 )+ \
			       (energy < 100 )*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 1e-3 )+ \
			       (energy < 100 )*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 1e-3 )+ \
			       (energy < 100 )*( abs(eta) <=0.09                     ) * ( 1e-3 )
}

# efficiency formula for c-jets (misidentification rate)
add EfficiencyFormula {4} {
    (energy >= 500 )*(abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 1e-1 )+ \
	(energy >= 500 )*(abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-2 )+ \
	(energy >= 500 )*(abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-2 )+ \
	(energy >= 500 )*(abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 3e-2 )+ \
	(energy >= 500 )*(abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 4e-2 )+ \
	(energy >= 500 )*(abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 4e-2 )+ \
	(energy >= 500 )*(abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 5e-2 )+ \
	(energy >= 500 )*(abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 5e-2 )+ \
	(energy >= 500 )*(abs(eta) <=0.09                     ) * ( 6e-2 )+ \
	(energy < 500 && energy >= 250)*(abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 1e-1 )+ \
	(energy < 500 && energy >= 250)*(abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-2 )+ \
	(energy < 500 && energy >= 250)*(abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-2 )+ \
	(energy < 500 && energy >= 250)*(abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 2e-2 )+ \
	(energy < 500 && energy >= 250)*(abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 1e-2 )+ \
	(energy < 500 && energy >= 250)*(abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 1e-2 )+ \
	(energy < 500 && energy >= 250)*(abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 1e-2 )+ \
	(energy < 500 && energy >= 250)*(abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 1e-2 )+ \
	(energy < 500 && energy >= 250)*(abs(eta) <=0.09                     ) * ( 1e-2 )+ \
	(energy < 250 && energy >= 100)*(abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 1e-1 )+ \
	(energy < 250 && energy >= 100)*(abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-2 )+ \
	(energy < 250 && energy >= 100)*(abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 2e-2 )+ \
	(energy < 250 && energy >= 100)*(abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 1e-2 )+ \
	(energy < 250 && energy >= 100)*(abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 1e-2 )+ \
	(energy < 250 && energy >= 100)*(abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 8e-3 )+ \
	(energy < 250 && energy >= 100)*(abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 7e-3 )+ \
	(energy < 250 && energy >= 100)*(abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 7e-3 )+ \
	(energy < 250 && energy >= 100)*(abs(eta) <=0.09                     ) * ( 7e-3 )+ \
	(energy < 100 )* (abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 1e-1 )+ \
	(energy < 100 )* (abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 7e-2 )+ \
	(energy < 100 )* (abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 4e-2 )+ \
	(energy < 100 )* (abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 2e-2 )+ \
	(energy < 100 )* (abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 1e-2 )+ \
	(energy < 100 )* (abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 1e-2 )+ \
	(energy < 100 )* (abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 1e-2 )+ \
	(energy < 100 )* (abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 9e-3 )+ \
	(energy < 100 )* (abs(eta) <=0.09                     ) * ( 9e-3 )
}
