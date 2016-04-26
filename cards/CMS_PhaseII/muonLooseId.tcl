 # tracking + LooseID efficiency formula for muons
    set EfficiencyFormula {
        
        (pt <= 5.0) * (0.00) + \
		(pt > 5 && pt < 10) * (abs(eta) > 0 && abs(eta) < 0.2) * 0.967 + \
		(pt > 5 && pt < 10) * (abs(eta) > 0.2 && abs(eta) < 0.4) * 0.968 + \
		(pt > 5 && pt < 10) * (abs(eta) > 0.4 && abs(eta) < 0.6) * 0.986 + \
		(pt > 5 && pt < 10) * (abs(eta) > 0.6 && abs(eta) < 0.8) * 0.982 + \
		(pt > 5 && pt < 10) * (abs(eta) > 0.8 && abs(eta) < 1) * 0.972 + \
		(pt > 5 && pt < 10) * (abs(eta) > 1 && abs(eta) < 1.2) * 0.974 + \
		(pt > 5 && pt < 10) * (abs(eta) > 1.2 && abs(eta) < 1.4) * 0.969 + \
		(pt > 5 && pt < 10) * (abs(eta) > 1.4 && abs(eta) < 1.6) * 0.988 + \
		(pt > 5 && pt < 10) * (abs(eta) > 1.6 && abs(eta) < 1.8) * 0.989 + \
		(pt > 5 && pt < 10) * (abs(eta) > 1.8 && abs(eta) < 2) * 0.992 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2 && abs(eta) < 2.2) * 0.979 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2.2 && abs(eta) < 2.4) * 0.976 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2.4 && abs(eta) < 2.5) * 0.882 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2.5 && abs(eta) < 2.6) * 0.883 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2.6 && abs(eta) < 2.7) * 0.851 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2.7 && abs(eta) < 2.8) * 0.828 + \

		(pt > 10 && pt < 20) * (abs(eta) > 0 && abs(eta) < 0.2) * 0.989 + \
		(pt > 10 && pt < 20) * (abs(eta) > 0.2 && abs(eta) < 0.4) * 0.989 + \
		(pt > 10 && pt < 20) * (abs(eta) > 0.4 && abs(eta) < 0.6) * 0.995 + \
		(pt > 10 && pt < 20) * (abs(eta) > 0.6 && abs(eta) < 0.8) * 0.995 + \
		(pt > 10 && pt < 20) * (abs(eta) > 0.8 && abs(eta) < 1) * 0.987 + \
		(pt > 10 && pt < 20) * (abs(eta) > 1 && abs(eta) < 1.2) * 0.989 + \
		(pt > 10 && pt < 20) * (abs(eta) > 1.2 && abs(eta) < 1.4) * 0.981 + \
		(pt > 10 && pt < 20) * (abs(eta) > 1.4 && abs(eta) < 1.6) * 0.99 + \
		(pt > 10 && pt < 20) * (abs(eta) > 1.6 && abs(eta) < 1.8) * 0.997 + \
		(pt > 10 && pt < 20) * (abs(eta) > 1.8 && abs(eta) < 2) * 0.997 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2 && abs(eta) < 2.2) * 0.989 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2.2 && abs(eta) < 2.4) * 0.995 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2.4 && abs(eta) < 2.5) * 0.888 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2.5 && abs(eta) < 2.6) * 0.894 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2.6 && abs(eta) < 2.7) * 0.886 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2.7 && abs(eta) < 2.8) * 0.883 + \

        (pt > 20 && pt < 40) * (abs(eta) > 0 && abs(eta) < 0.2) * 0.981 + \
		(pt > 20 && pt < 40) * (abs(eta) > 0.2 && abs(eta) < 0.4) * 0.994 + \
		(pt > 20 && pt < 40) * (abs(eta) > 0.4 && abs(eta) < 0.6) * 0.995 + \
		(pt > 20 && pt < 40) * (abs(eta) > 0.6 && abs(eta) < 0.8) * 0.997 + \
		(pt > 20 && pt < 40) * (abs(eta) > 0.8 && abs(eta) < 1) * 0.986 + \
		(pt > 20 && pt < 40) * (abs(eta) > 1 && abs(eta) < 1.2) * 0.991 + \
		(pt > 20 && pt < 40) * (abs(eta) > 1.2 && abs(eta) < 1.4) * 0.987 + \
		(pt > 20 && pt < 40) * (abs(eta) > 1.4 && abs(eta) < 1.6) * 0.992 + \
		(pt > 20 && pt < 40) * (abs(eta) > 1.6 && abs(eta) < 1.8) * 0.993 + \
		(pt > 20 && pt < 40) * (abs(eta) > 1.8 && abs(eta) < 2) * 0.999 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2 && abs(eta) < 2.2) * 0.994 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2.2 && abs(eta) < 2.4) * 0.997 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2.4 && abs(eta) < 2.5) * 0.894 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2.5 && abs(eta) < 2.6) * 0.879 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2.6 && abs(eta) < 2.7) * 0.871 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2.7 && abs(eta) < 2.8) * 0.864 + \

        (pt > 40) * (abs(eta) > 0 && abs(eta) < 0.2) * 0.986 + \
		(pt > 40) * (abs(eta) > 0.2 && abs(eta) < 0.4) * 0.989 + \
		(pt > 40) * (abs(eta) > 0.4 && abs(eta) < 0.6) * 0.996 + \
		(pt > 40) * (abs(eta) > 0.6 && abs(eta) < 0.8) * 0.994 + \
		(pt > 40) * (abs(eta) > 0.8 && abs(eta) < 1) * 0.991 + \
		(pt > 40) * (abs(eta) > 1 && abs(eta) < 1.2) * 0.989 + \
		(pt > 40) * (abs(eta) > 1.2 && abs(eta) < 1.4) * 0.985 + \
		(pt > 40) * (abs(eta) > 1.4 && abs(eta) < 1.6) * 0.993 + \
		(pt > 40) * (abs(eta) > 1.6 && abs(eta) < 1.8) * 0.996 + \
		(pt > 40) * (abs(eta) > 1.8 && abs(eta) < 2) * 0.997 + \
		(pt > 40) * (abs(eta) > 2 && abs(eta) < 2.2) * 0.99 + \
		(pt > 40) * (abs(eta) > 2.2 && abs(eta) < 2.4) * 0.993 + \
		(pt > 40) * (abs(eta) > 2.4 && abs(eta) < 2.5) * 0.905 + \
		(pt > 40) * (abs(eta) > 2.5 && abs(eta) < 2.6) * 0.86 + \
		(pt > 40) * (abs(eta) > 2.6 && abs(eta) < 2.7) * 0.845 + \
		(pt > 40) * (abs(eta) > 2.7 && abs(eta) < 2.8) * 0.855 + \

        (abs(eta) > 2.8) * (0.00)

    }

