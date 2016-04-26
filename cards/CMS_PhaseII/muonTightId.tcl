   # tracking + TightID efficiency formula for muons
    set EfficiencyFormula {

        (pt <= 5.0) * (0.00) + \
		(pt > 5 && pt < 10) * (abs(eta) > 0 && abs(eta) < 0.2) * 0.949 + \
		(pt > 5 && pt < 10) * (abs(eta) > 0.2 && abs(eta) < 0.4) * 0.894 + \
		(pt > 5 && pt < 10) * (abs(eta) > 0.4 && abs(eta) < 0.6) * 0.966 + \
		(pt > 5 && pt < 10) * (abs(eta) > 0.6 && abs(eta) < 0.8) * 0.964 + \
		(pt > 5 && pt < 10) * (abs(eta) > 0.8 && abs(eta) < 1) * 0.923 + \
		(pt > 5 && pt < 10) * (abs(eta) > 1 && abs(eta) < 1.2) * 0.94 + \
		(pt > 5 && pt < 10) * (abs(eta) > 1.2 && abs(eta) < 1.4) * 0.959 + \
		(pt > 5 && pt < 10) * (abs(eta) > 1.4 && abs(eta) < 1.6) * 0.972 + \
		(pt > 5 && pt < 10) * (abs(eta) > 1.6 && abs(eta) < 1.8) * 0.962 + \
		(pt > 5 && pt < 10) * (abs(eta) > 1.8 && abs(eta) < 2) * 0.972 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2 && abs(eta) < 2.2) * 0.947 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2.2 && abs(eta) < 2.4) * 0.911 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2.4 && abs(eta) < 2.5) * 0.825 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2.5 && abs(eta) < 2.6) * 0.821 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2.6 && abs(eta) < 2.7) * 0.783 + \
		(pt > 5 && pt < 10) * (abs(eta) > 2.7 && abs(eta) < 2.8) * 0.746 + \

		(pt > 10 && pt < 20) * (abs(eta) > 0 && abs(eta) < 0.2) * 0.975 + \
		(pt > 10 && pt < 20) * (abs(eta) > 0.2 && abs(eta) < 0.4) * 0.935 + \
		(pt > 10 && pt < 20) * (abs(eta) > 0.4 && abs(eta) < 0.6) * 0.986 + \
		(pt > 10 && pt < 20) * (abs(eta) > 0.6 && abs(eta) < 0.8) * 0.976 + \
		(pt > 10 && pt < 20) * (abs(eta) > 0.8 && abs(eta) < 1) * 0.954 + \
		(pt > 10 && pt < 20) * (abs(eta) > 1 && abs(eta) < 1.2) * 0.978 + \
		(pt > 10 && pt < 20) * (abs(eta) > 1.2 && abs(eta) < 1.4) * 0.975 + \
		(pt > 10 && pt < 20) * (abs(eta) > 1.4 && abs(eta) < 1.6) * 0.98 + \
		(pt > 10 && pt < 20) * (abs(eta) > 1.6 && abs(eta) < 1.8) * 0.982 + \
		(pt > 10 && pt < 20) * (abs(eta) > 1.8 && abs(eta) < 2) * 0.985 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2 && abs(eta) < 2.2) * 0.962 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2.2 && abs(eta) < 2.4) * 0.946 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2.4 && abs(eta) < 2.5) * 0.839 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2.5 && abs(eta) < 2.6) * 0.844 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2.6 && abs(eta) < 2.7) * 0.803 + \
		(pt > 10 && pt < 20) * (abs(eta) > 2.7 && abs(eta) < 2.8) * 0.791 + \

        (pt > 20 && pt < 40) * (abs(eta) > 0 && abs(eta) < 0.2) * 0.965 + \
		(pt > 20 && pt < 40) * (abs(eta) > 0.2 && abs(eta) < 0.4) * 0.938 + \
		(pt > 20 && pt < 40) * (abs(eta) > 0.4 && abs(eta) < 0.6) * 0.99 + \
		(pt > 20 && pt < 40) * (abs(eta) > 0.6 && abs(eta) < 0.8) * 0.977 + \
		(pt > 20 && pt < 40) * (abs(eta) > 0.8 && abs(eta) < 1) * 0.959 + \
		(pt > 20 && pt < 40) * (abs(eta) > 1 && abs(eta) < 1.2) * 0.976 + \
		(pt > 20 && pt < 40) * (abs(eta) > 1.2 && abs(eta) < 1.4) * 0.984 + \
		(pt > 20 && pt < 40) * (abs(eta) > 1.4 && abs(eta) < 1.6) * 0.981 + \
		(pt > 20 && pt < 40) * (abs(eta) > 1.6 && abs(eta) < 1.8) * 0.978 + \
		(pt > 20 && pt < 40) * (abs(eta) > 1.8 && abs(eta) < 2) * 0.985 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2 && abs(eta) < 2.2) * 0.969 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2.2 && abs(eta) < 2.4) * 0.964 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2.4 && abs(eta) < 2.5) * 0.863 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2.5 && abs(eta) < 2.6) * 0.828 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2.6 && abs(eta) < 2.7) * 0.807 + \
		(pt > 20 && pt < 40) * (abs(eta) > 2.7 && abs(eta) < 2.8) * 0.764 + \

        (pt > 40) * (abs(eta) > 0 && abs(eta) < 0.2) * 0.975 + \
		(pt > 40) * (abs(eta) > 0.2 && abs(eta) < 0.4) * 0.935 + \
		(pt > 40) * (abs(eta) > 0.4 && abs(eta) < 0.6) * 0.991 + \
		(pt > 40) * (abs(eta) > 0.6 && abs(eta) < 0.8) * 0.983 + \
		(pt > 40) * (abs(eta) > 0.8 && abs(eta) < 1) * 0.957 + \
		(pt > 40) * (abs(eta) > 1 && abs(eta) < 1.2) * 0.979 + \
		(pt > 40) * (abs(eta) > 1.2 && abs(eta) < 1.4) * 0.981 + \
		(pt > 40) * (abs(eta) > 1.4 && abs(eta) < 1.6) * 0.983 + \
		(pt > 40) * (abs(eta) > 1.6 && abs(eta) < 1.8) * 0.981 + \
		(pt > 40) * (abs(eta) > 1.8 && abs(eta) < 2) * 0.984 + \
		(pt > 40) * (abs(eta) > 2 && abs(eta) < 2.2) * 0.964 + \
		(pt > 40) * (abs(eta) > 2.2 && abs(eta) < 2.4) * 0.956 + \
		(pt > 40) * (abs(eta) > 2.4 && abs(eta) < 2.5) * 0.855 + \
		(pt > 40) * (abs(eta) > 2.5 && abs(eta) < 2.6) * 0.817 + \
		(pt > 40) * (abs(eta) > 2.6 && abs(eta) < 2.7) * 0.764 + \
		(pt > 40) * (abs(eta) > 2.7 && abs(eta) < 2.8) * 0.752 + \

        (abs(eta) > 2.8) * (0.00)

    }

