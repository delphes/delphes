Delphes Validation
==================

This directory contains a set of script that run on a given Delphes detector card and produce a report in the form of a pdf file.

The script requires a configuration file.
To run locally using python multiprocessing:

```
python3 submit.py  --config config/cfg_fcchh_I.py launch_local
```

To run on condor:

```
python3 submit.py --config config/cfg_fcchh_I.py launch_condor --queue longlunch --priority group_u_CMST3.all
```

When jobs (local or batch) are done, collect, hadd and produce final report:

```
python3 submit.py --config config/cfg_fcchh_I.py --collect
```

Configuration file
==================

Examples can be found in the ```config``` directory. The validation will run on the card specified under the variable ```card```.

In order to run properly, the following collections have to be defined in the ```TreeWriter``` section of the card:
```
"Particle", "Track", "Tower", "EFlowTrack", "EFlowPhoton", "EFlowNeutralHadron", "Electron", "Muon", "Photon", "PFJet", "Jet", "CaloJet", "GenJet"
```
