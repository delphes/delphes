# EIC DELPHES Instructions

## Introduction

The Electron-Ion Collider User Community has been studying the effects of detector technology choices on physics potential for the EIC. The Particle ID group has created standalone PID performance classes the assess the behavior of specific technologies, geometries, etc. To integrate this into DELPHES simulation, follow these instructions.

## Code

* First, make sure you use the DELPHES fork maintained by Stephen Sekula:

```
git clone git@github.com:stephensekula/delphes.git
```

* Next, you need to check out the PID code into the delphes/external/ folder. For now, take this code from the fork also maintained by Stephen Sekula (it's been updated to allow linking and compiling using C++ compilers, rather than just running as macros in ROOT)

```
cd delphes/external/
git clone git@gitlab.com:stephensekula/pid.git
cd -
```

* Now you can compile DELPHES as usual:

```
make -j
```

## EICPIDDetector Class

This class allows you to construct multiple PID detectors, each drawing its performance information from a class in the PID package. The setup is configured in TCL. You pass the module a list of smeared tracked, and it returns a list of tracks with the Candidate.PID data member altered. See below for usage information.

Here is an example of configuring the mRICH detector:

```
set ExecutionPath {
  # Declare all your core modules. Then:
  mRICHPID

  TreeWriter
}


module EICPIDDetector mRICHPID {
  set InputArray HCal/eflowTracks
  set OutputArray tracks

  ##
  ## mRICH Settings 
  ##
  set DetectorName mRICH
  # Pizel size in mm
  set PixelSize 3.0
  # Track resolution, in dp/p
  set TrackResolution 0.00175
  # Time resolution in ns
  set TimeResolution 1.0

  # K-Pi separation
  add Hypotheses {321} {211}
}

module TreeWriter TreeWriter {
  # Add all the core branches, then add this:
  add Branch mRICHPID/tracks mRICHTrack Track
}
```

You will now have in your output Delphes TTree a list of tracks called ```mRICHTrack``` whose PID data member has been modified from the original input track. Instead of storing one number (a PDG ID), it stores two 16-bit numbers concatenated together. The lowest 16 bits are the reconstructed identity of the track, determined from the PID detector. The highest 16 bits are the original identity of the track ("truth information"). You can retrieve the two numbers in C++ like so:

```
// True PID
(Track.PID & 0xffff0000) >> 16)
// Reco. PID
(Track.PID & 0xffff)
```
