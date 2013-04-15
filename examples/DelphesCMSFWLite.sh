
INC="`echo -I$CMSSW_FWLITE_INCLUDE_PATH | sed 's/:/ -I/g'` -I$CMSSW_RELEASE_BASE/src"

LIB="`echo -L$LD_LIBRARY_PATH | sed 's/include/lib/g; s/:/ -L/g'` `echo -L$CMSSW_FWLITE_INCLUDE_PATH | sed 's/include/lib/g; s/:/ -L/g'` -lFWCoreFWLite -lDataFormatsFWLite -lDataFormatsPatCandidates -lDataFormatsLuminosity -lCommonToolsUtils"

g++ -std=c++0x -I. -Iexternal $INC `root-config --cflags --ldflags --evelibs` -L. -lDelphes $LIB -o DelphesCMSFWLite examples/DelphesCMSFWLite.cpp

