#!/bin/sh

if [ -z "$PROMC" ]
then
  echo "** ERROR: PROMC variable is not set. Is ProMC installed?"
else
  echo ">> PROMC is set to $PROMC"
  echo ">> Building DelphesProMC";
  
  /bin/mkdir -p external/ProMC
  
  /bin/cp $PROMC/src/* external/ProMC/

  INC="-I$PROMC/include"

  LIB="external/ProMC/*.cc -L$PROMC/lib -lprotoc -lprotobuf -lprotobuf-lite -lcbook -lz"

  g++ -I. -Iexternal $INC `root-config --cflags --ldflags --libs` -lEG -L. -lDelphes $LIB -o DelphesProMC examples/DelphesProMC.cpp
fi
