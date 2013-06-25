
INC="-I$PROMC/include -I$PROMC/src"

LIB="-L$PROMC/lib -lprotoc -lprotobuf -lprotobuf-lite -lcbook -lz"

g++ -I. -Iexternal $INC `root-config --cflags --ldflags --libs` -lEG -L. -lDelphes $LIB -o DelphesProMC examples/DelphesProMC.cpp $PROMC/src/*.cc

