
INC="-I$PROMC/include -I$PROMC/src"

LIB="$PROMC/src/*.cc -L$PROMC/lib -lprotoc -lprotobuf -lprotobuf-lite -lcbook -lz"

g++ -I. -Iexternal $INC `root-config --cflags --ldflags --libs` -lEG -L. -lDelphes $LIB -o DelphesProMC examples/DelphesProMC.cpp

