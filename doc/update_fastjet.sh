version=3.0.3

#wget http://fastjet.fr/repo/fastjet-${version}.tar.gz

tar -zxf fastjet-${version}.tar.gz
mkdir -p fastjet/plugins
mkdir -p fastjet/internal
mkdir -p fastjet/tools

cd fastjet-${version}
./configure --enable-allcxxplugins
cd -

echo ${version} > fastjet/VERSION

cd fastjet
cp -a ../fastjet-${version}/AUTHORS .
cp -a ../fastjet-${version}/COPYING .
cp -a ../fastjet-${version}/src/*.cc .
cp -a ../fastjet-${version}/include/fastjet/*.h .
cp -a ../fastjet-${version}/include/fastjet/*.hh .
cp -a ../fastjet-${version}/include/fastjet/internal/*.hh internal
cp -a ../fastjet-${version}/include/fastjet/internal/*.icc internal
cp -a ../fastjet-${version}/tools/*.cc tools
cp -a ../fastjet-${version}/tools/fastjet/tools/*.hh tools
cd -

cd fastjet/plugins

mkdir -p ATLASCone/fastjet
cp -a ../../fastjet-${version}/plugins/ATLASCone/*.hh ATLASCone
cp -a ../../fastjet-${version}/plugins/ATLASCone/*.cc ATLASCone
cp -a ../../fastjet-${version}/plugins/ATLASCone/README ATLASCone
cp -a ../../fastjet-${version}/plugins/ATLASCone/fastjet/*.hh ATLASCone/fastjet

mkdir -p CDFCones/fastjet
cp -a ../../fastjet-${version}/plugins/CDFCones/*.cc CDFCones
cp -a ../../fastjet-${version}/plugins/CDFCones/CDFcode/*.hh CDFCones
cp -a ../../fastjet-${version}/plugins/CDFCones/CDFcode/*.cc CDFCones
cp -a ../../fastjet-${version}/plugins/CDFCones/fastjet/*.hh CDFCones/fastjet
 
mkdir -p CMSIterativeCone/fastjet
cp -a ../../fastjet-${version}/plugins/CMSIterativeCone/*.h CMSIterativeCone
cp -a ../../fastjet-${version}/plugins/CMSIterativeCone/*.cc CMSIterativeCone
cp -a ../../fastjet-${version}/plugins/CMSIterativeCone/README CMSIterativeCone
cp -a ../../fastjet-${version}/plugins/CMSIterativeCone/fastjet/*.hh CMSIterativeCone/fastjet

mkdir -p D0RunICone/fastjet
cp -a ../../fastjet-${version}/plugins/D0RunICone/*.h D0RunICone
cp -a ../../fastjet-${version}/plugins/D0RunICone/*.hpp D0RunICone
cp -a ../../fastjet-${version}/plugins/D0RunICone/*.cc D0RunICone
cp -a ../../fastjet-${version}/plugins/D0RunICone/fastjet/*.hh D0RunICone/fastjet

mkdir -p D0RunIICone/fastjet
cp -a ../../fastjet-${version}/plugins/D0RunIICone/*.h D0RunIICone
cp -a ../../fastjet-${version}/plugins/D0RunIICone/*.hpp D0RunIICone
cp -a ../../fastjet-${version}/plugins/D0RunIICone/*.cc D0RunIICone
cp -a ../../fastjet-${version}/plugins/D0RunIICone/fastjet/*.hh D0RunIICone/fastjet

mkdir -p EECambridge/fastjet
cp -a ../../fastjet-${version}/plugins/EECambridge/*.cc EECambridge
cp -a ../../fastjet-${version}/plugins/EECambridge/fastjet/*.hh EECambridge/fastjet

mkdir -p GridJet/fastjet
cp -a ../../fastjet-${version}/plugins/GridJet/*.cc GridJet
cp -a ../../fastjet-${version}/plugins/GridJet/fastjet/*.hh GridJet/fastjet

mkdir -p Jade/fastjet
cp -a ../../fastjet-${version}/plugins/Jade/*.cc Jade
cp -a ../../fastjet-${version}/plugins/Jade/fastjet/*.hh Jade/fastjet

mkdir -p NestedDefs/fastjet
cp -a ../../fastjet-${version}/plugins/NestedDefs/*.cc NestedDefs
cp -a ../../fastjet-${version}/plugins/NestedDefs/fastjet/*.hh NestedDefs/fastjet

cp -a ../../fastjet-${version}/plugins/README .

mkdir -p SISCone/fastjet
cp -a ../../fastjet-${version}/plugins/SISCone/SISConePlugin.cc SISCone
cp -a ../../fastjet-${version}/plugins/SISCone/siscone/AUTHORS SISCone
cp -a ../../fastjet-${version}/plugins/SISCone/siscone/COPYING SISCone
cp -a ../../fastjet-${version}/plugins/SISCone/siscone/siscone/*.h SISCone
cp -a ../../fastjet-${version}/plugins/SISCone/siscone/siscone/*.cpp SISCone
cp -a ../../fastjet-${version}/plugins/SISCone/fastjet/*.hh SISCone/fastjet
rename .cpp .cc SISCone/*.cpp


mkdir -p TrackJet/fastjet
cp -a ../../fastjet-${version}/plugins/TrackJet/*.cc TrackJet
cp -a ../../fastjet-${version}/plugins/TrackJet/fastjet/*.hh TrackJet/fastjet

cd -