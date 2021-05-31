#ifdef __CINT__

#pragma link C++ struct HepMC3::GenEventData+;
#pragma link C++ struct HepMC3::GenRunInfoData+;
#pragma link C++ struct HepMC3::GenParticleData+;
#pragma link C++ struct HepMC3::GenVertexData+;
#pragma link C++ class std::vector<HepMC3::GenParticleData>+;
#pragma link C++ class std::vector<HepMC3::GenVertexData>+;
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<std::string>+;
#pragma link C++ class HepMC3::FourVector+;
#pragma link C++ class HepMC3::Units+;
#pragma link C++ class HepMC3::Setup+;

/* To generate dictionaries for compatibility with HepMC3.0*/
#pragma link C++ typedef HepMC::GenEventData+;
#pragma link C++ typedef HepMC::GenRunInfoData+;
#pragma link C++ typedef HepMC::GenParticleData+;
#pragma link C++ typedef HepMC::GenVertexData+;
#pragma link C++ typedef std::vector<HepMC::GenParticleData>+;
#pragma link C++ typedef std::vector<HepMC::GenVertexData>+;
#pragma link C++ typedef HepMC::FourVector+;
#pragma link C++ typedef HepMC::Units+;
#pragma link C++ typedef HepMC::Setup+;


#endif
