/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
**  A library for HEP events storage and processing based on Google's ProtocolBuffers 
**
** Author: S.Chekanov (ANL). chekanov@anl.gov
** Copyright  2012 
** -------------------------------------------------------------------------*/


#ifndef ProMCBook_H 
#define ProMCBook_H

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "zipios++/zipios-config.h"
#include "zipios++/meta-iostreams.h"
#include "zipios++/zipoutputstream.h"
#include "zipios++/zipfile.h"
#include "zipios++/zipinputstream.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <exception>
#include "ProMC.pb.h"
#include "ProMCHeader.pb.h"
#include "ProMCStat.pb.h"
#include "ProMCDescription.pb.h"

using namespace promc;
using namespace std;
using namespace zipios ;



using std::cerr ;
using std::cout ;
using std::endl ;
using std::ios ;
using std::string ;
using std::exception ;
using std::istream;
using std::ofstream;
using std::ifstream;
using std::setw;


class ProMCBook{
public:

     ProMCBook();
     ~ProMCBook();
     ProMCBook(const char* filename, const char * option);
     void open(const char* filename, const char * option);
     ProMCEvent get();
     ProMCEvent event(long idx);
     int  next();
     int  getTimestamp() { return timestamp; };
     int  getVersion() { return version; }; 
     string getDescription() { return description; };
     void setDescription(int events, string s);
     int  getEvents() { return nev; };
     int  getEventsRequested() { return requested_nev; };
     ProMCHeader getHeader();
     void  setHeader(ProMCHeader h);
     void  setStatistics(ProMCStat h);
     void  close();
     void  write(ProMCEvent e ); // write histogram if any
     void  clear();
     unsigned int  size(); // all size
     static const int current_version=2;


private:

     int  timestamp;
     int  version;
     string description; 
     ProMCStat hStat;
     unsigned long nev;
     unsigned long requested_nev; 
     ZipOutputStream *outzip;
     ZipInputStream *inpzip;
     bool isToWrite;
     ZipFile *zfile;

};

#endif
