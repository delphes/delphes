/* ---------------------------------------------------------------------------
** This software is in the public domain, furnished "as is", without technical
** support, and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
**  A library for HEP events storage and processing based on Google's ProtocolBuffers 
**
** Author: S.Chekanov (ANL). chekanov@anl.gov
** Copyright  March, 2013 
** -------------------------------------------------------------------------*/



#include "ProMCBook.h"

bool replace(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = str.find(from);
	if(start_pos == std::string::npos)
		return false;
	str.replace(start_pos, from.length(), to);
	return true;
}



ProMCBook::ProMCBook(){ nev=0; isToWrite=false; }

ProMCBook::~ProMCBook(){ }

ProMCBook::ProMCBook(const char* filename, const char* option){

	open(filename,option);

}


/* This is a function to read or write ProMC file */
void ProMCBook::open(const char* filename, const char* option) {

	//cout << "PoMC Version: " << PROMCBOOK_VERSION << endl;


	GOOGLE_PROTOBUF_VERIFY_VERSION;
	nev=0;

	if (strcmp(option,"w")==0) {
		isToWrite=true;
		cout << "ProMCBook: Writing file = " << filename << endl;
		outzip = new ZipOutputStream(filename);
		description="none";

		outzip->putNextEntry( "version" ) ;
		std::stringstream out;
		out << current_version;
		(*outzip) << out.str();


	}  else if (strcmp(option,"r")==0) {
		isToWrite=false;
		ifstream file(filename);
		if (!file)
		{
			cout << "-> Error: File =" << filename << " does not exist " << endl;
			cout << "->        The program will exit.." << endl;
			exit(1);
		}
		file.close();
		cout << "ProMCBook: Reading file = " << filename << endl;


		zfile = new ZipFile( filename );
		unsigned nlength=zfile->size();
		// reverse iteration over all entries
		// we loop until the last event is found.
		// this is triggered by promc_description entry which follows
		// after the last event
		ConstEntries entries =zfile->entries() ;
		ConstEntries::iterator it ;
		vector<string> metatext;
		for( it = entries.end() ; it != entries.begin();) {
			--it;
			// cout << *(*it) << endl ;
			std::stringstream sout;
			sout << *(*it);
			string ss=sout.str();
			metatext.push_back(ss);
			std::size_t found = ss.find("promc_description");
			if (found != std::string::npos) {
				// istream * stt= zf.getInputStream(ss);
				//  cout << stt->rdbuf() << endl;
				// cout << "Contents of entry 1: " << cout << stt->rdbuf() << endl;
				// sentry=ss;
				//ConstEntryPointer ent = zf.getEntry(ss, FileCollection::MATCH );
				//if ( ent == 0 ) break;
				//auto_ptr< istream > is( zf.getInputStream( ent ) ) ;
				//cout << "Contents of entry 1: " << ent->getName() << " :" << endl ;
				//cout << is->rdbuf() ;
				break;
			}
		}


		// calculate the total number of entries from the file
		nev=nlength-metatext.size()-3; // 3 is the count for version,description,header
                // this is intresting: to get promc_nevents, you must call previous entry promc_description
                istream * stt = zfile->getInputStream( "promc_description", FileCollection::IGNORE ); 
		// ConstEntryPointer ent = zfile->getEntry( "promc_description", FileCollection::IGNORE ) ;
		if ( stt != 0 ) {
			// auto_ptr< istream > is( zfile->getInputStream( ent ) ) ;
			//  cout << is->rdbuf() ;
			 std::stringstream sout;
			 sout << stt->rdbuf();
			 string ss=sout.str();
                         cout << "DEBUG=" << ss << endl;
			// int nnev = atoi(ss.c_str());
		}


		// do not close to allow for random access. Close at the very end.
		// zfile->close();



		// read it
		inpzip = new ZipInputStream(filename);
		try {
			// version
			std::stringstream sout;
			sout << inpzip->rdbuf();
			string ss=sout.str();
			version = atoi(ss.c_str());


			next();




			try {

				ProMCDescription h;
				std::stringstream sout;
				sout << inpzip->rdbuf();
				h.ParseFromIstream(&sout);

				// get info
				version = h.version();
				description=h.description();
				timestamp=h.timestamp();
				requested_nev=h.events();
				cout << "ProMCBook: version     = " << version << endl;
				cout << "ProMCBook: Nr entries  = " << nev << endl;
				cout << "ProMCBook: description = " << description << endl;


			} catch( IOException &e ) {
				cerr << "IOException caught in fetching description:" << endl ;
				cerr << e.what() << endl ;
			}
			catch( ... ) {
				cerr << "Unspecified exception caught in fetching description:" << endl ;
			}








		} catch( IOException &e ) {
			cerr << "IOException caught in ProMCBook:" << endl ;
			cerr << e.what() << endl ;
		}
		catch( ... ) {
			cerr << "Unspecified exception caught in ProMCBook:" << endl ;
		}

	}



	/*
	  try {
	    zipios::ConstEntryPointer entry =inpzip->getNextEntry() ;
	    if(entry->isValid()) { 
	          std::stringstream sout;
	          sout << inpzip->rdbuf();
	          header.ParseFromIstream(&sout);
	    }
	  } catch( IOException &e ) {
	    cerr << "IOException caught in ProMCBook:" << endl ;
	    cerr << e.what() << endl ;
	   }
	   catch( ... ) {
	      cerr << "Unspecified exception caught in ProMCBook:" << endl ;
	  }
	*/


}




/**

Return the number of records.

@return Return the size of all records (excluding metadata)
**/
unsigned int  ProMCBook::size()
{
	return  nev;
}




/**
 Get the next record.
 
 @return 0 if the record was extracted OK. 6, or 7 if there is a problem

**/

int  ProMCBook::next(){


	try {
		if (inpzip == NULL) return 1;
		if (inpzip->eof()) return 2;
		if (inpzip->fail()) return 3;
		if (inpzip->bad()) return 4;
		zipios::ConstEntryPointer entry =inpzip->getNextEntry() ;
		if(!entry->isValid()) return 5;

	} catch( IOException &e ) {
		// cerr << "IOException caught in main:" << endl ;
		// cerr << e.what() << endl ;
		return 6;
	}
	catch( ... ) {
		//cerr << "Unspecified exception caught in main:" << endl ;
		return 7;
	}


	return 0;
}



/**
 Get the record with the header file.
 @param Header file record.
**/
ProMCHeader ProMCBook::getHeader(){

	// go to next event
	next();

	ProMCHeader h;

	try {
		std::stringstream sout;
		sout << inpzip->rdbuf();
		h.ParseFromIstream(&sout);
	} catch( IOException &e ) {
		cerr << "IOException caught in main:" << endl ;
		cerr << e.what() << endl ;
	}
	catch( ... ) {
		cerr << "Unspecified exception caught in main:" << endl ;
	}

	return h;

}



/**
 Get a data record (event) using a random access. 
 Use a key to extract the record. The key value
 runs from "0" to size()-1. 

 @param key (long) of the record
 @return ProMCEvent record corresponding to this key.
**/
ProMCEvent ProMCBook::event(long  idx){

  ProMCEvent eve;

  idx=idx-1;

  std::stringstream ss;
  ss << idx;
  string key=ss.str();

  if (idx==-1) key="header";

   

   try {

                istream * stt = zfile->getInputStream( key, FileCollection::IGNORE );
                std::stringstream sout;
                sout << stt->rdbuf();
                eve.ParseFromIstream(&sout);
                return eve;
              }  catch( ... ) {
                cerr << "Unspecified exception caught in main for the key:" << key << endl ;
                return eve; 
              }
}


/**
  Get the next record. Make sure you can next() first.
  @return ProMCEvent record.

**/
ProMCEvent ProMCBook::get(){

	ProMCEvent eve;

	try {

		std::stringstream sout;
		sout << inpzip->rdbuf();
		eve.ParseFromIstream(&sout);
		return eve;

		//   if (inpzip->good())
		//                    return eve;


	} catch( IOException &e ) {
		cerr << "IOException caught in main:" << endl ;
		cerr << e.what() << endl ;
	}
	catch( ... ) {
		cerr << "Unspecified exception caught in main:" << endl ;
	}


	/*
	   zipios::ConstEntryPointer entry =inpzip->getNextEntry() ;

	  if (entry->isValid()) {

	    ProMCEvent eve;

	    // entry->toString().c_str()
	    fstream input;
	    input<<entry;
	//    std::istream *str = &zipstream;

	    eve.ParseFromIstream(&input);
	    

	   } 
	*/

	return eve;

}



/**
 Clear all streams (dummmy)
**/
void ProMCBook::clear()
{

}

/**
 Set a header file.
 
 @param h Header file.

**/
void ProMCBook::setHeader( ProMCHeader h ) {

	std::string out;
	if (!h.SerializeToString(&out)) {
		cerr << "Failed to write header" << endl;
	}
	outzip->putNextEntry( "header" ) ;
	(*outzip) << out;
	h.Clear();
}


/**
  Set statistics information
  @param h Statistics info
**/
void ProMCBook::setStatistics( ProMCStat h ) {
	hStat=h;
}


/**

 Set the description information.
  @param requested_events Requested events (not the actual!)
  @param describe description

**/
void ProMCBook::setDescription( int requested_events, string  describe ) {


	description=describe;

	std::time_t t = std::time(0);  // t is an integer type

	ProMCDescription  eve;
	eve.set_version(current_version);
	eve.set_description(description);
	eve.set_events(requested_events);
	eve.set_timestamp((int)t);
	std::string out;
	if (!eve.SerializeToString(&out)) {
		cerr << "Failed to write description" << endl;
	}
	outzip->putNextEntry( "description" ) ;
	(*outzip) << out;
	eve.Clear();



}





/**
 Write an event.
 
 @param eve Event to be written.

**/
void ProMCBook::write( ProMCEvent eve ) {

	std::string out;
	if (!eve.SerializeToString(&out)) {
		cerr << "Failed to write event" << endl;
	}

	std::stringstream sout;
	sout << nev;
	outzip->putNextEntry( sout.str() ) ;
	(*outzip) << out;
	nev++;
	eve.Clear();
}



/**
 Close all files and write metadata.

**/
void ProMCBook::close() {


	if (isToWrite==true) {



		cout << " ##### Closing ProMC ##### " << endl;
		string filename;


		// description. always comes after the last event
		outzip->putNextEntry( "promc_description" ) ;
		std::stringstream out2;
		out2 << description;
		(*outzip) << out2.str();

		// write number of processed events
		outzip->putNextEntry( "promc_nevents" ) ;
		std::stringstream out1;
		out1 << nev;
		(*outzip) << out1.str();


		// write statistics
		std::string outS;
		if (!hStat.SerializeToString(&outS)) {
			cerr << "Failed to write statistics information" << endl;
			ProMCStat stat; // write dummy stattistics
			stat.SerializeToString(&outS);
			outzip->putNextEntry( "statistics" ) ;
			(*outzip) << outS;
			stat.Clear();
		}
		outzip->putNextEntry( "statistics" ) ;
		(*outzip) << outS;



		bool ierr=false;

		filename="proto/ProMCHeader.proto";
		ifstream ifs2( filename.c_str(), ios::in | ios::binary ) ;
		if(!ifs2.fail()) {
			replace(filename, "proto/", "");
			outzip->putNextEntry(filename);
			(*outzip) << ifs2.rdbuf() ;
			ifs2.close();
			if (outzip->fail())
				cout << " Problem with writing "<< filename << endl;
			//ierr=true;
		} else {
			cout << " Header file=" << filename << " not found" << endl;
			ierr=true;
		}

		filename="proto/ProMC.proto";
		ifstream ifs3( filename.c_str(), ios::in | ios::binary ) ;
		if(!ifs3.fail()) {
			replace(filename, "proto/", "");
			outzip->putNextEntry(filename);
			(*outzip) << ifs3.rdbuf() ;
			ifs3.close();
			if (outzip->fail())
				cout << " Problem with writing "<< filename << endl;
			//ierr=true;
		} else {
			cout << " Event record file=" << filename << " not found" << endl;
			ierr=true;
		}


		filename="proto/ProMCStat.proto";
		ifstream ifs4( filename.c_str(), ios::in | ios::binary ) ;
		if(!ifs4.fail()) {
			replace(filename, "proto/", "");
			outzip->putNextEntry(filename);
			(*outzip) << ifs4.rdbuf() ;
			ifs4.close();
			if (outzip->fail())
				cout << " Problem with writing "<< filename << endl;
			//ierr=true;
		} else {
			cout << " Statistics file=" << filename << " not found" << endl;
			ierr=true;
		}


		filename="proto/ProMCDescription.proto";
		ifstream ifs5( filename.c_str(), ios::in | ios::binary ) ;
		if(!ifs5.fail()) {
			replace(filename, "proto/", "");
			outzip->putNextEntry(filename);
			(*outzip) << ifs5.rdbuf() ;
			ifs5.close();
			if (outzip->fail())
				cout << " Problem with writing "<< filename << endl;
			//ierr=true;
		} else {
			cout << " Description file=" << filename << " not found" << endl;
			ierr=true;
		}


		if (ierr){
			cout << " -> Warning: Not self-describing file format." << endl;
			cout << "    To make it self-describing, put *.proto files to the directory proto/" << endl;
		}



		// write a log file if exists. Always the last one.
		filename="logfile.txt";
		ifstream ifs1( filename.c_str(), ios::in | ios::binary );
		if(!ifs1.fail()) {
			outzip->putNextEntry(filename);
			(*outzip) << ifs1.rdbuf() ;
			ifs1.close();
			cout << " Info: File=" << filename << " is attached" << endl;
			if (outzip->fail() || outzip->bad())
				cout << "Problem with writing "<< filename << endl;
		} else {
			cout << " Info: File=" << filename << " not found. No logfile info is written" << endl;
		}



		outzip->close();
		cout << " ##### ProMC file is closed #####" << endl;

	}


	// delete outzip;
	// delete inpzip;
	if (isToWrite==false) { zfile->close();
		                inpzip->close(); };
	clear();
}




