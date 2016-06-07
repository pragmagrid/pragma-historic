/******************************************************************************
** mpiformatdb.cpp
** Provides database segmentation and placement on shared storage.
** $Id: mpiformatdb.cpp,v 1.2 2005/02/02 15:52:42 cwwang Exp $
** $Revision: 1.2 $
** 
** This file is a part of mpiBLAST and is copyright (c) 2002, 2003
** by Aaron Darling of Los Alamos National Laboratory for the Regents
** of the University of California.
** 
** License:
** 
** Copyright 2002, 2003. The Regents of the University of California.
** 
** This material was produced under U.S. Government contract
** W-7405-ENG-36 for Los Alamos National Laboratory, which is operated
** by the University of California for the U.S. Department of Energy.
** The Government is granted for itself and other acting on its behalf
** a paid-up, non-exclusive, irrevocable worldwide license in the
** material to reproduce, prepare derivative works, and perform
** publicly and display publicly. Beginning five (5) years after
** November 6, 2002, subject to additional five-year worldwide
** renewals, the Government is granted for itself and others acting on
** its behalf a paid-up, nonexclusive, irrevocable worldwide license
** in this material to reproduce, prepare derivative works, distribute
** copies to the public, perform publicly and display publicly, and to
** permit others to do so. NEITHER THE UNITED STATES NOR THE UNITED
** STATES DEPARTMENT OF ENERGY, NOR THE UNIVERSITY OF CALIFORNIA, NOR
** ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
** ASSUMES ANY LEGAL LIABLITY OR RESPONSIBILITY FOR THE ACCURACY,
** COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
** OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
** PRIVATELY OWNED RIGHTS.
** 
** Additionally, this program is free software; you can distribute it
** and/or modify it under the terms of the GNU General Public License
** as published by the Free Software Foundation; either version 2 of
** the License, or any later version. Accordingly, this program is
** distributed in the hope that it will be useful, but WITHOUT ANY
** WARRANTY; without even the implied warranty of MERCHANTABILITY or
** FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
** License for more details. You should have received a copy of the
** GNU General Public License along with this program; if not, write
** to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
** Boston, MA 02111-1307, USA.
** 
** ------------------------------------------------------------
** 
** See the file 'COPYING' for the details of the GNU General Public License
** 
******************************************************************************/
/*
 * 	5/20/03 J. D. Gans Los Alamos National Lab
 *	- Fixed fragment size calculation to increase the size of the "runt"
 *	  fragment.
 *	- Changed the definition of the command line argument "-N X" so that,
 *	  when possible, X fragments are created and numbered from [0, ..., X - 1].
 *	- Added a command line argument "--decomp" that will display all
 *	  possible fragment decompositions and suggest decompositions that
 *	  maximize the relative size of the runt fragment (to improve load 
 *	  balancing).
 */
 
#include "mpiblast_config.hpp"
#include "mpiblast_util.hpp"
#include "mpiblast_types.h"
#include "fragment_list.hpp"
#include "file_util.hpp"
#include "db_spec.hpp"

#include <map>
#include <cstdio>
#include <sstream>
#include <math.h>
using namespace std;

extern "C"{
#include "ncbi.h"
//#include "ncbiwin.h"
}

#include "getopt.h"

void print_usage( const char* pname ){
	cout << "Usage: " << pname << " [--config-file=<CONFIGFILE>] <--nfrags=<NumberOfFrags>|--frag-size=<SizeOfFrags> (in MB)> <formatdb arguments>\n";
	cout << "Options:\n";
	cout << "\t-f  --config-file=<CONFIGFILE>  Specify the location of the mpiBLAST configuration file\n";
	cout << "\t-N  --nfrags=<NumberOfFrags>    Specifies the number of fragments to split the database into\n";
	cout << "\t-S  --frag-size=<SizeOfFrags>   Specifies the size for database fragments, can not be used with --nfrags\n";
	cout << "\t    --version                   Prints the version number\n";
	cout << "\t    --debug                     Prints out debugging information\n";
	cout << "\t    --decomp                    Enumerates all allowed database decompositions (these are approximate!)\n";
	cout << "\t    --update_db=<DBtoUpdate>    Update database by adding additional fragments\n";

}

// silly Windows defines _snprintf instead of snprintf
#if defined(WIN32) && !defined(snprintf)
#define snprintf _snprintf
#endif // WIN32

#define BUFFER_SIZE 10000
/**
 * count the number of seqence characters in seq_file (assuming it is FastA format)
 * @param seq_file	A C File object pointing to the open sequence file
 */
uint64 countChars( FILE *seq_file ){
	int i;
	uint64 total_chars = 0, counter = 0 ;
	char buffer[BUFFER_SIZE] ;
	char defline = 0 ; //defline is used as bool
	cout << "Reading input file  " ;
 	while ((fgets(buffer, BUFFER_SIZE, seq_file)) != 0) {
		defline = 0 ;
		for (i=0; buffer[i] != '\0'; i++) {
			if (buffer[i] == '>') {
				defline = 1 ;
			} else if ((buffer[i] == '\n') || (buffer[i] == '\r')) {
				defline = 0 ;
			} else if (defline == 0) {
				total_chars++ ;
			}
		}
		if ((counter%10000) == 2500){
			cout << "\b\\";
			cout.flush();
		}else if ((counter%10000) == 5000){
			cout << "\b-";
			cout.flush();
		}else if ((counter%10000) == 7500){
			cout << "\b/";
			cout.flush();
		}else if ((counter%10000) == 0){
			cout << "\b|";
			cout.flush();
		}

		counter++ ;
	}

	cout << "\nDone, read " << counter << " lines\n" ;	
	return total_chars;
}

/**
 * For a given database size (in number of bases) enumerate the allowed fragment decompositions
 * (i.e. number of fragements, fragment size and size of the last, or 'runt', fragment). This will
 * allow the user to make a sensible choice as to the number of fragments to use. For good load balancing,
 * one would like to make the runt fragment size as close to the fragment size as possible.
 * @param m_bases	The number of bases in the database
 */
void PrintAllowedDecomp(uint64 m_bases)
{
	int effective_size, remainder, i, num_frag, frag_size;
	double true_size, runt_size;
	map<int, double> runt_hash;
	
	true_size = m_bases / 1000000.0;
	
	// Convert bases to MB
	if(m_bases % 1000000 == 0){
		effective_size =  m_bases / 1000000;
	}
	else{
		effective_size = 1 + m_bases / 1000000;
	}
	
	cerr << "Actual database size is " << true_size << " mb" << endl;
	cerr << "Effective database size is " << effective_size << " mb" << endl;
	
	for(remainder = effective_size - 1;remainder >= 0;remainder--){
		for(i = 2;i < effective_size;i++){
			if((effective_size % i) == remainder){
				num_frag = effective_size / i;
				
				if(remainder > 0){
					num_frag ++;
				}
				
				runt_size = true_size - i*(num_frag - 1.0);
				
				if((runt_hash.find(num_frag) == runt_hash.end()) || 
				  (runt_size > runt_hash[num_frag])){
					
					runt_hash[num_frag] = runt_size;
				}
				
			}
		}
	}
	
	cerr << "Here are the approximate database fragmentations. Pick one that minimizes\n"
	     << "the difference between the fragment size and the runt fragment size.\n\n"
	     << "# of frags --> frag size : runt frag size \n";
	
	// Don't display decompositions with less than 3 fragments
	for(i = effective_size;i > 2;i --){
		if(runt_hash.find(i) != runt_hash.end()){
			frag_size = effective_size / i;
			
			if(effective_size % i != 0){
				frag_size ++;
			}
			
			cerr << i << " --> " << frag_size << " mb : " << runt_hash[i] << " mb";
			
			// What fragment combinations are most efficient? Here is an arbitrary standard.
			const double promising_cutoff = 0.03;
			
			double promising = (frag_size - runt_hash[i])/frag_size;
			
			if(promising <= promising_cutoff){
				cerr << "\t<-- Looks promising! (" << promising << ")";
			}
			
			cerr << endl;
		}
	}
}

/** 
 * parse an alias file and count the fragments in it
 */
uint getAliasFragmentCount( istream& alias_file ) {
	string line;
	uint fragment_count = 0;
	while( getline( alias_file, line ) ){
		stringstream line_str( line );
		// trim off comments
		getline( line_str, line, '#' );
		if( line.size() == 0 )
			continue;
		stringstream token_str( line );
		token_str >> line;
		// look for DBLIST token
		if( line != "DBLIST" )
			continue;
		while( token_str >> line ){
			fragment_count++;
		}
	}
	return fragment_count;
}

/**
 * add db frags to existing alias file for update_db
 */
void modifyAliasFile( const string& alias_orig, int high_frag_old, int high_frag_new ) {
	string alias_temp = alias_orig + "XXXXXX";
	getTempFileName( alias_temp );
	ifstream alias_in( alias_orig.c_str() );
	ofstream alias_out( alias_temp.c_str() );
	string line;
	string frag_addition = alias_orig.substr( 0, alias_orig.length() - 3 );
	char fragnum[4];
	try{
	while( getline( alias_in, line ) ){
		alias_out << line;
		// look for DBLIST token
		if( line.size() > 6 && line.substr(0,7) == "DBLIST " ) {
			for( int n = high_frag_old; n < high_frag_new; n++ ){
#ifdef MANY_FRAGMENTS_HACK				
				sprintf(fragnum,"%03d",n);
#else
				sprintf(fragnum,"%02d",n);
#endif
				alias_out << " " << frag_addition << fragnum;
			}
		}
		alias_out << endl;
	}
	rename( alias_temp.c_str(), alias_orig.c_str() );
	}catch( exception& e ){
		cerr << e.what() << endl;
	}
	return;
}

	
int main( int argc, char* argv[] ){
	MpiBlastConfig config;
	string localPath, blast_cl, blast_type = "p", input_filename;
	
	log_stream = &cout;

	if( argc < 2 ){
		if( argc == 0 )
			print_usage( "mpiformatdb" );
		else
			print_usage( argv[0] );
	}
	
	//
	// must read -i and -p options
	//
	const char* options = "-i:p:f:n:N:S:";
	int ac = argc;
	char** av = argv;
	int opt_code, opt_type = 0;
	opterr = 0;
	int config_opt = 0;
	int long_index = 0;
	vector< string > formatdb_opts;
	formatdb_opts.push_back( "formatdb" );
	
	string config_f_name = MpiBlastConfig::defaultConfigFileName();
	bool print_version = false;
	int nfrags = -1;
	int frag_size = -1;
	
	// Do we need to print a fragment decomposition?
	int print_decomp = false;

	// the full path name of the db to update
	string update_db;

	struct option long_opts[] = {
		{"config-file", required_argument, &config_opt, 'f' },
		{"nfrags", required_argument, &config_opt, 'N' },
		{"frag-size", required_argument, &config_opt, 'S' },
		{"version", no_argument, &config_opt, 6 },
		{"decomp", no_argument, &config_opt, 7 },
		{"debug", no_argument, &config_opt, 8 },
		{"update_db", required_argument, &config_opt, 9 },

		{0,0,0,0}	// signifies termination of option list
	};

	// get the mpiformatdb and interesting formatdb options with getopt_long
	while((opt_code = getopt_long( ac, av, options, long_opts, &long_index ) ) != EOF ){
		switch( opt_code ){
			case 0:
				if( config_opt == 'f' ){
					config_f_name = optarg;
				}
				if( config_opt == 'N' ){
					nfrags = atoi( optarg );
				}
				if( config_opt == 'S' ){
					frag_size = atoi( optarg );
				}
				if( config_opt == 6 ){
					print_version = true;
				}
				
				if( config_opt == 7 ){
					print_decomp = true;
				}
				if( config_opt == 8 ){
					debug_msg = true;
				}
	  			if( config_opt == 9 )
	   				update_db = optarg;
				break;
			case 'f':
				config_f_name = optarg;
				break;
			case 'N':
				nfrags = atoi( optarg );
				break;
			case 'S':
				frag_size = atoi( optarg );
				break;
			case 'i':
				input_filename = optarg;
				addOpt( formatdb_opts, opt_code, optarg );
				break;
			case 'p':
				addOpt( formatdb_opts, opt_code, optarg );
				if( strcmp( optarg, "F" ) == 0 )
					blast_type = "n";
				else
					blast_type = "p";
				break;
			case 'n':
				break;	// override the user specified -n
			case 1:
				formatdb_opts.push_back( optarg );
				break;
			case '?':
				addOpt( formatdb_opts, optopt, optarg );
				opt_type = optopt;
				break;
			default:
				addOpt( formatdb_opts, opt_code, optarg );
//				LOG_MSG << "Default " << (char)optopt << endl;
				opt_type = optopt;
				break;
		}
	}
	if( print_version ){
		cout << PACKAGE << " version " << VERSION << endl;
		cout << "Built " << __DATE__ << " at " << __TIME__ << endl;
		return 0;
	}

	// just in case an argument was left behind -- go get it
	for( int optI = optind; optI < ac; optI++ ){
		formatdb_opts.push_back( av[ optI ] );
	}
	
	//If we didn't specify any of -N, -S, --decomp, die
	if( nfrags <= 0 && frag_size <= 0  && print_decomp != true ){
		cerr << "You must specify one of: --nfrags, --frag-size, --decomp\n";
		return -1;
	}


	// read the configuration file
	try{
		config = MpiBlastConfig( config_f_name );
	}catch(const char *error){
		cerr << "Fatal Error:" << endl;
		cerr << error << endl;
	}
	localPath = config.localPath();

	// check for valid option combinations
	if( input_filename == "" ){
		print_usage( argv[0] );
		cerr << "You must specify an input file in FastA format\n";
		return -1;
	}

	//If update_db, check to make sure that we have a pre-formatted db
	string formatted_db = getFilename( input_filename );
	string dbspec_file ;
	uint64 ebases_old = 0;
	int high_frag_old = 0;

/*
 * Here we do checks for updating a database, then read in the necessary info regarding the old one.
 * First we check if the .dbs file exists, and whether we were given the full path or just the dbname
 * Next we make sure we're working with the same type of db (nucleotide or protein)
 * Finally we read in the info from the .dbs file
 */
	if ( update_db.size() > 0 ) {
		dbspec_file = update_db + ".dbs" ;
		if ( !doesFileExist(dbspec_file) ){
			update_db = config.sharedPath() + PATH_SEPARATOR + update_db;
			dbspec_file = update_db + ".dbs";
			if ( !doesFileExist(dbspec_file) ){
				cerr << dbspec_file << " does not exist, cannot update this db.\n";
				return -1;
			}
		}
		string al_orig = update_db + "." + blast_type + "al";
		if ( !doesFileExist(al_orig) ){
			cerr << al_orig << " does not exist.  Make sure the update and existing database are both the same type (NT or AA)\n";
			return -1;
		}
		
		try{
			DbSpecFile db_spec = DbSpecFile( dbspec_file );
			ebases_old = db_spec.dbSize();
			high_frag_old = db_spec.fragmentCount();
		}catch( const char* error ){
			cerr << error << endl;
			return -1;
		}
	}
	else {
		dbspec_file = config.sharedPath() + PATH_SEPARATOR + formatted_db + ".dbs";
	}

	// open the FastA sequence file
 	FILE *db_input_fp;

	if (!(db_input_fp = fopen( input_filename.c_str(), "r" ))) {
 		cerr << "Error: Unable to open " << input_filename << endl ;
 	}

	/*	
	* Here we read the database file to get the database size
	* We then approximate the effective size, which is really based on the relative entropy
	* of both the query and db sequence. Karlin and  Altschul  (1990)
	* This tends result in some E-values being .1 higher
	* (12/2230 matches on a 58MB database and an 85KB query file with 162 query sequences) E-values 
	*  We then add to this the total from the db we're updating, which is 0 if we're not updating
	*/
	uint64 bases = countChars( db_input_fp );
	uint64 ebases = (uint64)(bases * .998) ;
	
	if(print_decomp){
		// Display the allowed number of fragment/fragment size decompositions
		// to let the user select the decomposition that is best for them.
		PrintAllowedDecomp(bases);
		
		return 0;
	}
	
	// if we get here, we will be fragmenting the database.
	// check that the user specified a fragmentation
	if( nfrags <= 0 && frag_size <= 0 ){
		print_usage( argv[0] );
		cerr << "You must specify a positive value for either --nfrags or --frag-size\n";
		return -1;
	}
	
	if( nfrags > 0 && frag_size > 0 ){
		print_usage( argv[0] );
		cerr << "You may only specify one of --nfrags or --frag-size\n";
		return -1;
	}

	// formatdb takes its -v option in MB
	if(bases % 1000000 == 0){
		// It is pretty unlikely that bases % 1000000 == 0, but completeness demands
		// that we allow for the possiblilty ...
		bases /= 1000000;
	}
	else{
		bases = 1 + bases/1000000;
	}

	// If we didn't specify frag size on the cmd line, figure it out
	if( frag_size <= 0 ){
		if(bases % nfrags == 0){
			frag_size = bases / nfrags;
		}
		else{
			// Since the number of fragments does not divide the database size
			// exactly, we need to round the fragment size up to the nearest
			// integer.
			frag_size = 1 + bases / nfrags;
		}
	}
	if( frag_size == 0 ){
		cerr << "The database is too small for the requested number of fragments\n";
		return -1;
	}
	if( nfrags <= 0 ){
		nfrags = bases / frag_size;
		
		if( bases % frag_size != 0 )
			nfrags++;
	}
		
	if( nfrags == 0 ){
		cerr << "The database is too small for the requested fragment size\n";
		return -1;
	}
		
	formatdb_opts.push_back( "-v" );
	ostringstream oss;
	oss << frag_size;
	formatdb_opts.push_back( oss.str().c_str() );

	string base_name = config.sharedPath() + PATH_SEPARATOR + formatted_db;
	formatdb_opts.push_back( "-n" );
	formatdb_opts.push_back( base_name );
	
	cout << "Trying to break " << formatted_db << " (" << bases << " MB) into " << nfrags << 
		" fragments of " << frag_size << " MB\n";
	
	string formatdb_cl;
	for( uint optI = 0; optI < formatdb_opts.size(); optI++ )
		formatdb_cl += formatdb_opts[ optI ] + " ";
	LOG_MSG << "Executing: " << formatdb_cl << endl;
	// execute formatdb
	initNCBI( formatdb_opts );
	
	int retval = Main();

	cleanupNCBI();
	if( retval != 0 ){
		cerr << "There was an error executing formatdb.  Check formatdb.log\n";
		
		// Only return error values <= 0 (since this program now returns the number of
		// fragments actually created).
		return retval > 0 ? -retval : retval;
	}
	
	// Get the actual number of db fragments created from the alias file
	// If we can't open the alias file, check if we only generated 1 fragment
	// if only 1 frag, rename it unless we're not updating a database
	string al_dest = config.sharedPath() + PATH_SEPARATOR + formatted_db + "." + blast_type + "al";
	int fragmentI;
	if( debug_msg )
		LOG_MSG << "Opening " << al_dest << endl;

	ifstream alias_file( al_dest.c_str() );
	const vector< string >& extensions = FragmentExtensions( blast_type );
	const vector< string >& extensions_optional = FragmentExtensions_optional( blast_type );

	if( !alias_file.is_open() ){ 
		// can't open alias file, either 1 or 0 fragments were created
		string nfrag_check = config.sharedPath() + PATH_SEPARATOR + formatted_db + "." + blast_type + "sq";
		ifstream sq_file( nfrag_check.c_str() );
		if( !sq_file.is_open() ){  //doesn't look like we made anything. die
			cerr << "Error opening " << al_dest << endl;
			return -1;
		}

		//only created one frag, rename it and create an alias file for it
		sq_file.close();

		// first rename the database fragment
		string old_frag, new_frag, frag_number;

#ifdef MANY_FRAGMENTS_HACK
		frag_number = "000";
#else
		frag_number = "00";
#endif
		for( uint extI = 0; extI < extensions.size(); extI++ ){
			old_frag = base_name + extensions[extI];
			new_frag = base_name + "." + frag_number + extensions[extI];
			rename( old_frag.c_str(), new_frag.c_str() );
		}
		for( uint extI = 0; extI < extensions_optional.size(); extI++ ){
			old_frag = base_name + extensions_optional[extI];
			new_frag = base_name + "." + frag_number + extensions_optional[extI];
			rename( old_frag.c_str(), new_frag.c_str() );
		}
		
		// now make an alias file for it
		ofstream alias_out( al_dest.c_str() );
		alias_out << "DBLIST " << base_name << "." << frag_number << endl;
		alias_out.close();

		fragmentI = 1;
	}
	else {
		fragmentI = getAliasFragmentCount( alias_file );
	}

	cout << "Created " << fragmentI << " fragments.\n";
	// If a database is being updated rename the new fragments and add them
	// to the existing alias file.
	if (update_db.size() > 0) {
		char new_number[4];
		char old_number[4];

		for (int n = 0; n < fragmentI ; n++){
#ifdef MANY_FRAGMENTS_HACK
			snprintf( new_number, 4, "%03d", (n+high_frag_old) );
			snprintf( old_number, 4, "%03d", n );
#else
			snprintf( new_number, 3, "%02d", (n+high_frag_old) );
			snprintf( old_number, 3, "%02d", n );
#endif
			for( uint extI = 0; extI < extensions.size(); extI++ ){
				string old_frag = base_name + "." + old_number + extensions[extI];
				string new_frag = update_db + "." + new_number + extensions[extI];
				rename( old_frag.c_str(), new_frag.c_str() );
			}
			for( uint extI = 0; extI < extensions_optional.size(); extI++ ){
				string old_frag = base_name + "." + old_number + extensions_optional[extI];
				string new_frag = update_db + "." + new_number + extensions_optional[extI];
				rename( old_frag.c_str(), new_frag.c_str() );

			}
		}
		//Remove created [.nal|.pal] file and modify original one
		remove( al_dest.c_str() );
		string al_orig = update_db + "." + blast_type + "al";
		modifyAliasFile( al_orig, high_frag_old, (fragmentI + high_frag_old) );

	}

	// write the mpiBLAST database specification file
	ofstream dbspec_out( dbspec_file.c_str() );
	if( !dbspec_out.is_open() ){
		cerr << "Error opening " << dbspec_file << endl;
		return 0;
	}
	if( update_db.size() > 0 )
		DbSpecFile::write( dbspec_out, update_db, ebases + ebases_old, fragmentI + high_frag_old );
	else
		DbSpecFile::write( dbspec_out, formatted_db, ebases, fragmentI );
	dbspec_out.close();
	
	// mpiformatdb now returns the number of fragments actually created. This allows
	// automated scripts to record the number of fragments and go on to call mpiblast with
	// the correct number of processors. A return value of <= 0 indicates an error has
	// occured!
	//return 0;

// tar mpiblast db files to *.tgz for mpiblast-g2 use
        cout << "generate tar files of database for gasscopy" << endl;  
        int i,j,k,m,n,frag_order;
        double j1;
        string st_frag,tar_cmd,rm_cmd;

// if not configure with --enable-many-fragment
        frag_order = 2; 

        for (i=0 ; i<=nfrags-1 ; i++){
            m=i;
            if (m != 0) j1= log10(m*1.0);
            else j1= 0.0;
            j = int(j1)+1; 
            st_frag= ""; 
            n=1;
            for (k=0 ; k<frag_order-j; k++) st_frag+= "0";
            for (k=0 ; k<j-1  ; k++) n=n*10;
            for (k=0 ; k<j  ; k++){
                st_frag+= m/n+48;
                m=m%n;
                n=n/10;
            }
            rm_cmd = "/bin/rm -f "+input_filename+".";
            rm_cmd+= st_frag;
            rm_cmd+= ".tgz";
            LOG_MSG << "Executing: " << rm_cmd << endl;
            system(rm_cmd.c_str());
            tar_cmd = "/bin/tar -czf "+input_filename+".";
            tar_cmd+= st_frag;
            tar_cmd+= ".tgz ";
            tar_cmd+= input_filename;
            tar_cmd+= ".";
            tar_cmd+= st_frag;
            tar_cmd+= ".* ";
            LOG_MSG << "Executing: " << tar_cmd << endl;
            j=system(tar_cmd.c_str());
            if (j != 0){
               rm_cmd = "/bin/rm -f "+input_filename+".";
               rm_cmd+= st_frag;
               rm_cmd+= ".tgz";
               LOG_MSG << "Executing: " << rm_cmd << endl;
               system(rm_cmd.c_str());
               cout << "success to generate tar files for gasscopy" << endl;  
               break;
            }
            cout << "success to generate tar files for gasscopy" << endl;
        }
	return fragmentI;
}
