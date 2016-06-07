/******************************************************************************
** mpiblast.cpp
** Implements the MpiBlast class and associated helper functions
** $Id: mpiblast.cpp,v 1.15 2005/01/27 06:51:22 cwwang Exp $
** $Revision: 1.15 $
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
 *	- Improved memory management.
 *	- Added a new function, sendResults(), that returns the blast results
 *	  from a worker to the master using less memory than sendFile(). This
 *	  allows the use of larger queries (i.e. querying the entire y. pestis
 *	  genome against the nr database can now be completed on a master node
 *	  with 2 GB of ram memory). Rewrote the receiveResults() function to 
 *	  work with the new sendResults() function.
 *	7/7/03 J. D. Gans Los Alamos National Lab
 *	- Rewrote sendResults() [again!] -- This time to fix a problem that
 *	  caused tcp/ip timeouts for MPI. For more info, see the "Unsolved
 * 	  problems" section of http://www-unix.mcs.anl.gov/mpi/mpich/buglist-tbl.html
 *	  and look for the "errno=110 or connection timed out" bug. The problem occurs
 * 	  with both mpich (1.2.5) and LAM MPI (6.5.9). The change to
 *	  sendResults() amounts to writing the temp blast output file of the worker 
 *	  in binary ASN.1 format and then sending each record individually to the master.
 *	  This is not a hack as the binary files will give mpiblast a smaller disk
 * 	  space foot print on the worker nodes (and consume less bandwidth).
 *	9/7/03 J. D. Gans Los Alamos National Lab
 *	- Tried using ElectricFence (http://perens.com/FreeSoftware) to check for 
 *	  memory allocation errors but it fails with the error: "Electric Fence Exiting: 
 *	  mmap() failed: Cannot allocate memory". Need to increase the preprocess mmap memory
 *	  using sysctl vm.max_proc_mmap (will this work for Linux?). Update --
 *	  Electric Fence will work as long as only a handfull of sequences are
 *	  quried (i.e. two).
 *	- Modified code to allow distributed database access when compling blast output.
 */
 
#include "getopt.h"
#include <math.h>
#include "mpi.h"

#include <signal.h>
#include <algorithm>
#include <list>
#include "mpiblast.hpp"
#include "mpiblast_util.hpp"
#include "embed_rank.hpp"

#include "blast_hooks.h"
#include "ncbi_sizeof.h"
#include "distributed_bioseq.h"

using namespace std;

#ifdef MPE
#include <mpe.h>
int cpstart, cpend, blaststart, blastend, mergestart, mergeend; //MPE ids
int outputstart, outputend; //more MPE ids
#endif

string localPath;
static string alias_basename, alias_filename , local_query_tmpfile , local_results_tmpfile;

static MpiBlast* mpiblast_ptr = NULL;

// Local functions
static bool running_blast = false;
void checkBlastErrorExit();
void terminateProgram( int sig );


/**
 *  Scans through the list of SeqAlignPtrs separating it into several lists,
 *  each of which corresponds to a single BLAST result.
 */
void BreakUpResults( vector< SeqAlignPtr >& results_vec, SeqAlignPtr sap ){
	SeqAlignPtr cur_a = sap;
	results_vec.clear();
	results_vec.push_back( cur_a );
	int cur_count = 1;
	//	cerr << cur_count << "\t" << "cur_a: " << cur_a << endl;
	while( cur_a->next != NULL ){
		//		cerr << cur_count << "\t" << "cur_a: " << cur_a->next << endl;
		SeqIdPtr cur_id = m_TxGetSubjectIdFromSeqAlign( cur_a ); 
		SeqIdPtr next_id = m_TxGetSubjectIdFromSeqAlign( cur_a->next );
		if( m_SeqIdPtrCmp( cur_id, next_id ) == 0 ){
			// break the chain at cur_a, add cur_a->next
			SeqAlignPtr tmp_next = cur_a->next;
			cur_a->next = NULL;
			//			cerr << "tmp_next\t" << tmp_next << endl;
			results_vec.push_back( tmp_next );
			cur_a = tmp_next;
			//			cerr << "Segment count is: " << cur_count << endl;
			cur_count = 1;
		}else{
			cur_a = cur_a->next;
			cur_count++;
		}
	}
}

void mergeSeqAlignResults( const vector< SeqAlignPtr >& results, const vector< SeqAlignPtr >& new_results, 
			  vector< SeqAlignPtr >& merged_results ){
	uint resultI = 0;
	uint resultJ = 0;
	while( resultI != results.size() && resultJ != new_results.size() ){

		int a_score = 0, b_score = 0, a_number, b_number;
		Nlm_FloatHi a_bit_score = 0, b_bit_score = 0, a_evalue, b_evalue;

		m_GetScoreAndEvalue( results[ resultI ], &a_score, &a_bit_score, &a_evalue, &a_number );
		m_GetScoreAndEvalue( new_results[ resultJ ], &b_score, &b_bit_score, &b_evalue, &b_number );
//		cerr << "a_score: " << a_bit_score << "\t b_score: " << b_bit_score << endl;
		if( a_bit_score >= b_bit_score ){
			merged_results.push_back( results[ resultI ] );
			resultI++;
		}else{
			merged_results.push_back( new_results[ resultJ ] );
			resultJ++;
		}
	}
	for( ; resultI != results.size(); resultI++ ){
		merged_results.push_back( results[ resultI ] );
	}
	for( ; resultJ != new_results.size(); resultJ++ ){
		merged_results.push_back( new_results[ resultJ ] );
	}
}


void mergeResults( map< int, vector< SeqAlignPtr > >& query_map, vector< SeqAlignPtr >& new_results, int* result_list )
{
	for( uint resultI = 0; resultI < new_results.size(); resultI++ ){

		int result_id = result_list[ resultI ];
		vector< SeqAlignPtr > cur_results;
		BreakUpResults( cur_results, new_results[ resultI ] );

		map< int, vector< SeqAlignPtr > >::iterator query_iter = query_map.find( result_id );
		
		if( query_iter == query_map.end() ){
			// no results for this query yet.  add it.
			query_map.insert( map< int, vector< SeqAlignPtr > >::value_type( result_id, cur_results ) );
		}else{
			// merge these results with the existing query results.
			vector< SeqAlignPtr > merged_results;
			mergeSeqAlignResults( query_iter->second, cur_results, merged_results );
			query_map.erase( query_iter );
			query_map.insert( map< int, vector< SeqAlignPtr > >::value_type( result_id, merged_results ) );
		}
	}
}

// DEBUG --  Electric fence global variables
//extern int EF_ALLOW_MALLOC_0;
//extern int EF_PROTECT_BELOW;

int main( int argc, char* argv[] )
{
	// DEBUG --  Electric fence global variables
	//EF_ALLOW_MALLOC_0 = 1;
	//EF_PROTECT_BELOW = 1;
	
	// Perform a quick version check (and exit the program if
	// the version flag is found). If we don't check now, there's no way for
	// the user to query the version with out having mpi up and running.
	for(int i = 1;i < argc;i++){
		if(strcmp(argv[i], "--version") == 0){
			cerr << argv[0] << " version " 
			     << VERSION << endl;
			exit(1);
		}
	}
	
#ifdef MPE
	char mpelogfmt[20] = "MPE_LOG_FORMAT=SLOG";
	putenv(mpelogfmt);
#endif
	double realstarttime = clock();

	/* Start up MPI */
	
	// Original version (no MPI thread support)
	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &node_count);

	double progstarttime = clock();
#ifdef MPE
	cpstart = MPE_Log_get_event_number(); 
	cpend = MPE_Log_get_event_number(); 
	blaststart = MPE_Log_get_event_number(); 
	blastend = MPE_Log_get_event_number(); 
	mergestart = MPE_Log_get_event_number(); 
	mergeend = MPE_Log_get_event_number(); 

	MPE_Describe_state(cpstart, cpend, "CopyFagment","red");
	MPE_Describe_state(blaststart, blastend, "RunBlastall","DarkKhaki");
	MPE_Describe_state(mergestart, mergeend, "MergeResults","grey");
#endif

	prog_start = MPI_Wtime();
	log_stream = &cerr;
	MpiBlast mpb;
	
	if(node_count < 2){ 
		//sorry, we can't run on only 1 node.
		terminateProgram(0);
	}
	
	signal( SIGINT, terminateProgram );
	signal( SIGTERM, terminateProgram );
	signal( SIGSEGV, terminateProgram );

	// register an exit function to abort MPI if blast exits
	// without an error code
	atexit( &checkBlastErrorExit );
	
	mpiblast_ptr = &mpb;
	
	int retcode = mpb.main( argc, argv );
	
	prog_end = MPI_Wtime();
	
#ifdef GASSCOPY
        //MPI_Barrier(MPI_COMM_WORLD);

        if(rank == 0) {
               string out_file = mpb.getOutputFilename();
               string rslt_url = mpb.getResultURL();
 	       string curr_dir = get_current_dir_name();
 	       string output_fpath = "file://"+curr_dir+PATH_SEPARATOR+out_file;

               if( profile_msg ) {
                   PROFILE_NEW_MSG << "Fetch output took ";
               }

               double output_start = MPI_Wtime();

 	       if( gassCopyFile(output_fpath, rslt_url) != 0 ) {
 	       		LOG_MSG << "Unable to put " << output_fpath << " to "
 	       		        << rslt_url << endl;
 	       	     
 	       		throw __FILE__ "(main): Unable to put result to remote wdir";
 	       }

               double output_end = MPI_Wtime();

 	       if( profile_msg ){
                   PROFILE_MSG << (output_end - output_start) << endl;
                   PROFILE_NEW_MSG << "Total execution time is: " << (MPI_Wtime() - prog_start) << endl;
 	       }
        }

        // Returning profiles from each node
        if(profile_msg) {

               string prof_file = mpb.getProfileFilename();
               string prof_url  = mpb.getProfileURL();
 	       string curr_dir = get_current_dir_name();
 	       string output_fpath = "file://"+curr_dir+PATH_SEPARATOR+prof_file;
 	       if( debug_msg ){
 	           LOG_MSG << "profile fpath:" << output_fpath << endl;
 	           LOG_MSG << "profile file:" << prof_file << endl;
 	           LOG_MSG << "profile url:" << prof_url << endl;
 	       }


 	       if( gassCopyFile(output_fpath, prof_url) != 0 ) {
 	       		LOG_MSG << "Unable to put " << output_fpath << " to "
 	       		        << prof_url << endl;
 	       	     
 	       		throw __FILE__ "(main): Unable to put result to remote wdir";
 	       }
        }
#endif

        /* Move final profile_msg to the place before sending it back to client
	if( profile_msg ){
		PROFILE_NEW_MSG << "Total execution time is: " << (MPI_Wtime() - prog_start) <<endl;
	}
        */
	
	if( debug_msg ){
		LOG_MSG << "rank " << rank << " exiting successfully\n";
	}
	
	if( debug_msg ){
		LOG_MSG << "MPI startup  time is " << (double)((progstarttime - realstarttime)/CLOCKS_PER_SEC) << endl;
	}

	mpiblast_ptr = NULL;

//**************  start to remove scratch files ***********
        LOG_MSG << "remove scratch files" ;
	MPI_Barrier(MPI_COMM_WORLD);
	printf("remove %s %i\n",localPath.c_str(),rank);
	localPath="rm -rf " + localPath;
	int fss = system(localPath.c_str());
        if (fss == 0 )
           printf("successful to remove scratch files\n");
        else
	   printf("fail to remove scratch files\n");

	MPI_Finalize();
	
	return retcode;
}

void checkBlastErrorExit(){
	if( running_blast )
		terminateProgram( -1 );
}

void terminateProgram( int sig ){
	if(sig==0) 
		cerr << "Sorry, mpiBLAST must be run on 2 or more nodes\n";
	else {
		cerr << rank << "\t" << (MPI_Wtime()-prog_start) 
		     << "\tBailing out with signal "<< sig << endl;
		removeFile( alias_filename );
		removeFile( local_query_tmpfile );
		removeFile( local_results_tmpfile );
		if( mpiblast_ptr != NULL )
			if( mpiblast_ptr->remove_db )
				mpiblast_ptr->cleanupDB(1);
	}
	MPI_Abort( MPI_COMM_WORLD, 0 );
}

MpiBlast::~MpiBlast(){ 
	log_stream = &cerr;
	profile_stream = &cerr;
	for( uint cleanupI = 0; cleanupI < cleanup_files.size(); cleanupI++ ){
		removeFile( cleanup_files[ cleanupI ] );
	}
}

/* Add by H.C.Lee for getting parameters from MpiBlast object */
string MpiBlast::getOutputFilename() {
        return output_filename;        
}

string MpiBlast::getProfileFilename() {
	ostringstream pro_file_str;
	pro_file_str << pro_phile_name << "." << rank;
        return pro_file_str.str();
}

string MpiBlast::getResultURL() {
        return result_url;        
}

string MpiBlast::getProfileURL() {
        return profile_url;        
}

// Send the results of the blast search (stored in b_filename, the *BINARY* ASN.1 file of 
// sequence alignments). To keep the
// memory requirements to a minimum, the file is returned to the master in chunks corresponding to
// a number of Seq-annot records (10 for now). Note that if records are returned to
// the master one at a time, MPI becomes prone to the "errno=110" or "Connection timed
// out" error when running on linux (at least with kernel 2.4.18). The problem occurs
// with both mpich (1.2.5) and LAM MPI (6.5.9). This error is
// allegedly with the linux tcp/ip implementation. I have found that this problem is
// sensitive to the *NUMBER* of messages sent (in addition to the size of each
// message). The scheme of collecting results into "chunks" before sending to the
// master seems to correct the problem. More work may be required if the problem
// returns.
void MpiBlast::sendResults( string& b_filename, int dest, int tag )
{	
	
	// The results are assumed to be in ASN.1 binary format.  Convert back to SeqAnnot
	// records and send to the master.
	AsnIoPtr aip;
	if( use_binary_asn )
		aip = AsnIoOpen((char*)(b_filename.c_str()), "rb");
	else
		aip = AsnIoOpen((char*)(b_filename.c_str()), "r");

	if(aip == NULL){
		throw __FILE__ "(MpiBlast::sendResults): AsnIoOpen has failed!";
	}

	SeqAnnotPtr sap;
	list<SeqAnnotPtr> sap_list;
	long buf_size = 0;
	
	long last_file_pos = 0;
	
	// How many records should we send at once? This is the "chunk" size mentioned
	// above and was set by hand to prevent timeout errors running under linux.
	const int rec_chunk_size = 10;
	
	while(1){
		sap = (SeqAnnotPtr)SeqAnnotAsnRead(aip, NULL);
		
		if(sap != NULL){
			sap_list.push_back(sap);
		
			// Get the current file position after reading the SeqAnnot structure
			long curr_file_pos = AsnIoTell(aip);
			
			// Compute the size of the SeqAnnot structure currently in memory
			buf_size += curr_file_pos - last_file_pos;
			
			// DEBUG
			//if((curr_file_pos - last_file_pos) != SeqAnnotSize(sap)){
			//	cerr << "Structure size error!" << endl;
			//	cerr << "curr_file_pos - last_file_pos = " << curr_file_pos - last_file_pos << endl;
			//	cerr << "SeqAnnotSize(sap)" << SeqAnnotSize(sap) << endl;
			//	throw __FILE__ "(MpiBlast::sendResults): Memory buffer size mismatch";
			//}
			
			last_file_pos = curr_file_pos;
		}
		
		// Have we either 
		// 	(a) read the last record (and sap == NULL) or 
		//	(b) read rec_chunk_size number of records
		// and need to send the results to the master.
		if((sap_list.size() > 0) && ((sap == NULL) || (sap_list.size() % rec_chunk_size == 0))){
		  	if(buf_size <= 0){
				throw __FILE__ "(MpiBlast::sendResults): Memory buffer has zero size";
			}
			
			// Prepare a memory buffer to store all of the SeqAnnot structures
			//unsigned char *buf = new unsigned char [buf_size];
			unsigned char *buf = (unsigned char*)calloc(buf_size, sizeof(unsigned char));
			
			if(buf == NULL){
				throw __FILE__ "(MpiBlast::sendResults): Unable to prepare memory buffer";
			}
			AsnIoMemPtr aimp;
			if( use_binary_asn )
				aimp = AsnIoMemOpen("wb", buf, buf_size);
			else
				aimp = AsnIoMemOpen("w", buf, buf_size);
			
			if(aimp == NULL){
				throw __FILE__ "(MpiBlast::sendResults): AsnIoMemOpen has failed!";
			}
			
			list<SeqAnnotPtr>::iterator sap_iter;
			
			// Fill the buffer with the SeqAnnot structures
			for(sap_iter = sap_list.begin(); sap_iter != sap_list.end();sap_iter++){
				
				if( !SeqAnnotAsnWrite( (*sap_iter), aimp->aip, NULL ) ){
					throw __FILE__ "(MpiBlast::sendResults): SeqAnnotAsnWrite has failed!";
				}
				
				// Clean up as we go ...
				SeqAnnotFree( (*sap_iter) );
			}
			AsnIoFlush( aimp->aip );

			// We now have a buffer filled with valid SeqAnnot structures to send to the master
			MPI_Send(buf, aimp->count, MPI_BYTE, 0, BLAST_RESULTS_TAG, MPI_COMM_WORLD);

			AsnIoMemClose(aimp);
			
		  	// Flush the list
			sap_list.clear();
			
			// Reset the buffer size
			buf_size = 0;
			
			// Clean up buffer
			free(buf);
		}
		
		// Terminate the while loop if we're read our last record
		if(sap == NULL){
			break;
		}
	}
	
	AsnIoClose(aip);	
}

void MpiBlast::receiveResults( int src, int tag, vector<SeqAlignPtr>& results, int m_result_count )
{
	MPI_Status status;
		
	if( debug_msg) {
		LOG_MSG << "receiveResults(" << src << ") to read " << m_result_count << " records" << endl;
	}
	
	// Allocate space for the incoming data
	results.clear();
	results.reserve( m_result_count );
	
	int rec_read = 0;
	
	while(rec_read < m_result_count){
		int resultsize;
		
		MPI_Probe( src, tag, MPI_COMM_WORLD, &status );
		
		MPI_Get_count( &status, MPI_BYTE, &resultsize );
		
		if(resultsize <= 0){
			throw __FILE__ "(MpiBlast::receiveResults): Results buffer has zero size";
		}
		
		//unsigned char *buf = new unsigned char [resultsize];
		unsigned char *buf = (unsigned char*)calloc(resultsize, sizeof(unsigned char));
		
		if(buf == NULL){
			cerr << "Unable to allocate memory for "
			     << resultsize
			     << " byte result buffer (rank = "
			     << rank
			     << ")"
			     << endl;
			throw __FILE__ "(MpiBlast::receiveResults): Unable to allocate memory";
		}
		
		MPI_Recv( buf, resultsize, MPI_BYTE, src, tag, MPI_COMM_WORLD, &status );
		// we now have the results in ASN.1 binary format.  Convert back to seqalign
		//ASNIO_TEXT | ASNIO_IN
		//AsnIoMemPtr aimp = AsnIoMemOpen( "rb", buf, resultsize );
		AsnIoMemPtr aimp;
		if( use_binary_asn )
			aimp = AsnIoMemOpen( "rb", buf, resultsize);
		else
			aimp = AsnIoMemOpen( "r", buf, resultsize);

		if(aimp == NULL){
			throw __FILE__ "(MpiBlast::receiveResults): AsnIoMemOpen has failed!";
		}
		
		SeqAnnotPtr sap;
		
		while( (sap = (SeqAnnotPtr)SeqAnnotAsnRead(aimp->aip, NULL)) ){
			// This a potentially risky design choice, but let's pack the rank of this worker into the
			// extended field of each SeqId contianed in the SeqAnnot structures that we are sending:
			// SeqAnnot->{SeqAlign->{SeqId, SeqId, ...}, SeqAlign->{SeqId, SeqId, ... }, ... }
			// 
			// This will cause problems if NCBI eventually uses the extended data field. A quick grep
			// through the NCBI source (as of 9/15/03) indicates that the extended field is not
			// currently used.
			//
			// This will cause additional problems when we eventually make mpiBLAST fault tolerant (i.e.
			// if rank "x" dies, and another rank "y" takes its place, then we must change all of the
			// embedded ranks).
			//
			// Node that we only have 8 bits to store the rank in! If there are more than
			// 254 worker nodes (255 total nodes) then we must resort to the slower SeqId
			// broadcast method.
			setSeqAnnotRank(sap, src, node_count);
			
				
			// Save the data
			results.push_back( (SeqAlignPtr)sap->data );
			
			// Free the un-needed memory associated witht the SeqAnnot structure (but protect the
			// SeqAlign structure contained within).
			
			sap->data = NULL;
			sap = SeqAnnotFree(sap);
			
			rec_read ++;
		}
								
		if( debug_msg) {
			LOG_MSG << "(" << src << ") read " << rec_read << " records in " << resultsize
				<< " bytes" << endl;
		}
	
		// Clean up 
		AsnIoMemClose(aimp);
		
		// Free the result buffer!
		free(buf);
	}
}

void MpiBlast::cleanupDB(int reason) {
#ifndef WIN32
	unixCleanupDB( reason );
#else
	cerr << "Cleanup database not implemented for Win32\n";
#endif
}

void MpiBlast::unixCleanupDB(int reason) {
#ifndef WIN32
	string pscmd = "/bin/ps -x | grep mpiblast | egrep -v 'grep|mpirun' | wc -l";
	char buf[16];
	FILE *ptr;
	int nmpiblasts;
	if(reason == 1) //force cleanup, we don't care what else is running
		nmpiblasts = 2;
	else {
		if ((ptr = popen(pscmd.c_str(), "r")) != NULL)
		fgets(buf, 15, ptr);
		nmpiblasts = atoi(buf);
		fclose(ptr);
	}
	if( debug_msg) {
		LOG_MSG << "Entered cleanupDB, reason "<< reason << " nmpiblasts: "<< nmpiblasts << endl;
		system(pscmd.c_str());
	}
//	if(nmpiblasts == 1) { //we're the last one
		string clean_db = localPath + "/" + database_name + ".??*.???  " + localPath + "/" + database_name + ".mbf" ;
	if( debug_msg )
		LOG_MSG << "Cleaning up database files: " << clean_db << endl ;
	removeFile( clean_db, true );
//	}
#endif
}

void MpiBlast::broadcastFile( string& b_filename ){
	ifstream b_file;
	char* b_buf;
	
	b_file.open( b_filename.c_str(), ios::binary );
	
	if( !b_file.is_open() ){
		cerr << "Error opening \"" << b_filename << "\"\n";
		
		throw __FILE__ "(MpiBlast::broadcastFile): Unable to open file";
	}
	
	b_file.seekg( 0, ios::end );
	
	uint b_filesize = b_file.tellg();
	
	if( debug_msg ){
		LOG_MSG << "broadcasting file size of " << b_filesize << endl;
	}
	
	MPI_Bcast( &b_filesize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
	
	if( debug_msg ){
		LOG_MSG << "file size broadcasted\n";
	}
	
	//b_buf = new char[ b_filesize ];
	b_buf = (char*)calloc(b_filesize, sizeof(unsigned char));
	
	if(b_buf == NULL){
		throw __FILE__ "(MpiBlast::broadcastFile): Unable to allocate memory";
	}
	
	b_file.seekg( 0, ios::beg );
	
	if( b_file.read( b_buf, b_filesize ) ){
		uint b_read = b_file.gcount();
		if( b_read != b_filesize ){
		
			free(b_buf);
			
			cerr << "Only read " << b_read << " of " << b_filesize << " bytes\n";
			cerr << "Error completely reading \"" << b_filename << "\"\n";
			
			throw __FILE__ "(MpiBlast::broadcastFile): Error completely reading file";
		}
		
		if( debug_msg ){
			LOG_MSG << "broadcasting file\n";
		}
		
		MPI_Bcast( b_buf, b_read, MPI_BYTE, 0, MPI_COMM_WORLD );
		
		if( debug_msg ){
			LOG_MSG << "file broadcasted\n";
		}
	}else{
		free(b_buf);
		
		cerr << "Error reading\"" << b_filename << "\"\n";
		
		throw __FILE__ "(MpiBlast::broadcastFile):Error reading file";
	}
	
	b_file.close();
	
	free(b_buf);
}

void MpiBlast::recvBroadcastFile( string& b_filename ){
	ofstream b_file;
	char* b_buf;

	if( debug_msg ){
		LOG_MSG << "waiting for file size broadcast\n";
	}
	
	uint b_filesize;
	
	MPI_Bcast( &b_filesize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
	
	if( debug_msg ){
		LOG_MSG << "received file size broadcast of " << b_filesize << endl;
	}
	
	//b_buf = new char[ b_filesize ];
	b_buf = (char*)calloc(b_filesize, sizeof(unsigned char));
	
	if(b_buf == NULL){
		throw __FILE__ "(MpiBlast::recvBroadcastFile): Unable to allocate memory";
	}
	
	if( debug_msg ){
		LOG_MSG << "opening receive file " << b_filename << endl;
	}
	
	b_file.open( b_filename.c_str() );

	if(!b_file.is_open()){
		free(b_buf);
		
		cerr << "Error opening file: " << b_filename.c_str() << endl;
		
		throw __FILE__ "(MpiBlast::recvBroadcastFile): Unable to open file";
	}
	
	if( debug_msg ){
		LOG_MSG << "receiving file to " << b_filename << endl;
	}
	
	MPI_Bcast( b_buf, b_filesize, MPI_BYTE, 0, MPI_COMM_WORLD );
	
	if( debug_msg ){
		LOG_MSG << "received file broadcast\n";
	}
	
	b_file.write( b_buf, b_filesize );
	b_file.close();
	
	free(b_buf);
}

void MpiBlast::master()
{
	map< int, vector< SeqAlignPtr > > results_map;

	running_blast = true;

	initNCBI( master_opts );
	int rval = initBLAST();
	running_blast = false;
	
	if( rval != 0 ){
		// bail out on error for now.
		LOG_MSG << "Error: got initBLAST exit code " << rval << endl;
		throw __FILE__ "(MpiBlast::master): Received non-zero initBLAST exit code";
	}

	if( debug_msg ) {
		LOG_MSG << "Init blast error code " << rval << endl;
	}
	
	// load the queries just to count how many there are
	int query_count = loadQueries();
	
	outputHtmlHeader();
	MPI_Status status;
	string master_merge_filename = config.localPath() + PATH_SEPARATOR + getFilename( output_filename ) + ".merge";
	ifstream master_merge_file( master_merge_filename.c_str() );
	if( master_merge_file.is_open() ){
		master_merge_file.close();
		removeFile( master_merge_filename );
	}

	set< int > unsearched_fragments;
	int finished_nodes = 0;
	for( int fragI = 0; fragI < db_spec.fragmentCount(); fragI++ ){
		unsearched_fragments.insert( fragI );
	}
	
	// get the list of fragments that each worker already has
	vector< vector< int > > node_fragments;
	vector< int > fvec;
	node_fragments.push_back( fvec ); // <-- Do we need this????

	for( int nodeI = 1; nodeI < node_count; nodeI++ ){
	
		RecvIntVec(fvec, nodeI, DB_FRAGMENTS_TAG);
		
		// Add the vector to the array of vectors even if there are no fragments to download
		node_fragments.push_back(fvec);
		
		if( debug_msg ){
			LOG_MSG << "Received fragments from node " << nodeI << endl;
		}
	}

	double broadcast_start = MPI_Wtime();
	if( profile_msg ){
		PROFILE_NEW_MSG << "Broadcast query took ";
	}

	// broadcast the query to every node
	broadcastFile( query_filename );
	double broadcast_end = MPI_Wtime();
	
	if( profile_msg ){
		PROFILE_MSG << broadcast_end - broadcast_start << endl;
		//PROFILE_NEW_MSG << "Broadcast query took " << broadcast_end - broadcast_start << endl;
	}
	
	blast_job = BlastJob( db_spec.fragmentCount(), node_fragments );

	// tell the workers what fragments to query
	// ...
	while( true ){
		//
		// wait for a message from a worker
		//
		MPI_Status probe_status;
		MPI_Probe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &probe_status );
		if( debug_msg ){
			LOG_MSG << "Master got message with tag " << probe_status.MPI_TAG 
				<< " from node " << probe_status.MPI_SOURCE << endl;
		}
		
		if( probe_status.MPI_TAG == FRAGMENT_COPY_COMPLETE ){
			// copy completed
			int new_fragment;
			MPI_Recv( &new_fragment, 1, MPI_INT, probe_status.MPI_SOURCE, 
				  FRAGMENT_COPY_COMPLETE, MPI_COMM_WORLD, &status );
				  
			blast_job.CopyCompleted( probe_status.MPI_SOURCE, new_fragment );
			
			continue;
		}
		
		if( probe_status.MPI_TAG == BLAST_RETCODE_TAG ){
			// blast exited, get results from 'dis hear node.
			int node_ret_code;
			
			MPI_Recv( &node_ret_code, 1, MPI_INT, probe_status.MPI_SOURCE, 
				  BLAST_RETCODE_TAG, MPI_COMM_WORLD, &status );
				  
			if( node_ret_code != 0 ){
				// bail out on error for now."
				cerr << "An error has occured!" << endl;
				cerr << "Got blast exit code " << node_ret_code << " from node " << probe_status.MPI_SOURCE << endl;
				
				throw __FILE__ "(MpiBlast::master): Received non-zero exit code from worker";
			}


			double recv_start = MPI_Wtime();
			// receive the list of successful queries
			int* result_list;
			int result_count;
			MPI_Recv( &result_count, 1, MPI_INT, probe_status.MPI_SOURCE, 
				  SUCCESSFUL_ID_COUNT, MPI_COMM_WORLD, &status );
			
			if( result_count > 0 ){
				if( profile_msg ){
					PROFILE_NEW_MSG << "Receive results from ";
				}

				//result_list = new int[ result_count ];
				result_list = (int*)calloc(result_count, sizeof(int));
				
				if(result_list == NULL){
					throw __FILE__ "(MpiBlast::master): Unable to allocate memory for result_list";
				}
				
				if( debug_msg ){
					LOG_MSG << "Receiving " << result_count << " successful query results (" 
						<< probe_status.MPI_SOURCE << ")\n";
				}
				
				MPI_Recv( result_list, result_count, MPI_INT, probe_status.MPI_SOURCE, 
					  SUCCESSFUL_IDS_TAG, MPI_COMM_WORLD, &status );

				vector< SeqAlignPtr > cur_results;

				running_blast = true;
				receiveResults( probe_status.MPI_SOURCE, BLAST_RESULTS_TAG, cur_results, result_count);
				running_blast = false;

				double recv_end = MPI_Wtime();

				if( profile_msg ){
					//PROFILE_NEW_MSG << "Receive results from " << probe_status.MPI_SOURCE << " took " << recv_end - recv_start << endl;
					PROFILE_MSG << probe_status.MPI_SOURCE << " took " << recv_end - recv_start << endl;
				}
				
				if( debug_msg ){
					LOG_MSG << "Receive was successful -- about to merge (" 
						<< probe_status.MPI_SOURCE << ")\n";
				}
				
				double merge_start = MPI_Wtime();

				if( profile_msg ){
					PROFILE_NEW_MSG << "Merge fragments from ";
				}
				
				// aggregate results into single list
				mergeResults( results_map, cur_results, result_list );
				
				if( debug_msg ){
					LOG_MSG << "Query results have been merged\n";
				}
				
				free(result_list);
				
				double merge_end = MPI_Wtime();

				if( profile_msg ){
					PROFILE_MSG << probe_status.MPI_SOURCE << " took " << merge_end - merge_start << endl;
					//PROFILE_NEW_MSG << "Merge fragments from " << probe_status.MPI_SOURCE << " took " << merge_end - merge_start << endl;
				}
			}
			
			continue;

		}
		
		if( probe_status.MPI_TAG == WORKER_IDLE ){
			int idle_code;

			MPI_Recv( &idle_code, 1, MPI_INT, probe_status.MPI_SOURCE, 
				  WORKER_IDLE, MPI_COMM_WORLD, &status );

			// this worker's bored, give it something to do...
			int operation;
			int fragment_id;

			blast_job.getAssignment( probe_status.MPI_SOURCE, operation, fragment_id );

			if( debug_msg ){
				LOG_MSG << "Giving node " << probe_status.MPI_SOURCE << " assignment " << operation << " on fragment " << fragment_id << endl;
			}

			MPI_Send( &operation, 1, MPI_INT, probe_status.MPI_SOURCE, ASSIGNMENT_TYPE, MPI_COMM_WORLD );

			if( operation == SEARCH_COMPLETE ){
				finished_nodes++;
				if( finished_nodes == node_count - 1 )
					break;
				continue;
			}

			// send the fragment to work on.
			MPI_Send( &fragment_id, 1, MPI_INT, probe_status.MPI_SOURCE, 
				   DB_FRAGMENTS_TAG, MPI_COMM_WORLD );
					   
			continue;
		}
		
		// Catch any out of place messages
		cerr << "Error: Received msg " << probe_status.MPI_TAG << " from " << probe_status.MPI_SOURCE << endl;
		throw __FILE__ "Master(): Illegal message in main loop!";
	}

	// print out the results
#ifdef MPE
	MPE_Log_event(mergestart,0,"start merge");
#endif

	if( debug_msg ) {
		LOG_MSG << "<<<<<<<<<<< Start merge >>>>>>>>>>>" << endl;
	}
	
	// By default, do not use MPI to distribute sequence database look ups.
	if(use_mpi_db){
		if( debug_msg ) {
			LOG_MSG << "Enabling MPI Distributed DB queries" << endl;
		}
	
		if(!Enable_MPI_Distributed_DB()){
			cerr << "Error enabling MPI distributed DB" << endl;
			throw __FILE__ ":master: Error enabling MPI distributed DB";
		}
	}
	
	// DEBUG
	MpiBlastEnableBioseqFetch();
	
	double profile_time = MPI_Wtime();
	
	if(profile_msg){
		PROFILE_NEW_MSG << "Output results took ";
	}

	map< int, vector< SeqAlignPtr > >::iterator result_iter = results_map.begin();
	int cur_query_id = 0;
	for( ; result_iter != results_map.end(); result_iter++ ){

		// link up all the results for this query into a single list
		SeqAlignPtr cur_result_head = (*result_iter).second[0];
		SeqAlignPtr cur_result_tail = cur_result_head;
		
		while(cur_result_tail->next != NULL){
			cur_result_tail = cur_result_tail->next;
		}
			
		for( uint resultI = 1; resultI < (*result_iter).second.size(); resultI++ ){
			cur_result_tail->next = (*result_iter).second[ resultI ];
			
			while(cur_result_tail->next != NULL){
				cur_result_tail = cur_result_tail->next;
			}
		}
		
		running_blast = true;
		for( ; cur_query_id < (*result_iter).first; cur_query_id++ ){
			outputResults( NULL ); // there were no results for this query
		}
		
		// output the list of results
		outputResults( cur_result_head );
		running_blast = false;
		
		// After everything has been debugged, we need to free the memory associated with
		// cur_result_head:
		// cur_result_head = cur_result_head = SeqAlignFree(cur_result_head);
		cur_query_id++;
	}
	
	// output empty results for the remaining queries
	running_blast = true;
	for( ; cur_query_id < query_count; cur_query_id++ ){
		outputResults( NULL ); // there were no results for this query
	}
	running_blast = false;

	// DEBUG
	MpiBlastDisableBioseqFetch();
	
	profile_time = MPI_Wtime() - profile_time;
	
	if(profile_msg){
		PROFILE_MSG << profile_time << endl;
		//PROFILE_NEW_MSG << "Output results took " << profile_time << " sec" << endl;
	}
	
	outputHTMLfooter();
	
	if(use_mpi_db){
		if(!Disable_MPI_Distributed_DB()){
			cerr << "Error disabling MPI distributed DB" << endl;
		}
	}
	
	if( debug_msg ) {
		LOG_MSG << "*********** Finished merge ***********" << endl;
	}
	
	// Tell all of the workers to quit
	for(int i = 1;i < node_count;i++){
		int msg = WORKER_QUIT;

		MPI_Send( &msg, 1, MPI_INT, i, ASSIGNMENT_TYPE, MPI_COMM_WORLD );
	}
	
#ifdef MPE
	MPE_Log_event(mergeend,0,"end merge");
#endif
	cleanupBLAST();
	if( remove_db ){
		MPI_Barrier(MPI_COMM_WORLD);
		cleanupDB(1);
	}
}

void MpiBlast::worker(){

	double profile_time = MPI_Wtime();
		
	if( profile_msg ){
		PROFILE_NEW_MSG << "File setup took ";
	}

	//
	// create blast alias file
	//    
	alias_basename = config.localPath() + PATH_SEPARATOR + database_name + "XXXXXX";
	getTempFileName( alias_basename );
	alias_filename = alias_basename + "." + db_type + "al";
	moveFile( alias_basename, alias_filename );

	// create blast results output file
	stringstream rank_str;
	rank_str << rank;
	string local_results_file = config.localPath() + PATH_SEPARATOR + getFilename( output_filename ) + "." + rank_str.str();
	
	local_results_tmpfile = local_results_file + "XXXXXX";
	getTempFileName( local_results_tmpfile );

	// create local query tempfile
	local_query_tmpfile = local_query_filename + "XXXXXX";
	getTempFileName( local_query_tmpfile );

	// create blast command line
	addOpt( worker_opts, 'i', local_query_tmpfile.c_str() );
	if( use_binary_asn )
		addOpt( worker_opts, 'm', "11" );	/** ask for binary SeqAlign output */
	else
		addOpt( worker_opts, 'm', "10" );	/** ask for text SeqAlign output */
	addOpt( worker_opts, 'o', local_results_tmpfile.c_str() );
	addOpt( worker_opts, 'd', alias_basename.c_str() );
	addOpt( worker_opts, 'J', NULL );

	profile_time = MPI_Wtime() - profile_time;
	
	if( profile_msg ){
		PROFILE_MSG << profile_time << endl;
		//PROFILE_NEW_MSG << "File setup took " << profile_time << endl;
	}
		
	profile_time = MPI_Wtime();
	if( profile_msg ){
		PROFILE_NEW_MSG << "initNCBI took ";
	}
	
	// Init the NCBI library with the contrived blast command line
	running_blast = true;
	initNCBI( worker_opts );
	running_blast = false;
	
	profile_time = MPI_Wtime() - profile_time;
	
	if( profile_msg ){
		PROFILE_MSG << profile_time << endl;
		//PROFILE_NEW_MSG << "initNCBI took " << profile_time << endl;
	}
	
	profile_time = MPI_Wtime();

	if( profile_msg ){
		PROFILE_NEW_MSG << "Fragment list to master took ";
	}
	
	MPI_Status status;
	
	// These variables manage the database fragment list (for sending, receiving and writing
	// to disk).
	vector<int> fragment_list;
	int frag_id;
	
	// A list of ALL fragments processed by this worker during this run
	vector<int> grand_fragment_list;
	
	// tell the master node what fragments this node has
	string fragment_filename = config.localPath() + PATH_SEPARATOR + 
		database_name + FRAG_LIST_EXTENSION;
	
	frag_list_file = FragmentListFile( fragment_filename );
		
	frag_list_file.SendList( config, database_name, db_type );
		
	if( debug_msg ){
		LOG_MSG << "Fragment list sent." << endl;
	}
	
	profile_time = MPI_Wtime() - profile_time;
	
	if( profile_msg ){
		PROFILE_MSG << profile_time << endl;
		//PROFILE_NEW_MSG << "Fragment list to master took " << profile_time << endl;
	}
	
	profile_time = MPI_Wtime();

        if( profile_msg ) {
                PROFILE_NEW_MSG << "Receiving query took ";
        }	
	// receive the query file broadcast
	recvBroadcastFile( local_query_tmpfile );
	
	profile_time = MPI_Wtime() - profile_time;
	
	if( profile_msg ){
		PROFILE_MSG << profile_time << endl;
	}
	
	profile_time = MPI_Wtime();
	
	if( debug_msg ){
		LOG_MSG << "Query file received as " << local_query_filename << endl;
	}
	
	// loop until the master tells us we're finished.
	bool tell_master_idle = true;
	bool first_search = true;
	while( true ){
		MPI_Status probe_status;
		
		if(tell_master_idle){
			
			// tell master we are idle
			int idle_msg = WORKER_IDLE;
			MPI_Send( &idle_msg, 1, MPI_INT, 0, WORKER_IDLE, MPI_COMM_WORLD );
		
			if( debug_msg ){
				LOG_MSG << "Idle message sent" << endl;
			}
		}
		
		// get a fragment assignment or exit peacefully
		int assignment = 0;

		MPI_Probe( 0, MPI_ANY_TAG, MPI_COMM_WORLD, &probe_status );
		
		switch(probe_status.MPI_TAG){
			case ASSIGNMENT_TYPE:
				
				MPI_Recv( &assignment, 1, MPI_INT, 0, ASSIGNMENT_TYPE, MPI_COMM_WORLD, &status );
				break;
			case MPI_BLAST_FETCH_GI:
			case MPI_BLAST_FETCH_LOCAL:
				/* Only receive requests from the master */
				if( debug_msg )
					LOG_MSG << "Looking up bioseq\n";
				MPI_BSLookup(0, &probe_status);
				
				assignment = DO_NOTHING;
			
				break;
			default:
				cerr << "Unknown message tag (" << probe_status.MPI_TAG
				     << ") received by " << rank << endl;
				
				// An error has occured (i.e. unknwon message) so set the termination flag
				// to clean up as best we can.
				assignment = WORKER_QUIT;
				
				break;
			
		};
		
		if(assignment == DO_NOTHING){
			continue;
		}
		
		if( debug_msg ){
			LOG_MSG << "Assignment received is " << assignment << endl;
		}
		
		if( assignment == SEARCH_COMPLETE ){
			
			// If using MPI to distribute database lookups, workers should not quit after
			// receiveing a SEARCH_COMPLETE message.
			
			tell_master_idle = false;
			
			// This worker will not receive any more fragments. Write the final fragment file that
			// contains all of the fragments that have been processed to date
			if( grand_fragment_list.size() > 0 ){
				WriteAliasFile(alias_filename, grand_fragment_list);
			
				// DEBUG
				if( debug_msg )	
					LOG_MSG << "Enabling Bioseq Fetch\n";
				MpiBlastEnableBioseqFetch();
				if( debug_msg )	
					LOG_MSG << "Bioseq Fetch Enabled\n";
			}

			continue;
		}
		
		if( assignment == WORKER_QUIT ){
			// nothing left for this worker to do...
			if( remove_db ){
				// DEBUG
				MpiBlastDisableBioseqFetch();
				
				MPI_Barrier(MPI_COMM_WORLD);
				cleanupDB(0);
			}

			break;
		}
		
		
		MPI_Recv(&frag_id, 1, MPI_INT, 0, 
			DB_FRAGMENTS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		
		if( debug_msg ){
			LOG_MSG << "got fragment assignment " << frag_id << endl;
		}
		
		if( assignment == COPY_FRAGMENT ){
			
			double copy_start = MPI_Wtime();

			// copy the database fragments 
			// .nhr, .nin, .nsq, .nnd, .nni, .nsd, .nsi for nucleotide
			// .phr, .pin, .psq, .pnd, .pni, .psd, .psi for amino acids
			vector<string> extensions;
			vector<string> extensions_optional;
			extensions = FragmentExtensions(db_type);
			extensions_optional = FragmentExtensions_optional(db_type);
			if( profile_msg ){
				PROFILE_NEW_MSG << "Copy fragments";
			}
			
			
			char fragid[8];
			memset( fragid, 0, 8 );
#ifdef MANY_FRAGMENTS_HACK
			sprintf( fragid, "%03d", frag_id );
#else
			sprintf( fragid, "%02d", frag_id );
#endif
			if( debug_msg ){
				LOG_MSG << "copying fragment " << fragid << endl;
			}

			bool copy_success = true;
#ifdef MPE
			MPE_Log_event(cpstart,0,"start copy");
#endif

			if(compress_io) {

				string copy_fragment_name = database_name + "." + fragid + ".tgz";
				//LOG_MSG << copy_fragment_name << endl;

				string copy_src = config.sharedPath() + PATH_SEPARATOR + copy_fragment_name;
				string copy_dest = config.localPath() + PATH_SEPARATOR + copy_fragment_name;
				//LOG_MSG << copy_src << endl;
				//LOG_MSG << copy_dest << endl;
// Modified for switching gass copy
#ifdef GASSCOPY
				copy_src = config.sharedGassHost() + PATH_SEPARATOR + copy_src;
                
				if( gassCopyFile( copy_src, "file://"+copy_dest) != 0 )
#else                   
				if( copyFile( copy_src, copy_dest ) != 0 )
#endif                  
					copy_success = false;

				if(copy_success) { //untar db fragment only when copy success
					untarFile(copy_dest,config.localPath());
                                        removeFile(copy_dest);
				}
                        }
                        else {
				for( uint extI = 0; extI < extensions.size(); extI++ ){
					//string copy_fragment_name = database_name + ".";
					string copy_fragment_name = database_name + ".";
					copy_fragment_name += fragid + extensions[ extI ];
					string copy_src = config.sharedPath() + PATH_SEPARATOR + copy_fragment_name;
					string copy_dest = config.localPath() + PATH_SEPARATOR + copy_fragment_name;
// Modified for switching gass copy
#ifdef GASSCOPY
					copy_src = config.sharedGassHost() + PATH_SEPARATOR + copy_src;
					copy_dest = "file://" + copy_dest;
                        
					if( gassCopyFile( copy_src, copy_dest) != 0 )
#else                   
					if( copyFile( copy_src, copy_dest ) != 0 )
#endif                  
						copy_success = false;
				}
                        
				for( uint extI = 0; extI < extensions_optional.size(); extI++ ){
					string copy_fragment_name = database_name + ".";
					copy_fragment_name += fragid + extensions_optional[ extI ];
					string copy_src = config.sharedPath() + PATH_SEPARATOR + copy_fragment_name;
					string copy_dest = config.localPath() + PATH_SEPARATOR + copy_fragment_name;

// Modified for switching gass copy
#ifdef GASSCOPY
					copy_src = config.sharedGassHost() + PATH_SEPARATOR + copy_src;
					copy_dest = "file://" + copy_dest;
                                
					gassCopyFile( copy_src, copy_dest );
#else                           
					if( doesFileExist( copy_src ) )
						copyFile( copy_src, copy_dest );
#endif                          
				}
			}
#ifdef MPE
			MPE_Log_event(cpend,0,"end copy");
#endif
			if( copy_success ){
				frag_list_file.addFragment( frag_id );
			}
			else{
				cerr << "(" << rank << ") unable to copy fragment!" << endl;
				MPI_Abort( MPI_COMM_WORLD, -1 );
			}

			MPI_Send(&frag_id, 1, MPI_INT, 0, 
				   FRAGMENT_COPY_COMPLETE, MPI_COMM_WORLD );

			if( profile_msg ){
				PROFILE_MSG << " " << frag_id;
			}
			
			
			double copy_end = MPI_Wtime();

			if( profile_msg ){
				PROFILE_MSG << " took " << copy_end - copy_start << endl;
			}
			
			continue;
		}
		
		//profile_time = MPI_Wtime() - profile_time;
	
		//if( profile_msg ){
		//		PROFILE_MSG << "Copy frag time: " << profile_time << endl;
		//	}

		// else assignment is SEARCH_FRAGMENT
		// process fragment assignment
		profile_time = MPI_Wtime();
		//if( profile_msg ){
		//	PROFILE_NEW_MSG << "Alias file construction took ";
		//}

		//
		// fill alias file with the names of fragments to search.
		//    
		vector<int> single_frag_vector(1, frag_id);
		WriteAliasFile(alias_filename, single_frag_vector);
				
		// Save a record of each fragment processed 
		grand_fragment_list.push_back(frag_id);
				
		profile_time = MPI_Wtime() - profile_time;
		
		//if( profile_msg ){
		//	PROFILE_MSG << profile_time << endl;
		//}
		
		profile_time = MPI_Wtime();
		if( profile_msg ){
			PROFILE_NEW_MSG << "InitBlast() took ";
		}
		
		//
		// execute blast
		//
		if( debug_msg ){
			LOG_MSG << "node " << rank << " Executing BLAST search\n";
		}
		ValNodePtr results_id_list = NULL;
		
		// initialize BLAST with the current data set
		running_blast = true;
		int retcode = initBLAST();
		running_blast = false;
		
		profile_time = MPI_Wtime() - profile_time;
		
		if( profile_msg ){
			PROFILE_MSG << profile_time << endl;
			//PROFILE_NEW_MSG << "InitBlast() took " << profile_time << endl;
		}
		
		profile_time = MPI_Wtime();
		if( profile_msg ){
			PROFILE_NEW_MSG << "runBLAST() took ";
		}
		
		if( retcode == 0 ) {
			if( first_search ){
				loadQueries();
				first_search = false;
			}
#ifdef MPE
			MPE_Log_event(blaststart,0,"start blastall");
#endif
			running_blast = true;
			retcode = runBLAST( &results_id_list );
			running_blast = false;
#ifdef MPE
			MPE_Log_event(blastend,0,"end blastall");
#endif
		}
		else{
			// bail out on error for now.
			LOG_MSG << "Error: got initBLAST exit code " << retcode << endl;
			throw __FILE__ "(MpiBlast::worker): Received non-zero initBLAST exit code ";
		}
		
		profile_time = MPI_Wtime() - profile_time;
		
		if( profile_msg ){
			PROFILE_MSG << profile_time << endl;
			//PROFILE_NEW_MSG << "runBLAST() took " << profile_time << endl;
		}
		
		profile_time = MPI_Wtime();
		if( profile_msg ){
			PROFILE_NEW_MSG << "cleanupBLAST() took ";
		}
		
		running_blast = true;
		cleanupBLAST();
		running_blast = false;
		
		profile_time = MPI_Wtime() - profile_time;
		
		if( profile_msg ){
			PROFILE_MSG << profile_time << endl;
			//PROFILE_NEW_MSG << "cleanupBLAST() took " << profile_time << endl;
		}
		
		profile_time = MPI_Wtime();

		if( profile_msg ){
			PROFILE_NEW_MSG << "Return results to master took ";
		}
		
		if( debug_msg ){
			LOG_MSG << "node " << rank << " blast retcode= " << retcode << endl;
		}

		//
		// message master that blast has completed
		//
		MPI_Send( &retcode, 1, MPI_INT, 0, BLAST_RETCODE_TAG, MPI_COMM_WORLD );

		if( retcode != 0 ){
			LOG_MSG << "blastall exited with code: " << retcode << endl;
			MPI_Abort( MPI_COMM_WORLD, -1 );
		}

		// send list of queries that had hits
		int query_result_count = 0;
		for( ValNodePtr r_iter = results_id_list; r_iter != NULL; r_iter = r_iter->next ){
			query_result_count++;
		}
		
		if( debug_msg ){
			LOG_MSG << "Sending " << query_result_count << " successful query results\n";
		}
		
		// Don't send until the master is ready to receive this data
		MPI_Ssend( &query_result_count, 1, MPI_INT, 0, SUCCESSFUL_ID_COUNT, MPI_COMM_WORLD );
//		MPI_Send( &query_result_count, 1, MPI_INT, 0, SUCCESSFUL_ID_COUNT, MPI_COMM_WORLD );
		
		if( query_result_count > 0 ){
			//int* query_results = new int[ query_result_count ];
			int* query_results = (int*)calloc(query_result_count, sizeof(int));
			
			if(query_results == NULL){
				throw __FILE__ "(MpiBlast::worker): Unable to allocate memory for query_results";
			}
		
			int resultI = 0;
			for( ValNodePtr r_iter = results_id_list; r_iter != NULL; r_iter = r_iter->next ){
				//query_results[ resultI++ ] = r_iter->data.intvalue;
				query_results[ resultI ] = r_iter->data.intvalue;
				
				resultI++;
			}
			
			MPI_Send( query_results, query_result_count, MPI_INT, 0, 
				   SUCCESSFUL_IDS_TAG, MPI_COMM_WORLD );
			
			free(query_results);
			
			//
			// send master the blast results
			//
			sendResults( local_results_tmpfile, 0, BLAST_RESULTS_TAG );
		}
		profile_time = MPI_Wtime() - profile_time;

		if( profile_msg ){
			PROFILE_MSG << profile_time << endl;
			//PROFILE_NEW_MSG << "Return results to master took " << profile_time << endl;
		}
	}
	
	profile_time = MPI_Wtime();
	if( profile_msg ){
		PROFILE_NEW_MSG << "File clean up took ";
	}

	// free query sequence memory
	cleanupQueries();
		
	// clean up temporary files
	removeFile( alias_filename );
	removeFile( local_query_tmpfile );
	removeFile( local_results_tmpfile );
	
	profile_time = MPI_Wtime() - profile_time;
		
	if( profile_msg ){
		PROFILE_MSG << profile_time << endl;
		//PROFILE_NEW_MSG << "File clean up took " << profile_time << endl;
	}
}

void  MpiBlast::WriteAliasFile(string& alias_filename, vector< int >& fragment_list)
{
	if(fragment_list.size() == 0){
		return;	// don't write a file
		// this isn't an error condition when # workers > # db frags
//		throw __FILE__ "(WriteAliasFile) Empty fragment_list";
	}
	
	ofstream alias_file( alias_filename.c_str() );

	if( !alias_file.is_open() ){
		cerr << "Error opening " << alias_filename << endl;
		throw __FILE__ "(MpiBlast::worker): Unable to open alias file";
	}

	alias_file << "TITLE " << config.localPath() << PATH_SEPARATOR << database_name << endl;
	alias_file << "DBLIST";

	if( profile_msg ){
		PROFILE_NEW_MSG << "Search fragments";
	}
	
	list<int>::iterator iter;
	
	for( uint iter = 0; iter != fragment_list.size(); iter++ ){
		alias_file << " " << database_name << ".";
		char fragid[8];
		
		memset(fragid, 0, 8);
		
#ifdef MANY_FRAGMENTS_HACK
		sprintf( fragid, "%03d", fragment_list[ iter ] );
#else
		sprintf( fragid, "%02d", fragment_list[ iter ] );
#endif

		alias_file << fragid;
		if( debug_msg ){
			LOG_MSG << "node " << rank << " fragid: " << fragid << endl;
		}

		if( profile_msg ){
			PROFILE_MSG << " " << fragid;
		}
	}

	if( profile_msg ){
		PROFILE_MSG << endl;
	}

	alias_file << endl;
	alias_file.close();
		
}

// need getpid for debugging under win32
//#include <process.h>

int MpiBlast::main( int argc, char* argv[] )
{
	// for debugging under win32, wait so the process can be attached.
//	LOG_MSG << "Process id: " << _getpid() << endl;
//	Sleep( 20000 );

	if( argc == 0 ){
		cerr << "Incorrect usage\n";
		MPI_Abort( MPI_COMM_WORLD, -1 );
		throw __FILE__ "(MpiBlast::main): Incorrect usage";
	}
	exec_path = argv[0];
	exec_path = getPath( exec_path );

	//
	// must read -d, -i, and -o options
	//
	const char* options = "-d:i:o:p:m:J:O:";
	int ac = argc;
	char** av = argv;
	int opt_code, opt_type = 0;
	opterr = 0;
	remove_db = false;
	string log_file_name;
	//string pro_phile_name;
#ifdef GASSCOPY
	string wdir_url;
	string query_url;
	/* string result_url; */ //Should be global
#endif
	string config_f_name = MpiBlastConfig::defaultConfigFileName();

	int config_opt = 0;
	int long_index = 0;

	master_opts.push_back( "blastall" );
	worker_opts.push_back( "blastall" );
	struct option long_opts[] = {
		{"config-file", 1, &config_opt, 1 },
		{"log-file", 1, &config_opt, 2 },
		{"debug", 0, &config_opt, 3 },
		{"pro-phile", 1, &config_opt, 4 },
		{"removedb", 0, &config_opt, 5 },
		{"version", 0, &config_opt, 6 },
//		{"bin-seqalign", 0, &config_opt, 7 },	// don't need this, users can use -m 11 output option
		{"disable-mpi-db", 0, &config_opt, 7 },
		{"compress-io", 0, &config_opt, 8 },  // using compress-io or not
#ifdef GASSCOPY
		{"wdir-url", 1, &config_opt, 9 },  // the url for getting query sequence
#endif
		{0,0,0,0}	// signifies termination of option list
	};
	
	while((opt_code = getopt_long( ac, av, options, long_opts, &long_index ) ) != EOF ){
		switch( opt_code ){
			case 0:
				if( config_opt == 1 ){
					config_f_name = optarg;
				}
				
				if( config_opt == 2 ){
					log_file_name = optarg;
				}
				
				if( config_opt == 3 ){
					debug_msg = true;
					debug_bsfetch = true;
					debug_bslookup = true;
				}
				
				if( config_opt == 4 ){
					pro_phile_name = optarg;
					profile_msg = true;
				}
				
				if( config_opt == 5 ){
					remove_db = true;
				}
				
				if( config_opt == 6 ){
					cout << PACKAGE << " version " << VERSION << endl;
				}
								
				if( config_opt == 7 ){
					// Allow the user to disable MPI distributed database queries
					use_mpi_db = false;
				}
				if( config_opt == 8 ){
					// Using compress IO 
					compress_io = true;
				}
#ifdef GASSCOPY
				if( config_opt == 9 ){
					// Allow the user to disable MPI distributed database queries
					wdir_url = optarg;
				}
#endif
				
				break;

			case 'i':
				query_filename = optarg;
				addOpt( master_opts, opt_code, optarg );
				break;
			case 'o':
				output_filename = optarg;
				addOpt( master_opts, opt_code, optarg );
				break;
			case 'd':
				database_name = optarg;
				break;
			case 'p':
				addOpt( master_opts, opt_code, optarg );
				addOpt( worker_opts, opt_code, optarg );
				blast_type = optarg;
				break;
			case 'm':
				// worker always uses -m 11
				addOpt( master_opts, opt_code, optarg );
				break;
			case 'J':
				// worker already uses -m 11
				addOpt( master_opts, opt_code, optarg );
				break;
			case 'O':
				// Only the master is allowed to see this
				// flag since workers must use -m 11. Don't let the
				// -O force ASCII ASN.1 output on the workers! 
				addOpt( master_opts, opt_code, optarg );
				break;
			case 1:
				worker_opts.push_back( optarg );
				master_opts.push_back( optarg );
				break;
			case '?':
				addOpt( worker_opts, optopt, optarg );
				addOpt( master_opts, optopt, optarg );
				opt_type = optopt;
				break;
			default:
				addOpt( worker_opts, opt_type, optarg );
				addOpt( master_opts, opt_type, optarg );
				opt_type = optopt;
				break;
		}
	}

	try{
		config = MpiBlastConfig( config_f_name );		
		localPath = config.localPath();

		db_type = "n"; //Default to nucleotide
		if( blast_type == "blastp" || blast_type == "blastx" ){
			db_type = "p";
		}
				
		// set location of substitution matrices if not set by the user
		if( getenv( "BLASTMAT" ) == NULL ){
#ifdef WIN32
			string blastmat_env = "BLASTMAT=" + config.sharedPath();
			_putenv( blastmat_env.c_str() );
#else
			setenv("BLASTMAT", config.sharedPath().c_str(), 1);
#endif
		}
				
		//
		// open the log file
		//
		if( log_file_name != "" ){
			ostringstream log_file_str;
			log_file_str << log_file_name << "." << rank;
			log_file.open( log_file_str.str().c_str() );
			
			if( !log_file.is_open() ){
				cerr << "Error opening log file \"" << log_file_name << "\"\n";
				throw __FILE__ "(main): Unable to open log file";
			}
			
			if( debug_msg ){
				cerr << "logging to " << log_file_str.str() << endl;
			}
			
			log_stream = &log_file;
		}

		//
		// open the profile file
		//
		if( profile_msg ){
			ostringstream pro_file_str;
			pro_file_str << pro_phile_name << "." << rank;
			pro_phile.open( pro_file_str.str().c_str() );
			
			if( !pro_phile.is_open() ){
				cerr << "Error opening profile file \"" << pro_phile_name << "\"\n";
				throw __FILE__ "(main): Unable to open profile file";
			}
			
			profile_stream = &pro_phile;
		}else{
			profile_stream = &cerr;
		}

		query_file = getFilename( query_filename );
		local_query_filename = config.localPath() + PATH_SEPARATOR +  query_file;
#ifdef GASSCOPY
		string result_file = getFilename( output_filename );
		query_url   = wdir_url + PATH_SEPARATOR + query_file;
		result_url  = wdir_url + PATH_SEPARATOR + result_file;

	        stringstream rank_str;
	        rank_str << rank;
                profile_url = wdir_url + PATH_SEPARATOR + pro_phile_name + "." + rank_str.str();
#endif
		
		// remove the database if the user requested it...
		// doing it here allows the user to cleanup a database without running a query
		if( remove_db ){
			MPI_Barrier(MPI_COMM_WORLD);
			cleanupDB(1);
			remove_db = false;
		}

		//
		// Read database specification from shared storage
		//
		// Old version
		// db_spec_filename = config.sharedPath() + PATH_SEPARATOR + database_name + ".dbs";
		// 
		// Since the .dbs file may only be accessible via rcp, first copy the file to local storage and then
		// read its contents
		long double db_size;	// this is really an integer, but MPI doesn't explicitly support 64-bit ints
		
		if(rank == 0){
			// check if we need to copy the dbs file...
			db_spec_filename = config.localPath() + PATH_SEPARATOR + database_name + ".dbs";
			string::size_type src_pos = db_spec_filename.find(":", 0);

//Add for globus_gass_copy by H.C.Lee
#ifdef GASSCOPY
			db_spec_filename = config.localPath() + PATH_SEPARATOR + database_name + ".dbs";
				
			if( gassCopyFile( config.sharedGassHost() + PATH_SEPARATOR + config.sharedPath() + PATH_SEPARATOR + database_name + ".dbs", "file://"+db_spec_filename) != 0){
			     cerr << "(" << rank << ") Unable to gass copy " << database_name << ".dbs to local storage" << endl;
			     throw __FILE__ "(main): Unable to copy .dbs to local storage";
			}
			db_spec = DbSpecFile( db_spec_filename );

			// Clean up as we go ...
			unlink( db_spec_filename.c_str() );
#else
			if( src_pos == string::npos ){
				db_spec_filename = config.sharedPath() + PATH_SEPARATOR + database_name + ".dbs";
				db_spec = DbSpecFile( db_spec_filename );
			}else{
				db_spec_filename = config.localPath() + PATH_SEPARATOR + database_name + ".dbs";
			
				if( copyFile( config.sharedPath() + PATH_SEPARATOR + database_name + ".dbs", db_spec_filename) != 0){
				     cerr << "(" << rank << ") Unable to copy " << database_name << ".dbs to local storage" << endl;
				     throw __FILE__ "(main): Unable to copy .dbs to local storage";
				}
				db_spec = DbSpecFile( db_spec_filename );

				// Clean up as we go ...
				unlink( db_spec_filename.c_str() );
			}
#endif		
			db_size = db_spec.dbSize();
			
			// Send the db size to the worker
			MPI_Bcast( &db_size, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD );
		}
		else{
			// Get the db size from the master
			MPI_Bcast( &db_size, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD );
			
			stringstream db_size_ss;
			
			db_size_ss << (uint64)db_size;
			
			//
			// add blast db size to command line so e-values are correct
			//
			addOpt( worker_opts, 'z', db_size_ss.str().c_str() );
		}
		
		

		//
		// add database name for master
		//
		// The old version
		// string shared_db_name = config.sharedPath() + PATH_SEPARATOR + database_name;
		// Assume that the database resides on either the master, or a valid nfs mount
		string shared_db_name;
		string::size_type pos = config.sharedPath().find(":", 0);
		if( pos == string::npos){
			// NFS mount
			shared_db_name = config.sharedPath() + PATH_SEPARATOR + database_name;
		}
		else{
			// On master
			pos ++;
			shared_db_name = config.sharedPath().substr(pos, config.sharedPath().length() - pos)
				+ PATH_SEPARATOR + database_name;
		}
		
		addOpt( master_opts, 'd', shared_db_name.c_str() );

		if( rank == 0 ){
			// Where does the query reside?
			pos = query_filename.find(":", 0);
			string local_query;

			int clean_up_query = false;
#ifdef GASSCOPY
			if( gassCopyFile(query_url, "file://"+local_query_filename) != 0 ) {
				cerr << "Unable to get " << query_url << " to "
				     << local_query_filename << endl;
					     
				throw __FILE__ "(main): Unable to get query to master";
			}
			// Make sure the master now reads from the local copy
			query_filename = local_query_filename;

			// Overwrite the query value in the masters argument list
			addOpt( master_opts, 'i', local_query_filename.c_str() );
			clean_up_query = true;
#else
			if( pos == string::npos){
				// The query is either (a) on the masters local disk or (b) visible to the master
				// via an NFS mount. No path name adjustment is required.

			}
			else{
				// The query resides on another host (not local to the master). Use rcp to copy the query to the
				// masters local disk.
				if( copyFile( query_filename, local_query_filename ) != 0 ){
					cerr << "Unable to copy " << query_filename << " to "
					     << local_query_filename << endl;
					     
					throw __FILE__ "(main): Unable to copy query to master";
				}

				// Make sure the master now reads from the local copy
				query_filename = local_query_filename;

				// Overwrite the query value in the masters argument list
				addOpt( master_opts, 'i', local_query_filename.c_str() );
				clean_up_query = true;
			}
#endif				

			master();

			if(clean_up_query){
				// Remove the the copy of the query on the masters local disk
				unlink( local_query_filename.c_str() );
			}
		}
		else {
			// Clean up existing crap
			worker();
		}

#ifdef GASSCOPY
                //Make sure all workers and master are reach this point
                MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
	catch(const char *error){
		cerr << "Fatal Error:" << endl;
		cerr << error << endl;
		MPI_Abort( MPI_COMM_WORLD, -1 );
		
		return -1;
	}catch(exception& e ){
		cerr << "Fatal Exception!" << endl;
		cerr << e.what() << endl;
		MPI_Abort( MPI_COMM_WORLD, -1 );
		return -1;
	}
/*	catch(...){
		cerr << "Fatal Unhandled Exception!" << endl;
		MPI_Abort( MPI_COMM_WORLD, -1 );
		
		return -1;
	}
*/	
	cleanupNCBI();
	return 0;
}
