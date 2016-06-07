/******************************************************************************
** mpiblast_util.cpp
** Provides utility functions used by mpiBLAST.
** $Id: mpiblast_util.cpp,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
** $Revision: 1.1.1.2 $
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

#include "mpiblast_util.hpp"
#include "mpiblast_types.h"

extern "C"{
#include "ncbi.h"
#include "objseq.h"
#include "blast_hooks.h"
#include "ncbithr.h"
#include "ncbiwin.h"
#include "connect/ncbi_core_c.h"
}

using namespace std;

/* Logging for debugging and profiling */
bool debug_msg = false;		/* Default, set to true with --debug */
std::ostream* log_stream; 		/**< Output stream for logging information */
int rank;			/**< Rank of the current MPI process */
int node_count;		/**< Number of MPI processes allocated for the job */
bool profile_msg = false;	/* Default, set to true with --pro-phile=filename */
ostream* profile_stream;	/**< Output stream for profile information */
double realstarttime ;		/* Given value immediately on first line of main, before MPI_Init */
double prog_start;		/* Given value after MPI_Init and MPE defs */
double prog_end ;		/* Given value right before MPI_Finalize */


void initNCBI( vector< string >& ncbi_opts ){
	//
	// initialize ncbi library
	//
	int ncbi_argc = ncbi_opts.size();
	//char** ncbi_argv = new char*[ ncbi_opts.size() ];
	char** ncbi_argv = NULL;
	
	if( debug_msg ) {
		LOG_MSG << "initializing ncbi ...";
	}
	
	if(ncbi_argc > 0){
		ncbi_argv = (char**)calloc(ncbi_argc, sizeof(char*));
	}
	
	for(uint optI = 0; optI < (uint)ncbi_argc; optI++ ){
		if( debug_msg )
			CONT_LOG_MSG << ncbi_opts[ optI ] << " ";
		//ncbi_argv[ optI ] = new char[ ncbi_opts[ optI ].size() + 1 ];
		ncbi_argv[ optI ] =(char*)calloc((ncbi_opts[ optI ].size() + 1), sizeof(char));
		strcpy( ncbi_argv[ optI ], ncbi_opts[ optI ].c_str() );
	}

	if( debug_msg )
		CONT_LOG_MSG << endl;
		
	Nlm_SetupArguments( ncbi_argc, ncbi_argv );
	
	// Note! Do not clean up the memory allocated in this function!
	// (i.e. ncbi_argv and ncbi_argv[i]). It is now owned by the
	// ncbi library. It is not clear if it needs to be deleted with
	// a call to Nlm_FreeCmdLineArguments(char** argv) in cleanupNCBI().
	// Since this would require saving the pointer, ncbi_argv, for
	// future reference, just forget about it (not a lot of memory
	// is lost anyway ...).
	
#ifdef MSC_VIRT
	if ( !_vheapinit(0, 1250, _VM_ALLSWAP) )
	{
		ErrPost( CTX_NCBIOBJ, 1, "Can't open virtual memory" );
		return 1;
	}
#endif

	/* Initialize connection library's logger, registry and lock */
	CONNECT_Init( 0 );
	
	if( debug_msg ) {
		LOG_MSG << "\n(" << rank << ") done initializing ncbi." << endl;
	}
}


void cleanupNCBI(){
	//
	// cleanup ncbi library
	//
	NlmThreadJoinAll();

	Nlm_FreeConfigStruct();
	ErrSetLogfile( NULL, 0 );
	Nlm_ReleaseAppContext();

#ifdef MSC_VIRT
	_vheapterm();
#endif
	NlmThreadDestroyAll();
}


void addOpt(vector< string >& opt_vector, int opt_character, const char* opt_argument)
{
	char opt_char = opt_character;
	string opt_string = "-";
	opt_string += opt_char;
	opt_vector.push_back( opt_string );
	if( opt_argument != NULL ){
		opt_vector.push_back( opt_argument );
	}
}

void SendIntVec(vector<int> &vec, int dest, int tag)
{
	int count = vec.size();
	
	MPI_Send(&count, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);

	if(count > 0){
		int *tmp = (int*)calloc(count, sizeof(int));
		
		if(tmp == NULL){
			throw __FILE__ "SendIntVec: Unable to allocate memory for buffer";
		}
		
		for(int i = 0;i < count;i++){
			tmp[i] = vec[i];
		}

		MPI_Send(tmp, count, MPI_INT, dest, 
			tag, MPI_COMM_WORLD);
			
		free(tmp);
	}

}

void RecvIntVec(vector<int> &vec, int src, int tag)
{
	int count;
	
	MPI_Recv(&count, 1, MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	vec = vector<int>(count);
	
	if(count > 0){
		int *tmp = (int*)calloc(count, sizeof(int));
		
		if(tmp == NULL){
			throw __FILE__ "RecvIntVec: Unable to allocate memory for buffer";
		}
		
		MPI_Recv(tmp, count, MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		for(int i = 0;i < count;i++){
			vec[i] = tmp[i];
		}
			
		free(tmp);
	}
}
