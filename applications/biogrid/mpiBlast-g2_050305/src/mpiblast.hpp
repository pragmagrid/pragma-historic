/******************************************************************************
** mpiblast.hpp
** Interface and variable definitions for mpiBLAST
** $Id: mpiblast.hpp,v 1.6 2004/10/20 12:52:59 hclee Exp $
** $Revision: 1.6 $
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

#ifndef __mpiblast_hpp__
#define __mpiblast_hpp__

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <vector>
#ifdef USING_MPI
#include "mpi.h"
#endif
#ifdef __GNUG__
#include "unistd.h"
#endif
#include "stdio.h"

#include "file_util.hpp"
#include "mpiblast_util.hpp"
#include "mpiblast_config.hpp"
#include "mpiblast_types.h"
#include "fragment_list.hpp"
#include "db_spec.hpp"
#include "blastjob.hpp"

#include "ncbi.h"
#include "objseq.h"
#include "ncbithr.h"
#include "ncbiwin.h"
#include "connect/ncbi_core_c.h"

// Tags for each type of MPI communication
#define	DB_NAME_TAG		1	/* The database name */
#define	DB_FRAGMENTS_TAG	2	/* A list of database fragments */
#define	HOSTNAME_TAG		3	/* A hostname (obsolete) */
#define	DB_FRAGMENT_COUNT_TAG	4	/* The number of fragments */
#define	QUERY_NAME_LENGTH_TAG	5	/* The length of a query file name */
#define	QUERY_NAME_TAG		6	/* A query file name */
#define	BLAST_RETCODE_TAG	7	/* A BLAST return code */
#define	DB_NAME_LEN_TAG		8	/* The length of a database name */
#define	ASSIGNMENT_TYPE		9	/* The type of worker assignment given */
#define	WORKER_IDLE		10	/* That a worker is idle */
#define	FRAGMENT_COPY_COMPLETE	11	/* That a fragment copy has completed */
#define	BLAST_RESULTS_TAG	12	/* Actual BLAST results */
#define	SUCCESSFUL_IDS_TAG	13	/* A list of queries that had hits to the DB */
#define	SUCCESSFUL_ID_COUNT	14	/* The number of queries that hit the DB */


// assignment types
#define SEARCH_COMPLETE		0 /* Code sent to a worker when the search is complete */ 
#define COPY_FRAGMENT		1 /* Code indicating that the worker should copy a fragment */
#define SEARCH_FRAGMENT		2 /* Code telling a worker to search a fragment */
#define WORKER_QUIT		3 /* Code telling a worker to quit */
#define DO_NOTHING		4 /* place holder telling a worker to do nothing */



#ifdef USE_NCBI_ASCII
static bool use_binary_asn = false;
#else
static bool use_binary_asn = true;
#endif

/**
 * The MpiBlast class provides the functionality to perform BLAST
 * searches on a cluster of nodes connected with MPI.
 */
class MpiBlast{
	public:
		MpiBlast()
		{ 
			log_stream = &std::cerr;  
			profile_stream = &std::cerr; 
			remove_db = false;
			use_mpi_db = true;
                        compress_io = false;
		}
		
		~MpiBlast();

                std::string getOutputFilename();
                std::string getProfileFilename();
                std::string getResultURL();
                std::string getProfileURL();

		/**
		 * Send results from a workers ASN.1 format file
		 */
		void sendResults( std::string& b_filename, int src, int tag);
		
		/**
		 * Receives results in ASN.1 format from a worker
		 */
		void receiveResults( int src, int tag, std::vector< SeqAlignPtr >& results, 
				     int m_result_count );

		/**
		 * Reads the file named in b_filename and broadcasts it to all other nodes
		 */
		void broadcastFile( std::string& b_filename );
		/**
		 * Receives a broadcasted file and writes it to the file named in b_filename
		 */
		void recvBroadcastFile( std::string& b_filename );
		
		/**
		 * Removes all local files associated with the database
		 * reason allows us force cleanup when called from terminateProgram
		 */
		void cleanupDB( int reason );


		/**
		 * performs the master functions
		 */
		void master();
		/**
		 * Performs the worker functions
		 */
		void worker();
		
		/**
		 * fill alias file with the names of database fragments
		 */		
		void WriteAliasFile( std::string &alias_filename, std::vector< int > &fragment_list);
				    
		/**
		 * Call this with the program arguments to start mpiBLAST
		 */
		int main( int argc, char* argv[] );
		
		/** True if the user requested the database be removed */
		bool remove_db;	
		
		/** True if the user has enabled MPI distributed database queries */
		bool use_mpi_db;

		/** True if the user has compressed database fragments */
		bool compress_io;
	protected:
		BlastJob blast_job;		/**< Class to manage scheduling of work in each BLAST search job */
		MpiBlastConfig config;	/**< Class to read in configuration file */
		FragmentListFile frag_list_file;	/**< Class to read and modify local fragment list file */
		std::string exec_path;		/**< Directory that mpiblast was executed from */
		std::ofstream log_file, pro_phile;	/**< output file streams for logging and profiling information */
		std::string query_filename;	/**< Full path to query file */
		std::string output_filename;	/**< Full path to output file */
		std::string database_name;	/**< Name of database to query */
	        std::string pro_phile_name;     /**< Profile filename **/

#ifdef GASSCOPY
	        std::string result_url;
	        std::string profile_url;
#endif
		std::string blast_cl;		/**< Base BLAST command line including user specified options */
		std::string blast_type;		/**< Type of blast search ( blastn, blastp, etc. ) */
		std::string db_type;			/**< n for nucleotide, p for protein */
		std::string query_file;		/**< Just the filename of the query file (pathname stripped) */
		std::string local_query_filename;	/**< path to a local copy of the query for workers */
		std::string db_spec_filename;	/**< name of the database specification file */
		DbSpecFile db_spec;		/**< Class to read the database specification file */

		std::vector< std::string > master_opts;	/**< Command line options for the master node to generate output */
		std::vector< std::string > worker_opts;	/**< Command line opts for each worker node to search correctly */
		
		std::vector< std::string > cleanup_files;	/**< Files to remove upon termination */

		/**
		 * cleanupDB for unices
		 */
		void unixCleanupDB(int reason);
};


#endif  // __mpiblast_hpp__
