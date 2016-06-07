/******************************************************************************
** mpiblast_util.cpp
** Interface definitions of various mpiBLAST utility functions
** $Id: mpiblast_util.hpp,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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

#ifndef __MPI_BLAST_UTIL
#define __MPI_BLAST_UTIL

#include <vector>
#include <string>
#include <iostream>
#include "mpi.h"

// Logging for debugging and profiling: All variables are declared in mpiblast_util.cpp
extern bool debug_msg;		/* Default, set to true with --debug */
extern bool profile_msg;	/* Default, set to true with --pro-phile=filename */
extern std::ostream* log_stream; 		/**< Output stream for logging information */
extern std::ostream* profile_stream;	/**< Output stream for profile information */
extern double realstarttime ;		/* Given value immediately on first line of main, before MPI_Init */
extern double prog_start;		/* Given value after MPI_Init and MPE defs */
extern double prog_end ;		/* Given value right before MPI_Finalize */
extern int rank;			/**< Rank of the current MPI process */
extern int node_count;		/**< Number of MPI processes allocated for the job */

#ifdef USING_MPI
#define LOG_MSG (*log_stream) << "[" << rank << "]\t" << MPI_Wtime() - prog_start << '\t'
#define PROFILE_NEW_MSG (*profile_stream) << rank << '\t' << MPI_Wtime() - prog_start << '\t'
#else
#define LOG_MSG (*log_stream)
#define PROFILE_NEW_MSG (*profile_stream)
#endif // USING_MPI

#define CONT_LOG_MSG (*log_stream) 
#define PROFILE_MSG (*profile_stream) 

/**
 * Initializes the NCBI library with a particular vector of options
 */
void initNCBI( std::vector< std::string >& ncbi_opts );
void cleanupNCBI();

/**
 * Add a command line option to an option vector
 */
void addOpt( std::vector< std::string >& opt_vector, int opt_character, const char*
	opt_argument );

void SendIntVec( std::vector<int> &vec, int dest, int tag);

void RecvIntVec( std::vector<int> &vec, int src, int tag);

#endif // __MPI_BLAST_UTIL
