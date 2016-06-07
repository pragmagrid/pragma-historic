/******************************************************************************
** fragment_list.hpp
** Provides an API to the fragment list on worker nodes
** $Id: fragment_list.hpp,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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

#ifndef __MPI_BLAST_FRAGMENT_LIST
#define __MPI_BLAST_FRAGMENT_LIST

#include <string>
#include <vector>
#include <set>
#include "mpiblast_config.hpp"

/**
 * API for reading and modifying a database fragment list file.
 * The fragment list file tracks the database fragments currently on
 * local storage.
 */
class FragmentListFile {
public:
	FragmentListFile(){};
	/**
	 * Opens and reads a fragment list file, throws a (const char *) exception
	 * if something fails.
	 */
  	FragmentListFile( const std::string& filename );
 	FragmentListFile( const FragmentListFile& dbsf )
	{
		*this = dbsf;
	};
	
  	FragmentListFile& operator=( const FragmentListFile& mbc );
  	~FragmentListFile()
	{
		// Do nothing
	};
  
  	/** 
  	 * Adds the fragment ID given in <code>fragmentI</code>
  	 * to the fragment list file
  	 */
  	void addFragment( int fragmentI );

	#ifdef OLD_VERSION	
  	/**
  	 * calls <code>new</code> to allocate an array of type <code>uint</code>
  	 * and fills it with the current list of fragments on local storage
  	 */
	int* allocateFragmentList();
	
	#endif // OLD_VERSION
	
	void SendList( const MpiBlastConfig& config, const std::string& database_name, const std::string& db_type );
	
	/**
	 * Returns the number of fragments currently on local storage
	 */
	inline int fragmentCount()
	{
		return (int)(fragment_list.size());
	};
private:
  	std::vector< int > fragment_list;
 	std::set< int > fragment_set;
  	std::string frag_filename;
};


/**
 * Extensions for sequence database files
 */
std::vector< std::string > FragmentExtensions(std::string db_type);
std::vector< std::string > FragmentExtensions_optional(std::string db_type);
std::vector< std::string > makeFragmentExtensions(std::string db_type);
std::vector< std::string > makeOptionalFragmentExtensions(std::string db_type);

#endif
