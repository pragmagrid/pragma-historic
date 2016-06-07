/******************************************************************************
** blastjob.hpp
** Provides the programming interface for the BlastJob class
** $Id: blastjob.hpp,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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

#ifndef __BlastJob_hpp__
#define __BlastJob_hpp__

#include <string>
#include <set>
#include <vector>

/**
 * Manages each Blast job by allocating work to worker nodes
 */
class BlastJob{
	public:
		BlastJob(){};
		BlastJob( int fragment_count, std::vector< std::vector< int > >& node_fragments );
		BlastJob( const BlastJob& bj );
		BlastJob& operator=( const BlastJob& bj );
		
		/**
		 * Calculates the best assignment for a node based on the current
		 * distribution of database fragments
		 */
		void getAssignment( int nodeI, int& operation, int& fragment_id );
		/**
		 * Call this when a fragment copy completes
		 */
		void CopyCompleted( int nodeI, int fragment_id );
	protected:
		int fragment_count;	/**< The total number of fragments for this database */
		std::vector< std::vector< int > > node_fragments;	/**< The set of fragments residing on each node -- index[node][fragment] */
		std::vector< std::set< int > > fragment_nodes;	/**< The set of nodes each fragment resides on -- index[fragment][node] */
		std::set< int > unassigned_fragments;	/**< Fragments that haven't been assigned yet */
		std::vector< std::set< int > > copy_fragments; /**< Fragments that are currently being copied somewhere */
};

#endif	//__BlastJob_h__
