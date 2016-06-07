/******************************************************************************
** blastjob.cpp
** Implements the BlastJob class - manages resource allocation for BLAST jobs
** $Id: blastjob.cpp,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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

#include "blastjob.hpp"
#include "mpiblast.hpp"
using namespace std;

BlastJob::BlastJob( int fragment_count, vector< vector< int > >& node_fragments ){
	this->fragment_count = fragment_count;
	this->node_fragments = node_fragments;
	uint fragmentI;
	for( fragmentI = 0; fragmentI < (uint)fragment_count; fragmentI++ ){
		set< int > frag_set;
		fragment_nodes.push_back( frag_set );
		copy_fragments.push_back( frag_set );
		unassigned_fragments.insert( fragmentI );
	}

	for( uint nodeI = 0; nodeI < node_fragments.size(); nodeI++ ){
		if( node_fragments[ nodeI ].size() == 0 ){
			continue;
		}
		
		if(debug_msg){
			LOG_MSG << "node " << nodeI << " has fragments ";
		}
		
		for( fragmentI = 0; fragmentI < node_fragments[ nodeI ].size(); fragmentI++ ){
			if( debug_msg ){
				(*log_stream) << node_fragments[ nodeI ][ fragmentI ] << " ";
			}
			
			fragment_nodes[ node_fragments[ nodeI ][ fragmentI ] ].insert( nodeI );
		}
		
		if(debug_msg){
			(*log_stream) << endl;
		}
	}
}

BlastJob::BlastJob( const BlastJob& bj ){
	(*this)=bj;
}

BlastJob& BlastJob::operator=( const BlastJob& bj ){
	fragment_count = bj.fragment_count;
	node_fragments = bj.node_fragments;
	fragment_nodes = bj.fragment_nodes;
	unassigned_fragments = bj.unassigned_fragments;
	copy_fragments = bj.copy_fragments;
	return *this;
}

void BlastJob::getAssignment( int nodeI, int& operation, int& fragment_id ){
	// find the lowest replication fragment for this node

	int best_fragment;
	int lowest_replication = 0;
	int fragmentI = 0;
	set< int >::iterator frag_iter = unassigned_fragments.begin();
	
	for(; frag_iter != unassigned_fragments.end(); frag_iter++ ){
		fragmentI = *frag_iter;
		if( fragment_nodes[ fragmentI ].find( nodeI ) != fragment_nodes[ fragmentI ].end() ){
			if( (fragment_nodes[ fragmentI ].size() < (uint)lowest_replication) || 
			   (lowest_replication == 0) ){
				best_fragment = fragmentI;
				lowest_replication = fragment_nodes[ fragmentI ].size();
			}
		}
	}

	if( lowest_replication != 0 ){
		unassigned_fragments.erase( best_fragment );
		operation = SEARCH_FRAGMENT;
		fragment_id = best_fragment;
		return;
	}

	// if lowest_replication is still 0 then this node either processed all of its fragments or never
	// had any fragments to process.  have it copy something.

	frag_iter = unassigned_fragments.begin();
	for(; frag_iter != unassigned_fragments.end(); frag_iter++ ){
		fragmentI = *frag_iter;
		if( fragment_nodes[ fragmentI ].size() == 0 &&
			copy_fragments[ fragmentI ].size() == 0 ){
			// give it this one.
			fragment_id = fragmentI;
			operation = COPY_FRAGMENT;
			copy_fragments[ fragmentI ].insert( nodeI );
			return;
		}
		// only count this fragment if it isn't already being copied
		if( copy_fragments[ fragmentI ].size() == 0 &&  
			( fragment_nodes[ fragmentI ].size() < (uint)lowest_replication || 
			lowest_replication == 0 ) ){
			best_fragment = fragmentI;
			lowest_replication = fragment_nodes[ fragmentI ].size();
		}
	}
	if( lowest_replication != 0 ){
		operation = COPY_FRAGMENT;
		fragment_id = best_fragment;
		copy_fragments[ fragmentI ].insert( nodeI );
		return;
	}
	operation = SEARCH_COMPLETE;
	fragment_id = 0;

}

void BlastJob::CopyCompleted( int nodeI, int fragment_id ) {
	fragment_nodes[ fragment_id ].insert( nodeI );
	copy_fragments[ fragment_id ].erase( nodeI );
	node_fragments[ nodeI ].push_back( fragment_id );
	
	if(debug_msg){
		LOG_MSG << "node " << nodeI << " copied fragment " << fragment_id << endl;
	}
}
