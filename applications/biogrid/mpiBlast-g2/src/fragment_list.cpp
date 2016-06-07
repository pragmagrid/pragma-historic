/******************************************************************************
** fragment_list.cpp
** Parses and provides an API to the fragment list on worker nodes
** $Id: fragment_list.cpp,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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

#include "fragment_list.hpp"
#include "mpiblast.hpp"
#include "mpiblast_util.hpp"
using namespace std;

FragmentListFile::FragmentListFile( const string& filename ){
  frag_filename = filename;
  ifstream frag_file( filename.c_str() );

  if( !frag_file.is_open() ){
	ofstream new_frag_file( filename.c_str() );
	new_frag_file.close();
	return;
  }

  int fragmentI;
  while( frag_file >> fragmentI ){
  	fragment_list.push_back( fragmentI );
  	fragment_set.insert( fragmentI );
  }
  
  frag_file.close();
}

FragmentListFile& FragmentListFile::operator=( const FragmentListFile& mbc ){
  fragment_list = mbc.fragment_list;
  fragment_set = mbc.fragment_set;
  frag_filename = mbc.frag_filename;
  return *this;
}


void FragmentListFile::addFragment( int fragmentI ){
	if( debug_msg ){
		LOG_MSG << "trying to add fragment: " << fragmentI << endl;
	}
	
	set< int >::iterator frag_iter = fragment_set.find( fragmentI );
	
	if( frag_iter != fragment_set.end() ){
		return; 	// the fragment is already in the list
	}
	
	// FIXME: silly lockfile synchronization method.
	// fix this somehow (it has a race condition in unix)
	
	string lock_filename = frag_filename + ".lock";
	while( true ){
		ifstream lock_file( lock_filename.c_str(), ios::in );
		if( lock_file.is_open() ){
			lock_file.close();
			if( debug_msg ){
				LOG_MSG << "fragment list is locked.  waiting for 1 s.\n";
			}
#ifdef WIN32
			Sleep( 1000 );	// milliseconds
#else
			sleep( 1 );	// seconds
#endif
			continue;
		}
		
		if( debug_msg ){
			LOG_MSG << "locking fragment list for fragment: " << fragmentI << endl;
		}
		
		ofstream new_lock_file( lock_filename.c_str() );
		break;
	}
	
	ofstream frag_file( frag_filename.c_str(), ios::app );

	if( !frag_file.is_open() ){
		if( debug_msg ){
			cerr << "Error opening: " << frag_filename << endl;
		}
		
		removeFile( lock_filename );
		
		throw __FILE__ "(addFragment): Unable to open file" ;
	}
	
	frag_file << fragmentI << endl;
	fragment_list.push_back( fragmentI );
	fragment_set.insert( fragmentI );
	frag_file.close();
	
	removeFile( lock_filename );
}

#ifdef OLD_VERSION
int* FragmentListFile::allocateFragmentList()
{
	//int* flist = new int[ fragment_list.size() ];
	
	if(fragment_list.size() <= 0){
		return NULL;
	}
	
	int* flist = (int*)malloc(fragment_list.size()*sizeof(int));
	
	for( uint fragmentI = 0; fragmentI < fragment_list.size(); fragmentI++ ){
		flist[ fragmentI ] = fragment_list[ fragmentI ];
	}
	
	return flist;
}
#endif // OLD_VERSION

void FragmentListFile::SendList( const MpiBlastConfig& config, const string& database_name, const string& db_type )
{
	/** begin DB fragment update code */

	vector< int > up_to_date_fragments;
	char frag_num_buff[4];
	for(int fragI = 0 ; fragI < (int)fragment_list.size(); fragI++){
#ifdef MANY_FRAGMENTS_HACK
		sprintf(frag_num_buff,"%03d",fragment_list[ fragI ]);
#else
		sprintf(frag_num_buff,"%02d",fragment_list[ fragI ]);
#endif
//Add for globus_gass_copy by H.C.Lee
#ifdef GASSCOPY
//Ignore time checking, always update database fragments for gass copy
//Nothing to do
#else
		string local_fragname = config.localPath() + PATH_SEPARATOR + database_name + "." + frag_num_buff + "." + db_type + "hr";
		string remote_fragname = config.sharedPath() + PATH_SEPARATOR + database_name + "." + frag_num_buff + "." + db_type + "hr";
		//if the remote frag has a later modification time than the local one, we remove all refs to  the local one
		//statFileMTime returns 0 if the file couldn't be stated (doesn't exist)
		bool out_of_date = false;
		time_t local_filetime = statFileMTime( local_fragname.c_str() );
		time_t remote_filetime = statFileMTime( remote_fragname.c_str() );
		if ( (local_filetime <= remote_filetime) || (local_filetime == 0) )
			out_of_date = true;

		if( !out_of_date ){
			up_to_date_fragments.push_back( fragment_list[ fragI ] );
		}else if( debug_msg ){
			LOG_MSG << local_fragname.c_str() << " modified on " << local_filetime << endl ;
			LOG_MSG << remote_fragname.c_str() << " modified on " << remote_filetime << endl ;
		}
#endif	
	}
	
	/** end DB fragment update code */

	SendIntVec(fragment_list, 0, DB_FRAGMENTS_TAG);
}
	
vector<string> FragmentExtensions(string db_type)
{
  static vector<string> n_extensions = makeFragmentExtensions( "n" );
  static vector<string> p_extensions = makeFragmentExtensions( "p" );
  if( db_type == "n" ){
  	return n_extensions;
  }
  return p_extensions;
}
vector<string> FragmentExtensions_optional(string db_type)
{
	static vector<string> n_extensions = makeOptionalFragmentExtensions( "n" );
	static vector<string> p_extensions = makeOptionalFragmentExtensions( "p" );
	if( db_type == "n" ){
  		return n_extensions;
	}
	return p_extensions;
}
vector<string> makeFragmentExtensions(string db_type)
{
	vector< string > extensions;
	
	extensions.push_back( "." + db_type + "hr" );
	extensions.push_back( "." + db_type + "in" );
	extensions.push_back( "." + db_type + "sq" );
	
	return extensions;
}
vector<string> makeOptionalFragmentExtensions(string db_type)
{
	vector<string> extensions;
	
	extensions.push_back( "."+ db_type + "nd" );
	extensions.push_back( "."+ db_type + "ni" );
	extensions.push_back( "."+ db_type + "sd" );
	extensions.push_back( "."+ db_type + "si" );
	
	return extensions;
}
