/******************************************************************************
** mpiblast_config.cpp
** Parses the mpiBLAST configuration file
** $Id: mpiblast_config.cpp,v 1.2 2005/01/05 07:01:13 hclee Exp $
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

#include "mpiblast_config.hpp"
#include "file_util.hpp"
#include <iostream>
using namespace std;

MpiBlastConfig::MpiBlastConfig( const string& filename )
{
	// Do we need to make a local copy of the config file?
  	string::size_type pos = filename.find(":", 0);
  
  	if(pos == string::npos){
  		// The file is local or accessible via nfs
		config_filename = filename;
  	}
  	else{
  		// Make a copy of the file in the current directory
		// Define the local name by (a) stripping off the hostname and
		// (b) appending a random six character string to the end.
		pos ++;
		config_filename = filename.substr(pos, filename.length() - pos) + "XXXXXX";
  
  		try{
			getTempFileName(config_filename);
		}
		catch(const char* error){
			cerr << "Error creating temp config file:" << endl;
			throw error;
		}
	
		copyFile(filename, config_filename);
  	}
  
  
  	// read config file
  	ifstream config_file( config_filename.c_str() );
	
  	if( !config_file.is_open() ){
    		cerr << "Error opening configuration file: " << config_filename << endl;
		cerr << "Try specifying the configuration file location with the --config-file option\n";
    		throw __FILE__ "(MpiBlastConfig): Unable to open config file";
  	}

//Add for globus_gass_copy by H.C.Lee 
#ifdef GASSCOPY
  	// read shared gass host
  	getline( config_file, shared_gass_host );
#endif
  	// read shared path name
  	getline( config_file, shared_db_path );
	
  	// read local path name
  	getline( config_file, local_db_path );
	
  	// read BLAST binaries path name
  	getline( config_file, blast_path );
	
  	config_file.close();
  
  	// Delete the temporary local file name if we created one
  	if(pos != string::npos){
  		unlink( config_filename.c_str() );
  	}

        // Create temp local path under local_db_path
	local_db_path = local_db_path + PATH_SEPARATOR + "mpiblast_g2_XXXXXX";
  	try{
		getTempDirName(local_db_path);
	}
	catch(const char* error){
		cerr << "Error creating temp local db path:" << endl;
		throw error;
	}
}

MpiBlastConfig::MpiBlastConfig( const MpiBlastConfig& mbc )
{
  *this = mbc;
}

MpiBlastConfig& MpiBlastConfig::operator=( const MpiBlastConfig& mbc )
{
  config_filename = mbc.config_filename;
  local_db_path = mbc.local_db_path;

//Add for globus_gass_copy by H.C.Lee 
#ifdef GASSCOPY
  shared_gass_host = mbc.shared_gass_host;
#endif
  shared_db_path = mbc.shared_db_path;
  blast_path = mbc.blast_path;
  return *this;
}


string MpiBlastConfig::defaultConfigFileName()
{
	static string file_name;

// check for definition of the default config file location
#ifdef DEFAULT_CONFIG_FILENAME
	file_name = DEFAULT_CONFIG_FILENAME;
#endif

#ifdef WIN32
	if( file_name == "" ){
		const char* home_path = getenv( "USERPROFILE" );
		if( home_path != NULL ){
			file_name = home_path;
			file_name += PATH_SEPARATOR;
			file_name += ".mpiblastrc";
			if( doesFileExist( file_name ) )
				return file_name;
		}
		// try the windows system directory 
		const char* winnt_path = getenv( "windir" );
		if( winnt_path != NULL ){
			file_name = winnt_path;
			file_name += PATH_SEPARATOR;
			file_name += "mpiblast.ini";
		}else
			file_name = "mpiblast.conf";
	}
#else
	if( file_name == "" ){
		const char* home_path = getenv( "HOME" );
		if( home_path != NULL ){
			file_name = home_path;
			file_name += PATH_SEPARATOR;
			file_name += ".mpiblastrc";
			if( doesFileExist( file_name ) )
				return file_name;
		}
		// use the mpiblast.conf in the installation path
		file_name = INSTALL_PREFIX;
#ifdef GASSCOPY
		file_name += PATH_SEPARATOR + "etc" + PATH_SEPARATOR + "mpiblast-gasscopy.conf";

#else
		file_name += PATH_SEPARATOR + "etc" + PATH_SEPARATOR + "mpiblast.conf";
#endif
	}
#endif
	return file_name;
}
