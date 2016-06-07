/******************************************************************************
** file_util.cpp
** Implements a cross-platform interface to the filesystem.
** $Id: file_util.cpp,v 1.4 2005/01/05 07:01:13 hclee Exp $
** $Revision: 1.4 $
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

#include "file_util.hpp"
#include "mpiblast.hpp"
#include "sys/types.h"
#include "sys/stat.h"

#ifdef GASSCOPY
extern "C"{
#include "gasscopy.h"
}
#endif

using namespace std;

bool doesFileExist( const string& file_name )
{
	// this is a platform agnostic method of checking for existence,
	// but opening the file is somewhat undesirable.
	ifstream check_file( file_name.c_str(), ios::in );
	if( check_file.is_open() ){
		check_file.close();
		return true;
	}
	return false;
}

/** Utility function to get the path component of a file path */
string getPath( const string& file_path )
{
  string::size_type slashloc;
  string::size_type search_loc = string::npos;
  for(;;){
    slashloc = file_path.rfind( PATH_SEPARATOR, search_loc );
    if( slashloc > 0 && slashloc != string::npos )
      if( file_path[ slashloc - 1 ] == '\\' ){
	search_loc = slashloc - 1;
	continue;
      }
    break;
  }
  string path_str;
  if( slashloc != string::npos )
    path_str = file_path.substr( 0, slashloc + 1 );
  else
    path_str = "";

  return path_str;
}

/** Utility function to get the file name component of a file path */
string getFilename( const string& file_path )
{
  string::size_type slashloc;
  string::size_type search_loc = string::npos;
  for(;;){
    slashloc = file_path.rfind( PATH_SEPARATOR, search_loc );
    if( slashloc > 0 && slashloc != string::npos )
      if( file_path[ slashloc - 1 ] == '\\' ){
	search_loc = slashloc - 1;
	continue;
      }
    break;
  }
  string file_str;
  if( slashloc != string::npos )
    file_str = file_path.substr( slashloc + 1 );
  else
    file_str = file_path;

  return file_str;
}

/** Utility function to get a name for a temp directory */
void getTempDirName( string& tempname )
{

	//char* temp_char_name = new char[ tempname.size() + 1 ];
	char* temp_char_name = (char*)malloc((tempname.size() + 1)*sizeof(char));
	
	if(temp_char_name == NULL){
		throw __FILE__ "(getTempFileName): Unable to allocate memory for temp file" ;
	}
	
	strcpy( temp_char_name, tempname.c_str() );
	
	if( debug_msg ){
		LOG_MSG << "Temp name base: " << tempname << endl;
	}
	
#ifdef WIN32
	char* win_name = mktemp( temp_char_name );
	tempname = win_name;
	// mkstemp actually creates the dir, windows mktemp() doesn't.
        // Need to clarify the method for creating temp dir
	//ofstream tmp_dir( win_name );
#else
	char* fd = mkdtemp( temp_char_name );
	
	if(fd == NULL){
		cerr << "Unable to make temp dir name " << temp_char_name 
		     << " from " << tempname << endl;
		cerr << "fd = " << fd << endl;		
		cerr << "Error Msg = \"" << strerror(errno) << "\"" << endl;
		
		throw __FILE__ "(getTempFileName): Unable to make temp dir" ;
	}
	
	// Close the file descriptor (we only needed to create an empty file)
	//close(fd);
	
	tempname = temp_char_name;
#endif
		
	if( debug_msg ){

		LOG_MSG << "Got temp name: " << tempname << endl;
	}
	
	//delete[] temp_char_name;
	free(temp_char_name);
}

/** Utility function to get a name for a temp file */
void getTempFileName( string& tempname )
{

	//char* temp_char_name = new char[ tempname.size() + 1 ];
	char* temp_char_name = (char*)malloc((tempname.size() + 1)*sizeof(char));
	
	if(temp_char_name == NULL){
		throw __FILE__ "(getTempFileName): Unable to allocate memory for temp file" ;
	}
	
	strcpy( temp_char_name, tempname.c_str() );
	
	if( debug_msg ){
		LOG_MSG << "Temp name base: " << tempname << endl;
	}
	
#ifdef WIN32
	char* win_name = mktemp( temp_char_name );
	tempname = win_name;
	// mkstemp actually creates the file, windows mktemp() doesn't.
	ofstream tmp_file( win_name );
#else
	int fd = mkstemp( temp_char_name );
	
	if(fd < 0){
		cerr << "Unable to make temp file name " << temp_char_name 
		     << " from " << tempname << endl;
		cerr << "fd = " << fd << endl;		
		cerr << "Error Msg = \"" << strerror(errno) << "\"" << endl;
		
		throw __FILE__ "(getTempFileName): Unable to make temp file" ;
	}
	
	// Close the file descriptor (we only needed to create an empty file)
	close(fd);
	
	tempname = temp_char_name;
#endif
		
	if( debug_msg ){

		LOG_MSG << "Got temp name: " << tempname << endl;
	}
	
	//delete[] temp_char_name;
	free(temp_char_name);
}

/** Utility function to move a file */
int moveFile( const string& src, const string& dest )
{
#ifdef WIN32
	return !MoveFileEx( src.c_str(), dest.c_str(), MOVEFILE_COPY_ALLOWED | MOVEFILE_REPLACE_EXISTING );
//	string mv_cmd = "move /Y " + src + " " + dest;
//	if( debug_msg )
//		LOG_MSG << mv_cmd << endl;
#else
	// The original version
	// string mv_cmd = "mv " + src + " " + dest;
	// return system( mv_cmd.c_str() );
	
	// Updated version that will use rcp if requested.
	// Note that we are expecting only the destination string to 
	// possibly have a hostname prepended to it.
	string::size_type dest_pos = dest.find(":", 0);
	
	// Do we need to use rcp?
	if( dest_pos == string::npos ){
		// Don't use rcp, use mv
		// This is the original behavior
		
		// First, check to make sure that src != dest
		if(src == dest){
			return EXIT_SUCCESS;
		}
		
		string mv_cmd;

		mv_cmd = "mv " + src + " " + dest;
		
		return system( mv_cmd.c_str() );
	}
	else{
		// Don't use cp, use rcp
		
		// First, check to make sure that src != dest
		// Note that we trust the HOSTNAME environment variable
		// Build the fully qualified src path, i.e.:
		//	hostname:/local/path/to/src
		string full_src_name = getenv( "HOSTNAME" );
		
		// Do we need to append the full path?
		if(src.find(PATH_SEPARATOR, 0) == string::npos){
			// Yes, append the path to src
			char dir_buf[1024];
			
			string curr_dir = getcwd(dir_buf, 1024);
			
			full_src_name = full_src_name + ":" + curr_dir + PATH_SEPARATOR + src;
			
		}
		else{
			// src already contains a path. Lets hope its an absolute path!
			full_src_name = full_src_name + ":" + src;
		}
		
		if(full_src_name == dest){
			return EXIT_SUCCESS;
		}
		
		string rcp_cmd;
		int ret_value;
		
		rcp_cmd = "rcp " + src + " " + dest;
		
		ret_value = system( rcp_cmd.c_str() );
		
		if(ret_value != EXIT_SUCCESS){
			// DEBUG
			cerr << "rcp command failed!" << endl;
			cerr << "source = " << src << endl;
			cerr << "dest = " << dest << endl;
			
			return ret_value;
		}
		
		// Unlink the src file to complete the move
		ret_value = unlink( src.c_str() );
		
		if(ret_value != EXIT_SUCCESS){
			// DEBUG
			cerr << "unlink command failed!" << endl;
			cerr << "source = " << src << endl;
			cerr << "dest = " << dest << endl;
			
			return ret_value;
		}
		
		return ret_value;
	}
	
#endif
}

#ifdef GASSCOPY
/** Utility function to copy a file with gasscopy */
int gassCopyFile( const string& src, const string& dest )
{
	// The original version
	// string cp_cmd = "cp " + src + " " + dest;
	// return system( cp_cmd.c_str() );
	
	// Updated version that will use rcp if requested
	// Note that we are expecting only the source string to 
	// possibly have a hostname.
	
	//string::size_type src_pos = src.find(":", 0);
	
	int ret_value;

        // Add by H.C.Lee for -
        // Using gasscopy instead of normal command line cp/rcp
        // User credential (i.e. grid-proxy) is required
        char* gass_src  = const_cast<char*>(src.c_str());
        char* gass_dest = const_cast<char*>(dest.c_str());
	ret_value = gasscopy(gass_src, gass_dest);
	if(ret_value != EXIT_SUCCESS){
		cerr << "gasscopy failed!" << endl;
		cerr << "source = " << src << endl;
		cerr << "dest = " << dest << endl;
		cerr << "ret_value = " << ret_value << endl;
	}
	return ret_value;
}
#endif

/** Utility function to copy a file */
int copyFile( const string& src, const string& dest )
{
#ifdef WIN32
	return !CopyFile( src.c_str(), dest.c_str(), false );
#else
	// The original version
	// string cp_cmd = "cp " + src + " " + dest;
	// return system( cp_cmd.c_str() );
	
	// Updated version that will use rcp if requested
	// Note that we are expecting only the source string to 
	// possibly have a hostname.
	
	string::size_type src_pos = src.find(":", 0);
	
	int ret_value;

	// Do we need to use rcp?
	if( src_pos == string::npos ){
		// This is the original behavior
		string cp_cmd = "cp " + src + " " + dest;
		
		ret_value = system( cp_cmd.c_str() );
		
		if(ret_value != EXIT_SUCCESS){
			cerr << "cp command failed!" << endl;
			cerr << "command: " << cp_cmd << endl;
			cerr << "source = " << src << endl;
			cerr << "dest = " << dest << endl;
			cerr << "ret_value = " << ret_value << endl;
		}
		
		return ret_value;
	}
	else{
		// Don't use cp, use rcp
		
		string rcp_cmd;
		
		rcp_cmd = "rcp " + src + " " + dest;
		
		ret_value = system( rcp_cmd.c_str() );
		
		if(ret_value != EXIT_SUCCESS){
			// DEBUG
			cerr << "rcp command failed!" << endl;
			cerr << "source = " << src << endl;
			cerr << "dest = " << dest << endl;
			
			return ret_value;
		}
		
		return ret_value;
	}
	// Keep the compiler happy -- should never get here
	return EXIT_FAILURE;
#endif
}

/** Utility function to untar a file in specific directory */
int untarFile( const string& filename, const string& dir, bool verbose)
{
	string tar_cmd;
#ifdef WIN32
	//return !DeleteFile( filename.c_str() );
#else
	if( verbose )
		tar_cmd = "/bin/tar -C " + dir +" -xzvf " + filename;
	else
		tar_cmd = "/bin/tar -C " + dir +" -xzf " + filename;
	return system( tar_cmd.c_str() );
#endif
}

/** Utility function to delete a file */
int removeFile( const string& filename, bool verbose )
{
	string rm_cmd;
#ifdef WIN32
	return !DeleteFile( filename.c_str() );
#else
	if( verbose )
		rm_cmd = "/bin/rm -fv " + filename;
	else
		rm_cmd = "/bin/rm -f " + filename;
	return system( rm_cmd.c_str() );
#endif
}

uint64 statFileSize( const char * path ) {
#ifdef WIN32
	WIN32_FILE_ATTRIBUTE_DATA file_data;
	uint64 f_size;
	GetFileAttributesEx( path, GetFileExInfoStandard, (void*)&file_data );
	f_size = file_data.nFileSizeHigh;
	f_size <<= 32;
	f_size += file_data.nFileSizeLow;
	return f_size;
#else
	struct stat stat_data;
	if( stat( path , &stat_data) ){
		perror(path);
	}
	return stat_data.st_size;
#endif
}

time_t statFileMTime( const char * path ) {
#ifdef WIN32
	struct __stat64 stat_data;
	if( _stat64( path, &stat_data ) ){
		perror(path);
	}
	return stat_data.st_mtime;
#else
	struct stat stat_data;
	if( stat( path , &stat_data) ){
		perror(path);
	}
	return stat_data.st_mtime;
#endif
}
