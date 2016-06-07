/******************************************************************************
** file_util.hpp
** Provides a cross-platform interface to the filesystem.
** $Id: file_util.hpp,v 1.3 2005/01/05 07:01:13 hclee Exp $
** $Revision: 1.3 $
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

#ifndef __MPIBLAST_FILE_UTIL
#define __MPIBLAST_FILE_UTIL

#include <string>
#include <fstream>
#include "mpiblast_types.h"

#ifdef WIN32
const std::string PATH_SEPARATOR = "\\";
#else
#include "unistd.h"
const std::string PATH_SEPARATOR = "/";
#endif

/** Extension for a database's local fragment list file */
#define	FRAG_LIST_EXTENSION 	".mbf"

bool doesFileExist(const std::string& file_name );
std::string getPath(const std::string& file_path );
std::string getFilename(const std::string& file_path);
void getTempFileName(std::string& tempname);
void getTempDirName(std::string& tempname);
int moveFile(const std::string& src, const std::string& dest);
/** copy file by globus_gass_copy - add by H.C.Lee */
int gassCopyFile(const std::string& src, const std::string& dest);
int copyFile(const std::string& src, const std::string& dest);
int removeFile(const std::string& filename, bool verbose = false);
int untarFile( const std::string& filename, const std::string& dir, bool verbose = false);

/** get the size in bytes of a particular file */
uint64 statFileSize( const char* path );

/** get the modification time of a file */
time_t statFileMTime( const char* path );

#endif // __MPIBLAST_FILE_UTIL
