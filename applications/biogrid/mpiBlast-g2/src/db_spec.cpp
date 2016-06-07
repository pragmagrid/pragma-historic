/******************************************************************************
** db_spec.cpp
** Implements an interface to the database specification file
** $Id: db_spec.cpp,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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

#include "db_spec.hpp"
#include "mpiblast.hpp"
using namespace std;

DbSpecFile::DbSpecFile( const string& filename )
{
  ifstream spec_file( filename.c_str() );

  if( !spec_file.is_open() ){
    cerr << "(" << rank << ") Error opening: " << filename << endl;
    throw __FILE__ "(DbSpecFile): Unable to open file";
  }
  
  string db_name;
  string nc;

  getline( spec_file, db_name );
  string db_size_line;
  getline( spec_file, db_size_line );
  stringstream db_sss( db_size_line );
  db_sss >> db_size;
  spec_file >> fragment_count;
  spec_file.close();
}

DbSpecFile::DbSpecFile( const DbSpecFile& dbsf )
{
  *this = dbsf;
}

DbSpecFile& DbSpecFile::operator=( const DbSpecFile& mbc )
{
  
  db_size = mbc.db_size;
  fragment_count = mbc.fragment_count;
  return *this;
}

void DbSpecFile::write( ostream& os, const string& db_name, uint64 db_size, uint fragment_count ) {
	os << db_name << endl;
	os << db_size << endl;
	os << fragment_count << endl;
}
