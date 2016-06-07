/******************************************************************************
** mpiblast_config.hpp
** Provides an interface to the mpiBLAST configuration file
** $Id: mpiblast_config.hpp,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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

#ifndef __MPI_BLAST_CONFIG
#define __MPI_BLAST_CONFIG

#include <string>

/**
 * API for reading the mpiBLAST configuration file.
 */
class MpiBlastConfig {
public:
  MpiBlastConfig(){};
  /**
   * Opens and reads a config file, throws a (const char*) exception if
   * something fails.
   */
  MpiBlastConfig( const std::string& filename );
  MpiBlastConfig( const MpiBlastConfig& mbc );
  MpiBlastConfig& operator=( const MpiBlastConfig& mbc );

// Add for globus_gass_copy by H.C.Lee
#ifdef GASSCOPY
  /** Returns the path to shared storage */
  const std::string& sharedGassHost() const { return shared_gass_host; }
#endif

  /** Returns the path to shared storage */
  const std::string& sharedPath() const { return shared_db_path; }

  /** Returns the path to local storage */
  const std::string& localPath() const { return local_db_path; }

  /** Returns the path to the NCBI BLAST binaries */
  const std::string& blastPath() const { return blast_path; }

	/**
	 * Returns the default path to the configuration file
	 */
	static std::string defaultConfigFileName();

private:
  std::string config_filename;	/**< The path to the config file */
  std::string local_db_path;		/**< The path to the local database storage */

// Add for globus_gass_copy by H.C.Lee
#ifdef GASSCOPY
  std::string shared_gass_host;	/**< The gass host to shared database storage*/
#endif
  std::string shared_db_path;	/**< The path to the shared database storage */
  std::string blast_path;		/**< The path to the NCBI BLAST binaries (obsolete) */
};

#endif
