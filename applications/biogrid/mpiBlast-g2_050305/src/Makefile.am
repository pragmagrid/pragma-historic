###############################################################################
# This file is a part of mpiBLAST and is copyright (c) 2002, 2003
# by Aaron Darling of Los Alamos National Laboratory for the Regents
# of the University of California.
# 
# License:
# 
# Copyright 2002, 2003. The Regents of the University of California.
# 
# This material was produced under U.S. Government contract
# W-7405-ENG-36 for Los Alamos National Laboratory, which is operated
# by the University of California for the U.S. Department of Energy.
# The Government is granted for itself and other acting on its behalf
# a paid-up, non-exclusive, irrevocable worldwide license in the
# material to reproduce, prepare derivative works, and perform
# publicly and display publicly. Beginning five (5) years after
# November 6, 2002, subject to additional five-year worldwide
# renewals, the Government is granted for itself and others acting on
# its behalf a paid-up, nonexclusive, irrevocable worldwide license
# in this material to reproduce, prepare derivative works, distribute
# copies to the public, perform publicly and display publicly, and to
# permit others to do so. NEITHER THE UNITED STATES NOR THE UNITED
# STATES DEPARTMENT OF ENERGY, NOR THE UNIVERSITY OF CALIFORNIA, NOR
# ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
# ASSUMES ANY LEGAL LIABLITY OR RESPONSIBILITY FOR THE ACCURACY,
# COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
# OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
# PRIVATELY OWNED RIGHTS.
# 
# Additionally, this program is free software; you can distribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of
# the License, or any later version. Accordingly, this program is
# distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
# License for more details. You should have received a copy of the
# GNU General Public License along with this program; if not, write
# to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
# 
# ------------------------------------------------------------
# 
# See the file 'COPYING' for the details of the GNU General Public License
# 
###############################################################################

CXX=mpiCC
CC=mpicc

if DEBUG
D_CXXFLAGS = -Wall -g
else
P_CXXFLAGS = -O2
endif
AM_CXXFLAGS = $(D_CXXFLAGS) $(P_CXXFLAGS) -DINSTALL_PREFIX="\"@prefix@\""
AM_CFLAGS = $(D_CXXFLAGS) $(P_CXXFLAGS) -DINSTALL_PREFIX="\"@prefix@\""

# EXTRA_DIST = mpiformatdb
# bin_SCRIPTS = mpiformatdb

bin_PROGRAMS = mpiblast-g2 mpiformatdb
mpiblastsrc= blastjob.cpp db_spec.cpp db_spec.hpp blast_hooks.c blast_hooks.h \
			embed_rank.cpp embed_rank.hpp file_util.cpp file_util.hpp fragment_list.cpp \
			fragment_list.hpp mpiblast_util.cpp mpiblast_util.hpp mpiblast.cpp blast_hooks.h \
			blastjob.hpp mpiblast.hpp mpiblast_config.cpp mpiblast_config.hpp \
			mpiblast_types.h ncbi_sizeof.c ncbi_sizeof.h distributed_bioseq.h \
			distributed_bioseq.c
mpiblast_SOURCES = $(mpiblastsrc)
EXTRA_mpiblast_SOURCES = getopt.c getopt.h getopt1.c

mpiblast_LDADD = @GETOPT_LONG_SYSTEM@ @CONFIG_NCBIDIR@/build/ncbithr.o  -lncbitool -lncbiobj -lncbi -L@CONFIG_NCBIDIR@/lib/  -lm -lpthread
mpiblast_DEPENDENCIES = @GETOPT_LONG_SYSTEM@
mpiblast_CPPFLAGS = -DUSING_MPI

mpiformatdbsrc= mpiformatdb.cpp mpiblast_config.cpp mpiblast_config.hpp \
				file_util.cpp file_util.hpp mpiblast_util.cpp mpiblast_util.hpp \
				fragment_list.cpp fragment_list.hpp db_spec.hpp db_spec.cpp
mpiformatdb_SOURCES = $(mpiformatdbsrc)

mpiformatdb_LDADD = @GETOPT_LONG_SYSTEM@ formatdb.o @CONFIG_NCBIDIR@/build/ncbithr.o  -lncbitool -lncbiobj -lncbi -L@CONFIG_NCBIDIR@/lib/  -lm -lpthread
mpiformatdb_DEPENDENCIES = @GETOPT_LONG_SYSTEM@ formatdb

mpiformatdb_INCLUDES = -I$(top_srcdir)/include/ -I@CONFIG_NCBIDIR@/include/ 

INCLUDES = -I$(top_srcdir)/include/ -I@CONFIG_NCBIDIR@/include/ 

formatdb:
	cp @CONFIG_NCBIDIR@/demo/formatdb.c ./
	$(CC) -c formatdb.c -I@CONFIG_NCBIDIR@/include/

MOSTLYCLEANFILES=formatdb.c
