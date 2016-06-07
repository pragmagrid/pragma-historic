/******************************************************************************
** distributed_bioseq.h
** Interface for functions to query remote bioseq databases.
** $Id: distributed_bioseq.h,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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
#ifndef __distributed_bioseq_h__
#define __distributed_bioseq_h__

#ifdef __cplusplus
extern "C"{
#endif

#include "ncbi.h"
#include "objseq.h"
#include "ncbithr.h"
#include "ncbiwin.h"
#include "connect/ncbi_core_c.h"

#include "mpi.h"
#include "blast_hooks.h"

/*
 * A set of unique tags for labeling MPI messages. These may need to be changes
 * to prevent collision with other pre-exisiting or user-defined tags.
 */
#define	MPI_BLAST_FETCH		101
#define	MPI_BLAST_FETCH_GI	102
#define	MPI_BLAST_FETCH_LOCAL	103

/*
 * MPI Blast distributed database query functions
 */
Boolean Enable_MPI_Distributed_DB();
Boolean Disable_MPI_Distributed_DB();
Boolean MPI_BSLookup(int dest, MPI_Status *status_ptr);

/*
 * Pack/Unpack an NCBI structure into/outof a memory buffer
 */
unsigned char* PackSeqId(SeqIdPtr sid, int *buf_size);
SeqIdPtr UnpackSeqId(unsigned char* buffer, int buf_size);

unsigned char* PackBioseq(BioseqPtr bsp, int *buf_size);
BioseqPtr UnpackBioseq(unsigned char* buffer, int buf_size);

extern Boolean debug_bsfetch;
extern Boolean debug_bslookup;

#ifdef __cplusplus
}
#endif

#endif // __distributed_bioseq_h__
