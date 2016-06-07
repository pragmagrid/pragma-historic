/******************************************************************************
** ncbi_sizeof.c
** Implementation of functions to calculate sizes of various NCBI toolkit structs
** $Id: ncbi_sizeof.c,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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
/*
 * This file written by J.D. Gans
 */
#include "ncbi_sizeof.h"

/* File scope variables */
long byte_count;

#ifdef WIN32
char* bitbucket = "NUL";
#else
char* bitbucket = "/dev/null";
#endif

/* Local functions */
int
LIBCALLBACK
CountBytes(Pointer iostruct, CharPtr buff, Uint2 bytes_to_write);

/* 
 * This is a dummy write function. It doesn't write anything but counts what would have been
 * written instead.
 */
int 
LIBCALLBACK
CountBytes(Pointer iostruct, CharPtr buff, Uint2 bytes_to_write)
{
	byte_count += bytes_to_write;
	
	return bytes_to_write;
}

int SeqIdSize(SeqIdPtr sid)
{
	if(sid == NULL){
		return 0;
	}
	
	/*
	     * Compute the size of the SeqId structure pointed to by sid
	     */
	return NCBISizeOf(sid, (AsnWriteFunc)SeqIdAsnWrite);
}

int BioseqSize(BioseqPtr bsp)
{
	if(bsp == NULL){
		return 0;
	}
	
	/*
	     * Compute the size of the Bioseq structure pointed to by bsp
	     */ 
	return NCBISizeOf(bsp, (AsnWriteFunc)BioseqAsnWrite);
}

int BioseqDataSize(BioseqPtr bsp)
{
	if(bsp == NULL){
		return 0;
	}
	
	/*
	     * Compute the size of the Bioseq structure pointed to by bsp
	     */ 
	return NCBISizeOf(bsp, (AsnWriteFunc)BioseqInstAsnWrite);

}

int SeqAnnotSize(SeqAnnotPtr sap)
{
	if(sap == NULL){
		return 0;
	}
	
	/*
	     * Compute the size of the Bioseq structure pointed to by bsp
	     */ 
	return NCBISizeOf(sap, (AsnWriteFunc)SeqAnnotAsnWrite);
}

/*
 * Compute the size of the NCBI structure pointed to by ptr. Return the size on
 * success and 0 on failure.
 */
int NCBISizeOf(Pointer ptr, AsnWriteFunc writefunc)
{
	static int mutex = FALSE;
	AsnIoPtr aip;
	
	while(mutex){
		/* Another thread is using this function! */
		fprintf(stderr, "NCBISizeOf: Unable to lock thread\n");
	}
	
	/* Give this thread a lock */
	mutex = TRUE;
	
	/*
     * Open a dummy AsnIoPtr (we don't actually want to write any data).
     */
	#ifdef USE_NCBI_ASCII
	aip = AsnIoOpen( bitbucket, "w");
	#else
	aip = AsnIoOpen( bitbucket, "wb");
	#endif
	
	if(aip == NULL){
		fprintf(stderr, "NCBISizeOf: Unable to get valid AsnIoPtr\n");
		return -1;
	}
	
	/* Zero the count buffer (a file scope variable) */
	byte_count = 0L;
	
	/* Insert our byte-counting function */
	aip->writefunc = (AsnIoFunc)CountBytes;
	
	if( !(*writefunc)(ptr, aip, NULL) ){
		return -1;
	}
	
	/* 
	 * Close the ASN file BEFORE checking byte_count to insure that all data
	 * has been flushed! 
	 */
	AsnIoClose(aip);
	
	/* free the lock */
	mutex = FALSE;
	
	return (int)(byte_count);
}
