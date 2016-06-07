/******************************************************************************
** distributed_bioseq.c
** Interface for functions to query remote bioseq databases.
** $Id: distributed_bioseq.c,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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
 * Written by J.D. Gans
 */
#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blast.h>
#include <blastpri.h>
#include <simutil.h>
#include <txalign.h>
#include <gapxdrop.h>
#include <sqnutils.h>
#include <xmlblast.h>
#include <mblast.h>
#ifdef BLAST_CS_API
#include <objblst3.h>
#include <netblap3.h>
#endif

#include "distributed_bioseq.h"
#include "blast_hooks.h"
#include "ncbi_sizeof.h"

BSFetchTop bsfetch_orig = NULL;
Boolean debug_bsfetch = FALSE;
Boolean debug_bslookup = FALSE;
const unsigned int timeout_secs = 5;	/**< The number of seconds to wait for workers to respond to a sequence request before timing out */

/* 
 * Note: Half hearted benchmarks with a partially loaded NFS server
 * show that NOT using PACK_FULL_BIOSEQ actually slows the code down!
 * More testing is needed to see what is going on here ...
 */
#define	PACK_FULL_BIOSEQ

/*
 * Fetch the Bioseq that matches the given SeqIdPtr using MPI. The request is sent to 
 * either the node known to possess the bioseq or all worker nodes, which query their 
 * local database fragments. The worker that has the matching Bioseq returns
 * it to the master node (which is doing the requesting).
 *
 * Optimization Possibility; Is it possible to "prefetch" the needed sequences? The identity of each
 * sequence is known from the blast results that each woker collects. To save time, these sequences could
 * be pre-packed into memory buffers, ready to be sent to the master node. More profiling is requried to
 * determine if this would result in an appreciable speed up.
 */
static BioseqPtr LIBCALLBACK MPI_BSFetchFunc(SeqIdPtr sid, Uint1 ld_type)
{
	BioseqPtr bsp = NULL;
	ObjMgrPtr omp = NULL;
	unsigned char* buffer = NULL;
	int buf_size = 0;
	time_t start_time;
	int found_bioseq;

	static int node_count = -1;
	int i = 0;
	MPI_Status status;
	
	if( debug_bsfetch )
		fprintf( stderr, "bsfetching\n" );

	if(sid == NULL){
		return NULL;
	}
		
	/* Assume, for now, that only the master node (rank == 0) will be calling this function */
	/*int rank;*/
	
	/*
	 * Here are the required steps to fetch a bioseq from a remote
	 * process:
	 * 1) Pack the SeqIdPtr into a memory buffer
	 * 2) Send the SeqIdPtr buffer to a node or all nodes (see below).
	 * 3) Wait for a response from the node that has the requested bioseq
	 * 4) Unpack the response into a Bioseq
	 * 5) Return a pointer to the Bioseq
	 */
	
	/* we only need to determine the node_count the first time this function is called. */
	if(node_count < 0){
		MPI_Comm_size(MPI_COMM_WORLD, &node_count);
	}

	switch(sid->choice){
		case SEQID_GI:
			if(sid->extended > 0){
				if( debug_bsfetch )
					fprintf( stderr, "Looking for bioseq GI %d from node %d\n", sid->data.intvalue, sid->extended );
				/* 
			 	 * The extended field contains the rank of the 
				 * worker that has the bioseq we need!
				 */
				if(MPI_Send(&(sid->data.intvalue), 1, MPI_INT, sid->extended, 
					MPI_BLAST_FETCH_GI, MPI_COMM_WORLD) != MPI_SUCCESS){

					 fprintf(stderr, "MPI_BSFetchFunc: Error sending SeqId to %d\n", i);
					 return NULL;
				}
			}
			else{
				if( debug_bsfetch )
					fprintf( stderr, "Looking for bioseq GI %d from everybody\n", sid->data.intvalue );
				/* 
				 * We have no idea which worker has the Bioseq we need,
				 * so we had better ask them all.
			  	 * Send to all ranks except ourselves 
				 * (assume that our rank is 0).
			  	 */
				for(i = 1;i < node_count;i++){
					/* 
					 * We need to send this int as fast as possible.
					 * Don't use MPI_Send (syncronous) as this will block until
					 * a matching receive has been posted). Don't bother with
					 * MPI_Bcast() here (profiling reveals that this for loop 
					 * accounts for < 1% of the time gathering Bioseqs).
					 */
					if(MPI_Send(&(sid->data.intvalue), 1, MPI_INT, i, 
						MPI_BLAST_FETCH_GI, MPI_COMM_WORLD) != MPI_SUCCESS){

						 fprintf(stderr, "MPI_BSFetchFunc: Error sending SeqId to %d\n", i);
						 return NULL;
					}
				}		
			}
			break;
		case SEQID_LOCAL:
			if(sid->extended > 0){
				buffer = PackSeqId(sid, &buf_size);
				if( debug_bsfetch )
					fprintf( stderr, "Looking for local %s bioseq from %d\n",
							((ObjectIdPtr)(sid->data.ptrvalue))->str,  sid->extended );

				/*
				 * The extended field contains the rank of the
				 * worker that has the bioseq we need!
				 */
				if(MPI_Send(buffer, buf_size, MPI_BYTE, sid->extended,
                        MPI_BLAST_FETCH_LOCAL, MPI_COMM_WORLD) != MPI_SUCCESS){

					fprintf(stderr, "MPI_BSFetchFunc: Error sending SeqId to %d\n", sid->extended);
					return NULL;
                }

                /*
                 * Free the send buffer so we can use it to receive the incoming Bioseq
                 * data.
                 */
                free(buffer);
			}else{
				/*
				 * Don't know where this Bioseq came from, try looking it up locally
				 */
				bsp = bsfetch_orig( sid, ld_type );
				if( debug_bsfetch ){
					if( bsp )
						fprintf( stderr, "Found local bioseq %s on this node\n",
                                        ((ObjectIdPtr)(sid->data.ptrvalue))->str );
					else
						fprintf( stderr, "Couldn't find local bioseq %s on this node\n",
                                        ((ObjectIdPtr)(sid->data.ptrvalue))->str );
				}
				return bsp;
			}

	        break;
		case SEQID_GENERAL:
			/* 
			 * The database has not been created with the -O T flag! In principle 
			 * it should be possible to convert the indicies stored in a 
			 * SEQID_GENERAL SeqId to match the indicies of a given fragment, but it
			 * hardly seems worth it when the problem can be fixed by formatting the 
			 * database with the -O T flag.
			 * The index currently stored referes to the entire database, not the 
			 * individual fragment that contains the sequence. We have a "global"
			 * address (to the entire database) not a "local" address (to a
			 * particular fragment).
			 */
			bsp = bsfetch_orig( sid, ld_type );
			
			return bsp;
			
			break;
		default:
			buffer = PackSeqId(sid, &buf_size);

			if(sid->extended > 0){
				if( debug_bsfetch )
					fprintf( stderr, "Default: looking for packed bioseq from %d\n", sid->extended );
				/* 
				 * The extended field contains the rank of the 
				 * worker that has the bioseq we need!
				 */
				if(MPI_Send(buffer, buf_size, MPI_BYTE, sid->extended, 
					MPI_BLAST_FETCH_LOCAL, MPI_COMM_WORLD) != MPI_SUCCESS){

					fprintf(stderr, "MPI_BSFetchFunc: Error sending SeqId to %d\n", i);
					return NULL;
				}
			}
			else{
		        /* 
	        	 * Send this buffer to all ranks except ourselves 
	        	 * (assume that our rank is 0).
	       		 */
				if( debug_bsfetch )
					fprintf( stderr, "Default: looking for packed bioseq from everybody\n" );
				for(i = 1;i < node_count;i++){
				    /* 
				 	 * We need to send this buffer as fast as possible.
					 * Don't use MPI_Send (synchronous) as this will block until
					 * a matching receive has been posted). Don't bother with
					 * MPI_Bcast() here (profiling reveals that this for loop 
					 * accounts for < 1% of the time gathering Bioseqs).
					 */
					if(MPI_Send(buffer, buf_size, MPI_BYTE, i, 
						MPI_BLAST_FETCH_LOCAL, MPI_COMM_WORLD) != MPI_SUCCESS){

						 fprintf(stderr, "MPI_BSFetchFunc: Error sending SeqId to %d\n", i);
						 return NULL;
					}
				}
			}
		
			/*
    	 	 * Free the send buffer so we can use it to receive the incoming Bioseq
    		 * data.
	   		 */
			free(buffer);
		
			break;
	};
	
	/**
	 * For some reason, a bioseq may not be found on the workers.  Rather than let
	 * mpiBLAST stall out, use a spin-wait timeout loop until this problem can be debugged.
	 */
	start_time = time(NULL);
	found_bioseq = 0;
	while( start_time + timeout_secs > time(NULL) ){
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_BLAST_FETCH, MPI_COMM_WORLD, &found_bioseq, &status);
		if( found_bioseq )
			break;
	}
	if( !found_bioseq ){
		// timed out trying to get the bioseq from a worker, try to find it in the shared
		// database
		fprintf( stderr, "Error: Timed out waiting for biosequence from workers\n" );
		bsp = bsfetch_orig( sid, ld_type );
		return bsp;
	}
	MPI_Get_count(&status, MPI_BYTE, &buf_size);
	
	/*
	 * An optimization idea -- can we use a persistent bioseq buffer? Reallocation would
	 * only be required when the incoming message is larger than the last message. This
	 * could cut down on the overhead of memory management.
	 */
	buffer = (unsigned char*)calloc(buf_size, sizeof(unsigned char));
	
	if(buffer == NULL){
		fprintf(stderr, "MPI_BSFetchFunc: Error allocating memory for Bioseq\n");
		return NULL;
	}
	
	/* #define	USE_NON_BLOCKING_DB_LOOKUP */
	
	#ifdef USE_NON_BLOCKING_DB_LOOKUP
	{
	MPI_Request hack_req;
	int rec_status;
	double wait;
	const double timeout = 0.5; // sec
	
	MPI_Irecv(buffer, buf_size, MPI_BYTE, status.MPI_SOURCE, MPI_BLAST_FETCH,
		MPI_COMM_WORLD, &hack_req);
	
	wait = MPI_Wtime();
	
	// Now test for the completion of this receive
	while(1){
		MPI_Test(&hack_req, &rec_status, MPI_STATUS_IGNORE);
		
		if(rec_status){
			// Success!
			break;
		}
		
		if(MPI_Wtime() - wait > timeout){
			FILE *fout;
			
			// Tell the user we've timed out ...
			fprintf(stderr, "MPI_BSFetchFunc: Timeout (%.3f sec)\n", 
				timeout);
			
			fprintf(stderr, "Waiting for %d\n",  status.MPI_SOURCE);
			fprintf(stderr, "uchar size is %d\n",  sizeof(unsigned char));
			fprintf(stderr, "buff size = %d\n", buf_size);
			fprintf(stderr, "sid->choice = %d\n", sid->choice);
			fprintf(stderr, "sid->data = %d\n", sid->data.intvalue);
			fprintf(stderr, "sid->next = %p\n", sid->next);
			
			MPI_Cancel(&hack_req);
			
			fprintf(stderr, "Canceled MPI_Irecv\n");
			
			// Free the receive buffer
			free(buffer);
			
			fout = fopen("/tmp/packetdump.out", "w");
			
			if(fout == NULL){
				fprintf(stderr, "Unable to write buffer info\n");
			}
			else{
				fwrite(buffer, sizeof(unsigned char), buf_size, fout);
				
				fclose(fout);
			}
			
			// Retry
			return MPI_BSFetchFunc(sid, ld_type);
		}
	}
	}				
	#else	
	if(MPI_Recv(buffer, buf_size, MPI_BYTE, status.MPI_SOURCE, 
		   MPI_BLAST_FETCH, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){
		   
		fprintf(stderr, "MPI_BSFetchFunc: Error receiving Bioseq\n");
		return NULL;
	}
	
	#endif /* USE_NON_BLOCKING_DB_LOOKUP */
	
	bsp = UnpackBioseq(buffer, buf_size);
	
	#ifndef PACK_FULL_BIOSEQ
	bsp->id = sid;
	#endif /* PACK_FULL_BIOSEQ */
	
	free(buffer);
	
	if(bsp == NULL){
		fprintf(stderr, "MPI_BSFetchFunc: Error unpacking Bioseq\n");
		return NULL;
	}
	
	/* 
	 * Tell the object manager that this is only a temporary structure and can be
	 * "reaped" as needed. Note that we can do better than letting NCBI reap
	 * old structures. Use a circular buffer to periodically delete Bioseqs that
	 * we have provided (using BioseqFree perhaps?).
	 */
	if(ld_type == BSFETCH_TEMP){
		omp = ObjMgrWriteLock();
		ObjMgrSetTempLoad(omp, (Pointer)bsp);	
		ObjMgrUnlock();
	}
	
	//#define	SET_TRIPWIRE
	
	#ifdef SET_TRIPWIRE
	{
	int tripwire = 0; // false
	
	MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, 
		&tripwire, MPI_STATUS_IGNORE);
	
	if(tripwire){
		BioseqPtr bsp2;
		AsnIoPtr aip;
		
		fprintf(stderr, "****** Double tap detected! ******\n");
		fprintf(stderr, "sid->choice = %d\n", sid->choice);
		fprintf(stderr, "sid->extended = %d\n", sid->extended);
		fprintf(stderr, "sid->data = %d\n", sid->data.intvalue);
		if(sid->choice == SEQID_LOCAL){
			fprintf(stderr, "sid->data (str) = %s\n",
				(char*)(((ObjectIdPtr)(sid->data.ptrvalue))->str));
			fprintf(stderr, "sid->data (id) = %d\n",
				(((ObjectIdPtr)(sid->data.ptrvalue))->id));
		}
		
		fprintf(stderr, "sid->next = %p\n", sid->next);
		
		fprintf(stderr, "Loading second seq ...\n");
		
		MPI_Probe(MPI_ANY_SOURCE, MPI_BLAST_FETCH, MPI_COMM_WORLD, &status);	
		MPI_Get_count(&status, MPI_BYTE, &buf_size);
		
		buffer = (unsigned char*)calloc(buf_size, sizeof(unsigned char));
	
		if(buffer == NULL){
			fprintf(stderr, "MPI_BSFetchFunc: Error allocating memory for Bioseq\n");
			return NULL;
		}
	
		if(MPI_Recv(buffer, buf_size, MPI_BYTE, status.MPI_SOURCE, 
		   MPI_BLAST_FETCH, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){
		   
			fprintf(stderr, "MPI_BSFetchFunc: Error receiving Bioseq\n");
			return NULL;
		}
		
		bsp2 = UnpackBioseq(buffer, buf_size);
		
		aip = AsnIoOpen("/tmp/bioseq_doubletap.txt", "w");
		
		if(!aip){
			fprintf(stderr, "Can't open /tmp/bioseq_doubletap.txt\n");
		}
		else{
			BioseqAsnWrite(bsp, aip, NULL);
			fprintf(stderr, "Wrote bsp\n");
			BioseqAsnWrite(bsp2, aip, NULL);
			fprintf(stderr, "Wrote bsp2\n");
			
			AsnIoFlush( aimp->aip );
			AsnIoClose(aip);
		}
		
		// Wait for help!
		MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	}
	#endif //SET_TRIPWIRE
	
	return bsp;
}

/* Use this function to enable distributed database searching with MPI */
Boolean Enable_MPI_Distributed_DB()
{
	SeqMgrPtr smp;
	
	smp = SeqMgrWriteLock();
	
	if(smp == NULL){
		return FALSE;
	}
	
	/* Save the original lookup function */
	bsfetch_orig = smp->bsfetch;
	smp->bsfetch = MPI_BSFetchFunc;
	
	return SeqMgrUnlock();
}

/* Use this function to disable distributed database searching with MPI */
Boolean Disable_MPI_Distributed_DB()
{
	SeqMgrPtr smp;
	
	smp = SeqMgrWriteLock();
	
	if(smp == NULL){
		return FALSE;
	}
	
	if(bsfetch_orig == NULL){
		return FALSE;
	}
	
	/* Restore the default fetch function -- see NCBI's seqmgr.c for details */
	smp->bsfetch = bsfetch_orig;
	
	return SeqMgrUnlock();
}


/*
 * Respond to a remote database query and return the requested Bioseq (if found).
 * Return True if a bioseq is found, false otherwise.
 */
Boolean MPI_BSLookup(int dest, MPI_Status *status_ptr)
{
	BioseqPtr bsp = NULL;
	unsigned char* buffer = NULL;
	int buf_size = 0;
	SeqId si;
	SeqIdPtr sid;
		
	if( debug_bslookup )
		fprintf( stderr, "bslooking_up\n" );
	/*
	 * Here are the required steps:
	 * 1) Unpack the SeqId that identifies the desired Bioseq
	 * 2) Lookup the Bioseq referred to by the SeqId
	 * 3) If found, return the Bioseq
	 */
	if(status_ptr->MPI_TAG == MPI_BLAST_FETCH_GI){
		memset(&si, 0, sizeof(SeqId));
		
		si.choice = SEQID_GI;

		if( debug_bslookup )
			fprintf(stderr, "Receiving GI request...\n");
		MPI_Recv(&(si.data.intvalue), 1, MPI_INT, dest, 
			MPI_BLAST_FETCH_GI, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
		sid = &si;
	}
	else{ // msg == MPI_BLAST_FETCH_LOCAL
	
		if( debug_bslookup )
			fprintf(stderr, "Receiving Local request...\n");
		MPI_Get_count(status_ptr, MPI_BYTE, &buf_size);
	
		buffer = (unsigned char*)calloc(buf_size, sizeof(unsigned char));

		if(buffer == NULL){
			fprintf(stderr, "MPI_BSLookup: Error allocating memory for SeqId\n");
			return FALSE;
		}
		
		MPI_Recv(buffer, buf_size, MPI_BYTE, dest, 
			MPI_BLAST_FETCH_LOCAL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		sid = UnpackSeqId(buffer, buf_size);

		free(buffer);

		if(sid == NULL){
			fprintf(stderr, "MPI_BSLookup: Error unpacking SeqId\n");
			return FALSE;
		}
	}
			
	/* 
	 * Look up the Bioseq
	 */
	if( debug_bslookup )
		fprintf( stderr, "retrieving Bioseq\n" );
	bsp = BioseqLockById(sid);
	
	if(bsp == NULL){		
		/* Free the SeqIdPtr if we need to */
		if(sid != &si){
			SeqIdFree(sid);
		}
		
		/* The Bioseq was not found, there is nothing further to do ... */
		return FALSE;
	}
		
	/*
	 * We found the Bioseq -- pack it up and prepare for sending
	 */ 
	if( debug_bslookup )
		fprintf( stderr, "packing Bioseq\n" );
	buffer = PackBioseq(bsp, &buf_size);
	
	if(buffer == NULL){
		fprintf(stderr, "MPI_BSLookup: Error packing Bioseq\n");
		return FALSE;
	}
		
	#ifdef USE_NON_BLOCKING_DB_LOOKUP
	{
	MPI_Request send_req;
	double wait;
	const double timeout = 0.5;
	int send_status;
	int rank;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	MPI_Irsend(buffer, buf_size, MPI_BYTE, dest, 
		  MPI_BLAST_FETCH, MPI_COMM_WORLD, &send_req);
	
	BioseqUnlock(bsp);
	
	/* Free the SeqIdPtr if we allocated it */
	
	if(sid != &si){
		SeqIdFree(sid);
	}
	
	wait = MPI_Wtime();
	
	// Now test for the completion of this receive
	while(1){
		MPI_Test(&send_req, &send_status, MPI_STATUS_IGNORE);
		
		if(send_status){
			// Success!
			break;
		}
		
		if(MPI_Wtime() - wait > timeout){
			// Tell the user we've timed out ...
			fprintf(stderr, "[%d] MPI_BSLookup: Timeout (%.3f sec)\n", 
				rank, timeout);
						
			MPI_Cancel(&send_req);
			
			fprintf(stderr, "Canceled MPI_Irsend\n");
			
			// Free the receive buffer
			free(buffer);
			
			// Return to the main wait loop and listen for the resend msg
			// from the master
			return FALSE;
		}
	}
	}
	#else
	if( debug_bslookup )
		fprintf( stderr, "Sending Bioseq\n" );
	MPI_Send(buffer, buf_size, MPI_BYTE, dest, 
		MPI_BLAST_FETCH, MPI_COMM_WORLD);
		
	BioseqUnlock(bsp);
	
	/* Free the SeqIdPtr if we allocated it */
	if(sid != &si){
		SeqIdFree(sid);
	}
	
	#endif //USE_NON_BLOCKING_DB_LOOKUP

	free(buffer);
		
	return TRUE;
}

/*
 * Write the contents of the SeqId structure to a memory buffer. The buffer
 * is returned and the buffer size is updated. The user is responsible for
 * freeing the memory returned by this function. 
 */
unsigned char* PackSeqId(SeqIdPtr sid, int *buf_size)
{
	unsigned char* buffer = NULL;
	AsnIoMemPtr aimp;
	
	if(sid == NULL){
		fprintf(stderr, "PackSeqId:NULL SeqId pointer\n");
	}
	
	*buf_size = SeqIdSize(sid);
	
	if(*buf_size <= 0){
		fprintf(stderr, "PackSeqId: Error computing SeqId Size\n");
		return NULL;
	}else if( debug_bsfetch || debug_bslookup ){
		fprintf(stderr, "Got seq id size: %d\n", *buf_size );
	}
	
	buffer = (unsigned char*)calloc(*buf_size, sizeof(unsigned char));
	
	if(buffer == NULL){
		fprintf(stderr, "PackSeqId: Error allocating memory for SeqId\n");
		return NULL;
	}
		
	#ifdef USE_NCBI_ASCII
	/*
	 * Use ascii ASN encoding for debugging
	 */
	aimp = AsnIoMemOpen("w", buffer, *buf_size);
	#else
	/*
	 * Use binary ASN encoding to save space
	 */
	aimp = AsnIoMemOpen("wb", buffer, *buf_size);
	#endif // USE_NCBI_ASCII
	
	if(aimp == NULL){
		fprintf(stderr, "PackSeqId: AsnIoMemOpen Error\n");
		return NULL;
	}
	
	if( !SeqIdAsnWrite(sid, aimp->aip, NULL) ){
		fprintf(stderr, "PackSeqId: Error writing SeqId to buffer\n");
		return NULL;
	}
	AsnIoFlush( aimp->aip );
	AsnIoMemClose(aimp);
		
	return buffer;
}

SeqIdPtr UnpackSeqId(unsigned char* buffer, int buf_size)
{
	AsnIoMemPtr aimp;
	SeqIdPtr sid = NULL;
	
	#ifdef USE_NCBI_ASCII
	/*
	 * Use ascii ASN.1 encoding for debugging
	 */
	aimp = AsnIoMemOpen("r", buffer, buf_size);
	#else
	/*
	 * Use binary ASN.1 encoding to save space
	 */
	aimp = AsnIoMemOpen("rb", buffer, buf_size);
	#endif // USE_NCBI_ASCII
	
	
	if(aimp == NULL){
		fprintf(stderr, "UnpackSeqId: AsnIoMemOpen Error\n");
		return NULL;
	}
	
	sid = SeqIdAsnRead(aimp->aip, NULL);
	
	if(sid == NULL){
		fprintf(stderr, "UnpackSeqId: Error reading SeqId from buffer\n");
		return NULL;
	}
	
	AsnIoMemClose(aimp);
	
	return sid;
}

unsigned char* PackBioseq(BioseqPtr bsp, int *buf_size)
{
	
	unsigned char* buffer = NULL;
	AsnIoMemPtr aimp;
	
	#ifdef PACK_FULL_BIOSEQ
	*buf_size = BioseqSize(bsp);
	#else
	*buf_size = BioseqDataSize(bsp);
	#endif /*PACK_FULL_BIOSEQ*/
	
	if(*buf_size <= 0){
		fprintf(stderr, "PackBioseq: Error computing Bioseq Size\n");
		return NULL;
	}
	
	buffer = (unsigned char*)calloc(*buf_size, sizeof(unsigned char));
	
	if(buffer == NULL){
		fprintf(stderr, "PackBioseq: Error allocating memory for Bioseq\n");
		return NULL;
	}
	
	#ifdef USE_NCBI_ASCII
	/*
	 * Use ascii ASN.1 encoding for debugging
	 */
	aimp = AsnIoMemOpen("w", buffer, *buf_size);
	#else
	/*
	 * Use binary ASN.1 encoding to save space
	 */
	aimp = AsnIoMemOpen("wb", buffer, *buf_size);
	#endif
	
	if(aimp == NULL){
		fprintf(stderr, "PackBioseq: AsnIoMemOpen Error\n");
		return NULL;
	}
	
	#ifdef PACK_FULL_BIOSEQ
	if( !BioseqAsnWrite(bsp, aimp->aip, NULL) ){
		fprintf(stderr, "PackBioseq: Error writing Bioseq to buffer\n");
		return NULL;
	}
	#else
	if( !BioseqInstAsnWrite(bsp, aimp->aip, NULL) ){
		fprintf(stderr, "PackBioseq: Error writing Bioseq to buffer\n");
		return NULL;
	}
	#endif /*PACK_FULL_BIOSEQ*/
	AsnIoFlush( aimp->aip );
	AsnIoMemClose(aimp);
	
	return buffer;
}

BioseqPtr UnpackBioseq(unsigned char* buffer, int buf_size)
{
	AsnIoMemPtr aimp;
	BioseqPtr bsp = NULL;
	
	#ifdef USE_NCBI_ASCII
	/*
	 * Use ascii ASN.1 encoding for debugging
	 */
	aimp = AsnIoMemOpen("r", buffer, buf_size);
	#else
	/*
	 * Use binary ASN.1 encoding to save space
	 */
	aimp = AsnIoMemOpen("rb", buffer, buf_size);
	#endif
	
	if(aimp == NULL){
		fprintf(stderr, "UnpackBioseq: AsnIoMemOpen Error\n");
		return NULL;
	}
	
	#ifdef PACK_FULL_BIOSEQ
	bsp = BioseqAsnRead(aimp->aip, NULL);
	#else
	bsp = BioseqNew();
	BioseqInstAsnRead(bsp, aimp->aip, NULL);
	#endif /* PACK_FULL_BIOSEQ */
	
	if(bsp == NULL){
		fprintf(stderr, "UnpackBioseq: Error reading Bioseq from buffer\n");
		return NULL;
	}
	
	AsnIoMemClose(aimp);
	
	return bsp;

}
//#endif /* USING_MPI */
