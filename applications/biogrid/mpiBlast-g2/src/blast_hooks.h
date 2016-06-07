/******************************************************************************
** blast_hooks.h
** Interface to the NCBI blastall.c functionality
** $Id: blast_hooks.h,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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
#ifndef __mpiblastall_h__
#define __mpiblastall_h__

#ifdef __cplusplus
extern "C"{
#endif

#include <stdlib.h>

#ifndef Int2
#define Int2 short
#endif

#ifndef Int4
#define Int4 int
#endif

#ifndef Nlm_FloatHi
#define	Nlm_FloatHi double
#endif

#ifndef Boolean
#define	Boolean  unsigned char
#endif

#ifndef Uint1
#define	Uint1  unsigned char
#endif

/** Reads command line and initializes BLAST data structures */
Int2 initBLAST ();
/** Executes BLAST and stores successful result ids in result_id_list */
Int2 runBLAST( ValNodePtr* result_id_list );
/** Outputs an HTML header IF the HTML command line option was set */
void outputHtmlHeader();
/** Outputs the results in <code>sappy</code> */
Int2 outputResults( SeqAlignPtr sappy );
/** Outputs an HTML footer IF the HTML command line option was set */
void outputHTMLfooter();
/** Cleans up BLAST data structures */
void cleanupBLAST();


/** Load query sequences from the multi-FastA file */
int loadQueries( void );
/** Free memory used to store the query sequences */
void cleanupQueries( void );

void MpiBlastEnableBioseqFetch();
void MpiBlastDisableBioseqFetch();

/** Stub for the NCBI SeqIdComp function */
int m_SeqIdPtrCmp( SeqIdPtr a, SeqIdPtr b );
/** Stub for the NCBI TxGetSubjectIdFromSeqAlign function */
SeqIdPtr  m_TxGetSubjectIdFromSeqAlign( SeqAlignPtr seqalign );
/** Stub for the NCBI GetScoreAndEvalue function */
Boolean m_GetScoreAndEvalue( SeqAlignPtr seqalign, Int4 *score, Nlm_FloatHi *bit_score, Nlm_FloatHi *evalue, Int4 *number );


#ifdef __cplusplus
}
#endif

#endif // __mpiblastall_h__
