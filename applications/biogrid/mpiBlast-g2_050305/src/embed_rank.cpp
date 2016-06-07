/******************************************************************************
** embed_rank.cpp
** Implements rank embedding into NCBI SeqAlign structs
** $Id: embed_rank.cpp,v 1.1.1.2 2005/01/27 03:02:21 cwwang Exp $
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

/* Implemented by J.D. Gans */
#include "embed_rank.hpp"
#include "mpiblast_util.hpp"
#include <iostream>
using namespace std;

extern "C" {
#include "sequtil.h"
extern SeqIdPtr SeqLocId (SeqLocPtr anp);
}

void setSeqIdRank( SeqIdPtr sip_p, int rank, int node_count );
void setSeqIdRank( SeqIdPtr sip_p, int rank, int node_count ){
	while(sip_p != NULL){
		// Note that the extended member variable is an unsigned short, so
		// this trick of embeding the worker rank will only work if the
		// node_count is less than 256 (which corresponds to less than 255
		// workers).
		if(node_count < 256){
			sip_p->extended = rank;
		}
		
		sip_p = sip_p->next;
	}
}

// Store the "rank" in the "extended" field of all of the SeqId structures contained
// in the given SeqAlignPtr. This will allow a direct query of the worker that
// generated a given blast hit and disambiguate blast hits with the same local
// id string.
void setSeqAlignRank( SeqAlignPtr sap_p, int rank, int node_count )
{
	
    SeqIdPtr sip;
    DenseSegPtr dsp;
    DenseDiagPtr ddp;
    StdSegPtr ssp;
    SeqLocPtr slp;
    SeqAlignPtr sap;

    while(sap_p !=NULL) {
        switch(sap_p->segtype) {
        case SAS_DENDIAG:
            ddp = (DenseDiagPtr)(sap_p->segs);
            while(ddp) {
				setSeqIdRank( ddp->id, rank, node_count );
                ddp = ddp->next;
            }
            break;
            
        case SAS_DENSEG:
            dsp = (DenseSegPtr)(sap_p->segs);
			setSeqIdRank( dsp->ids, rank, node_count );
            break;
            
        case SAS_STD:
            ssp = (StdSegPtr)(sap_p->segs);
			setSeqIdRank( ssp->ids, rank, node_count );
            for(slp = ssp->loc; slp != NULL; slp = slp->next) {
                sip = SeqLocId(slp);
				setSeqIdRank( sip, rank, node_count );
            }
            break;
		case SAS_PACKED:
			setSeqIdRank( ((PackSegPtr)sap_p->segs)->ids, rank, node_count );
			break;
        case SAS_DISC:
            for(sap = (SeqAlignPtr)sap_p->segs; sap != NULL; sap = sap->next) {
                setSeqAlignRank( sap, rank, node_count );
            }
            break;

        default:
			LOG_MSG << "Error:  unable to embed rank, segtype is " << sap_p->segtype << endl;
            break;
        }
        sap_p = sap_p->next;
    }
}

void setSeqAnnotRank( SeqAnnotPtr sap, int rank, int node_count ){
	setSeqAlignRank( (SeqAlignPtr)sap->data, rank, node_count );
}
