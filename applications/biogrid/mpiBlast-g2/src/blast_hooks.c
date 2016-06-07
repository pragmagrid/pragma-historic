/* $Id: blast_hooks.c,v 1.2 2004/10/19 09:07:11 hclee Exp $
**************************************************************************
*                                                                         *
*                             COPYRIGHT NOTICE                            *
*                                                                         *
* This software/database is categorized as "United States Government      *
* Work" under the terms of the United States Copyright Act.  It was       *
* produced as part of the author's official duties as a Government        *
* employee and thus can not be copyrighted.  This software/database is    *
* freely available to the public for use without a copyright notice.      *
* Restrictions can not be placed on its present or future use.            *
*                                                                         *
* Although all reasonable efforts have been taken to ensure the accuracy  *
* and reliability of the software and data, the National Library of       *
* Medicine (NLM) and the U.S. Government do not and can not warrant the   *
* performance or results that may be obtained by using this software,     *
* data, or derivative works thereof.  The NLM and the U.S. Government     *
* disclaim any and all warranties, expressed or implied, as to the        *
* performance, merchantability or fitness for any particular purpose or   *
* use.                                                                    *
*                                                                         *
* In any work or product derived from this material, proper attribution   *
* of the author(s) as the source of the software or data would be         *
* appreciated.                                                            *
*                                                                         *
************************************************************************** 
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

#include "blast_hooks.h"

/* Local function definition */
void SeqIdDataFree(SeqIdPtr sip);

/* Used by the callback function. */
FILE *global_fp=NULL;
/*
	Callback to print out ticks, in UNIX only due to file systems
	portability issues.
*/

static int LIBCALLBACK
tick_callback(Int4 sequence_number, Int4 number_of_positive_hits)

{
#ifdef OS_UNIX
    /* #ifndef BLAST_CS_API */
    fprintf(global_fp, "%s", ".");
    fflush(global_fp);
    /* #endif */
#endif
    return 0;
}

static Int2
BlastGetMaskingLoc(FILE *infp, FILE *outfp, CharPtr instructions)
{
	BioseqPtr bsp;
	Char buffer[50];
	SeqEntryPtr sep;
	SeqLocPtr slp, slp_start, tmp_slp;

	if (infp == NULL || outfp == NULL || instructions == NULL)
		return 1;

	while ((sep=FastaToSeqEntryEx(infp, TRUE, NULL, TRUE)) != NULL) 
	{
		bsp = NULL;
		SeqEntryExplore(sep, &bsp, FindNuc);

		if (bsp == NULL)
		{
	  	 	ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
	   		return 2;
		}
		SeqIdWrite(bsp->id, buffer, PRINTID_FASTA_LONG, 50);
		fprintf(outfp, ">%s\n", buffer);
		slp_start = slp = BlastBioseqFilter(bsp, instructions);
        	while (slp)
        	{
               		tmp_slp=NULL;
               		while((tmp_slp = SeqLocFindNext(slp, tmp_slp))!=NULL)
               	 	{
				fprintf(outfp, "%ld %ld\n", (long) (1+SeqLocStart(tmp_slp)), (long) (1+SeqLocStop(tmp_slp)));
                 	}
                	slp = slp->next;
        	}

		slp_start = SeqLocSetFree(slp_start);
		sep = SeqEntryFree(sep);
	}

	return 0;
}

#define DO_NOT_SUPPRESS_BLAST_OP

#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs[] = {
    { "Program Name",           /* 0 */
      NULL, NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
    { "Database",               /* 1 */
      "nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},
    { "Query File",             /* 2 */
      "stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
    { "Expectation value (E)",  /* 3 */
      "10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},
	{ "alignment view options:\n0 = pairwise,\n1 = query-anchored showing identities,\n2 = query-anchored no identities,\n3 = flat query-anchored, show identities,\n4 = flat query-anchored, no identities,\n5 = query-anchored no identities and blunt ends,\n6 = flat query-anchored, no identities and blunt ends,\n7 = XML Blast output,\n8 = tabular, \n9 tabular with comment lines\n10 ASN, text\n11 ASN, binary", /* 4 */
      "0", NULL, NULL, FALSE, 'm', ARG_INT, 0.0, 0, NULL},
    { "BLAST report Output File", /* 5 */
      "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
    { "Filter query sequence (DUST with blastn, SEG with others)", /* 6 */
      "T", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},
    { "Cost to open a gap (zero invokes default behavior)", /* 7 */
      "0", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL},
    { "Cost to extend a gap (zero invokes default behavior)", /* 8 */
      "0", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL},
    { "X dropoff value for gapped alignment (in bits) (zero invokes default "
      "behavior)\n      blastn 30, megablast 20, tblastx 0, all others 15", /* 9 */
      "0", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},
    { "Show GI's in deflines",  /* 10 */
      "F", NULL, NULL, FALSE, 'I', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Penalty for a nucleotide mismatch (blastn only)", /* 11 */
      "-3", NULL, NULL, FALSE, 'q', ARG_INT, 0.0, 0, NULL},
    { "Reward for a nucleotide match (blastn only)", /* 12 */
      "1", NULL, NULL, FALSE, 'r', ARG_INT, 0.0, 0, NULL},
    { "Number of database sequences to show one-line descriptions for (V)", /* 13 */
      "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},
    { "Number of database sequence to show alignments for (B)", /* 14 */
      "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},
    { "Threshold for extending hits, default if zero\n" /* 15 */
      "      blastp 11, blastn 0, blastx 12, tblastn 13\n"
      "      tblastx 13, megablast 0",
      "0", NULL, NULL, FALSE, 'f', ARG_INT, 0.0, 0, NULL},
    { "Perfom gapped alignment (not available with tblastx)", /* 16 */
        "T", NULL, NULL, FALSE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Query Genetic code to use", /* 17 */
      "1", NULL, NULL, FALSE, 'Q', ARG_INT, 0.0, 0, NULL},
    { "DB Genetic code (for tblast[nx] only)", /* 18 */
      "1", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL},
    { "Number of processors to use", /* 19 */
      "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},
    { "SeqAlign file",          /* 20 */
      NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},
    { "Believe the query defline", /* 21 */
      "F", NULL, NULL, FALSE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Matrix",                 /* 22 */
      "BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},
    { "Word size, default if zero (blastn 11, megablast 28, "
        "all others 3)", /* 23 */
      "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},
    { "Effective length of the database (use zero for the real size)", /* 24 */
      "0", NULL, NULL, FALSE, 'z', ARG_FLOAT, 0.0, 0, NULL},
    { "Number of best hits from a region to keep (off by default, if used a value of 100 is recommended)", /* 25 */
      "0", NULL, NULL, FALSE, 'K', ARG_INT, 0.0, 0, NULL},
    { "0 for multiple hit, 1 for single hit",/* 26 */
       "0",  NULL, NULL, FALSE, 'P', ARG_INT, 0.0, 0, NULL},
    { "Effective length of the search space (use zero for the real size)", /* 27 */
      "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},
    { "Query strands to search against database (for blast[nx], and tblastx)\n"
      "       3 is both, 1 is top, 2 is bottom", /* 28 */
      "3", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},
    { "Produce HTML output",    /* 29 */
      "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Restrict search of database to list of GI's",              /* 30 */
      NULL, NULL, NULL, TRUE, 'l', ARG_STRING, 0.0, 0, NULL},
    {"Use lower case filtering of FASTA sequence", /* 31 */
     "F", NULL,NULL,TRUE,'U',ARG_BOOLEAN, 0.0,0,NULL},
    { "X dropoff value for ungapped extensions in bits (0.0 invokes default "
      "behavior)\n      blastn 20, megablast 10, all others 7", /* 32 */
      "0.0", NULL, NULL, FALSE, 'y', ARG_FLOAT, 0.0, 0, NULL},
    { "X dropoff value for final gapped alignment in bits " /* 33 */
      "(0.0 invokes default behavior)\n"
      "      blastn/megablast 50, tblastx 0, all others 25",
      "0", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},
    { "PSI-TBLASTN checkpoint file", /* 34 */
      NULL, NULL, NULL, TRUE, 'R', ARG_FILE_IN, 0.0, 0, NULL},
    { "MegaBlast search",       /* 35 */
      "F", NULL, NULL, FALSE, 'n', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Location on query sequence",/* 36 */
      NULL, NULL, NULL, TRUE, 'L', ARG_STRING, 0.0, 0, NULL},
    { "Multiple Hits window size, default if zero (blastn/megablast 0, "
        "all others 40", /* 37 */
      "0", NULL, NULL, FALSE, 'A', ARG_INT, 0.0, 0, NULL},
#ifdef DO_NOT_SUPPRESS_BLAST_OP
    { "Frame shift penalty (OOF algorithm for blastx)", /* 38 */
      "0", NULL, NULL, FALSE, 'w', ARG_INT, 0.0, 0, NULL},
    { "Length of the largest intron allowed in tblastn for linking HSPs (0 disables linking)", /* 39 */
      "0", NULL, NULL, FALSE, 't', ARG_INT, 0.0, 0, NULL}, 
#endif
/*--KM
   seems ok to add another param b/c NUMARG is defined based on
    sizeof(myargs) itself
   made optional=TRUE but this may change?
*/
    { "Number of concatenated queries, for blastn and tblastn", /* 40 */
      "0", NULL, NULL, TRUE, 'B', ARG_INT, 0.0, 0, NULL}
};




/* Needed for Mega BLAST only */
#define MAX_NUM_QUERIES 16383 /* == 1/2 INT2_MAX */

/*
 * these variables were local to blastall's main function,
 * they were made global so the main function can be split up.
 */
AsnIoPtr aip, xml_aip;
BioseqPtr fake_bsp = NULL, query_bsp, bsp;
BioSourcePtr source;
BLAST_MatrixPtr matrix;
Int4Ptr PNTR txmatrix;
BLAST_OptionsBlkPtr options;
BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
BlastPruneSapStructPtr prune;
Boolean db_is_na, query_is_na, show_gi, believe_query=FALSE;
Boolean html = FALSE;
CharPtr params_buffer=NULL;
Int4 number_of_descriptions, number_of_alignments;
SeqAlignPtr  seqalign;
SeqAnnotPtr seqannot = NULL;
SeqEntryPtr sep;
TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
Uint1 align_type, align_view, err_ticket;
Uint4 align_options, print_options;
ValNodePtr mask_loc, mask_loc_start = NULL, vnp, next_mask_loc = NULL;
ValNodePtr other_returns, error_returns;
CharPtr blast_program, blast_database, blast_inputfile, blast_outputfile;
FILE *infp, *outfp;
Char buf[256] = { '\0' } ;
/* Mega BLAST related variables */
SeqAlignPtr sap, next_seqalign, PNTR seqalignp;
Int4 num_bsps, m_index;
SeqLocPtr last_mask, mask_slp, slp = NULL, tmp_slp;
Int2 ctr = 1;
Char prefix[2];
Boolean done = TRUE;
int (LIBCALLBACK *handle_results)(VoidPtr srch);       
Int4 from = 0, to = -1;
const char *dummystr;
Uint1 num_queries;         /*--KM for concatenated queries in blastn, tblastn */
Uint1 num_iters;
Uint1 sap_iter;
SeqAlignPtr curr_seqalign;
SeqAlignPtrArray sap_array;                /*--KM for separating seqaligns to test concat printing, temporary?*/
SeqAlignPtrArrayPtr sap_arr_ptr;
SeqAnnotPtr curr_seqannot;
SeqAnnotPtrArray seq_annot_arr;
Uint1 bsp_iter;
BspArray fake_bsp_arr;     /*--KM the array of fake_bsps for indiv. queries */
Boolean concat_done, nuc_concat;
QueriesPtr mult_queries = NULL;    /*--KM, AM: stores information related to
                                                query multipolexing, to put in search */
BioseqPtr curr_bsp;

/* AM: Support for query multiplexing. */
Uint4 num_spacers;

ValNodePtr sep_list;	/**< a list of all sep data structures */
ValNodePtr cur_sep = NULL;		/**< tracks the current position in sep_list */

#ifdef MPE
#include <mpi.h>
#include <mpe.h>
	int fakebioseq_start = 0;
	int fakebioseq_end = 0;
	int blastengine_start = 0;
	int blastengine_end = 0;
	int addaligninfo_start = 0;
	int addaligninfo_end = 0;
	int asnwrite_start = 0;
	int asnwrite_end = 0;
#endif

Int2 initBLAST ()
 
{

#ifdef MPE
	fakebioseq_start = MPE_Log_get_event_number(); 
	fakebioseq_end = MPE_Log_get_event_number(); 
	blastengine_start = MPE_Log_get_event_number(); 
	blastengine_end = MPE_Log_get_event_number(); 
	addaligninfo_start = MPE_Log_get_event_number(); 
	addaligninfo_end = MPE_Log_get_event_number(); 
	asnwrite_start = MPE_Log_get_event_number(); 
	asnwrite_end = MPE_Log_get_event_number(); 

	MPE_Describe_state(fakebioseq_start, fakebioseq_end, "GetFakeBioseq","purple");
	MPE_Describe_state(blastengine_start, blastengine_end, "blast_engine","yellow");
	MPE_Describe_state(addaligninfo_start, addaligninfo_end, "AddAlignInfo","blue");
	MPE_Describe_state(asnwrite_start, asnwrite_end, "AsnWrite","orange");
#endif


    StringCpy(buf, "blastall ");
    StringNCat(buf, BlastGetVersionNumber(), sizeof(buf)-StringLen(buf));
    if (! GetArgs (buf, NUMARG, myargs)) {
        return (1);
    }
    
    UseLocalAsnloadDataAndErrMsg ();
    
    if (! SeqEntryLoad())
        return 1;
    
    ErrSetMessageLevel(SEV_WARNING);
    
    blast_program = myargs[0].strvalue;

    blast_database = myargs[1].strvalue;
    blast_inputfile = myargs[2].strvalue;
    blast_outputfile = myargs[5].strvalue;

    if (myargs[29].intvalue)
        html = TRUE;
    
    if ((infp = FileOpen(blast_inputfile, "r")) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open input file %s\n", blast_inputfile);
        return (1);
    }
	
	align_view = (Int1) myargs[4].intvalue;
    outfp = NULL;
	if (align_view != 7 && align_view != 10 && align_view != 11 && blast_outputfile != NULL) {
        if ((outfp = FileOpen(blast_outputfile, "w")) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", blast_outputfile);
            return (1);
        }
    }
    
    if (StringCmp("filter", blast_program) == 0) {
        BlastGetMaskingLoc(infp, outfp, myargs[6].strvalue);
        FileClose(outfp);
        FileClose(infp);	
        return 0;
    }
    
    align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);

    if(align_view < 7) {
        if (StringICmp("blastx", blast_program) == 0) {
            if (align_view != 0) {
                ErrPostEx(SEV_FATAL, 0, 0, "This option is not available with blastx");
                return 1;
            }
        } else if (StringICmp("tblastx", blast_program) == 0) {
            if (align_view != 0) {
                ErrPostEx(SEV_FATAL, 0, 0, "This option is not available with tblastx");
                return 1;
            }
        }
    }
    
    believe_query = FALSE;
    if (myargs[21].intvalue != 0)
        believe_query = TRUE;

    if (believe_query == FALSE && (myargs[20].strvalue || align_view == 10 || align_view ==11)) {
        ErrPostEx(SEV_FATAL, 0, 0, "-J option must be TRUE to produce a SeqAlign file");
    }
    
    options = BLASTOptionNewEx(blast_program, (Boolean) myargs[16].intvalue, (Boolean) myargs[35].intvalue);
    if (options == NULL)
        return 3;
    
    if (align_view == 8 && options->is_megablast_search) {
       options->output = (VoidPtr) outfp;
       handle_results = MegaBlastPrintAlignInfo;
    } else 
       handle_results = NULL;

    BLASTOptionSetGapParams(options, myargs[22].strvalue, 0, 0); 
    options->expect_value  = (Nlm_FloatHi) myargs[3].floatvalue;
    number_of_descriptions = myargs[13].intvalue;	
    number_of_alignments = myargs[14].intvalue;	
    options->hitlist_size = MAX(number_of_descriptions, number_of_alignments);
	if (StringICmp("blastn", blast_program) == 0) {
		options->penalty = myargs[11].intvalue;
		options->reward = myargs[12].intvalue;
		if (options->reward > 1) {
			/* Scale the default values for gap costs; will be overridden
			later, if command line values are non-zero */
			options->gap_open *= options->reward;
			options->gap_extend *= options->reward;
		}
	} else {
		if (myargs[15].intvalue != 0) {
			options->threshold_second = myargs[15].intvalue;
		}
	}

    if (myargs[7].intvalue != 0)
        options->gap_open = myargs[7].intvalue;
    if (myargs[8].intvalue != 0)
        options->gap_extend = myargs[8].intvalue;
    if (myargs[9].intvalue != 0)
        options->gap_x_dropoff = myargs[9].intvalue;

	/* Multiple hits does not apply to blastn. The two-pass method is defunct */
    if (StringICmp("blastn", blast_program)) {
        if (myargs[26].intvalue == 1) {
            options->two_pass_method  = FALSE;
            options->multiple_hits_only  = FALSE;
        } else {
            /* all other inputs, including the default 0 use 2-hit method */
           options->two_pass_method  = FALSE;
           options->multiple_hits_only  = TRUE;
		}
    }
    else
    { /* Reverse these for blastn for now. */
        if (myargs[26].intvalue == 1) {
            options->two_pass_method  = FALSE;
            options->multiple_hits_only  = TRUE;
        } else {
            options->two_pass_method  = FALSE;
            options->multiple_hits_only  = FALSE;
	    }
    }

    if(myargs[33].intvalue != 0) 
        options->gap_x_dropoff_final = myargs[33].intvalue;

    if (StringICmp(myargs[6].strvalue, "T") == 0) {
        if (StringICmp("blastn", blast_program) == 0)
            options->filter_string = StringSave("D");
        else
            options->filter_string = StringSave("S");
    } else {
        options->filter_string = StringSave(myargs[6].strvalue);
    }
    
    show_gi = (Boolean) myargs[10].intvalue;
    
    options->genetic_code = myargs[17].intvalue;
    options->db_genetic_code = myargs[18].intvalue;
    options->number_of_cpus = myargs[19].intvalue;
    if (myargs[23].intvalue != 0) {
        options->wordsize = myargs[23].intvalue;
    }
    
    if (options->is_megablast_search) {
       options->cutoff_s2 = options->wordsize*options->reward;
       options->cutoff_s = (options->wordsize + 4)*options->reward;
    }

    options->db_length = (Int8) myargs[24].floatvalue;
    
    options->hsp_range_max  = myargs[25].intvalue;
    if (options->hsp_range_max != 0)
        options->perform_culling = TRUE;
    if (myargs[27].floatvalue)
        options->searchsp_eff = (Nlm_FloatHi) myargs[27].floatvalue;
    
    options->strand_option = myargs[28].intvalue;

    if(myargs[32].floatvalue != 0.0) {
        options->dropoff_2nd_pass  = (Nlm_Int4)myargs[32].floatvalue;
        if(options->dropoff_1st_pass > options->dropoff_2nd_pass)
            options->dropoff_1st_pass = options->dropoff_2nd_pass;
    }

    if (myargs[37].intvalue != 0)
        options->window_size = myargs[37].intvalue;

    print_options = 0;
    align_options = 0;
    align_options += TXALIGN_COMPRESS;
    align_options += TXALIGN_END_NUM;
    if (StringICmp("blastx", blast_program) == 0) {
        align_options += TXALIGN_BLASTX_SPECIAL;
    }
    if (show_gi) {
        align_options += TXALIGN_SHOW_GI;
        print_options += TXALIGN_SHOW_GI;
    }
    if (myargs[16].intvalue == 0)
        print_options += TXALIGN_SHOW_NO_OF_SEGS;
    
    if (align_view) {
        align_options += TXALIGN_MASTER;
        if (align_view == 1 || align_view == 3)
            align_options += TXALIGN_MISMATCH;
        if (align_view == 3 || align_view == 4 || align_view == 6)
            align_options += TXALIGN_FLAT_INS;
        if (align_view == 5 || align_view == 6)
            align_options += TXALIGN_BLUNT_END;
    } else {
        align_options += TXALIGN_MATRIX_VAL;
        align_options += TXALIGN_SHOW_QS;
    }
    
    if (html) {
        align_options += TXALIGN_HTML;
        print_options += TXALIGN_HTML;
    }

    if (myargs[30].strvalue) {
        options->gifile = StringSave(myargs[30].strvalue);
    }
    
    /* 
       Out-of-frame option is valid only for blastx, tblastn and 
       psitblastnsearches
    */

#ifdef DO_NOT_SUPPRESS_BLAST_OP
    if(myargs[38].intvalue > 0) {
        if (!StringICmp("blastx", blast_program) || 
            !StringICmp("tblastn", blast_program)||
	    !StringICmp("psitblastn", blast_program)) {
           if (!StringICmp("blastx", blast_program)) {
              options->is_ooframe = TRUE;
              options->shift_pen = myargs[38].intvalue;
           }
        }
    }
#endif
        
#ifdef DO_NOT_SUPPRESS_BLAST_OP
    /* Input longest intron length is in nucleotide scale; in the lower level
       code it will be used in protein scale */
    if (myargs[39].intvalue > 0) 
       options->longest_intron = MAX(myargs[39].intvalue, MAX_INTRON_LENGTH);
#endif

    aip = NULL;
    if (myargs[20].strvalue != NULL) {
    		
		aip = AsnIoOpen (myargs[20].strvalue,"w");
	
        if (aip == NULL) {
                ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", myargs[20].strvalue);
                return 1;
        }
    }
	else if (align_view == 10 || align_view == 11)
	{
		const char* mode = (align_view == 10) ? "w" : "wb";
		
		if ((aip = AsnIoOpen (blast_outputfile, (char*) mode)) == NULL) {
			ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", myargs[20].strvalue);
			return 1;
		}
	}


                  /* Futamura: Setting up the psitblastn options */
    if (NULL != myargs[34].strvalue) {
          options->recoverCheckpoint = TRUE;
          options->freqCheckpoint = TRUE;
    }
    options->CheckpointFileName=myargs[34].strvalue;

	/*--KM get number of queries for concatenated blastn/tblastn queries */
	options->NumQueries=myargs[40].intvalue;
	num_queries = options->NumQueries;
	if (num_queries>0 &&
		!( (StringICmp("blastn",  blast_program) == 0) ||
		(StringICmp("tblastn", blast_program) == 0)   ) ) {
		
		ErrPostEx(SEV_FATAL, 0, 0, "blast: Can't concat with program %s\n", myargs[0].strvalue);
		return 1;
	}
	
	/* AM: Query concatenation is not consistent with ungapped search */
	if( num_queries > 0 && !myargs[16].intvalue )
	{
		ErrPostEx( SEV_FATAL, 0, 0,
			"blast: Query concatenation is inconsistent with ungapped search\n" );
		return 1;
	}
	
	/* --KM set bool value if DNA and concat needed, need for Fasta->seq functions */
	if (num_queries>0 && query_is_na == TRUE) {
		nuc_concat = TRUE;
	} else {
		nuc_concat = FALSE;
	}

	concat_done = FALSE;       /*--KM */

	return 0;
}

/**
 * function to parse the FastA format query file, load the queries, and count the number loaded
 */
int loadQueries( void ) 
{
	int query_count = 0;
	while (TRUE) {
		if (!options->is_megablast_search) {
			if(myargs[31].intvalue) {
				sep = FastaToSeqEntryForDb (infp, query_is_na, NULL, 0/*believe_query*/, 
					NULL, NULL, &options->query_lcase_mask);
			} else {
				sep = FastaToSeqEntryEx(infp, query_is_na, NULL, 0/*believe_query*/ );
			}
			
			if(sep == NULL){
				break; /* no more queries, can go to finish with next break */
          		}
			
	  		ValNodeAddPointer( &sep_list, 0, sep );
	  		query_count++;
       		}
	}
	return query_count;
}

void cleanupQueries( void ) {
	ValNodePtr cur_sep = sep_list;
	for( ; cur_sep != NULL; cur_sep = cur_sep->next )
		cur_sep->data.ptrvalue = SeqEntryFree( (SeqEntryPtr)cur_sep->data.ptrvalue );
	sep_list = ValNodeFree( sep_list );
}

/**
 * function to load values into fake_bsp for the current bioseq
 * returns 0 on success
 * returns > 0 on error condition
 * returns -1 if there are no more bioseqs to load
 */
int getFakeBioseq( Boolean preloaded ){
		if (options->is_megablast_search) {
			StrCpy(prefix, "");
			slp = NULL;
			num_bsps = 0;
			done = TRUE;
			SeqMgrHoldIndexing(TRUE);
			mask_slp = last_mask = NULL;
			while ((sep=FastaToSeqEntryForDb(infp, query_is_na, NULL,
					   believe_query, prefix, &ctr, 
					   &mask_slp)) != NULL) {
				if ((Boolean)myargs[31].intvalue) {
					if (mask_slp) {
						if (!last_mask)
							options->query_lcase_mask = last_mask = mask_slp;
						else {
							last_mask->next = mask_slp;
							last_mask = last_mask->next;
						}
						mask_slp = NULL;
					}
				} else {
					mask_slp = SeqLocSetFree(mask_slp);
				}
				query_bsp = NULL;
				if (query_is_na) 
					SeqEntryExplore(sep, &query_bsp, FindNuc);
				else
					SeqEntryExplore(sep, &query_bsp, FindProt);
	     
				if (query_bsp == NULL) {
					ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
					return 2;
				}
	     
				/* Only for the first query */
				if (num_bsps == 0) {
					to = MIN(to, query_bsp->length - 1);
                 
					/* -1 means end of sequence */
					if (to < 0)
						to = query_bsp->length - 1;
					if (from >= query_bsp->length || to < 0) {
						ErrPostEx(SEV_FATAL, 0, 0, 
                               "Location outside of the query sequence range\n");
						return 3;
					}
                 slp = SeqLocIntNew(from, to, options->strand_option, 
                                    SeqIdFindBest(query_bsp->id, SEQID_GI));
				} else 
					ValNodeAddPointer(&slp, SEQLOC_WHOLE,
                                   SeqIdDup(SeqIdFindBest(query_bsp->id,
                                                          SEQID_GI)));
				num_bsps++;
				if (num_bsps >= MAX_NUM_QUERIES) {
					done = FALSE;
					break;
				}
				/*sep = MemFree(sep);*/ /* Do not free the underlying Bioseq */
			}
			SeqMgrHoldIndexing(FALSE);
			if (num_bsps == 0) 
				return -1;
		} else {

			/* not megablast */
			
			/*--KM make array of fake_bsp's if concat. query */
			if (concat_done)
			   return -1;
			if (num_queries > 0)  {
			   fake_bsp_arr = (BspArray) MemNew(sizeof(BioseqPtr)*num_queries);
			}
			num_iters = (num_queries>0) ? num_queries : 1;
			for (bsp_iter=0; bsp_iter<num_iters; bsp_iter++) {

				if( preloaded ){
					if(cur_sep == NULL){
						sep = NULL;
					}else{
						sep = (SeqEntryPtr)cur_sep->data.ptrvalue;
						cur_sep = cur_sep->next;
					}
				}else{
					if(myargs[31].intvalue) {	  
						sep = FastaToSeqEntryForDb (infp, query_is_na, NULL, believe_query, 
							NULL, NULL, &options->query_lcase_mask);
					} else {	  
						sep = FastaToSeqEntryEx(infp, query_is_na, NULL, believe_query);
					}
				}

				/* if concat and num_queries has not been reached and sep is NULL, crap out */
				if (sep == NULL && bsp_iter < num_queries) {   /* implies num_queries>0 */
					ErrPostEx(SEV_FATAL, 0, 0, "blast: Only %n queries found!\n", bsp_iter);
					return (1);
				}
				if(sep == NULL)
					return -1;	/* go to finish, all bioseqs have been parsed */

				query_bsp = NULL;
				if (query_is_na) {
					SeqEntryExplore(sep, &query_bsp, FindNuc);
				} else {
					SeqEntryExplore(sep, &query_bsp, FindProt);
				}
   	       
				if (query_bsp == NULL) {
					ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
					return 2;
				}
	
			   if (num_queries>0) {
			      *(fake_bsp_arr + bsp_iter) = query_bsp;
			   }
			}
			if ( (sep == NULL && num_queries ==0) || (num_queries>0 && concat_done) )
			   return -1;  /* go to finish */

			/* --KM */
			if (num_queries>0) {
				concat_done = TRUE;   /* --KM to prevent futher looping */

				/* AM: Determine the number of query separators. */
				num_spacers = GetNumSpacers( options, believe_query, fake_bsp_arr );

				if( num_spacers%2 ) ++num_spacers;

				/* --KM make the concatenated fake_bsp */
				/* AM: Added num_spacers. */
				if( query_is_na )
					fake_bsp = (BioseqPtr)
			                BlastMakeFakeBspConcat(fake_bsp_arr, num_queries, query_is_na, num_spacers);
				else
					fake_bsp = (BioseqPtr)
			                BlastMakeFakeBspConcat(fake_bsp_arr, num_queries, query_is_na, num_spacers);

				/* construct the MultQueries struct here*/
				mult_queries = (QueriesPtr) BlastMakeMultQueries(fake_bsp_arr, num_queries, query_is_na, num_spacers);
			} else {

				if(believe_query){
					fake_bsp = query_bsp;
				}
				else{ 
					fake_bsp = BlastMakeFakeBioseq(query_bsp, NULL);
				}
	  		}
		  err_ticket = BlastSetUserErrorString(NULL, query_bsp->id, believe_query);
        
          /* If fake_bsp created mask should be updated to use it's id */
          BLASTUpdateSeqIdInSeqInt(options->query_lcase_mask, fake_bsp->id);
        
          source = BioSourceNew();
          source->org = OrgRefNew();
          source->org->orgname = OrgNameNew();
          source->org->orgname->gcode = options->genetic_code;
          ValNodeAddPointer(&(query_bsp->descr), Seq_descr_source, source);
       }
	   return 0;
}

Int2 runBLAST( ValNodePtr* result_id_list ){
    int queryI = -1;
	ValNodePtr new_node;
	cur_sep = sep_list;
    /* --- Main loop over all FASTA entries in the input file ---- */

    if (myargs[36].strvalue) {       
        CharPtr delimiters = " ,;";
        CharPtr location;
        location = myargs[36].strvalue;
        from = atoi(StringTokMT(location, delimiters, &location)) - 1;
        to = atoi(location) - 1;
        from = MAX(from, 0);
    }

    while (TRUE) {
		int fake_bsp_rval = 0;
		queryI++;
#ifdef MPE
			MPE_Log_event(fakebioseq_start,0,"get fakebioseq start");
#endif
		fake_bsp_rval = getFakeBioseq( TRUE );
#ifdef MPE
			MPE_Log_event(fakebioseq_end,0,"get fakebioseq end");
#endif
		if( fake_bsp_rval > 0 )
			return fake_bsp_rval;
		else if( fake_bsp_rval == -1 )
			break;

       global_fp = outfp;
       
        other_returns = NULL;
        error_returns = NULL;

        if (options->is_megablast_search) {
           seqalignp = BioseqMegaBlastEngineByLoc(slp, blast_program,
			           blast_database, options, &other_returns, 
                                   &error_returns, 
                                   align_view < 7 ? tick_callback : NULL,
                                   NULL, NULL, 0, handle_results);
           seqalign = NULL;
	   for (m_index=0; m_index<num_bsps; m_index++) { 
	      if (seqalignp[m_index]) {
                 if (seqalign == NULL) 
                    sap = seqalign = seqalignp[m_index];
                 else
                    sap->next = seqalignp[m_index];
                 while (sap->next != NULL)
                    sap = sap->next;
              }
           }
           seqalignp = (SeqAlign**)MemFree(seqalignp);
        } else if (!myargs[36].strvalue) {       
           /* KM added mult_queries param */
#ifdef MPE
			MPE_Log_event(blastengine_start,0,"blastengine start");
#endif
           seqalign = BioseqBlastEngineWithCallbackMult(fake_bsp, blast_program, blast_database, options, &other_returns, &error_returns, align_view < 7 ? tick_callback : NULL, handle_results, mult_queries);
#ifdef MPE
			MPE_Log_event(blastengine_end,0,"blastengine end");
#endif
        } else { /* Location on query provided */
           to = MIN(to, fake_bsp->length - 1);
           
           /* -1 means end of sequence */
           if (to < 0)
              to = fake_bsp->length - 1;
           if (from >= fake_bsp->length || to < 0) {
              ErrPostEx(SEV_FATAL, 0, 0, 
                        "Location outside of the query sequence range\n");
              return 3;
           }
           slp = SeqLocIntNew(from, to, options->strand_option, 
                              fake_bsp->id);
           seqalign = BioseqBlastEngineByLocWithCallbackMult(slp, blast_program, blast_database, options, &other_returns, &error_returns, align_view < 7 ? tick_callback : NULL, NULL, NULL, 0, handle_results, mult_queries);
           
        }
/** Added by AED */

		if( seqalign != NULL ){
			new_node = ValNodeAdd( result_id_list );
			new_node->choice = 0;
			new_node->data.intvalue = queryI;

			seqannot = SeqAnnotNew();
			seqannot->type = 2;
#ifdef MPE
			MPE_Log_event(addaligninfo_start,0,"addaligninfo start");
#endif
			AddAlignInfoToSeqAnnot(seqannot, align_type);
#ifdef MPE
			MPE_Log_event(addaligninfo_end,0,"addaligninfo end");
#endif
			seqannot->data = seqalign;
#ifdef MPE
			MPE_Log_event(asnwrite_start,0,"asnwrite start");
#endif
			if (aip) {
				SeqAnnotAsnWrite((SeqAnnotPtr) seqannot, aip, NULL);
				AsnIoReset(aip);
			}
#ifdef MPE
			MPE_Log_event(asnwrite_end,0,"asnwrite end");
#endif
		}
		seqannot = SeqAnnotFree(seqannot);

/** End added by AED */

        BlastErrorPrint(error_returns);

        dbinfo = NULL;
        ka_params = NULL;
        ka_params_gap = NULL;
        params_buffer = NULL;
        mask_loc = NULL;
        matrix = NULL;
        txmatrix = NULL;
        for (vnp=other_returns; vnp; vnp = vnp->next) {
            switch (vnp->choice) {
            case TXDBINFO:
                //dbinfo = (TxDfDbInfoPtr)(vnp->data.ptrvalue);
                break;
            case TXKABLK_NOGAP:
                ka_params = (BLAST_KarlinBlk *)vnp->data.ptrvalue;
                break;
            case TXKABLK_GAP:
                ka_params_gap = (BLAST_KarlinBlk *)vnp->data.ptrvalue;
                break;
            case TXPARAMETERS:
                params_buffer = (char*)vnp->data.ptrvalue;
                break;
            case TXMATRIX:
                matrix = (BLAST_MatrixPtr)(vnp->data.ptrvalue);
                if (matrix)
                   txmatrix = BlastMatrixToTxMatrix(matrix);
                break;
            case SEQLOC_MASKING_NOTSET:
            case SEQLOC_MASKING_PLUS1:
            case SEQLOC_MASKING_PLUS2:
            case SEQLOC_MASKING_PLUS3:
            case SEQLOC_MASKING_MINUS1:
            case SEQLOC_MASKING_MINUS2:
            case SEQLOC_MASKING_MINUS3:
                ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
                break;
            default:
                break;
            }
        }	

/** start 'cleanup the mess' code */

        slp = SeqLocSetFree(slp);
        matrix = BLAST_MatrixDestruct(matrix);
        if (txmatrix)
           txmatrix = TxMatrixDestruct(txmatrix);

        init_buff_ex(85);
        dbinfo_head = dbinfo;
        
        dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
        
        if (ka_params) {
           MemFree(ka_params);
        }
        
        if (ka_params_gap) {
           MemFree(ka_params_gap);
        }

        MemFree(params_buffer);
        free_buff();
        mask_loc = mask_loc_start;
        while (mask_loc) {
           SeqLocSetFree((ValNode*)(mask_loc->data.ptrvalue));
           mask_loc = mask_loc->next;
        }
        ValNodeFree(mask_loc_start);
        
        if(!believe_query)
           fake_bsp = BlastDeleteFakeBioseq(fake_bsp);
        
        other_returns = ValNodeFree(other_returns);

#ifndef BLAST_CS_API
        /* This is freed earlier in client-server case */
        options->query_lcase_mask = SeqLocSetFree(options->query_lcase_mask);
        ReadDBBioseqFetchDisable();
#endif
        
        if (!options->is_megablast_search) 
           BlastDeleteUserErrorString(err_ticket);

/** end cleanup code */

    } /* while(TRUE)  - main loop of the program over all FASTA entries */
	return 0;
}

/*
 * Printing out html header here 
 */
void outputHtmlHeader(){
    if(align_view < 7) {
       if (html) {
          fprintf(outfp, "<HTML>\n<TITLE>BLAST Search Results</TITLE>\n");
          fprintf(outfp, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
                  "VLINK=\"#660099\" ALINK=\"#660099\">\n");
          fprintf(outfp, "<PRE>\n");
       }
    } else if (align_view == 7) { 
       xml_aip = AsnIoOpen(blast_outputfile, "wx");
	}
	if(align_view >= 7 && myargs[40].intvalue > 1)
	{
		ErrPostEx( SEV_FATAL, 0, 0,
			"blast: Query concatenation is currently not supported with -m > 7");
//		return 1;
	}
}

/*
 * This function is a modified version of :
 * 	static SeqIdPtr SeqIdFree(SeqIdPtr sip)
 * from the ncbi toolkit (see objloc.c). It is responsible for 
 * freeing all of the memory associated with the first data.ptrvalue variable.
 */
void SeqIdDataFree(SeqIdPtr sip)
{
        switch(sip->choice) {
         	case SEQID_LOCAL:	/* local */
            		sip->data.ptrvalue = 
				ObjectIdFree(sip->data.ptrvalue);
            		break;
         	case SEQID_GIBBSQ:      /* gibbseq */
         	case SEQID_GIBBMT:      /* gibbmt */
            		break;
         	case SEQID_GIIM:      	/* giimid */
            		GiimFree(sip->data.ptrvalue);
            		break;
        	case SEQID_GENBANK:	/* genbank */
        	case SEQID_EMBL:      	/* embl */
        	case SEQID_PIR:		/* pir   */
        	case SEQID_SWISSPROT:	/* swissprot */
        	case SEQID_OTHER:	/* other */
        	case SEQID_DDBJ:
                case SEQID_TPG:		/* Third Party Annot/Seq Genbank */
             	case SEQID_TPE:		/* Third Party Annot/Seq EMBL */
             	case SEQID_TPD:		/* Third Party Annot/Seq DDBJ */
         	case SEQID_PRF:
            		sip->data.ptrvalue = 
				TextSeqIdFree(sip->data.ptrvalue);
            		break;
         	case SEQID_PATENT:      /* patent seq id */
            		sip->data.ptrvalue = 
				PatentSeqIdFree(sip->data.ptrvalue);
            		break;
         	case SEQID_GENERAL:     /* general */
            		sip->data.ptrvalue = DbtagFree(sip->data.ptrvalue);
            		break;
         	case SEQID_GI:		/* gi */
            		break;
         	case SEQID_PDB:
         		sip->data.ptrvalue = PDBSeqIdFree(sip->data.ptrvalue);
                	break;
        };
}
/**
 * set the data inside of id to the data inside of new_id, deleting
 * space allocated for id->data
 */
void setQueryId( SeqIdPtr id, SeqIdPtr new_id, Boolean free_data ){
	if(new_id == NULL){
		id->data.ptrvalue = NULL;
	}
	else{
//		fprintf( stderr, "Changing id choice %d ext %d to choice %d ext %d\n", id->choice, id->extended , new_id->choice, new_id->extended );
//		fprintf( stderr, "old name: %s, new name %s\n", ((ObjectIdPtr)(id->data.ptrvalue))->str, ((ObjectIdPtr)(new_id->data.ptrvalue))->str );
		/* Delete the original contents */
//		if( free_data )
//			SeqIdDataFree(id);

		id->choice = new_id->choice;
		id->extended = new_id->extended;
		id->data = new_id->data;
	}
}

/**
 * Change the query SeqIdPtrs inside the SeqAlign to store the SeqId
 * contained by new_id
 */
void setSeqAlignIds( SeqAlignPtr sappy, SeqIdPtr new_id )
{
	SeqAlignPtr sap_iter;
	SeqIdPtr tmp, sip;
    SeqLocPtr slp;
	
	for( sap_iter = sappy; sap_iter != NULL; sap_iter = sap_iter->next ){
		tmp = NULL;
		
		switch( sap_iter->segtype ){
			case SAS_DENDIAG:
				tmp = ((DenseDiagPtr)sap_iter->segs)->id;
			break;
			case SAS_DENSEG:
				tmp = ((DenseSegPtr)sap_iter->segs)->ids;
			break;
			case SAS_STD:
				tmp = ((StdSegPtr)sap_iter->segs)->ids;
        	    for(slp = ((StdSegPtr)sap_iter->segs)->loc; slp != NULL; slp = slp->next) {
        	        sip = SeqLocId(slp);
					if( SeqIdComp( sip, tmp ) == 0 )
						setQueryId( sip, new_id, FALSE );
	            }
			break;
			case SAS_PACKED:
				tmp = ((PackSegPtr)sap_iter->segs)->ids;
			break;
			case SAS_DISC:
			break;
			case SAS_COMPSEQ:
			break;
			default:
				fprintf(stderr, "Unknown seq type!\n");
				exit(1);
			break;
		}
		if(tmp != NULL){
			setQueryId( tmp, new_id, TRUE );
		}
	}
}

void MpiBlastPrintReference( int html, Int4 line_length, FILE* outfp ) {
	if (outfp == NULL)
		return;
	
	if (html) {
		fprintf( outfp, "<b><a href=\"http://mpiblast.lanl.gov\">Reference</a>:</b>" );
		fprintf( outfp, "Aaron E. Darling, Lucas Carey, and Wu-chun Feng, " );
		fprintf( outfp, "\"The design, implementation, and evaluation of mpiBLAST.\"" );
		fprintf( outfp, "In Proceedings of <a href=\"http://clusterworld.com\">ClusterWorld 2003</a>, June 24-26 2003, San Jose, CA" );
	} else {
		fprintf( outfp, "Reference: Aaron E. Darling, Lucas Carey, and Wu-chun Feng,\n" );
		fprintf( outfp, "\"The design, implementation, and evaluation of mpiBLAST.\"\n" );
		fprintf( outfp, "In Proceedings of ClusterWorld 2003, June 24-26 2003, San Jose, CA\n\n" );
	}
}

/**
 * output is done here
 */
Int2 outputResults( SeqAlignPtr sappy )
{
	int fake_bsp_rval = 0;
	seqalign = sappy;
	global_fp = outfp;
	// NOTE: if this function is called more times than the number of queries
	// then it will re-cycle through the query list!
	if( cur_sep == NULL )
		cur_sep = sep_list;
	
	/* bioseqs are now pre-loaded here */
	fake_bsp_rval = getFakeBioseq( TRUE );
	
	if( fake_bsp_rval > 0 )
		return fake_bsp_rval;
	else if( fake_bsp_rval == -1 )
		return 0;

	
       if(align_view < 7) {
           init_buff_ex(90);
           BlastPrintVersionInfo(blast_program, html, outfp);
           fprintf(outfp, "\n");
           MpiBlastPrintReference(html, 90, outfp);
           fprintf(outfp, "\n");
           if (!options->is_megablast_search) {
              /* KM added loop here for concat case */
              num_iters = (num_queries>0) ? num_queries : 1;
              for (bsp_iter=0; bsp_iter<num_iters; bsp_iter++) {
                 curr_bsp = (num_queries>0) ? *(fake_bsp_arr + bsp_iter) : query_bsp;
                 AcknowledgeBlastQuery(curr_bsp, 70, outfp, believe_query, html);
              }
           }

            /* Here we first check, that database do no exists */

           if(!PrintDbInformation(blast_database, !db_is_na, 70, outfp, html))
                return 1;
            free_buff();
		if (options->is_ooframe)
        		ErrPostEx(SEV_WARNING, 0, 0, "Out-of-frame option selected, Expect values are only approximate and calculated not assuming out-of-frame alignments");
        }
#ifdef OS_UNIX
/*        if(align_view < 7) {
            fprintf(global_fp, "%s", "Searching");
        }
*/
#endif		
	if( sappy != NULL ){	  
//		setSeqAlignIds( sappy, fake_bsp->id );
	}
	
#ifndef BLAST_CS_API
// The function call to ReadDBBioseqFetchEnable has been moved into 
// mpiblast.cpp (in the outer loop).
//ReadDBBioseqFetchEnable ("blastall", blast_database, db_is_na, TRUE);
#endif

        ReadDBBioseqSetDbGeneticCode(options->db_genetic_code);

        tmp_slp = slp;
        if (slp)
           query_bsp = NULL;

        if (getenv("POST_BLAST_CLUSTER_HITS") != NULL)
           BlastClusterHitsFromSeqAlign(seqalign, blast_program, blast_database, 
                                        options, 0.9, 1.6, 0.5, TRUE);

        if (mask_loc) {
           mask_loc_start = mask_loc;
        }	
		else
	{	/* Could have become non-NUll for last query. */
           mask_loc_start = NULL;
	}
	  
        /* Print header in any case */
        if (align_view == 9) {
              	PrintTabularOutputHeader(blast_database, query_bsp, slp, 
                                       blast_program, 0, believe_query,
                                       global_fp);
		}

        if (seqalign) {
          if (num_queries > 0) { /* AM: Support for query multiplexing. */
             sap_array = mult_queries->sap_array_data->sap_array;
          }

          if (align_view == 8 || align_view == 9) {
/* --KM need to put a loop around this. seqaligns already broken up
   note the method for looping if num_aligns > 0 - reuse this method everywhere */
             num_iters = (num_queries>0) ? num_queries : 1;
             for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
                curr_seqalign = (num_queries>0) ? *(sap_array + sap_iter) : seqalign;
                BlastPrintTabulatedResults(curr_seqalign, query_bsp, slp,
                                         number_of_alignments,
                                         blast_program, 
                                         !options->gapped_calculation,
                                         believe_query, from, 0, global_fp,
                                         (align_view == 9));
					 
/*				seqalign = NULL;
*/              SeqAlignSetFree(curr_seqalign);
             }
           } else {
           while (seqalign) {
	   
              if (!options->is_megablast_search){
                 next_seqalign = NULL;
              } else {
                 SeqIdPtr sip, next_sip = NULL;
                 
                 sap = seqalign;
                 sip = TxGetQueryIdFromSeqAlign(seqalign);
		 
                 while (sap != NULL) { 
                    if (sap->next != NULL) {
                       next_sip = TxGetQueryIdFromSeqAlign(sap->next);

                       if (SeqIdComp(sip, next_sip) != SIC_YES) {
                          next_seqalign = sap->next;
                          sap->next = NULL;
                       }
                    } else{
                       next_seqalign = NULL;
		    }
                    sap = sap->next;
                 }
                 
                 while (tmp_slp && SeqIdComp(sip, SeqLocId(tmp_slp)) != SIC_YES)
                    tmp_slp = tmp_slp->next;
                 if (tmp_slp == NULL) /* Should never happen */
                    break;
                 /* Separate the mask locations list for this query */
                 if (!mask_loc && next_mask_loc) {
                    mask_loc = next_mask_loc;
                    next_mask_loc = NULL;
                 }
                 if (mask_loc) {
                    if (next_mask_loc) {
                       mask_loc->next = next_mask_loc;
                       mask_loc = next_mask_loc;
                    }
                    mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
                    next_mask_loc = mask_loc;
                    while (SeqIdComp(SeqLocId(mask_slp), sip) != SIC_YES) {
                       mask_loc = mask_loc->next;
                       if (!mask_loc)
                          break;
                       mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
                    }
                    if (mask_loc) {
                       next_mask_loc = mask_loc->next;
                       mask_loc->next = NULL;
                    }
                 }
                 
                 bsp = BioseqLockById(SeqLocId(tmp_slp));
                 init_buff_ex(85);
                 fprintf(outfp, "\n");
                 AcknowledgeBlastQuery(bsp, 70, outfp, believe_query, html);
                 free_buff();
                 BioseqUnlock(bsp);
              }
	      
              if(align_view == 7 && !options->is_ooframe) {
                 if (options->is_megablast_search) {
                    bsp = BioseqLockById(SeqLocId(tmp_slp));
                    BXMLPrintOutput(xml_aip, seqalign, 
                                    options, blast_program, blast_database, 
                                    bsp, other_returns, 0, NULL);
                    BioseqUnlock(bsp);
                 AsnIoReset(xml_aip);
                 seqalign = SeqAlignSetFree(seqalign);

                 } else {
                    num_iters = (num_queries>0) ? num_queries : 1;
                    for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
                       curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
                       BXMLPrintOutput(xml_aip, curr_seqalign,
                                    options, blast_program, blast_database, 
                                    fake_bsp, other_returns, 0, NULL);
		 
                       AsnIoReset(xml_aip);
                     curr_seqalign = SeqAlignSetFree(curr_seqalign);
					 if( num_queries == 0 )
					 	seqalign = curr_seqalign;

                     }  /* for loop over sap-array (concat) */
			     }
              } else {
		
                /* create the array of SeqAnnotPtrs, if necessary */

                num_iters = (num_queries > 0) ? num_queries : 1;
                 for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
                    curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
                    if ( (num_queries > 0) && (sap_iter == 0) ) {
                       seq_annot_arr = (SeqAnnotPtrArray) MemNew(sizeof(SeqAnnotPtr)*num_queries);
                    }
                 seqannot = SeqAnnotNew();
                 seqannot->type = 2;
                 AddAlignInfoToSeqAnnot(seqannot, align_type);
                 seqannot->data = curr_seqalign;
		 
                 if (num_queries > 0) {
                 	*(seq_annot_arr + sap_iter) = seqannot;
                 }
                 } /* make seqannots over the sap_iters from concat, or the single seqalign */

		 		
                 if (outfp) { /* Uncacheing causes problems with ordinal nos. vs. gi's. */
                    ObjMgrSetHold();
                    /* print deflines */
		    
                    for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
                       curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
			
                       init_buff_ex(85);
			
                       PrintDefLinesFromSeqAlignEx2(curr_seqalign, 80, outfp,
                                       print_options, FIRST_PASS, NULL,
                                       number_of_descriptions, NULL, NULL);
			
                       free_buff();
                    } /* print deflines, looped if concat */
 
                    for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
                       curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
                       curr_seqannot = (num_queries > 0) ? *(seq_annot_arr + sap_iter) : seqannot;
			
                       prune = BlastPruneHitsFromSeqAlign(curr_seqalign,
                                       number_of_alignments, NULL);
				       
                       curr_seqannot->data = prune->sap;
			
                       if(options->is_ooframe) {
                           OOFShowBlastAlignment(curr_seqalign, /*mask*/ NULL,
                                              outfp, align_options, txmatrix);
                       } else {
                           if (align_view != 0)
                             ShowTextAlignFromAnnot(curr_seqannot, 60, outfp, NULL, NULL,
                                               align_options, txmatrix, mask_loc, NULL);
                           else
                             ShowTextAlignFromAnnot(curr_seqannot, 60, outfp, NULL, NULL,
                                       align_options, txmatrix, mask_loc,
                                               FormatScoreFunc);
                       }
                       
			if(aip == NULL){
				curr_seqannot->data = NULL;  
			}
			else{
			        /* 
		        	 * Restore the SeqAlign data that was temporarily replaced
				 * by the prune->sap data. The SeqAlign data will be needed
				 * since the user is writing a ASN.1 of SeqAnnot data.
		        	 */
				curr_seqannot->data = curr_seqalign;  
			}
			
		       prune = BlastPruneSapStructDestruct(prune);
                    } /* show text align, loop over seqalign/seqannots for concat */
		    
                    ObjMgrClearHold();
                    
                    ObjMgrFreeCache(0);

                 } /* if outfp */
		 
		
		 if (aip) {
		 	/* Substitue the "fake" bioseq ID for the query bioseq ID.
			 * This insures that query file SeqId's end up in the ASN.1
			 * output file (as opposed to the local SeqId's that are used
			 * internally -- i.e. loadQueries always sets belive_query
			 * to false!). 
			 */		
			for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
                       		curr_seqannot = 
					(num_queries > 0) ? *(seq_annot_arr + sap_iter) : seqannot;
				
		   		setSeqAlignIds( curr_seqannot->data, fake_bsp->id );
                   		SeqAnnotAsnWrite((SeqAnnotPtr) curr_seqannot, aip, NULL);
				setSeqAlignIds( curr_seqannot->data, NULL );
				AsnIoReset(aip);
                 	}
		 } /* if aip */
			
                 for (sap_iter=0; sap_iter < num_queries; sap_iter++) {
                    /* upper bound is num_queries, take care not to do this unless concat */
                    *(seq_annot_arr + sap_iter) = SeqAnnotFree(*(seq_annot_arr + sap_iter));
                 }
       /*--KM free seqalign array and all seqaligns?? */

              } /* end of else (not XML Printing) */
		
              if (options->is_megablast_search)
                 tmp_slp = tmp_slp->next;
		 
        /* --KM watch for memory leaks */
              if (seqannot && num_queries == 0)
                 seqannot = SeqAnnotFree(seqannot);
		 
		if( aip == NULL ){
		      seqalign = SeqAlignSetFree(seqalign);
	      	}
		
              seqalign = next_seqalign;
           } /* End of loop on all seqaligns */
           if (mask_loc && next_mask_loc)
              mask_loc->next = next_mask_loc;

           } /* end of align_view not tabular case */
        } else {         /* seqalign is NULL */
           if(align_view == 7 && !options->is_ooframe) {
              BlastErrorMsgPtr error_msg;
              CharPtr message;
              
              if (error_returns == NULL) {
                 message = "No hits found";
              } else {
                 error_msg = error_returns->data.ptrvalue;
                 message = error_msg->msg;
              }
              if (options->is_megablast_search) {
                 bsp = BioseqLockById(SeqLocId(tmp_slp));
                 BXMLPrintOutput(xml_aip, seqalign, 
                                 options, blast_program, blast_database, 
                                 bsp, other_returns, 0, NULL);
                 BioseqUnlock(bsp);
              } else {
                 BXMLPrintOutput(xml_aip, NULL, 
                                 options, blast_program, blast_database, 
                                 fake_bsp, other_returns, 0, message);
              }
              AsnIoReset(xml_aip);
           } else if (align_view < 8) {
              fprintf(outfp, "\n\n ***** No hits found ******\n\n");
           }
           if (error_returns != NULL) {
              for (vnp = error_returns; vnp; vnp = vnp->next) {
                 BlastDestroyErrorMessage((BlastErrorMsgPtr)vnp->data.ptrvalue);
              }
              ValNodeFree(error_returns);
           }
        }
        
        slp = SeqLocSetFree(slp);
        matrix = BLAST_MatrixDestruct(matrix);
        if (txmatrix)
           txmatrix = TxMatrixDestruct(txmatrix);
        
        if(html) {
           fprintf(outfp, "<PRE>\n");
        }
        
        init_buff_ex(85);
        dbinfo_head = dbinfo;
        
        if(align_view < 7 && done) {
           while (dbinfo) {
              PrintDbReport(dbinfo, 70, outfp);
              dbinfo = dbinfo->next;
           }
        }
        dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
        
        if (ka_params) {
           if(align_view < 7 && done) {
              PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, outfp, FALSE);
           }
           MemFree(ka_params);
        }
        
        if (ka_params_gap) {
           if(align_view < 7 && done) {
              PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, outfp, TRUE);
           }
           MemFree(ka_params_gap);
        }
        
        if(align_view < 7 && done) {
           PrintTildeSepLines(params_buffer, 70, outfp);
        }
        
        MemFree(params_buffer);
        free_buff();
        mask_loc = mask_loc_start;
        while (mask_loc) {
           SeqLocSetFree((ValNode*)mask_loc->data.ptrvalue);
           mask_loc = mask_loc->next;
        }
        ValNodeFree(mask_loc_start);
        
        if(!believe_query)
           fake_bsp = BlastDeleteFakeBioseq(fake_bsp);
        
        other_returns = ValNodeFree(other_returns);
        if (done) 
           sep = SeqEntryFree(sep);
#ifndef BLAST_CS_API
        /* This is freed earlier in client-server case */
        options->query_lcase_mask = SeqLocSetFree(options->query_lcase_mask);
	
	// The function call to ReadDBBioseqFetchDisable has been moved into 
	// mpiblast.cpp (in the outer loop).
        //ReadDBBioseqFetchDisable();
#endif
        if (html)
           fprintf(outfp, "</PRE>\n<P><HR><BR>\n<PRE>");
        
        if (!options->is_megablast_search) 
           BlastDeleteUserErrorString(err_ticket);
	return 0;
}

void outputHTMLfooter(){
    if(align_view < 7) {
        if (html) {
            fprintf(outfp, "</PRE>\n</BODY>\n</HTML>\n");
        }
	}
}

void cleanupBLAST(){
    aip = AsnIoClose(aip);
    
    if (align_view == 7)
        xml_aip = AsnIoClose(xml_aip);
    
    options = BLASTOptionDelete(options);

    // Add by H.C. for flushing outputs
    FileClose(outfp);

    FileClose(infp);
    
}

void MpiBlastEnableBioseqFetch()
{
	#ifndef BLAST_CS_API
    	ReadDBBioseqFetchEnable ("blastall", blast_database, db_is_na, TRUE);
	#endif
}

void MpiBlastDisableBioseqFetch()
{
	#ifndef BLAST_CS_API
        ReadDBBioseqFetchDisable();
	#endif
}

int m_SeqIdPtrCmp( SeqIdPtr a, SeqIdPtr b ){
	return SeqIdComp( a, b );
}

SeqIdPtr  m_TxGetSubjectIdFromSeqAlign( SeqAlignPtr seqalign ){
	return TxGetSubjectIdFromSeqAlign( seqalign );
}

Boolean m_GetScoreAndEvalue( SeqAlignPtr seqalign, Int4 *score, Nlm_FloatHi *bit_score, Nlm_FloatHi *evalue, Int4 *number ){
	return GetScoreAndEvalue( seqalign, score, bit_score, evalue, number );
}
