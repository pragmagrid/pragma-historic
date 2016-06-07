/*

  GASSCOPY.C
  ==========

  A simple example how to use GASS to copy a file
  from one url to another. 

  Arto Teräs <arto.teras@hip.fi>

*/

#include <stdlib.h>
#include <stdio.h>

#include <globus_common.h>
#include <globus_gass_copy.h>

int gasscopy(char* srcurl, char *desturl) {
   
    /* Maked by H.C.Lee 
    char*  srcurl;
    char*  desturl;
    */
    int             rval;            /*  return value (old APIs)     */
    globus_result_t rptr;            /*  return value (new APIs)     */
    globus_gass_copy_attr_t srcattr;
    globus_gass_copy_attr_t destattr;

    globus_gass_copy_handle_t handle;
    globus_gass_copy_handleattr_t handleattr;
 
    /* Maked by H.C.Lee 
    if (argc != 3 ) {
	fprintf(stderr, "GASSCOPY: Wrong number of arguments.\n");
	fprintf(stderr, "Usage: gasscopy <source url> <destination url>\n");
	exit(EXIT_FAILURE);
    }
    else {
	srcurl = argv[1];
	desturl = argv[2];
    } 
    */  

    rval = globus_module_activate(GLOBUS_GASS_COPY_MODULE);
    if (rval != GLOBUS_SUCCESS) {
	printf("GASSCOPY: failed\n");
	exit(EXIT_FAILURE);
    }

    rptr = globus_gass_copy_attr_init(&srcattr);
    if (rptr != GLOBUS_SUCCESS) {
	printf("GASSCOPY: Initializing srcattr failed\n");
	exit(EXIT_FAILURE);
    }

    rptr = globus_gass_copy_attr_init(&destattr);
    if (rptr != GLOBUS_SUCCESS) {
	printf("GASSCOPY: Initializing destattr failed\n");
	exit(EXIT_FAILURE);
    }

    rptr = globus_gass_copy_handleattr_init(&handleattr);
    if (rptr != GLOBUS_SUCCESS) {
	printf("GASSCOPY: Initializing handleattr failed\n");
	exit(EXIT_FAILURE);
    }

    rptr = globus_gass_copy_handle_init(&handle, &handleattr);
    if (rptr != GLOBUS_SUCCESS) {
	printf("GASSCOPY: Initializing handle failed\n");
	exit(EXIT_FAILURE);
    }

    rptr = globus_gass_copy_url_to_url(&handle, srcurl, &srcattr,
				       desturl, &destattr);
    if (rptr != GLOBUS_SUCCESS) {
	printf("GASSCOPY: Copy operation failed\n");
	exit(EXIT_FAILURE);
    }

    rptr = globus_gass_copy_handle_destroy(&handle);
    if (rptr != GLOBUS_SUCCESS) {
	printf("GASSCOPY: Destroying handle failed\n");
	exit(EXIT_FAILURE);
    }

    rval = globus_module_deactivate(GLOBUS_GASS_COPY_MODULE);
    if (rval != GLOBUS_SUCCESS) {
	printf("GASSCOPY: failed\n");
	exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
