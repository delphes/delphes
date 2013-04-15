/* 
 * tclPreserve.c --
 *
 *	This file contains a collection of procedures that are used
 *	to make sure that widget records and other data structures
 *	aren't reallocated when there are nested procedures that
 *	depend on their existence.
 *
 * Copyright (c) 1991-1994 The Regents of the University of California.
 * Copyright (c) 1994-1995 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * RCS: @(#) $Id: tclPreserve.c,v 1.1 2008-06-04 13:58:10 demin Exp $
 */

#include "tclInt.h"

/*
 * The following data structure is used to keep track of all the
 * Tcl_Preserve calls that are still in effect.  It grows as needed
 * to accommodate any number of calls in effect.
 */

typedef struct {
    ClientData clientData;	/* Address of preserved block. */
    int refCount;		/* Number of Tcl_Preserve calls in effect
				 * for block. */
    int mustFree;		/* Non-zero means Tcl_EventuallyFree was
				 * called while a Tcl_Preserve call was in
				 * effect, so the structure must be freed
				 * when refCount becomes zero. */
    Tcl_FreeProc *freeProc;	/* Procedure to call to free. */
} Reference;

static Reference *refArray;	/* First in array of references. */
static int inUse = 0;		/* Count of structures currently in use
				 * in refArray. */
#define INITIAL_SIZE 2

/*
 *----------------------------------------------------------------------
 *
 * Tcl_EventuallyFree --
 *
 *	Free up a block of memory, unless a call to Tcl_Preserve is in
 *	effect for that block.  In this case, defer the free until all
 *	calls to Tcl_Preserve have been undone by matching calls to
 *	Tcl_Release.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Ptr may be released by calling free().
 *
 *----------------------------------------------------------------------
 */

void
Tcl_EventuallyFree(clientData, freeProc)
    ClientData clientData;	/* Pointer to malloc'ed block of memory. */
    Tcl_FreeProc *freeProc;	/* Procedure to actually do free. */
{
    Reference *refPtr;
    int i;

    /*
     * See if there is a reference for this pointer.  If so, set its
     * "mustFree" flag (the flag had better not be set already!).
     */

    for (i = 0, refPtr = refArray; i < inUse; i++, refPtr++) {
	if (refPtr->clientData != clientData) {
	    continue;
	}
	if (refPtr->mustFree) {
	    panic("Tcl_EventuallyFree called twice for 0x%x\n", clientData);
        }
        refPtr->mustFree = 1;
	refPtr->freeProc = freeProc;
        return;
    }

    /*
     * No reference for this block.  Free it now.
     */

    if ((freeProc == TCL_DYNAMIC)
	    || (freeProc == (Tcl_FreeProc *) free)) {
	ckfree((char *) clientData);
    } else {
	(*freeProc)((char *)clientData);
    }
}
