/* 
 * tclAlloc.c --
 *
 *	This is a very fast storage allocator.  It allocates blocks of a
 *	small number of different sizes, and keeps free lists of each size.
 *	Blocks that don't exactly fit are passed up to the next larger size.
 *	Blocks over a certain size are directly allocated from the system.
 *
 * Copyright (c) 1983 Regents of the University of California.
 * Copyright (c) 1996-1997 Sun Microsystems, Inc.
 *
 * Portions contributed by Chris Kingsley, Jack Jansen and Ray Johnson.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * RCS: @(#) $Id: tclAlloc.c,v 1.1 2008-06-04 13:58:03 demin Exp $
 */

#include "tclInt.h"
#include "tclPort.h"

/*
 *----------------------------------------------------------------------
 *
 * TclpAlloc --
 *
 *	Allocate more memory.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

char *
TclpAlloc(
    unsigned int nbytes)	/* Number of bytes to allocate. */
{
    return (char*) malloc(nbytes);
}

/*
 *----------------------------------------------------------------------
 *
 * TclpFree --
 *
 *	Free memory.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

void
TclpFree(
    char *cp)		/* Pointer to memory to free. */
{   
    free(cp);
    return;
}

/*
 *----------------------------------------------------------------------
 *
 * TclpRealloc --
 *
 *	Reallocate memory.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

char *
TclpRealloc(
    char *cp,			/* Pointer to alloced block. */
    unsigned int nbytes)	/* New size of memory. */
{   
    return (char*) realloc(cp, nbytes);
}

