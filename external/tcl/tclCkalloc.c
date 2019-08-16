/* 
 * tclCkalloc.c --
 *
 *    Interface to malloc and free that provides support for debugging problems
 *    involving overwritten, double freeing memory and loss of memory.
 *
 * Copyright (c) 1991-1994 The Regents of the University of California.
 * Copyright (c) 1994-1996 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * This code contributed by Karl Lehenbauer and Mark Diekhans
 *
 * RCS: @(#) $Id: tclCkalloc.c,v 1.1 2008-06-04 13:58:04 demin Exp $
 */

#include "tclInt.h"
#include "tclPort.h"

#define FALSE	0
#define TRUE	1

/*
 *----------------------------------------------------------------------
 *
 * Tcl_Alloc --
 *	Interface to TclpAlloc.
 *	It does check that memory was actually allocated.
 *
 *----------------------------------------------------------------------
 */

char *
Tcl_Alloc (size)
    unsigned int size;
{
        char *result;

        result = TclpAlloc(size);
        if (result == NULL) 
                panic("unable to alloc %d bytes", size);
        return result;
}

char *
Tcl_DbCkalloc(size, file, line)
    unsigned int size;
    char        *file;
    int          line;
{
    char *result;

    result = (char *) TclpAlloc(size);

    if (result == NULL) {
        fflush(stdout);
        panic("unable to alloc %d bytes, %s line %d", size, file, 
              line);
    }
    return result;
}


/*
 *----------------------------------------------------------------------
 *
 * Tcl_Realloc --
 *	Interface to TclpRealloc.
 *	It does check that memory was actually allocated.
 *
 *----------------------------------------------------------------------
 */

char *
Tcl_Realloc(ptr, size)
    char *ptr;
    unsigned int size;
{
    char *result;

    result = TclpRealloc(ptr, size);
    if (result == NULL) 
	panic("unable to realloc %d bytes", size);
    return result;
}

char *
Tcl_DbCkrealloc(ptr, size, file, line)
    char *ptr;
    unsigned int size;
    char *file;
    int line;
{
    char *result;

    result = (char *) TclpRealloc(ptr, size);

    if (result == NULL) {
        fflush(stdout);
        panic("unable to realloc %d bytes, %s line %d", size, file, 
              line);
    }
    return result;
}

/*
 *----------------------------------------------------------------------
 *
 * Tcl_Free --
 *	Interface to TclpFree.
 *
 *----------------------------------------------------------------------
 */

void
Tcl_Free (ptr)
    char *ptr;
{
        TclpFree(ptr);
}
