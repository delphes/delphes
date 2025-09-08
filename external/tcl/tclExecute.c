/* 
 * tclExecute.c --
 *
 *	This file contains procedures that execute byte-compiled Tcl
 *	commands.
 *
 * Copyright (c) 1996-1997 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * RCS: @(#) $Id: tclExecute.c,v 1.1 2008-06-04 13:58:06 demin Exp $
 */

#include "tclInt.h"
#include "tclCompile.h"

#ifdef NO_FLOAT_H
#   include "../compat/float.h"
#else
#   include <float.h>
#endif
#ifndef TCL_NO_MATH
#include "tclMath.h"
#endif

/*
 * The stuff below is a bit of a hack so that this file can be used
 * in environments that include no UNIX, i.e. no errno.  Just define
 * errno here.
 */

#ifndef TCL_GENERIC_ONLY
#include "tclPort.h"
#else
#define NO_ERRNO_H
#endif

#ifdef NO_ERRNO_H
int errno;
#define EDOM 33
#define ERANGE 34
#endif

/*
 * Boolean flag indicating whether the Tcl bytecode interpreter has been
 * initialized.
 */

static int execInitialized = 0;

/*
 * The following global variable is use to signal matherr that Tcl
 * is responsible for the arithmetic, so errors can be handled in a
 * fashion appropriate for Tcl.  Zero means no Tcl math is in
 * progress;  non-zero means Tcl is doing math.
 */

int tcl_MathInProgress = 0;

/*
 * The variable below serves no useful purpose except to generate
 * a reference to matherr, so that the Tcl version of matherr is
 * linked in rather than the system version. Without this reference
 * the need for matherr won't be discovered during linking until after
 * libtcl.a has been processed, so Tcl's version won't be used.
 */

#ifdef NEED_MATHERR
extern int matherr();
int (*tclMatherrPtr)() = matherr;
#endif

/*
 * Array of instruction names.
 */

static char *opName[256];

/*
 * Mapping from expression instruction opcodes to strings; used for error
 * messages. Note that these entries must match the order and number of the
 * expression opcodes (e.g., INST_LOR) in tclCompile.h.
 */

static char *operatorStrings[] = {
    "||", "&&", "|", "^", "&", "==", "!=", "<", ">", "<=", ">=", "<<", ">>",
    "+", "-", "*", "/", "%", "+", "-", "~", "!",
    "BUILTIN FUNCTION", "FUNCTION"
};
    
/*
 * Macros for testing floating-point values for certain special cases. Test
 * for not-a-number by comparing a value against itself; test for infinity
 * by comparing against the largest floating-point value.
 */

#define IS_NAN(v) ((v) != (v))
#ifdef DBL_MAX
#   define IS_INF(v) (((v) > DBL_MAX) || ((v) < -DBL_MAX))
#else
#   define IS_INF(v) 0
#endif

/*
 * Macro to adjust the program counter and restart the instruction execution
 * loop after each instruction is executed.
 */

#define ADJUST_PC(instBytes) \
    pc += instBytes;  continue

/*
 * Macros used to cache often-referenced Tcl evaluation stack information
 * in local variables. Note that a DECACHE_STACK_INFO()-CACHE_STACK_INFO()
 * pair must surround any call inside TclExecuteByteCode (and a few other
 * procedures that use this scheme) that could result in a recursive call
 * to TclExecuteByteCode.
 */

#define CACHE_STACK_INFO() \
    stackPtr = eePtr->stackPtr; \
    stackTop = eePtr->stackTop

#define DECACHE_STACK_INFO() \
    eePtr->stackTop = stackTop

/*
 * Macros used to access items on the Tcl evaluation stack. PUSH_OBJECT
 * increments the object's ref count since it makes the stack have another
 * reference pointing to the object. However, POP_OBJECT does not decrement
 * the ref count. This is because the stack may hold the only reference to
 * the object, so the object would be destroyed if its ref count were
 * decremented before the caller had a chance to, e.g., store it in a
 * variable. It is the caller's responsibility to decrement the ref count
 * when it is finished with an object.
 */

#define STK_ITEM(offset)    (stackPtr[stackTop + (offset)])
#define STK_OBJECT(offset)  (STK_ITEM(offset).o)
#define STK_INT(offset)     (STK_ITEM(offset).i)
#define STK_POINTER(offset) (STK_ITEM(offset).p)

/*
 * WARNING! It is essential that objPtr only appear once in the PUSH_OBJECT
 * macro. The actual parameter might be an expression with side effects,
 * and this ensures that it will be executed only once. 
 */
    
#define PUSH_OBJECT(objPtr) \
    Tcl_IncrRefCount(stackPtr[++stackTop].o = (objPtr))
    
#define POP_OBJECT() \
    (stackPtr[stackTop--].o)

/*
 * Declarations for local procedures to this file:
 */

static void		CallTraceProcedure _ANSI_ARGS_((Tcl_Interp *interp,
			    Trace *tracePtr, Command *cmdPtr,
			    char *command, int numChars,
			    int objc, Tcl_Obj *objv[]));
static void		DupCmdNameInternalRep _ANSI_ARGS_((Tcl_Obj *objPtr,
			    Tcl_Obj *copyPtr));
static int		ExprAbsFunc _ANSI_ARGS_((Tcl_Interp *interp,
			    ExecEnv *eePtr, ClientData clientData));
static int		ExprBinaryFunc _ANSI_ARGS_((Tcl_Interp *interp,
			    ExecEnv *eePtr, ClientData clientData));
static int		ExprCallMathFunc _ANSI_ARGS_((Tcl_Interp *interp,
			    ExecEnv *eePtr, int objc, Tcl_Obj **objv));
static int		ExprDoubleFunc _ANSI_ARGS_((Tcl_Interp *interp,
			    ExecEnv *eePtr, ClientData clientData));
static int		ExprIntFunc _ANSI_ARGS_((Tcl_Interp *interp,
			    ExecEnv *eePtr, ClientData clientData));
static int		ExprRoundFunc _ANSI_ARGS_((Tcl_Interp *interp,
			    ExecEnv *eePtr, ClientData clientData));
static int		ExprUnaryFunc _ANSI_ARGS_((Tcl_Interp *interp,
			    ExecEnv *eePtr, ClientData clientData));
static void		FreeCmdNameInternalRep _ANSI_ARGS_((
    			    Tcl_Obj *objPtr));
static char *		GetSrcInfoForPc _ANSI_ARGS_((unsigned char *pc,
        		    ByteCode* codePtr, int *lengthPtr));
static void		GrowEvaluationStack _ANSI_ARGS_((ExecEnv *eePtr));
static void		IllegalExprOperandType _ANSI_ARGS_((
			    Tcl_Interp *interp, unsigned int opCode,
			    Tcl_Obj *opndPtr));
static void		InitByteCodeExecution _ANSI_ARGS_((
			    Tcl_Interp *interp));
static void		PrintByteCodeInfo _ANSI_ARGS_((ByteCode *codePtr));
static void		RecordTracebackInfo _ANSI_ARGS_((Tcl_Interp *interp,
			    unsigned char *pc, ByteCode *codePtr));
static int		SetCmdNameFromAny _ANSI_ARGS_((Tcl_Interp *interp,
			    Tcl_Obj *objPtr));
static void		UpdateStringOfCmdName _ANSI_ARGS_((Tcl_Obj *objPtr));

/*
 * Table describing the built-in math functions. Entries in this table are
 * indexed by the values of the INST_CALL_BUILTIN_FUNC instruction's
 * operand byte.
 */

BuiltinFunc builtinFuncTable[] = {
#ifndef TCL_NO_MATH
    {"acos", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) acos},
    {"asin", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) asin},
    {"atan", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) atan},
    {"atan2", 2, {TCL_DOUBLE, TCL_DOUBLE}, ExprBinaryFunc, (ClientData) atan2},
    {"ceil", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) ceil},
    {"cos", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) cos},
    {"cosh", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) cosh},
    {"exp", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) exp},
    {"floor", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) floor},
    {"fmod", 2, {TCL_DOUBLE, TCL_DOUBLE}, ExprBinaryFunc, (ClientData) fmod},
    {"hypot", 2, {TCL_DOUBLE, TCL_DOUBLE}, ExprBinaryFunc, (ClientData) hypot},
    {"log", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) log},
    {"log10", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) log10},
    {"pow", 2, {TCL_DOUBLE, TCL_DOUBLE}, ExprBinaryFunc, (ClientData) pow},
    {"sin", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) sin},
    {"sinh", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) sinh},
    {"sqrt", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) sqrt},
    {"tan", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) tan},
    {"tanh", 1, {TCL_DOUBLE}, ExprUnaryFunc, (ClientData) tanh},
#endif
    {"abs", 1, {TCL_EITHER}, ExprAbsFunc, 0},
    {"double", 1, {TCL_EITHER}, ExprDoubleFunc, 0},
    {"int", 1, {TCL_EITHER}, ExprIntFunc, 0},
    {"round", 1, {TCL_EITHER}, ExprRoundFunc, 0},
    {0},
};

/*
 * The structure below defines the command name Tcl object type by means of
 * procedures that can be invoked by generic object code. Objects of this
 * type cache the Command pointer that results from looking up command names
 * in the command hashtable. Such objects appear as the zeroth ("command
 * name") argument in a Tcl command.
 */

Tcl_ObjType tclCmdNameType = {
    "cmdName",				/* name */
    FreeCmdNameInternalRep,		/* freeIntRepProc */
    DupCmdNameInternalRep,		/* dupIntRepProc */
    UpdateStringOfCmdName,		/* updateStringProc */
    SetCmdNameFromAny			/* setFromAnyProc */
};

/*
 *----------------------------------------------------------------------
 *
 * InitByteCodeExecution --
 *
 *	This procedure is called once to initialize the Tcl bytecode
 *	interpreter.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	This procedure initializes the array of instruction names.
 *
 *----------------------------------------------------------------------
 */

static void
InitByteCodeExecution(interp)
    Tcl_Interp *interp;		/* Interpreter for which the Tcl variable
				 * "tcl_traceExec" is linked to control
				 * instruction tracing. */
{
    int i;
    
    Tcl_RegisterObjType(&tclCmdNameType);

    (VOID *) memset(opName, 0, sizeof(opName));
    for (i = 0;  instructionTable[i].name != NULL;  i++) {
	opName[i] = instructionTable[i].name;
    }
}

/*
 *----------------------------------------------------------------------
 *
 * TclCreateExecEnv --
 *
 *	This procedure creates a new execution environment for Tcl bytecode
 *	execution. An ExecEnv points to a Tcl evaluation stack. An ExecEnv
 *	is typically created once for each Tcl interpreter (Interp
 *	structure) and recursively passed to TclExecuteByteCode to execute
 *	ByteCode sequences for nested commands.
 *
 * Results:
 *	A newly allocated ExecEnv is returned. This points to an empty
 *	evaluation stack of the standard initial size.
 *
 * Side effects:
 *	The bytecode interpreter is also initialized here, as this
 *	procedure will be called before any call to TclExecuteByteCode.
 *
 *----------------------------------------------------------------------
 */

#define TCL_STACK_INITIAL_SIZE 2000

ExecEnv *
TclCreateExecEnv(interp)
    Tcl_Interp *interp;		/* Interpreter for which the execution
				 * environment is being created. */
{
    ExecEnv *eePtr = (ExecEnv *) ckalloc(sizeof(ExecEnv));

    eePtr->stackPtr = (StackItem *)
	ckalloc((unsigned) (TCL_STACK_INITIAL_SIZE * sizeof(StackItem)));
    eePtr->stackTop = -1;
    eePtr->stackEnd = (TCL_STACK_INITIAL_SIZE - 1);

    if (!execInitialized) {
        TclInitAuxDataTypeTable();
        InitByteCodeExecution(interp);
        execInitialized = 1;
    }

    return eePtr;
}
#undef TCL_STACK_INITIAL_SIZE

/*
 *----------------------------------------------------------------------
 *
 * TclDeleteExecEnv --
 *
 *	Frees the storage for an ExecEnv.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Storage for an ExecEnv and its contained storage (e.g. the
 *	evaluation stack) is freed.
 *
 *----------------------------------------------------------------------
 */

void
TclDeleteExecEnv(eePtr)
    ExecEnv *eePtr;		/* Execution environment to free. */
{
    ckfree((char *) eePtr->stackPtr);
    ckfree((char *) eePtr);
}

/*
 *----------------------------------------------------------------------
 *
 * TclFinalizeExecEnv --
 *
 *	Finalizes the execution environment setup so that it can be
 *	later reinitialized.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	After this call, the next time TclCreateExecEnv will be called
 *	it will call InitByteCodeExecution.
 *
 *----------------------------------------------------------------------
 */

void
TclFinalizeExecEnv()
{
    execInitialized = 0;
    TclFinalizeAuxDataTypeTable();
}

/*
 *----------------------------------------------------------------------
 *
 * GrowEvaluationStack --
 *
 *	This procedure grows a Tcl evaluation stack stored in an ExecEnv.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The size of the evaluation stack is doubled.
 *
 *----------------------------------------------------------------------
 */

static void
GrowEvaluationStack(eePtr)
    register ExecEnv *eePtr; /* Points to the ExecEnv with an evaluation
			      * stack to enlarge. */
{
    /*
     * The current Tcl stack elements are stored from eePtr->stackPtr[0]
     * to eePtr->stackPtr[eePtr->stackEnd] (inclusive).
     */

    int currElems = (eePtr->stackEnd + 1);
    int newElems  = 2*currElems;
    int currBytes = currElems * sizeof(StackItem);
    int newBytes  = 2*currBytes;
    StackItem *newStackPtr = (StackItem *) ckalloc((unsigned) newBytes);

    /*
     * Copy the existing stack items to the new stack space, free the old
     * storage if appropriate, and mark new space as malloc'ed.
     */
 
    memcpy((VOID *) newStackPtr, (VOID *) eePtr->stackPtr,
	   (size_t) currBytes);
    ckfree((char *) eePtr->stackPtr);
    eePtr->stackPtr = newStackPtr;
    eePtr->stackEnd = (newElems - 1); /* i.e. index of last usable item */
}

/*
 *----------------------------------------------------------------------
 *
 * TclExecuteByteCode --
 *
 *	This procedure executes the instructions of a ByteCode structure.
 *	It returns when a "done" instruction is executed or an error occurs.
 *
 * Results:
 *	The return value is one of the return codes defined in tcl.h
 *	(such as TCL_OK), and interp->objResultPtr refers to a Tcl object
 *	that either contains the result of executing the code or an
 *	error message.
 *
 * Side effects:
 *	Almost certainly, depending on the ByteCode's instructions.
 *
 *----------------------------------------------------------------------
 */

int
TclExecuteByteCode(interp, codePtr)
    Tcl_Interp *interp;		/* Token for command interpreter. */
    ByteCode *codePtr;		/* The bytecode sequence to interpret. */
{
    Interp *iPtr = (Interp *) interp;
    ExecEnv *eePtr = iPtr->execEnvPtr;
    				/* Points to the execution environment. */
    register StackItem *stackPtr = eePtr->stackPtr;
    				/* Cached evaluation stack base pointer. */
    register int stackTop = eePtr->stackTop;
    				/* Cached top index of evaluation stack. */
    Tcl_Obj **objArrayPtr = codePtr->objArrayPtr;
    				/* Points to the ByteCode's object array. */
    unsigned char *pc = codePtr->codeStart;
				/* The current program counter. */
    unsigned char opCode;	/* The current instruction code. */
    int opnd;			/* Current instruction's operand byte. */
    int pcAdjustment;		/* Hold pc adjustment after instruction. */
    int initStackTop = stackTop;/* Stack top at start of execution. */
    ExceptionRange *rangePtr;	/* Points to closest loop or catch exception
				 * range enclosing the pc. Used by various
				 * instructions and processCatch to
				 * process break, continue, and errors. */
    int result = TCL_OK;	/* Return code returned after execution. */
    Tcl_Obj *valuePtr, *value2Ptr, *namePtr, *objPtr;
    char *bytes;
    int length;
    long i;

    /*
     * This procedure uses a stack to hold information about catch commands.
     * This information is the current operand stack top when starting to
     * execute the code for each catch command. It starts out with stack-
     * allocated space but uses dynamically-allocated storage if needed.
     */

#define STATIC_CATCH_STACK_SIZE 5
    int (catchStackStorage[STATIC_CATCH_STACK_SIZE]);
    int *catchStackPtr = catchStackStorage;
    int catchTop = -1;

    /*
     * Make sure the catch stack is large enough to hold the maximum number
     * of catch commands that could ever be executing at the same time. This
     * will be no more than the exception range array's depth.
     */

    if (codePtr->maxExcRangeDepth > STATIC_CATCH_STACK_SIZE) {
	catchStackPtr = (int *)
	        ckalloc(codePtr->maxExcRangeDepth * sizeof(int));
    }

    /*
     * Make sure the stack has enough room to execute this ByteCode.
     */

    while ((stackTop + codePtr->maxStackDepth) > eePtr->stackEnd) {
        GrowEvaluationStack(eePtr); 
        stackPtr = eePtr->stackPtr;
    }

    /*
     * Loop executing instructions until a "done" instruction, a TCL_RETURN,
     * or some error.
     */

    for (;;) {
	opCode = *pc;

        switch (opCode) {
	case INST_DONE:
	    /*
	     * Pop the topmost object from the stack, set the interpreter's
	     * object result to point to it, and return.
	     */
	    valuePtr = POP_OBJECT();
	    Tcl_SetObjResult(interp, valuePtr);
	    TclDecrRefCount(valuePtr);
	    if (stackTop != initStackTop) {
		fprintf(stderr, "\nTclExecuteByteCode: done instruction at pc %u: stack top %d != entry stack top %d\n",
			(unsigned int)(pc - codePtr->codeStart),
			(unsigned int) stackTop,
			(unsigned int) initStackTop);
		fprintf(stderr, "  Source: ");
		TclPrintSource(stderr, codePtr->source, 150);
		panic("TclExecuteByteCode execution failure: end stack top != start stack top");
	    }
	    goto done;
	    
	case INST_PUSH1:
	    valuePtr = objArrayPtr[TclGetUInt1AtPtr(pc+1)];
	    PUSH_OBJECT(valuePtr);
	    ADJUST_PC(2);
	    
	case INST_PUSH4:
	    valuePtr = objArrayPtr[TclGetUInt4AtPtr(pc+1)];
	    PUSH_OBJECT(valuePtr);
	    ADJUST_PC(5);
	    
	case INST_POP:
	    valuePtr = POP_OBJECT();
	    TclDecrRefCount(valuePtr); /* finished with pop'ed object. */
	    ADJUST_PC(1);

	case INST_DUP:
	    valuePtr = stackPtr[stackTop].o;
	    PUSH_OBJECT(Tcl_DuplicateObj(valuePtr));
	    ADJUST_PC(1);

	case INST_CONCAT1:
	    opnd = TclGetUInt1AtPtr(pc+1);
	    {
		Tcl_Obj *concatObjPtr;
		int totalLen = 0;

		/*
		 * Concatenate strings (with no separators) from the top
		 * opnd items on the stack starting with the deepest item.
		 * First, determine how many characters are needed.
		 */

		for (i = (stackTop - (opnd-1));  i <= stackTop;  i++) {
		    valuePtr = stackPtr[i].o;
		    bytes = TclGetStringFromObj(valuePtr, &length);
		    if (bytes != NULL) {
			totalLen += length;
		    }
                }

		/*
		 * Initialize the new append string object by appending the
		 * strings of the opnd stack objects. Also pop the objects. 
		 */

		TclNewObj(concatObjPtr);
		if (totalLen > 0) {
		    char *p = (char *) ckalloc((unsigned) (totalLen + 1));
		    concatObjPtr->bytes = p;
		    concatObjPtr->length = totalLen;
		    for (i = (stackTop - (opnd-1));  i <= stackTop;  i++) {
			valuePtr = stackPtr[i].o;
			bytes = TclGetStringFromObj(valuePtr, &length);
			if (bytes != NULL) {
			    memcpy((VOID *) p, (VOID *) bytes,
			            (size_t) length);
			    p += length;
			}
			TclDecrRefCount(valuePtr);
		    }
		    *p = '\0';
		} else {
		    for (i = (stackTop - (opnd-1));  i <= stackTop;  i++) {
			valuePtr = stackPtr[i].o;
			Tcl_DecrRefCount(valuePtr);
		    }
		}
		stackTop -= opnd;
		
		PUSH_OBJECT(concatObjPtr);
		ADJUST_PC(2);
            }
	    
	case INST_INVOKE_STK4:
	    opnd = TclGetUInt4AtPtr(pc+1);
	    pcAdjustment = 5;
	    goto doInvocation;

	case INST_INVOKE_STK1:
	    opnd = TclGetUInt1AtPtr(pc+1);
	    pcAdjustment = 2;
	    
	    doInvocation:
	    {
		char *cmdName;
		Command *cmdPtr;   /* Points to command's Command struct. */
		int objc = opnd;   /* The number of arguments. */
		Tcl_Obj **objv;	   /* The array of argument objects. */
		Tcl_Obj *objv0Ptr; /* Holds objv[0], the command name. */
		int newPcOffset = 0;
				   /* Instruction offset computed during
				    * break, continue, error processing.
				    * Init. to avoid compiler warning. */
		Tcl_Command cmd;
		
		/*
		 * If the interpreter was deleted, return an error.
		 */
		
		if (iPtr->flags & DELETED) {
		    Tcl_ResetResult(interp);
		    Tcl_AppendToObj(Tcl_GetObjResult(interp),
		            "attempt to call eval in deleted interpreter", -1);
		    Tcl_SetErrorCode(interp, "CORE", "IDELETE",
			    "attempt to call eval in deleted interpreter",
			    (char *) NULL);
		    result = TCL_ERROR;
		    goto checkForCatch;
		}
    
		objv = &(stackPtr[stackTop - (objc-1)].o);
		objv0Ptr = objv[0];
		cmdName = TclGetStringFromObj(objv0Ptr, (int *) NULL);
		
		/*
		 * Find the procedure to execute this command. If there
		 * isn't one, then see if there is a command "unknown". If
		 * so, invoke it, passing it the original command words as
		 * arguments.
		 *
		 * We convert the objv[0] object to be a CmdName object.
		 * This caches a pointer to the Command structure for the
		 * command; this pointer is held in a ResolvedCmdName
		 * structure the object's internal rep. points to.
		 */

		cmd = Tcl_GetCommandFromObj(interp, objv0Ptr);
		cmdPtr = (Command *) cmd;
		
		/*
		 * If the command is still not found, handle it with the
		 * "unknown" proc.
		 */

		if (cmdPtr == NULL) {
		    cmd = Tcl_FindCommand(interp, "unknown",
                            (Tcl_Namespace *) NULL, /*flags*/ TCL_GLOBAL_ONLY);
                    if (cmd == (Tcl_Command) NULL) {
			Tcl_ResetResult(interp);
			Tcl_AppendStringsToObj(Tcl_GetObjResult(interp),
			        "invalid command name \"", cmdName, "\"",
				(char *) NULL);
			result = TCL_ERROR;
			goto checkForCatch;
		    }
		    cmdPtr = (Command *) cmd;
		    stackTop++; /* need room for new inserted objv[0] */
		    for (i = objc;  i >= 0;  i--) {
			objv[i+1] = objv[i];
		    }
		    objc++;
		    objv[0] = Tcl_NewStringObj("unknown", -1);
		    Tcl_IncrRefCount(objv[0]);
		}
		
		/*
		 * Call any trace procedures.
		 */

		if (iPtr->tracePtr != NULL) {
		    Trace *tracePtr, *nextTracePtr;

		    for (tracePtr = iPtr->tracePtr;  tracePtr != NULL;
		            tracePtr = nextTracePtr) {
			nextTracePtr = tracePtr->nextPtr;
			if (iPtr->numLevels <= tracePtr->level) {
			    int numChars;
			    char *cmd = GetSrcInfoForPc(pc, codePtr,
				    &numChars);
			    if (cmd != NULL) {
				DECACHE_STACK_INFO();
				CallTraceProcedure(interp, tracePtr, cmdPtr,
				        cmd, numChars, objc, objv);
				CACHE_STACK_INFO();
			    }
			}
		    }
		}
		
		/*
		 * Finally, invoke the command's Tcl_ObjCmdProc. First reset
		 * the interpreter's string and object results to their
		 * default empty values since they could have gotten changed
		 * by earlier invocations.
		 */
		
		Tcl_ResetResult(interp);

		iPtr->cmdCount++;
		DECACHE_STACK_INFO();
		result = (*cmdPtr->objProc)(cmdPtr->objClientData, interp,
					    objc, objv);
		CACHE_STACK_INFO();

		/*
		 * If the interpreter has a non-empty string result, the
		 * result object is either empty or stale because some
		 * procedure set interp->result directly. If so, move the
		 * string result to the result object, then reset the
		 * string result.
		 */

		if (*(iPtr->result) != 0) {
		    (void) Tcl_GetObjResult(interp);
		}
		
		/*
		 * Pop the objc top stack elements and decrement their ref
		 * counts. 
		 */
		
		i = (stackTop - (objc-1));
		while (i <= stackTop) {
		    valuePtr = stackPtr[i].o;
		    TclDecrRefCount(valuePtr);
		    i++;
		}
		stackTop -= objc;

		/*
		 * Process the result of the Tcl_ObjCmdProc call.
		 */
		
		switch (result) {
		case TCL_OK:
		    /*
		     * Push the call's object result and continue execution
		     * with the next instruction.
		     */
		    PUSH_OBJECT(Tcl_GetObjResult(interp));
		    ADJUST_PC(pcAdjustment);
		    
		case TCL_BREAK:
		case TCL_CONTINUE:
		    /*
		     * The invoked command requested a break or continue.
		     * Find the closest enclosing loop or catch exception
		     * range, if any. If a loop is found, terminate its
		     * execution or skip to its next iteration. If the
		     * closest is a catch exception range, jump to its
		     * catchOffset. If no enclosing range is found, stop
		     * execution and return the TCL_BREAK or TCL_CONTINUE.
		     */
		    rangePtr = TclGetExceptionRangeForPc(pc,
                            /*catchOnly*/ 0, codePtr);
		    if (rangePtr == NULL) {
			goto abnormalReturn; /* no catch exists to check */
		    }
		    switch (rangePtr->type) {
		    case LOOP_EXCEPTION_RANGE:
			if (result == TCL_BREAK) {
			    newPcOffset = rangePtr->breakOffset;
			} else if (rangePtr->continueOffset == -1) {
			    goto checkForCatch;
			} else {
			    newPcOffset = rangePtr->continueOffset;
			}
			break;
		    case CATCH_EXCEPTION_RANGE:
			goto processCatch; /* it will use rangePtr */
		    default:
			panic("TclExecuteByteCode: unrecognized ExceptionRange type %d\n", rangePtr->type);
		    }
		    result = TCL_OK;
		    pc = (codePtr->codeStart + newPcOffset);
		    continue;	/* restart outer instruction loop at pc */
		    
		case TCL_ERROR:
		    /*
		     * The invoked command returned an error. Look for an
		     * enclosing catch exception range, if any.
		     */
		    goto checkForCatch;

		case TCL_RETURN:
		    /*
		     * The invoked command requested that the current
		     * procedure stop execution and return. First check
		     * for an enclosing catch exception range, if any.
		     */
		    goto checkForCatch;

		default:
		    goto checkForCatch;
		} /* end of switch on result from invoke instruction */
	    }
	    
	case INST_EVAL_STK:
	    objPtr = POP_OBJECT();
	    DECACHE_STACK_INFO();
	    result = Tcl_EvalObj(interp, objPtr);
	    CACHE_STACK_INFO();
	    if (result == TCL_OK) {
		/*
		 * Normal return; push the eval's object result.
		 */
		
		PUSH_OBJECT(Tcl_GetObjResult(interp));
		TclDecrRefCount(objPtr);
		ADJUST_PC(1);
	    } else if ((result == TCL_BREAK) || (result == TCL_CONTINUE)) {
		/*
		 * Find the closest enclosing loop or catch exception range,
		 * if any. If a loop is found, terminate its execution or
		 * skip to its next iteration. If the closest is a catch
		 * exception range, jump to its catchOffset. If no enclosing
		 * range is found, stop execution and return that same
		 * TCL_BREAK or TCL_CONTINUE.
		 */

		int newPcOffset = 0; /* Pc offset computed during break,
				      * continue, error processing. Init.
				      * to avoid compiler warning. */

		rangePtr = TclGetExceptionRangeForPc(pc, /*catchOnly*/ 0,
			codePtr);
		if (rangePtr == NULL) {
		    Tcl_DecrRefCount(objPtr);
		    goto abnormalReturn;    /* no catch exists to check */
		}
		switch (rangePtr->type) {
		case LOOP_EXCEPTION_RANGE:
		    if (result == TCL_BREAK) {
			newPcOffset = rangePtr->breakOffset;
		    } else if (rangePtr->continueOffset == -1) {
			Tcl_DecrRefCount(objPtr);
			goto checkForCatch;
		    } else {
			newPcOffset = rangePtr->continueOffset;
		    }
		    result = TCL_OK;
		    break;
		case CATCH_EXCEPTION_RANGE:
		    Tcl_DecrRefCount(objPtr);
		    goto processCatch;  /* it will use rangePtr */
		default:
		    panic("TclExecuteByteCode: unrecognized ExceptionRange type %d\n", rangePtr->type);
		}
		Tcl_DecrRefCount(objPtr);
		pc = (codePtr->codeStart + newPcOffset);
		continue;	/* restart outer instruction loop at pc */
	    } else { /* eval returned TCL_ERROR, TCL_RETURN, unknown code */
		Tcl_DecrRefCount(objPtr);
		goto checkForCatch;
	    }

	case INST_EXPR_STK:
	    objPtr = POP_OBJECT();
	    Tcl_ResetResult(interp);
	    DECACHE_STACK_INFO();
	    result = Tcl_ExprObj(interp, objPtr, &valuePtr);
	    CACHE_STACK_INFO();
	    if (result != TCL_OK) {
		Tcl_DecrRefCount(objPtr);
		goto checkForCatch;
	    }
	    stackPtr[++stackTop].o = valuePtr; /* already has right refct */
	    TclDecrRefCount(objPtr);
	    ADJUST_PC(1);

	case INST_LOAD_SCALAR4:
	    opnd = TclGetInt4AtPtr(pc+1);
	    pcAdjustment = 5;
	    goto doLoadScalar;

	case INST_LOAD_SCALAR1:
	    opnd = TclGetUInt1AtPtr(pc+1);
	    pcAdjustment = 2;
	    
	    doLoadScalar:
	    DECACHE_STACK_INFO();
	    valuePtr = TclGetIndexedScalar(interp, opnd,
					   /*leaveErrorMsg*/ 1);
	    CACHE_STACK_INFO();
	    if (valuePtr == NULL) {
		result = TCL_ERROR;
		goto checkForCatch;
            }
	    PUSH_OBJECT(valuePtr);
	    ADJUST_PC(pcAdjustment);

	case INST_LOAD_SCALAR_STK:
	    namePtr = POP_OBJECT();
	    DECACHE_STACK_INFO();
	    valuePtr = Tcl_ObjGetVar2(interp, namePtr, (Tcl_Obj *) NULL, 
				      TCL_LEAVE_ERR_MSG);
	    CACHE_STACK_INFO();
	    if (valuePtr == NULL) {
		Tcl_DecrRefCount(namePtr);
		result = TCL_ERROR;
		goto checkForCatch;
            }
	    PUSH_OBJECT(valuePtr);
	    TclDecrRefCount(namePtr);
	    ADJUST_PC(1);

	case INST_LOAD_ARRAY4:
	    opnd = TclGetUInt4AtPtr(pc+1);
	    pcAdjustment = 5;
	    goto doLoadArray;

	case INST_LOAD_ARRAY1:
	    opnd = TclGetUInt1AtPtr(pc+1);
	    pcAdjustment = 2;
	    
	    doLoadArray:
	    {
		Tcl_Obj *elemPtr = POP_OBJECT();
		
		DECACHE_STACK_INFO();
		valuePtr = TclGetElementOfIndexedArray(interp, opnd,
	                elemPtr, /*leaveErrorMsg*/ 1);
		CACHE_STACK_INFO();
		if (valuePtr == NULL) {
		    Tcl_DecrRefCount(elemPtr);
		    result = TCL_ERROR;
		    goto checkForCatch;
		}
		PUSH_OBJECT(valuePtr);
		TclDecrRefCount(elemPtr);
	    }
	    ADJUST_PC(pcAdjustment);

	case INST_LOAD_ARRAY_STK:
	    {
		Tcl_Obj *elemPtr = POP_OBJECT();
		
		namePtr = POP_OBJECT();
		DECACHE_STACK_INFO();
		valuePtr = Tcl_ObjGetVar2(interp, namePtr, elemPtr,
		        TCL_LEAVE_ERR_MSG);
		CACHE_STACK_INFO();
		if (valuePtr == NULL) {
		    Tcl_DecrRefCount(namePtr);
		    Tcl_DecrRefCount(elemPtr);
		    result = TCL_ERROR;
		    goto checkForCatch;
		}
		PUSH_OBJECT(valuePtr);
		TclDecrRefCount(namePtr);
		TclDecrRefCount(elemPtr);
	    }
	    ADJUST_PC(1);

	case INST_LOAD_STK:
	    namePtr = POP_OBJECT();
	    DECACHE_STACK_INFO();
	    valuePtr = Tcl_ObjGetVar2(interp, namePtr, NULL,
		    TCL_PARSE_PART1|TCL_LEAVE_ERR_MSG);
	    CACHE_STACK_INFO();
	    if (valuePtr == NULL) {
		Tcl_DecrRefCount(namePtr);
		result = TCL_ERROR;
		goto checkForCatch;
	    }
	    PUSH_OBJECT(valuePtr);
	    TclDecrRefCount(namePtr);
	    ADJUST_PC(1);
	    
	case INST_STORE_SCALAR4:
	    opnd = TclGetUInt4AtPtr(pc+1);
	    pcAdjustment = 5;
	    goto doStoreScalar;

	case INST_STORE_SCALAR1:
	    opnd = TclGetUInt1AtPtr(pc+1);
	    pcAdjustment = 2;
	    
	  doStoreScalar:
	    valuePtr = POP_OBJECT();
	    DECACHE_STACK_INFO();
	    value2Ptr = TclSetIndexedScalar(interp, opnd, valuePtr,
					      /*leaveErrorMsg*/ 1);
	    CACHE_STACK_INFO();
	    if (value2Ptr == NULL) {
		Tcl_DecrRefCount(valuePtr);
		result = TCL_ERROR;
		goto checkForCatch;
	    }
	    PUSH_OBJECT(value2Ptr);
	    TclDecrRefCount(valuePtr);
	    ADJUST_PC(pcAdjustment);

	case INST_STORE_SCALAR_STK:
	    valuePtr = POP_OBJECT();
	    namePtr = POP_OBJECT();
	    DECACHE_STACK_INFO();
	    value2Ptr = Tcl_ObjSetVar2(interp, namePtr, NULL, valuePtr,
	            TCL_LEAVE_ERR_MSG);
	    CACHE_STACK_INFO();
	    if (value2Ptr == NULL) {
		Tcl_DecrRefCount(namePtr);
		Tcl_DecrRefCount(valuePtr);
		result = TCL_ERROR;
		goto checkForCatch;
	    }
	    PUSH_OBJECT(value2Ptr);
	    TclDecrRefCount(namePtr);
	    TclDecrRefCount(valuePtr);
	    ADJUST_PC(1);

	case INST_STORE_ARRAY4:
	    opnd = TclGetUInt4AtPtr(pc+1);
	    pcAdjustment = 5;
	    goto doStoreArray;

	case INST_STORE_ARRAY1:
	    opnd = TclGetUInt1AtPtr(pc+1);
	    pcAdjustment = 2;
	    
	    doStoreArray:
	    {
		Tcl_Obj *elemPtr;

		valuePtr = POP_OBJECT();
		elemPtr = POP_OBJECT();
		DECACHE_STACK_INFO();
		value2Ptr = TclSetElementOfIndexedArray(interp, opnd,
		        elemPtr, valuePtr, TCL_LEAVE_ERR_MSG);
		CACHE_STACK_INFO();
		if (value2Ptr == NULL) {
		    Tcl_DecrRefCount(elemPtr);
		    Tcl_DecrRefCount(valuePtr);
		    result = TCL_ERROR;
		    goto checkForCatch;
		}
		PUSH_OBJECT(value2Ptr);
		TclDecrRefCount(elemPtr);
		TclDecrRefCount(valuePtr);
	    }
	    ADJUST_PC(pcAdjustment);

	case INST_STORE_ARRAY_STK:
	    {
		Tcl_Obj *elemPtr;

		valuePtr = POP_OBJECT();
		elemPtr = POP_OBJECT();
		namePtr = POP_OBJECT();
		DECACHE_STACK_INFO();
		value2Ptr = Tcl_ObjSetVar2(interp, namePtr, elemPtr,
		        valuePtr, TCL_LEAVE_ERR_MSG);
		CACHE_STACK_INFO();
		if (value2Ptr == NULL) {
		    Tcl_DecrRefCount(namePtr);
		    Tcl_DecrRefCount(elemPtr);
		    Tcl_DecrRefCount(valuePtr);
		    result = TCL_ERROR;
		    goto checkForCatch;
		}
		PUSH_OBJECT(value2Ptr);
		TclDecrRefCount(namePtr);
		TclDecrRefCount(elemPtr);
		TclDecrRefCount(valuePtr);
	    }
	    ADJUST_PC(1);

	case INST_STORE_STK:
	    valuePtr = POP_OBJECT();
	    namePtr = POP_OBJECT();
	    DECACHE_STACK_INFO();
	    value2Ptr = Tcl_ObjSetVar2(interp, namePtr, NULL, valuePtr,
		    TCL_PARSE_PART1|TCL_LEAVE_ERR_MSG);
	    CACHE_STACK_INFO();
	    if (value2Ptr == NULL) {
		Tcl_DecrRefCount(namePtr);
		Tcl_DecrRefCount(valuePtr);
		result = TCL_ERROR;
		goto checkForCatch;
	    }
	    PUSH_OBJECT(value2Ptr);
	    TclDecrRefCount(namePtr);
	    TclDecrRefCount(valuePtr);
	    ADJUST_PC(1);

	case INST_INCR_SCALAR1:
	    opnd = TclGetUInt1AtPtr(pc+1);
	    valuePtr = POP_OBJECT(); 
	    if (valuePtr->typePtr != &tclIntType) {
		result = tclIntType.setFromAnyProc(interp, valuePtr);
		if (result != TCL_OK) {
		    Tcl_DecrRefCount(valuePtr);
		    goto checkForCatch;
		}
	    }
	    i = valuePtr->internalRep.longValue;
	    DECACHE_STACK_INFO();
	    value2Ptr = TclIncrIndexedScalar(interp, opnd, i);
	    CACHE_STACK_INFO();
	    if (value2Ptr == NULL) {
		Tcl_DecrRefCount(valuePtr);
		result = TCL_ERROR;
		goto checkForCatch;
	    }
	    PUSH_OBJECT(value2Ptr);
	    TclDecrRefCount(valuePtr);
	    ADJUST_PC(2);

	case INST_INCR_SCALAR_STK:
	case INST_INCR_STK:
	    valuePtr = POP_OBJECT();
	    namePtr = POP_OBJECT();
	    if (valuePtr->typePtr != &tclIntType) {
		result = tclIntType.setFromAnyProc(interp, valuePtr);
		if (result != TCL_OK) {
		    Tcl_DecrRefCount(namePtr);
		    Tcl_DecrRefCount(valuePtr);
		    goto checkForCatch;
		}
	    }
	    i = valuePtr->internalRep.longValue;
	    DECACHE_STACK_INFO();
	    value2Ptr = TclIncrVar2(interp, namePtr, (Tcl_Obj *) NULL, i,
	        /*part1NotParsed*/ (opCode == INST_INCR_STK));
	    CACHE_STACK_INFO();
	    if (value2Ptr == NULL) {
		Tcl_DecrRefCount(namePtr);
		Tcl_DecrRefCount(valuePtr);
		result = TCL_ERROR;
		goto checkForCatch;
	    }
	    PUSH_OBJECT(value2Ptr);
	    Tcl_DecrRefCount(namePtr);
	    Tcl_DecrRefCount(valuePtr);
	    ADJUST_PC(1);

	case INST_INCR_ARRAY1:
	    {
		Tcl_Obj *elemPtr;

		opnd = TclGetUInt1AtPtr(pc+1);
		valuePtr = POP_OBJECT();
		elemPtr = POP_OBJECT();
		if (valuePtr->typePtr != &tclIntType) {
		    result = tclIntType.setFromAnyProc(interp, valuePtr);
		    if (result != TCL_OK) {
			Tcl_DecrRefCount(elemPtr);
			Tcl_DecrRefCount(valuePtr);
			goto checkForCatch;
		    }
		}
		i = valuePtr->internalRep.longValue;
		DECACHE_STACK_INFO();
		value2Ptr = TclIncrElementOfIndexedArray(interp, opnd,
		        elemPtr, i);
		CACHE_STACK_INFO();
		if (value2Ptr == NULL) {
		    Tcl_DecrRefCount(elemPtr);
		    Tcl_DecrRefCount(valuePtr);
		    result = TCL_ERROR;
		    goto checkForCatch;
		}
		PUSH_OBJECT(value2Ptr);
		Tcl_DecrRefCount(elemPtr);
		Tcl_DecrRefCount(valuePtr);
	    }
	    ADJUST_PC(2);
	    
	case INST_INCR_ARRAY_STK:
	    {
		Tcl_Obj *elemPtr;

		valuePtr = POP_OBJECT();
		elemPtr = POP_OBJECT();
		namePtr = POP_OBJECT();
		if (valuePtr->typePtr != &tclIntType) {
		    result = tclIntType.setFromAnyProc(interp, valuePtr);
		    if (result != TCL_OK) {
			Tcl_DecrRefCount(namePtr);
			Tcl_DecrRefCount(elemPtr);
			Tcl_DecrRefCount(valuePtr);
			goto checkForCatch;
		    }
		}
		i = valuePtr->internalRep.longValue;
		DECACHE_STACK_INFO();
		value2Ptr = TclIncrVar2(interp, namePtr, elemPtr, i,
					/*part1NotParsed*/ 0);
		CACHE_STACK_INFO();
		if (value2Ptr == NULL) {
		    Tcl_DecrRefCount(namePtr);
		    Tcl_DecrRefCount(elemPtr);
		    Tcl_DecrRefCount(valuePtr);
		    result = TCL_ERROR;
		    goto checkForCatch;
		}
		PUSH_OBJECT(value2Ptr);
		Tcl_DecrRefCount(namePtr);
		Tcl_DecrRefCount(elemPtr);
		Tcl_DecrRefCount(valuePtr);
	    }
	    ADJUST_PC(1);
	    
	case INST_INCR_SCALAR1_IMM:
	    opnd = TclGetUInt1AtPtr(pc+1);
	    i = TclGetInt1AtPtr(pc+2);
	    DECACHE_STACK_INFO();
	    value2Ptr = TclIncrIndexedScalar(interp, opnd, i);
	    CACHE_STACK_INFO();
	    if (value2Ptr == NULL) {
		result = TCL_ERROR;
		goto checkForCatch;
	    }
	    PUSH_OBJECT(value2Ptr);
	    ADJUST_PC(3);

	case INST_INCR_SCALAR_STK_IMM:
	case INST_INCR_STK_IMM:
	    namePtr = POP_OBJECT();
	    i = TclGetInt1AtPtr(pc+1);
	    DECACHE_STACK_INFO();
	    value2Ptr = TclIncrVar2(interp, namePtr, (Tcl_Obj *) NULL, i,
		    /*part1NotParsed*/ (opCode == INST_INCR_STK_IMM));
	    CACHE_STACK_INFO();
	    if (value2Ptr == NULL) {
		result = TCL_ERROR;
		Tcl_DecrRefCount(namePtr);
		goto checkForCatch;
	    }
	    PUSH_OBJECT(value2Ptr);
	    TclDecrRefCount(namePtr);
	    ADJUST_PC(2);

	case INST_INCR_ARRAY1_IMM:
	    {
		Tcl_Obj *elemPtr;

		opnd = TclGetUInt1AtPtr(pc+1);
		i = TclGetInt1AtPtr(pc+2);
		elemPtr = POP_OBJECT();
		DECACHE_STACK_INFO();
		value2Ptr = TclIncrElementOfIndexedArray(interp, opnd,
		        elemPtr, i);
		CACHE_STACK_INFO();
		if (value2Ptr == NULL) {
		    Tcl_DecrRefCount(elemPtr);
		    result = TCL_ERROR;
		    goto checkForCatch;
		}
		PUSH_OBJECT(value2Ptr);
		Tcl_DecrRefCount(elemPtr);
	    }
	    ADJUST_PC(3);
	    
	case INST_INCR_ARRAY_STK_IMM:
	    {
		Tcl_Obj *elemPtr;

		i = TclGetInt1AtPtr(pc+1);
		elemPtr = POP_OBJECT();
		namePtr = POP_OBJECT();
		DECACHE_STACK_INFO();
		value2Ptr = TclIncrVar2(interp, namePtr, elemPtr, i,
		        /*part1NotParsed*/ 0);
		CACHE_STACK_INFO();
		if (value2Ptr == NULL) {
		    Tcl_DecrRefCount(namePtr);
		    Tcl_DecrRefCount(elemPtr);
		    result = TCL_ERROR;
		    goto checkForCatch;
		}
		PUSH_OBJECT(value2Ptr);
		Tcl_DecrRefCount(namePtr);
		Tcl_DecrRefCount(elemPtr);
	    }
	    ADJUST_PC(2);

	case INST_JUMP1:
	    opnd = TclGetInt1AtPtr(pc+1);
	    ADJUST_PC(opnd);

	case INST_JUMP4:
	    opnd = TclGetInt4AtPtr(pc+1);
	    ADJUST_PC(opnd);

	case INST_JUMP_TRUE4:
	    opnd = TclGetInt4AtPtr(pc+1);
	    pcAdjustment = 5;
	    goto doJumpTrue;

	case INST_JUMP_TRUE1:
	    opnd = TclGetInt1AtPtr(pc+1);
	    pcAdjustment = 2;
	    
	    doJumpTrue:
	    {
		int b;
		
		valuePtr = POP_OBJECT();
		if (valuePtr->typePtr == &tclIntType) {
		    b = (valuePtr->internalRep.longValue != 0);
		} else if (valuePtr->typePtr == &tclDoubleType) {
		    b = (valuePtr->internalRep.doubleValue != 0.0);
		} else {
		    result = Tcl_GetBooleanFromObj(interp, valuePtr, &b);
		    if (result != TCL_OK) {
			Tcl_DecrRefCount(valuePtr);
			goto checkForCatch;
		    }
		}
		if (b) {
		    TclDecrRefCount(valuePtr);
		    ADJUST_PC(opnd);
		} else {
		    TclDecrRefCount(valuePtr);
		    ADJUST_PC(pcAdjustment);
		}
	    }
	    
	case INST_JUMP_FALSE4:
	    opnd = TclGetInt4AtPtr(pc+1);
	    pcAdjustment = 5;
	    goto doJumpFalse;

	case INST_JUMP_FALSE1:
	    opnd = TclGetInt1AtPtr(pc+1);
	    pcAdjustment = 2;
	    
	    doJumpFalse:
	    {
		int b;
		
		valuePtr = POP_OBJECT();
		if (valuePtr->typePtr == &tclIntType) {
		    b = (valuePtr->internalRep.longValue != 0);
		} else if (valuePtr->typePtr == &tclDoubleType) {
		    b = (valuePtr->internalRep.doubleValue != 0.0);
		} else {
		    result = Tcl_GetBooleanFromObj(interp, valuePtr, &b);
		    if (result != TCL_OK) {
			Tcl_DecrRefCount(valuePtr);
			goto checkForCatch;
		    }
		}
		if (b) {
		    TclDecrRefCount(valuePtr);
		    ADJUST_PC(pcAdjustment);
		} else {
		    TclDecrRefCount(valuePtr);
		    ADJUST_PC(opnd);
		}
	    }
	    
	case INST_LOR:
	case INST_LAND:
	    {
		/*
		 * Operands must be boolean or numeric. No int->double
		 * conversions are performed.
		 */
		
		int i1, i2;
		int iResult;
		char *s;
		Tcl_ObjType *t1Ptr, *t2Ptr;
		
		value2Ptr = POP_OBJECT();
		valuePtr  = POP_OBJECT();
		t1Ptr = valuePtr->typePtr;
		t2Ptr = value2Ptr->typePtr;
		
		if ((t1Ptr == &tclIntType) || (t1Ptr == &tclBooleanType)) {
		    i1 = (valuePtr->internalRep.longValue != 0);
		} else if (t1Ptr == &tclDoubleType) {
		    i1 = (valuePtr->internalRep.doubleValue != 0.0);
		} else {	/* FAILS IF NULL STRING REP */
		    s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
		    if (TclLooksLikeInt(s)) {
			result = Tcl_GetLongFromObj((Tcl_Interp *) NULL,
				valuePtr, &i);
			i1 = (i != 0);
		    } else {
			result = Tcl_GetBooleanFromObj((Tcl_Interp *) NULL,
				valuePtr, &i1);
			i1 = (i1 != 0);
		    }
		    if (result != TCL_OK) {
			IllegalExprOperandType(interp, opCode, valuePtr);
			Tcl_DecrRefCount(valuePtr);
			Tcl_DecrRefCount(value2Ptr);
			goto checkForCatch;
		    }
		}
		
		if ((t2Ptr == &tclIntType) || (t2Ptr == &tclBooleanType)) {
		    i2 = (value2Ptr->internalRep.longValue != 0);
		} else if (t2Ptr == &tclDoubleType) {
		    i2 = (value2Ptr->internalRep.doubleValue != 0.0);
		} else {	/* FAILS IF NULL STRING REP */
		    s = Tcl_GetStringFromObj(value2Ptr, (int *) NULL);
		    if (TclLooksLikeInt(s)) {
			result = Tcl_GetLongFromObj((Tcl_Interp *) NULL,
				value2Ptr, &i);
			i2 = (i != 0);
		    } else {
			result = Tcl_GetBooleanFromObj((Tcl_Interp *) NULL,
				value2Ptr, &i2);
			i2 = (i2 != 0);
		    }
		    if (result != TCL_OK) {
			IllegalExprOperandType(interp, opCode, value2Ptr);
			Tcl_DecrRefCount(valuePtr);
			Tcl_DecrRefCount(value2Ptr);
			goto checkForCatch;
		    }
		}
		
		/*
		 * Reuse the valuePtr object already on stack if possible.
		 */

		if (opCode == INST_LOR) {
		    iResult = (i1 || i2);
		} else {
		    iResult = (i1 && i2);
		}
		if (Tcl_IsShared(valuePtr)) {
		    PUSH_OBJECT(Tcl_NewLongObj(iResult));
		    TclDecrRefCount(valuePtr);
		} else {	/* reuse the valuePtr object */
		    Tcl_SetLongObj(valuePtr, iResult);
		    ++stackTop; /* valuePtr now on stk top has right r.c. */
		}
		TclDecrRefCount(value2Ptr);
	    }
	    ADJUST_PC(1);

	case INST_EQ:
	case INST_NEQ:
	case INST_LT:
	case INST_GT:
	case INST_LE:
	case INST_GE:
	    {
		/*
		 * Any type is allowed but the two operands must have the
	         * same type. We will compute value op value2.
		 */

		Tcl_ObjType *t1Ptr, *t2Ptr;
		char *s1 = NULL;   /* Init. avoids compiler warning. */
		char *s2 = NULL;   /* Init. avoids compiler warning. */
		long i2 = 0;	   /* Init. avoids compiler warning. */
		double d1 = 0.0;   /* Init. avoids compiler warning. */
		double d2 = 0.0;   /* Init. avoids compiler warning. */
		long iResult = 0;  /* Init. avoids compiler warning. */

		value2Ptr = POP_OBJECT();
		valuePtr  = POP_OBJECT();
		t1Ptr = valuePtr->typePtr;
		t2Ptr = value2Ptr->typePtr;
		
		if ((t1Ptr != &tclIntType) && (t1Ptr != &tclDoubleType)) {
		    s1 = Tcl_GetStringFromObj(valuePtr, &length);
		    if (TclLooksLikeInt(s1)) { /* FAILS IF NULLS */
			(void) Tcl_GetLongFromObj((Tcl_Interp *) NULL,
				valuePtr, &i);
		    } else {
			(void) Tcl_GetDoubleFromObj((Tcl_Interp *) NULL,
				valuePtr, &d1);
		    }
		    t1Ptr = valuePtr->typePtr;
		}
		if ((t2Ptr != &tclIntType) && (t2Ptr != &tclDoubleType)) {
		    s2 = Tcl_GetStringFromObj(value2Ptr, &length);
		    if (TclLooksLikeInt(s2)) { /* FAILS IF NULLS */
			(void) Tcl_GetLongFromObj((Tcl_Interp *) NULL,
				value2Ptr, &i2);
		    } else {
			(void) Tcl_GetDoubleFromObj((Tcl_Interp *) NULL,
				value2Ptr, &d2);
		    }
		    t2Ptr = value2Ptr->typePtr;
		}

		if (((t1Ptr != &tclIntType) && (t1Ptr != &tclDoubleType))
		        || ((t2Ptr != &tclIntType) && (t2Ptr != &tclDoubleType))) {
		    /*
		     * One operand is not numeric. Compare as strings.
		     * THIS FAILS IF AN OBJECT'S STRING REP CONTAINS NULLS.
		     */
		    int cmpValue;
		    s1 = TclGetStringFromObj(valuePtr, &length);
		    s2 = TclGetStringFromObj(value2Ptr, &length);
		    cmpValue = strcmp(s1, s2);
		    switch (opCode) {
		    case INST_EQ:
			iResult = (cmpValue == 0);
			break;
		    case INST_NEQ:
			iResult = (cmpValue != 0);
			break;
		    case INST_LT:
			iResult = (cmpValue < 0);
			break;
		    case INST_GT:
			iResult = (cmpValue > 0);
			break;
		    case INST_LE:
			iResult = (cmpValue <= 0);
			break;
		    case INST_GE:
			iResult = (cmpValue >= 0);
			break;
		    }
		} else if ((t1Ptr == &tclDoubleType)
		        || (t2Ptr == &tclDoubleType)) {
		    /*
		     * Compare as doubles.
		     */
		    if (t1Ptr == &tclDoubleType) {
			d1 = valuePtr->internalRep.doubleValue;
			if (t2Ptr == &tclIntType) {
			    d2 = value2Ptr->internalRep.longValue;
			} else {
			    d2 = value2Ptr->internalRep.doubleValue;
			}
		    } else {	/* t1Ptr is int, t2Ptr is double */
			d1 = valuePtr->internalRep.longValue;
			d2 = value2Ptr->internalRep.doubleValue;
		    }
		    switch (opCode) {
		    case INST_EQ:
			iResult = d1 == d2;
			break;
		    case INST_NEQ:
			iResult = d1 != d2;
			break;
		    case INST_LT:
			iResult = d1 < d2;
			break;
		    case INST_GT:
			iResult = d1 > d2;
			break;
		    case INST_LE:
			iResult = d1 <= d2;
			break;
		    case INST_GE:
			iResult = d1 >= d2;
			break;
		    }
		} else {
		    /*
		     * Compare as ints.
		     */
		    i  = valuePtr->internalRep.longValue;
		    i2 = value2Ptr->internalRep.longValue;
		    switch (opCode) {
		    case INST_EQ:
			iResult = i == i2;
			break;
		    case INST_NEQ:
			iResult = i != i2;
			break;
		    case INST_LT:
			iResult = i < i2;
			break;
		    case INST_GT:
			iResult = i > i2;
			break;
		    case INST_LE:
			iResult = i <= i2;
			break;
		    case INST_GE:
			iResult = i >= i2;
			break;
		    }
		}

		/*
		 * Reuse the valuePtr object already on stack if possible.
		 */
		
		if (Tcl_IsShared(valuePtr)) {
		    PUSH_OBJECT(Tcl_NewLongObj(iResult));
		    TclDecrRefCount(valuePtr);
		} else {	/* reuse the valuePtr object */
		    Tcl_SetLongObj(valuePtr, iResult);
		    ++stackTop; /* valuePtr now on stk top has right r.c. */
		}
		TclDecrRefCount(value2Ptr);
	    }
	    ADJUST_PC(1);
	    
	case INST_MOD:
	case INST_LSHIFT:
	case INST_RSHIFT:
	case INST_BITOR:
	case INST_BITXOR:
	case INST_BITAND:
	    {
		/*
		 * Only integers are allowed. We compute value op value2.
		 */

		long i2, rem, negative;
		long iResult = 0; /* Init. avoids compiler warning. */
		
		value2Ptr = POP_OBJECT();
		valuePtr  = POP_OBJECT(); 
		if (valuePtr->typePtr == &tclIntType) {
		    i = valuePtr->internalRep.longValue;
		} else {	/* try to convert to int */
		    result = Tcl_GetLongFromObj((Tcl_Interp *) NULL,
			    valuePtr, &i);
		    if (result != TCL_OK) {
			IllegalExprOperandType(interp, opCode, valuePtr);
			Tcl_DecrRefCount(valuePtr);
			Tcl_DecrRefCount(value2Ptr);
			goto checkForCatch;
		    }
		}
		if (value2Ptr->typePtr == &tclIntType) {
		    i2 = value2Ptr->internalRep.longValue;
		} else {
		    result = Tcl_GetLongFromObj((Tcl_Interp *) NULL,
			    value2Ptr, &i2);
		    if (result != TCL_OK) {
			IllegalExprOperandType(interp, opCode, value2Ptr);
			Tcl_DecrRefCount(valuePtr);
			Tcl_DecrRefCount(value2Ptr);
			goto checkForCatch;
		    }
		}

		switch (opCode) {
		case INST_MOD:
		    /*
		     * This code is tricky: C doesn't guarantee much about
		     * the quotient or remainder, but Tcl does. The
		     * remainder always has the same sign as the divisor and
		     * a smaller absolute value.
		     */
		    if (i2 == 0) {
			Tcl_DecrRefCount(valuePtr);
			Tcl_DecrRefCount(value2Ptr);
			goto divideByZero;
		    }
		    negative = 0;
		    if (i2 < 0) {
			i2 = -i2;
			i = -i;
			negative = 1;
		    }
		    rem  = i % i2;
		    if (rem < 0) {
			rem += i2;
		    }
		    if (negative) {
			rem = -rem;
		    }
		    iResult = rem;
		    break;
		case INST_LSHIFT:
		    iResult = i << i2;
		    break;
		case INST_RSHIFT:
		    /*
		     * The following code is a bit tricky: it ensures that
		     * right shifts propagate the sign bit even on machines
		     * where ">>" won't do it by default.
		     */
		    if (i < 0) {
			iResult = ~((~i) >> i2);
		    } else {
			iResult = i >> i2;
		    }
		    break;
		case INST_BITOR:
		    iResult = i | i2;
		    break;
		case INST_BITXOR:
		    iResult = i ^ i2;
		    break;
		case INST_BITAND:
		    iResult = i & i2;
		    break;
		}

		/*
		 * Reuse the valuePtr object already on stack if possible.
		 */
		
		if (Tcl_IsShared(valuePtr)) {
		    PUSH_OBJECT(Tcl_NewLongObj(iResult));
		    TclDecrRefCount(valuePtr);
		} else {	/* reuse the valuePtr object */
		    Tcl_SetLongObj(valuePtr, iResult);
		    ++stackTop; /* valuePtr now on stk top has right r.c. */
		}
		TclDecrRefCount(value2Ptr);
	    }
	    ADJUST_PC(1);
	    
	case INST_ADD:
	case INST_SUB:
	case INST_MULT:
	case INST_DIV:
	    {
		/*
		 * Operands must be numeric and ints get converted to floats
		 * if necessary. We compute value op value2.
		 */

		Tcl_ObjType *t1Ptr, *t2Ptr;
		long i2, quot, rem;
		double d1, d2;
		long iResult = 0;     /* Init. avoids compiler warning. */
		double dResult = 0.0; /* Init. avoids compiler warning. */
		int doDouble = 0;     /* 1 if doing floating arithmetic */
		
		value2Ptr = POP_OBJECT();
		valuePtr  = POP_OBJECT();
		t1Ptr = valuePtr->typePtr;
		t2Ptr = value2Ptr->typePtr;
		
		if (t1Ptr == &tclIntType) {
		    i  = valuePtr->internalRep.longValue;
		} else if (t1Ptr == &tclDoubleType) {
		    d1 = valuePtr->internalRep.doubleValue;
		} else {	     /* try to convert; FAILS IF NULLS */
		    char *s = Tcl_GetStringFromObj(valuePtr, &length);
		    if (TclLooksLikeInt(s)) {
			result = Tcl_GetLongFromObj((Tcl_Interp *) NULL,
				valuePtr, &i);
		    } else {
			result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL,
				valuePtr, &d1);
		    }
		    if (result != TCL_OK) {
			IllegalExprOperandType(interp, opCode, valuePtr);
			Tcl_DecrRefCount(valuePtr);
			Tcl_DecrRefCount(value2Ptr);
			goto checkForCatch;
		    }
		    t1Ptr = valuePtr->typePtr;
		}
		
		if (t2Ptr == &tclIntType) {
		    i2 = value2Ptr->internalRep.longValue;
		} else if (t2Ptr == &tclDoubleType) {
		    d2 = value2Ptr->internalRep.doubleValue;
		} else {	     /* try to convert; FAILS IF NULLS */
		    char *s = Tcl_GetStringFromObj(value2Ptr, &length);
		    if (TclLooksLikeInt(s)) {
			result = Tcl_GetLongFromObj((Tcl_Interp *) NULL,
				value2Ptr, &i2);
		    } else {
			result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL,
				value2Ptr, &d2);
		    }
		    if (result != TCL_OK) {
			IllegalExprOperandType(interp, opCode, value2Ptr);
			Tcl_DecrRefCount(valuePtr);
			Tcl_DecrRefCount(value2Ptr);
			goto checkForCatch;
		    }
		    t2Ptr = value2Ptr->typePtr;
		}

		if ((t1Ptr == &tclDoubleType) || (t2Ptr == &tclDoubleType)) {
		    /*
		     * Do double arithmetic.
		     */
		    doDouble = 1;
		    if (t1Ptr == &tclIntType) {
			d1 = i;       /* promote value 1 to double */
		    } else if (t2Ptr == &tclIntType) {
			d2 = i2;      /* promote value 2 to double */
		    }
		    switch (opCode) {
		    case INST_ADD:
			dResult = d1 + d2;
			break;
		    case INST_SUB:
			dResult = d1 - d2;
			break;
		    case INST_MULT:
			dResult = d1 * d2;
			break;
		    case INST_DIV:
			if (d2 == 0.0) {
			    Tcl_DecrRefCount(valuePtr);
			    Tcl_DecrRefCount(value2Ptr);
			    goto divideByZero;
			}
			dResult = d1 / d2;
			break;
		    }
		    
		    /*
		     * Check now for IEEE floating-point error.
		     */
		    
		    if (IS_NAN(dResult) || IS_INF(dResult)) {
			TclExprFloatError(interp, dResult);
			result = TCL_ERROR;
			Tcl_DecrRefCount(valuePtr);
			Tcl_DecrRefCount(value2Ptr);
			goto checkForCatch;
		    }
		} else {
		    /*
		     * Do integer arithmetic.
		     */
		    switch (opCode) {
		    case INST_ADD:
			iResult = i + i2;
			break;
		    case INST_SUB:
			iResult = i - i2;
			break;
		    case INST_MULT:
			iResult = i * i2;
			break;
		    case INST_DIV:
			/*
			 * This code is tricky: C doesn't guarantee much
			 * about the quotient or remainder, but Tcl does.
			 * The remainder always has the same sign as the
			 * divisor and a smaller absolute value.
			 */
			if (i2 == 0) {
			    Tcl_DecrRefCount(valuePtr);
			    Tcl_DecrRefCount(value2Ptr);
			    goto divideByZero;
			}
			if (i2 < 0) {
			    i2 = -i2;
			    i = -i;
			}
			quot = i / i2;
			rem  = i % i2;
			if (rem < 0) {
			    quot -= 1;
			}
			iResult = quot;
			break;
		    }
		}

		/*
		 * Reuse the valuePtr object already on stack if possible.
		 */
		
		if (Tcl_IsShared(valuePtr)) {
		    if (doDouble) {
			PUSH_OBJECT(Tcl_NewDoubleObj(dResult));
		    } else {
			PUSH_OBJECT(Tcl_NewLongObj(iResult));
		    } 
		    TclDecrRefCount(valuePtr);
		} else {	    /* reuse the valuePtr object */
		    if (doDouble) { /* NB: stack top is off by 1 */
			Tcl_SetDoubleObj(valuePtr, dResult);
		    } else {
			Tcl_SetLongObj(valuePtr, iResult);
		    }
		    ++stackTop; /* valuePtr now on stk top has right r.c. */
		}
		TclDecrRefCount(value2Ptr);
	    }
	    ADJUST_PC(1);
	    
	case INST_UPLUS:
	    {
	        /*
	         * Operand must be numeric.
	         */

		double d;
		Tcl_ObjType *tPtr;
		
		valuePtr = stackPtr[stackTop].o;
		tPtr = valuePtr->typePtr;
		if ((tPtr != &tclIntType) && (tPtr != &tclDoubleType)) {
		    char *s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
		    if (TclLooksLikeInt(s)) { /* FAILS IF NULLS */
			result = Tcl_GetLongFromObj((Tcl_Interp *) NULL,
				valuePtr, &i);
		    } else {
			result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL,
				valuePtr, &d);
		    }
		    if (result != TCL_OK) { 
			IllegalExprOperandType(interp, opCode, valuePtr);
			goto checkForCatch;
		    }
		}
	    }
	    ADJUST_PC(1);
	    
	case INST_UMINUS:
	case INST_LNOT:
	    {
		/*
		 * The operand must be numeric. If the operand object is
		 * unshared modify it directly, otherwise create a copy to
		 * modify: this is "copy on write". free any old string
		 * representation since it is now invalid.
		 */
		
		double d;
		Tcl_ObjType *tPtr;
		
		valuePtr = POP_OBJECT();
		tPtr = valuePtr->typePtr;
		if ((tPtr != &tclIntType) && (tPtr != &tclDoubleType)) {
		    char *s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
		    if (TclLooksLikeInt(s)) { /* FAILS IF NULLS */
			result = Tcl_GetLongFromObj((Tcl_Interp *) NULL,
				valuePtr, &i);
		    } else {
			result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL,
				valuePtr, &d);
		    }
		    if (result != TCL_OK) {
			IllegalExprOperandType(interp, opCode, valuePtr);
			Tcl_DecrRefCount(valuePtr);
			goto checkForCatch;
		    }
		    tPtr = valuePtr->typePtr;
		}
		
		if (Tcl_IsShared(valuePtr)) {
		    /*
		     * Create a new object.
		     */
		    if (tPtr == &tclIntType) {
			i = valuePtr->internalRep.longValue;
			objPtr = Tcl_NewLongObj(
			        (opCode == INST_UMINUS)? -i : !i);
		    } else {
			d = valuePtr->internalRep.doubleValue;
			if (opCode == INST_UMINUS) {
			    objPtr = Tcl_NewDoubleObj(-d);
			} else {
			    /*
			     * Should be able to use "!d", but apparently
			     * some compilers can't handle it.
			     */
			    objPtr = Tcl_NewLongObj((d==0.0)? 1 : 0);
			}
		    }
		    PUSH_OBJECT(objPtr);
		    TclDecrRefCount(valuePtr);
		} else {
		    /*
		     * valuePtr is unshared. Modify it directly.
		     */
		    if (tPtr == &tclIntType) {
			i = valuePtr->internalRep.longValue;
			Tcl_SetLongObj(valuePtr,
			        (opCode == INST_UMINUS)? -i : !i);
		    } else {
			d = valuePtr->internalRep.doubleValue;
			if (opCode == INST_UMINUS) {
			    Tcl_SetDoubleObj(valuePtr, -d);
			} else {
			    /*
			     * Should be able to use "!d", but apparently
			     * some compilers can't handle it.
			     */
			    Tcl_SetLongObj(valuePtr, (d==0.0)? 1 : 0);
			}
		    }
		    ++stackTop; /* valuePtr now on stk top has right r.c. */
		}
	    }
	    ADJUST_PC(1);
	    
	case INST_BITNOT:
	    {
		/*
		 * The operand must be an integer. If the operand object is
		 * unshared modify it directly, otherwise modify a copy. 
		 * Free any old string representation since it is now
		 * invalid.
		 */
		
		Tcl_ObjType *tPtr;
		
		valuePtr = POP_OBJECT();
		tPtr = valuePtr->typePtr;
		if (tPtr != &tclIntType) {
		    result = Tcl_GetLongFromObj((Tcl_Interp *) NULL,
			    valuePtr, &i);
		    if (result != TCL_OK) {   /* try to convert to double */
			IllegalExprOperandType(interp, opCode, valuePtr);
			Tcl_DecrRefCount(valuePtr);
			goto checkForCatch;
		    }
		}
		
		i = valuePtr->internalRep.longValue;
		if (Tcl_IsShared(valuePtr)) {
		    PUSH_OBJECT(Tcl_NewLongObj(~i));
		    TclDecrRefCount(valuePtr);
		} else {
		    /*
		     * valuePtr is unshared. Modify it directly.
		     */
		    Tcl_SetLongObj(valuePtr, ~i);
		    ++stackTop; /* valuePtr now on stk top has right r.c. */
		}
	    }
	    ADJUST_PC(1);
	    
	case INST_CALL_BUILTIN_FUNC1:
	    opnd = TclGetUInt1AtPtr(pc+1);
	    {
		/*
		 * Call one of the built-in Tcl math functions.
		 */

		BuiltinFunc *mathFuncPtr;

		if ((opnd < 0) || (opnd > LAST_BUILTIN_FUNC)) {
		    panic("TclExecuteByteCode: unrecognized builtin function code %d", opnd);
		}
		mathFuncPtr = &(builtinFuncTable[opnd]);
		DECACHE_STACK_INFO();
		tcl_MathInProgress++;
		result = (*mathFuncPtr->proc)(interp, eePtr,
		        mathFuncPtr->clientData);
		tcl_MathInProgress--;
		CACHE_STACK_INFO();
		if (result != TCL_OK) {
		    goto checkForCatch;
		}
	    }
	    ADJUST_PC(2);
		    
	case INST_CALL_FUNC1:
	    opnd = TclGetUInt1AtPtr(pc+1);
	    {
		/*
		 * Call a non-builtin Tcl math function previously
		 * registered by a call to Tcl_CreateMathFunc.
		 */
		
		int objc = opnd;   /* Number of arguments. The function name
				    * is the 0-th argument. */
		Tcl_Obj **objv;	   /* The array of arguments. The function
				    * name is objv[0]. */
		
		objv = &(stackPtr[stackTop - (objc-1)].o); /* "objv[0]" */
		DECACHE_STACK_INFO();
		tcl_MathInProgress++;
		result = ExprCallMathFunc(interp, eePtr, objc, objv);
		tcl_MathInProgress--;
		CACHE_STACK_INFO();
		if (result != TCL_OK) {
		    goto checkForCatch;
		}
		ADJUST_PC(2);
	    }

	case INST_TRY_CVT_TO_NUMERIC:
	    {
		/*
		 * Try to convert the topmost stack object to an int or
		 * double object. This is done in order to support Tcl's
		 * policy of interpreting operands if at all possible as
		 * first integers, else floating-point numbers.
		 */
		
		double d;
		char *s;
		Tcl_ObjType *tPtr;
		int converted, shared;

		valuePtr = stackPtr[stackTop].o;
		tPtr = valuePtr->typePtr;
		converted = 0;
		if ((tPtr != &tclIntType) && (tPtr != &tclDoubleType)) {
		    s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
		    if (TclLooksLikeInt(s)) { /* FAILS IF NULLS */
			result = Tcl_GetLongFromObj((Tcl_Interp *) NULL,
				valuePtr, &i);
		    } else {
			result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL,
				valuePtr, &d);
		    }
		    if (result == TCL_OK) {
			converted = 1;
		    }
		    result = TCL_OK; /* reset the result variable */
		    tPtr = valuePtr->typePtr;
		}

		/*
		 * Ensure that the topmost stack object, if numeric, has a
		 * string rep the same as the formatted version of its
		 * internal rep. This is used, e.g., to make sure that "expr
		 * {0001}" yields "1", not "0001". We implement this by
		 * _discarding_ the string rep since we know it will be
		 * regenerated, if needed later, by formatting the internal
		 * rep's value. Also check if there has been an IEEE
		 * floating point error.
		 */

		if ((tPtr == &tclIntType) || (tPtr == &tclDoubleType)) {
		    shared = 0;
		    if (Tcl_IsShared(valuePtr)) {
			shared = 1;
			if (tPtr == &tclIntType) {
			    i = valuePtr->internalRep.longValue;
			    objPtr = Tcl_NewLongObj(i);
			} else {
			    d = valuePtr->internalRep.doubleValue;
			    objPtr = Tcl_NewDoubleObj(d);
			}
			Tcl_IncrRefCount(objPtr);
			TclDecrRefCount(valuePtr);
			valuePtr = objPtr;
			tPtr = valuePtr->typePtr;
		    } else {
			Tcl_InvalidateStringRep(valuePtr);
		    }
		    stackPtr[stackTop].o = valuePtr;
		
		    if (tPtr == &tclDoubleType) {
			d = valuePtr->internalRep.doubleValue;
			if (IS_NAN(d) || IS_INF(d)) {
			    TclExprFloatError(interp, d);
			    result = TCL_ERROR;
			    goto checkForCatch;
			}
		    }
		    shared = shared;		/* lint, shared not used. */
		    converted = converted;	/* lint, converted not used. */
		}
	    }
	    ADJUST_PC(1);

	case INST_BREAK:
	    /*
	     * First reset the interpreter's result. Then find the closest
	     * enclosing loop or catch exception range, if any. If a loop is
	     * found, terminate its execution. If the closest is a catch
	     * exception range, jump to its catchOffset. If no enclosing
	     * range is found, stop execution and return TCL_BREAK.
	     */

	    Tcl_ResetResult(interp);
	    rangePtr = TclGetExceptionRangeForPc(pc, /*catchOnly*/ 0,
		    codePtr);
	    if (rangePtr == NULL) {
		result = TCL_BREAK;
		goto abnormalReturn; /* no catch exists to check */
	    }
	    switch (rangePtr->type) {
	    case LOOP_EXCEPTION_RANGE:
		result = TCL_OK;
		break;
	    case CATCH_EXCEPTION_RANGE:
		result = TCL_BREAK;
		goto processCatch; /* it will use rangePtr */
	    default:
		panic("TclExecuteByteCode: unrecognized ExceptionRange type %d\n", rangePtr->type);
	    }
	    pc = (codePtr->codeStart + rangePtr->breakOffset);
	    continue;	/* restart outer instruction loop at pc */

	case INST_CONTINUE:
            /*
	     * Find the closest enclosing loop or catch exception range,
	     * if any. If a loop is found, skip to its next iteration.
	     * If the closest is a catch exception range, jump to its
	     * catchOffset. If no enclosing range is found, stop
	     * execution and return TCL_CONTINUE.
	     */

	    Tcl_ResetResult(interp);
	    rangePtr = TclGetExceptionRangeForPc(pc, /*catchOnly*/ 0,
		    codePtr);
	    if (rangePtr == NULL) {
		result = TCL_CONTINUE;
		goto abnormalReturn;
	    }
	    switch (rangePtr->type) {
	    case LOOP_EXCEPTION_RANGE:
		if (rangePtr->continueOffset == -1) {
		    goto checkForCatch;
		} else {
		    result = TCL_OK;
		}
		break;
	    case CATCH_EXCEPTION_RANGE:
		result = TCL_CONTINUE;
		goto processCatch; /* it will use rangePtr */
	    default:
		panic("TclExecuteByteCode: unrecognized ExceptionRange type %d\n", rangePtr->type);
	    }
	    pc = (codePtr->codeStart + rangePtr->continueOffset);
	    continue;	/* restart outer instruction loop at pc */

	case INST_FOREACH_START4:
	    opnd = TclGetUInt4AtPtr(pc+1);
	    {
	        /*
		 * Initialize the temporary local var that holds the count
		 * of the number of iterations of the loop body to -1.
		 */

		ForeachInfo *infoPtr = (ForeachInfo *)
		    codePtr->auxDataArrayPtr[opnd].clientData;
		int iterTmpIndex = infoPtr->loopIterNumTmp;
		CallFrame *varFramePtr = iPtr->varFramePtr;
		Var *compiledLocals = varFramePtr->compiledLocals;
		Var *iterVarPtr;
		Tcl_Obj *oldValuePtr;

		iterVarPtr = &(compiledLocals[iterTmpIndex]);
		oldValuePtr = iterVarPtr->value.objPtr;
		if (oldValuePtr == NULL) {
		    iterVarPtr->value.objPtr = Tcl_NewLongObj(-1);
		    Tcl_IncrRefCount(iterVarPtr->value.objPtr);
		} else {
		    Tcl_SetLongObj(oldValuePtr, -1);
		}
		TclSetVarScalar(iterVarPtr);
		TclClearVarUndefined(iterVarPtr);
	    }
	    ADJUST_PC(5);
	
	case INST_FOREACH_STEP4:
	    opnd = TclGetUInt4AtPtr(pc+1);
	    {
	        /*
		 * "Step" a foreach loop (i.e., begin its next iteration) by
		 * assigning the next value list element to each loop var.
		 */

		ForeachInfo *infoPtr = (ForeachInfo *)
		    codePtr->auxDataArrayPtr[opnd].clientData;
		ForeachVarList *varListPtr;
		int numLists = infoPtr->numLists;
		int iterTmpIndex = infoPtr->loopIterNumTmp;
		CallFrame *varFramePtr = iPtr->varFramePtr;
		Var *compiledLocals = varFramePtr->compiledLocals;
		int iterNum, listTmpIndex, listLen, numVars;
		int varIndex, valIndex, j;
		Tcl_Obj *listPtr, *elemPtr, *oldValuePtr;
		List *listRepPtr;
		Var *iterVarPtr, *listVarPtr;
		int continueLoop = 0;

		/*
		 * Increment the temp holding the loop iteration number.
		 */

		iterVarPtr = &(compiledLocals[iterTmpIndex]);
		oldValuePtr = iterVarPtr->value.objPtr;
		iterNum = (oldValuePtr->internalRep.longValue + 1);
		Tcl_SetLongObj(oldValuePtr, iterNum);
		
		/*
		 * Check whether all value lists are exhausted and we should
		 * stop the loop.
		 */

		listTmpIndex = infoPtr->firstListTmp;
		for (i = 0;  i < numLists;  i++) {
		    varListPtr = infoPtr->varLists[i];
		    numVars = varListPtr->numVars;

		    listVarPtr = &(compiledLocals[listTmpIndex]);
		    listPtr = listVarPtr->value.objPtr;
		    result = Tcl_ListObjLength(interp, listPtr, &listLen);
		    if (result != TCL_OK) {
			goto checkForCatch;
		    }
		    if (listLen > (iterNum * numVars)) {
			continueLoop = 1;
		    }
		    listTmpIndex++;
		}

		/*
		 * If some var in some var list still has a remaining list
		 * element iterate one more time. Assign to var the next
		 * element from its value list. We already checked above
		 * that each list temp holds a valid list object.
		 */
		
		if (continueLoop) {
		    listTmpIndex = infoPtr->firstListTmp;
		    for (i = 0;  i < numLists;  i++) {
			varListPtr = infoPtr->varLists[i];
			numVars = varListPtr->numVars;

			listVarPtr = &(compiledLocals[listTmpIndex]);
			listPtr = listVarPtr->value.objPtr;
			listRepPtr = (List *)
			        listPtr->internalRep.otherValuePtr;
			listLen = listRepPtr->elemCount;
			
			valIndex = (iterNum * numVars);
			for (j = 0;  j < numVars;  j++) {
			    int setEmptyStr = 0;
			    if (valIndex >= listLen) {
				setEmptyStr = 1;
				elemPtr = Tcl_NewObj();
			    } else {
				elemPtr = listRepPtr->elements[valIndex];
			    }
			    
			    varIndex = varListPtr->varIndexes[j];
			    DECACHE_STACK_INFO();
			    value2Ptr = TclSetIndexedScalar(interp,
			           varIndex, elemPtr, /*leaveErrorMsg*/ 1);
			    CACHE_STACK_INFO();
			    if (value2Ptr == NULL) {
				if (setEmptyStr) {
				    Tcl_DecrRefCount(elemPtr); /* unneeded */
				}
				result = TCL_ERROR;
				goto checkForCatch;
			    }
			    valIndex++;
			}
			listTmpIndex++;
		    }
		}
		
		/*
		 * Now push a "1" object if at least one value list had a
		 * remaining element and the loop should continue.
		 * Otherwise push "0".
		 */

		PUSH_OBJECT(Tcl_NewLongObj(continueLoop));
	    }
	    ADJUST_PC(5);

	case INST_BEGIN_CATCH4:
	    /*
	     * Record start of the catch command with exception range index
	     * equal to the operand. Push the current stack depth onto the
	     * special catch stack.
	     */
	    catchStackPtr[++catchTop] = stackTop;
	    ADJUST_PC(5);

	case INST_END_CATCH:
	    catchTop--;
	    result = TCL_OK;
	    ADJUST_PC(1);

	case INST_PUSH_RESULT:
	    PUSH_OBJECT(Tcl_GetObjResult(interp));
	    ADJUST_PC(1);

	case INST_PUSH_RETURN_CODE:
	    PUSH_OBJECT(Tcl_NewLongObj(result));
	    ADJUST_PC(1);

	default:
	    panic("TclExecuteByteCode: unrecognized opCode %u", opCode);
	} /* end of switch on opCode */

	/*
	 * Division by zero in an expression. Control only reaches this
	 * point by "goto divideByZero".
	 */
	
        divideByZero:
	Tcl_ResetResult(interp);
	Tcl_AppendToObj(Tcl_GetObjResult(interp), "divide by zero", -1);
	Tcl_SetErrorCode(interp, "ARITH", "DIVZERO", "divide by zero",
			 (char *) NULL);
	result = TCL_ERROR;
	
	/*
	 * Execution has generated an "exception" such as TCL_ERROR. If the
	 * exception is an error, record information about what was being
	 * executed when the error occurred. Find the closest enclosing
	 * catch range, if any. If no enclosing catch range is found, stop
	 * execution and return the "exception" code.
	 */
	
        checkForCatch:
	if ((result == TCL_ERROR) && !(iPtr->flags & ERR_ALREADY_LOGGED)) {
	    RecordTracebackInfo(interp, pc, codePtr);
        }
	rangePtr = TclGetExceptionRangeForPc(pc, /*catchOnly*/ 1, codePtr);
	if (rangePtr == NULL) {
	    goto abnormalReturn;
	}

	/*
	 * A catch exception range (rangePtr) was found to handle an
	 * "exception". It was found either by checkForCatch just above or
	 * by an instruction during break, continue, or error processing.
	 * Jump to its catchOffset after unwinding the operand stack to
	 * the depth it had when starting to execute the range's catch
	 * command.
	 */

        processCatch:
	while (stackTop > catchStackPtr[catchTop]) {
	    valuePtr = POP_OBJECT();
	    TclDecrRefCount(valuePtr);
	}
	pc = (codePtr->codeStart + rangePtr->catchOffset);
	continue;		/* restart the execution loop at pc */
    } /* end of infinite loop dispatching on instructions */

    /*
     * Abnormal return code. Restore the stack to state it had when starting
     * to execute the ByteCode.
     */

    abnormalReturn:
    while (stackTop > initStackTop) {
	valuePtr = POP_OBJECT();
	Tcl_DecrRefCount(valuePtr);
    }

    /*
     * Free the catch stack array if malloc'ed storage was used.
     */

    done:
    if (catchStackPtr != catchStackStorage) {
	ckfree((char *) catchStackPtr);
    }
    eePtr->stackTop = initStackTop;
    return result;
#undef STATIC_CATCH_STACK_SIZE
}

/*
 *----------------------------------------------------------------------
 *
 * IllegalExprOperandType --
 *
 *	Used by TclExecuteByteCode to add an error message to errorInfo
 *	when an illegal operand type is detected by an expression
 *	instruction. The argument opCode holds the failing instruction's
 *	opcode and opndPtr holds the operand object in error.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	An error message is appended to errorInfo.
 *
 *----------------------------------------------------------------------
 */

static void
IllegalExprOperandType(interp, opCode, opndPtr)
    Tcl_Interp *interp;		/* Interpreter to which error information
				 * pertains. */
    unsigned int opCode;	/* The instruction opcode being executed
				 * when the illegal type was found. */
    Tcl_Obj *opndPtr;		/* Points to the operand holding the value
				 * with the illegal type. */
{
    Tcl_ResetResult(interp);
    if ((opndPtr->bytes == NULL) || (opndPtr->length == 0)) {
	Tcl_AppendStringsToObj(Tcl_GetObjResult(interp),
		"can't use empty string as operand of \"",
		operatorStrings[opCode - INST_LOR], "\"", (char *) NULL);
    } else {
	Tcl_AppendStringsToObj(Tcl_GetObjResult(interp), "can't use ",
		((opndPtr->typePtr == &tclDoubleType) ?
		    "floating-point value" : "non-numeric string"),
		" as operand of \"", operatorStrings[opCode - INST_LOR],
		"\"", (char *) NULL);
    }
}

/*
 *----------------------------------------------------------------------
 *
 * CallTraceProcedure --
 *
 *	Invokes a trace procedure registered with an interpreter. These
 *	procedures trace command execution. Currently this trace procedure
 *	is called with the address of the string-based Tcl_CmdProc for the
 *	command, not the Tcl_ObjCmdProc.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Those side effects made by the trace procedure.
 *
 *----------------------------------------------------------------------
 */

static void
CallTraceProcedure(interp, tracePtr, cmdPtr, command, numChars, objc, objv)
    Tcl_Interp *interp;		/* The current interpreter. */
    register Trace *tracePtr;	/* Describes the trace procedure to call. */
    Command *cmdPtr;		/* Points to command's Command struct. */
    char *command;		/* Points to the first character of the
				 * command's source before substitutions. */
    int numChars;		/* The number of characters in the
				 * command's source. */
    register int objc;		/* Number of arguments for the command. */
    Tcl_Obj *objv[];		/* Pointers to Tcl_Obj of each argument. */
{
    Interp *iPtr = (Interp *) interp;
    register char **argv;
    register int i;
    int length;
    char *p;

    /*
     * Get the string rep from the objv argument objects and place their
     * pointers in argv. First make sure argv is large enough to hold the
     * objc args plus 1 extra word for the zero end-of-argv word.
     * THIS FAILS IF AN OBJECT'S STRING REP CONTAINS NULLS.
     */
    
    argv = (char **) ckalloc((unsigned)(objc + 1) * sizeof(char *));
    for (i = 0;  i < objc;  i++) {
	argv[i] = Tcl_GetStringFromObj(objv[i], &length);
    }
    argv[objc] = 0;

    /*
     * Copy the command characters into a new string.
     */

    p = (char *) ckalloc((unsigned) (numChars + 1));
    memcpy((VOID *) p, (VOID *) command, (size_t) numChars);
    p[numChars] = '\0';
    
    /*
     * Call the trace procedure then free allocated storage.
     */
    
    (*tracePtr->proc)(tracePtr->clientData, interp, iPtr->numLevels,
                      p, cmdPtr->proc, cmdPtr->clientData, objc, argv);

    ckfree((char *) argv);
    ckfree((char *) p);
}

/*
 *----------------------------------------------------------------------
 *
 * RecordTracebackInfo --
 *
 *      Procedure called by TclExecuteByteCode to record information
 *      about what was being executed when the error occurred.
 *
 * Results:
 *      None.
 *
 * Side effects:
 *      Appends information about the command being executed to the
 *      "errorInfo" variable. Sets the errorLine field in the interpreter
 *      to the line number of that command. Sets the ERR_ALREADY_LOGGED
 *      bit in the interpreter's execution flags.
 *
 *----------------------------------------------------------------------
 */

static void
RecordTracebackInfo(interp, pc, codePtr)
    Tcl_Interp *interp;         /* The interpreter in which the error
                                 * occurred. */
    unsigned char *pc;          /* The program counter value where the error                                 * occurred. This points to a bytecode
                                 * instruction in codePtr's code. */
    ByteCode *codePtr;          /* The bytecode sequence being executed. */
{
    register Interp *iPtr = (Interp *) interp;
    char *cmd, *ellipsis;
    char buf[200];
    register char *p;
    int numChars;
    
    /*
     * Record the command in errorInfo (up to a certain number of
     * characters, or up to the first newline).
     */
    
    iPtr->errorLine = 1;
    cmd = GetSrcInfoForPc(pc, codePtr, &numChars);
    if (cmd != NULL) {
        for (p = codePtr->source;  p != cmd;  p++) {
            if (*p == '\n') {
                iPtr->errorLine++;
            }
        }
        for ( ;  (isspace(UCHAR(*p)) || (*p == ';'));  p++) {
            if (*p == '\n') {
                iPtr->errorLine++;
            }
        }
	
        ellipsis = "";
        if (numChars > 150) {
            numChars = 150;
            ellipsis = "...";
        }
        if (!(iPtr->flags & ERR_IN_PROGRESS)) {
            sprintf(buf, "\n    while executing\n\"%.*s%s\"",
                    numChars, cmd, ellipsis);
        } else {
            sprintf(buf, "\n    invoked from within\n\"%.*s%s\"",
                    numChars, cmd, ellipsis);
        }
        Tcl_AddObjErrorInfo(interp, buf, -1);
        iPtr->flags |= ERR_ALREADY_LOGGED;
    }
}

/*
 *----------------------------------------------------------------------
 *
 * GetSrcInfoForPc --
 *
 *	Given a program counter value, finds the closest command in the
 *	bytecode code unit's CmdLocation array and returns information about
 *	that command's source: a pointer to its first byte and the number of
 *	characters.
 *
 * Results:
 *	If a command is found that encloses the program counter value, a
 *	pointer to the command's source is returned and the length of the
 *	source is stored at *lengthPtr. If multiple commands resulted in
 *	code at pc, information about the closest enclosing command is
 *	returned. If no matching command is found, NULL is returned and
 *	*lengthPtr is unchanged.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

static char *
GetSrcInfoForPc(pc, codePtr, lengthPtr)
    unsigned char *pc;		/* The program counter value for which to
				 * return the closest command's source info.
				 * This points to a bytecode instruction
				 * in codePtr's code. */
    ByteCode *codePtr;		/* The bytecode sequence in which to look
				 * up the command source for the pc. */
    int *lengthPtr;		/* If non-NULL, the location where the
				 * length of the command's source should be
				 * stored. If NULL, no length is stored. */
{
    register int pcOffset = (pc - codePtr->codeStart);
    int numCmds = codePtr->numCommands;
    unsigned char *codeDeltaNext, *codeLengthNext;
    unsigned char *srcDeltaNext, *srcLengthNext;
    int codeOffset, codeLen, codeEnd, srcOffset, srcLen, delta, i;
    int bestDist = INT_MAX;	/* Distance of pc to best cmd's start pc. */
    int bestSrcOffset = -1;	/* Initialized to avoid compiler warning. */
    int bestSrcLength = -1;	/* Initialized to avoid compiler warning. */

    if ((pcOffset < 0) || (pcOffset >= codePtr->numCodeBytes)) {
	return NULL;
    }

    /*
     * Decode the code and source offset and length for each command. The
     * closest enclosing command is the last one whose code started before
     * pcOffset.
     */

    codeDeltaNext = codePtr->codeDeltaStart;
    codeLengthNext = codePtr->codeLengthStart;
    srcDeltaNext  = codePtr->srcDeltaStart;
    srcLengthNext = codePtr->srcLengthStart;
    codeOffset = srcOffset = 0;
    for (i = 0;  i < numCmds;  i++) {
	if ((unsigned int) (*codeDeltaNext) == (unsigned int) 0xFF) {
	    codeDeltaNext++;
	    delta = TclGetInt4AtPtr(codeDeltaNext);
	    codeDeltaNext += 4;
	} else {
	    delta = TclGetInt1AtPtr(codeDeltaNext);
	    codeDeltaNext++;
	}
	codeOffset += delta;

	if ((unsigned int) (*codeLengthNext) == (unsigned int) 0xFF) {
	    codeLengthNext++;
	    codeLen = TclGetInt4AtPtr(codeLengthNext);
	    codeLengthNext += 4;
	} else {
	    codeLen = TclGetInt1AtPtr(codeLengthNext);
	    codeLengthNext++;
	}
	codeEnd = (codeOffset + codeLen - 1);

	if ((unsigned int) (*srcDeltaNext) == (unsigned int) 0xFF) {
	    srcDeltaNext++;
	    delta = TclGetInt4AtPtr(srcDeltaNext);
	    srcDeltaNext += 4;
	} else {
	    delta = TclGetInt1AtPtr(srcDeltaNext);
	    srcDeltaNext++;
	}
	srcOffset += delta;

	if ((unsigned int) (*srcLengthNext) == (unsigned int) 0xFF) {
	    srcLengthNext++;
	    srcLen = TclGetInt4AtPtr(srcLengthNext);
	    srcLengthNext += 4;
	} else {
	    srcLen = TclGetInt1AtPtr(srcLengthNext);
	    srcLengthNext++;
	}
	
	if (codeOffset > pcOffset) {      /* best cmd already found */
	    break;
	} else if (pcOffset <= codeEnd) { /* this cmd's code encloses pc */
	    int dist = (pcOffset - codeOffset);
	    if (dist <= bestDist) {
		bestDist = dist;
		bestSrcOffset = srcOffset;
		bestSrcLength = srcLen;
	    }
	}
    }

    if (bestDist == INT_MAX) {
	return NULL;
    }
    
    if (lengthPtr != NULL) {
	*lengthPtr = bestSrcLength;
    }
    return (codePtr->source + bestSrcOffset);
}

/*
 *----------------------------------------------------------------------
 *
 * TclGetExceptionRangeForPc --
 *
 *	Procedure that given a program counter value, returns the closest
 *	enclosing ExceptionRange that matches the kind requested.
 *
 * Results:
 *	In the normal case, catchOnly is 0 (false) and this procedure
 *	returns a pointer to the most closely enclosing ExceptionRange
 *	structure regardless of whether it is a loop or catch exception
 *	range. This is appropriate when processing a TCL_BREAK or
 *	TCL_CONTINUE, which will be "handled" either by a loop exception
 *	range or a closer catch range. If catchOnly is nonzero (true), this
 *	procedure ignores loop exception ranges and returns a pointer to the
 *	closest catch range. If no matching ExceptionRange is found that
 *	encloses pc, a NULL is returned.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

ExceptionRange *
TclGetExceptionRangeForPc(pc, catchOnly, codePtr)
    unsigned char *pc;		/* The program counter value for which to
				 * search for a closest enclosing exception
				 * range. This points to a bytecode
				 * instruction in codePtr's code. */
    int catchOnly;		/* If 0, consider either loop or catch
				 * ExceptionRanges in search. Otherwise
				 * consider only catch ranges (and ignore
				 * any closer loop ranges). */
    ByteCode* codePtr;		/* Points to the ByteCode in which to search
				 * for the enclosing ExceptionRange. */
{
    ExceptionRange *rangeArrayPtr;
    int numRanges = codePtr->numExcRanges;
    register ExceptionRange *rangePtr;
    int codeOffset = (pc - codePtr->codeStart);
    register int i, level;

    if (numRanges == 0) {
	return NULL;
    }
    rangeArrayPtr = codePtr->excRangeArrayPtr;

    for (level = codePtr->maxExcRangeDepth;  level >= 0;  level--) {
	for (i = 0;  i < numRanges;  i++) {
	    rangePtr = &(rangeArrayPtr[i]);
	    if (rangePtr->nestingLevel == level) {
		int start = rangePtr->codeOffset;
		int end   = (start + rangePtr->numCodeBytes);
		if ((start <= codeOffset) && (codeOffset < end)) {
		    if ((!catchOnly)
			    || (rangePtr->type == CATCH_EXCEPTION_RANGE)) {
			return rangePtr;
		    }
		}
	    }
	}
    }
    return NULL;
}

/*
 *----------------------------------------------------------------------
 *
 * Math Functions --
 *
 *	This page contains the procedures that implement all of the
 *	built-in math functions for expressions.
 *
 * Results:
 *	Each procedure returns TCL_OK if it succeeds and pushes an
 *	Tcl object holding the result. If it fails it returns TCL_ERROR
 *	and leaves an error message in the interpreter's result.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

static int
ExprUnaryFunc(interp, eePtr, clientData)
    Tcl_Interp *interp;		/* The interpreter in which to execute the
				 * function. */
    ExecEnv *eePtr;		/* Points to the environment for executing
				 * the function. */
    ClientData clientData;	/* Contains the address of a procedure that
				 * takes one double argument and returns a
				 * double result. */
{
    StackItem *stackPtr;        /* Cached evaluation stack base pointer. */
    register int stackTop;	/* Cached top index of evaluation stack. */
    register Tcl_Obj *valuePtr;
    Tcl_ObjType *tPtr;
    double d, dResult;
    long i;
    int result = TCL_OK;
    
    double (*func) _ANSI_ARGS_((double)) =
	(double (*)_ANSI_ARGS_((double))) clientData;

    /*
     * Set stackPtr and stackTop from eePtr.
     */
    
    CACHE_STACK_INFO();

    /*
     * Pop the function's argument from the evaluation stack. Convert it
     * to a double if necessary.
     */

    valuePtr = POP_OBJECT();
    tPtr = valuePtr->typePtr;
    
    if (tPtr == &tclIntType) {
	d = (double) valuePtr->internalRep.longValue;
    } else if (tPtr == &tclDoubleType) {
	d = valuePtr->internalRep.doubleValue;
    } else {			/* FAILS IF STRING REP HAS NULLS */
	char *s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
	
	if (TclLooksLikeInt(s)) {
	    result = Tcl_GetLongFromObj((Tcl_Interp *) NULL, valuePtr, &i);
	    d = (double) valuePtr->internalRep.longValue;
	} else {
	    result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL, valuePtr, &d);
	}
	if (result != TCL_OK) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendToObj(Tcl_GetObjResult(interp),
	            "argument to math function didn't have numeric value", -1);
	    goto done;
	}
    }

    errno = 0;
    dResult = (*func)(d);
    if ((errno != 0) || IS_NAN(dResult) || IS_INF(dResult)) {
	TclExprFloatError(interp, dResult);
	result = TCL_ERROR;
	goto done;
    }
    
    /*
     * Push a Tcl object holding the result.
     */

    PUSH_OBJECT(Tcl_NewDoubleObj(dResult));
    
    /*
     * Reflect the change to stackTop back in eePtr.
     */

    done:
    Tcl_DecrRefCount(valuePtr);
    DECACHE_STACK_INFO();
    return result;
}

static int
ExprBinaryFunc(interp, eePtr, clientData)
    Tcl_Interp *interp;		/* The interpreter in which to execute the
				 * function. */
    ExecEnv *eePtr;		/* Points to the environment for executing
				 * the function. */
    ClientData clientData;	/* Contains the address of a procedure that
				 * takes two double arguments and
				 * returns a double result. */
{
    StackItem *stackPtr;        /* Cached evaluation stack base pointer. */
    register int stackTop;	/* Cached top index of evaluation stack. */
    register Tcl_Obj *valuePtr, *value2Ptr;
    Tcl_ObjType *tPtr;
    double d1, d2, dResult;
    long i;
    char *s;
    int result = TCL_OK;
    
    double (*func) _ANSI_ARGS_((double, double))
	= (double (*)_ANSI_ARGS_((double, double))) clientData;

    /*
     * Set stackPtr and stackTop from eePtr.
     */
    
    CACHE_STACK_INFO();

    /*
     * Pop the function's two arguments from the evaluation stack. Convert
     * them to doubles if necessary.
     */

    value2Ptr = POP_OBJECT();
    valuePtr  = POP_OBJECT();

    tPtr = valuePtr->typePtr;
    if (tPtr == &tclIntType) {
	d1 = (double) valuePtr->internalRep.longValue;
    } else if (tPtr == &tclDoubleType) {
	d1 = valuePtr->internalRep.doubleValue;
    } else {			/* FAILS IF STRING REP HAS NULLS */
	s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
	if (TclLooksLikeInt(s)) {
	    result = Tcl_GetLongFromObj((Tcl_Interp *) NULL, valuePtr, &i);
	    d1 = (double) valuePtr->internalRep.longValue;
	} else {
	    result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL, valuePtr, &d1);
	}
	if (result != TCL_OK) {
            badArg:
	    Tcl_ResetResult(interp);
	    Tcl_AppendToObj(Tcl_GetObjResult(interp),
	            "argument to math function didn't have numeric value", -1);
	    goto done;
	}
    }

    tPtr = value2Ptr->typePtr;
    if (tPtr == &tclIntType) {
	d2 = value2Ptr->internalRep.longValue;
    } else if (tPtr == &tclDoubleType) {
	d2 = value2Ptr->internalRep.doubleValue;
    } else {			/* FAILS IF STRING REP HAS NULLS */
	s = Tcl_GetStringFromObj(value2Ptr, (int *) NULL);
	if (TclLooksLikeInt(s)) {
	    result = Tcl_GetLongFromObj((Tcl_Interp *) NULL, value2Ptr, &i);
	    d2 = (double) value2Ptr->internalRep.longValue;
	} else {
	    result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL, value2Ptr, &d2);
	}
	if (result != TCL_OK) {
	    goto badArg;
	}
    }

    errno = 0;
    dResult = (*func)(d1, d2);
    if ((errno != 0) || IS_NAN(dResult) || IS_INF(dResult)) {
	TclExprFloatError(interp, dResult);
	result = TCL_ERROR;
	goto done;
    }

    /*
     * Push a Tcl object holding the result.
     */

    PUSH_OBJECT(Tcl_NewDoubleObj(dResult));
    
    /*
     * Reflect the change to stackTop back in eePtr.
     */

    done:
    Tcl_DecrRefCount(valuePtr);
    Tcl_DecrRefCount(value2Ptr);
    DECACHE_STACK_INFO();
    return result;
}

static int
ExprAbsFunc(interp, eePtr, clientData)
    Tcl_Interp *interp;		/* The interpreter in which to execute the
				 * function. */
    ExecEnv *eePtr;		/* Points to the environment for executing
				 * the function. */
    ClientData clientData;	/* Ignored. */
{
    StackItem *stackPtr;        /* Cached evaluation stack base pointer. */
    register int stackTop;	/* Cached top index of evaluation stack. */
    register Tcl_Obj *valuePtr;
    Tcl_ObjType *tPtr;
    long i, iResult;
    double d, dResult;
    int result = TCL_OK;

    /*
     * Set stackPtr and stackTop from eePtr.
     */
    
    CACHE_STACK_INFO();

    /*
     * Pop the argument from the evaluation stack.
     */

    valuePtr = POP_OBJECT();
    tPtr = valuePtr->typePtr;
    
    if (tPtr == &tclIntType) {
	i = valuePtr->internalRep.longValue;
    } else if (tPtr == &tclDoubleType) {
	d = valuePtr->internalRep.doubleValue;
    } else {			/* FAILS IF STRING REP HAS NULLS */
	char *s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
	
	if (TclLooksLikeInt(s)) {
	    result = Tcl_GetLongFromObj((Tcl_Interp *) NULL, valuePtr, &i);
	} else {
	    result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL, valuePtr, &d);
	}
	if (result != TCL_OK) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendToObj(Tcl_GetObjResult(interp),
	            "argument to math function didn't have numeric value", -1);
	    goto done;
	}
	tPtr = valuePtr->typePtr;
    }

    /*
     * Push a Tcl object with the result.
     */
    
    if (tPtr == &tclIntType) {
	if (i < 0) {
	    iResult = -i;
	    if (iResult < 0) {
		Tcl_ResetResult(interp);
		Tcl_AppendToObj(Tcl_GetObjResult(interp),
		        "integer value too large to represent", -1);
		Tcl_SetErrorCode(interp, "ARITH", "IOVERFLOW",
			"integer value too large to represent", (char *) NULL);
		result = TCL_ERROR;
		goto done;
	    }
	} else {
	    iResult = i;
	}	    
	PUSH_OBJECT(Tcl_NewLongObj(iResult));
    } else {
	if (d < 0.0) {
	    dResult = -d;
	} else {
	    dResult = d;
	}
	if (IS_NAN(dResult) || IS_INF(dResult)) {
	    TclExprFloatError(interp, dResult);
	    result = TCL_ERROR;
	    goto done;
	}
	PUSH_OBJECT(Tcl_NewDoubleObj(dResult));
    }
    
    /*
     * Reflect the change to stackTop back in eePtr.
     */

    done:
    Tcl_DecrRefCount(valuePtr);
    DECACHE_STACK_INFO();
    return result;
}

static int
ExprDoubleFunc(interp, eePtr, clientData)
    Tcl_Interp *interp;		/* The interpreter in which to execute the
				 * function. */
    ExecEnv *eePtr;		/* Points to the environment for executing
				 * the function. */
    ClientData clientData;	/* Ignored. */
{
    StackItem *stackPtr;        /* Cached evaluation stack base pointer. */
    register int stackTop;	/* Cached top index of evaluation stack. */
    register Tcl_Obj *valuePtr;
    double dResult;
    long i;
    int result = TCL_OK;

    /*
     * Set stackPtr and stackTop from eePtr.
     */
    
    CACHE_STACK_INFO();

    /*
     * Pop the argument from the evaluation stack.
     */

    valuePtr = POP_OBJECT();
    if (valuePtr->typePtr == &tclIntType) {
	dResult = (double) valuePtr->internalRep.longValue;
    } else if (valuePtr->typePtr == &tclDoubleType) {
	dResult = valuePtr->internalRep.doubleValue;
    } else {			/* FAILS IF STRING REP HAS NULLS */
	char *s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
	
	if (TclLooksLikeInt(s)) {
	    result = Tcl_GetLongFromObj((Tcl_Interp *) NULL, valuePtr, &i);
	    dResult = (double) valuePtr->internalRep.longValue;
	} else {
	    result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL, valuePtr,
		    &dResult);
	}
	if (result != TCL_OK) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendToObj(Tcl_GetObjResult(interp),
		    "argument to math function didn't have numeric value", -1);
	    goto done;
	}
    }

    /*
     * Push a Tcl object with the result.
     */

    PUSH_OBJECT(Tcl_NewDoubleObj(dResult));

    /*
     * Reflect the change to stackTop back in eePtr.
     */

    done:
    Tcl_DecrRefCount(valuePtr);
    DECACHE_STACK_INFO();
    return result;
}

static int
ExprIntFunc(interp, eePtr, clientData)
    Tcl_Interp *interp;		/* The interpreter in which to execute the
				 * function. */
    ExecEnv *eePtr;		/* Points to the environment for executing
				 * the function. */
    ClientData clientData;	/* Ignored. */
{
    StackItem *stackPtr;        /* Cached evaluation stack base pointer. */
    register int stackTop;	/* Cached top index of evaluation stack. */
    register Tcl_Obj *valuePtr;
    Tcl_ObjType *tPtr;
    long i = 0;			/* Initialized to avoid compiler warning. */
    long iResult;
    double d;
    int result = TCL_OK;

    /*
     * Set stackPtr and stackTop from eePtr.
     */
    
    CACHE_STACK_INFO();

    /*
     * Pop the argument from the evaluation stack.
     */

    valuePtr = POP_OBJECT();
    tPtr = valuePtr->typePtr;
    
    if (tPtr == &tclIntType) {
	i = valuePtr->internalRep.longValue;
    } else if (tPtr == &tclDoubleType) {
	d = valuePtr->internalRep.doubleValue;
    } else {			/* FAILS IF STRING REP HAS NULLS */
	char *s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
	
	if (TclLooksLikeInt(s)) {
	    result = Tcl_GetLongFromObj((Tcl_Interp *) NULL, valuePtr, &i);
	} else {
	    result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL, valuePtr, &d);
	}
	if (result != TCL_OK) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendToObj(Tcl_GetObjResult(interp),
		    "argument to math function didn't have numeric value", -1);
	    goto done;
	}
	tPtr = valuePtr->typePtr;
    }

    /*
     * Push a Tcl object with the result.
     */
    
    if (tPtr == &tclIntType) {
	iResult = i;
    } else {
	if (d < 0.0) {
	    if (d < (double) (long) LONG_MIN) {
		tooLarge:
		Tcl_ResetResult(interp);
		Tcl_AppendToObj(Tcl_GetObjResult(interp),
		        "integer value too large to represent", -1);
		Tcl_SetErrorCode(interp, "ARITH", "IOVERFLOW",
			"integer value too large to represent", (char *) NULL);
		result = TCL_ERROR;
		goto done;
	    }
	} else {
	    if (d > (double) LONG_MAX) {
		goto tooLarge;
	    }
	}
	if (IS_NAN(d) || IS_INF(d)) {
	    TclExprFloatError(interp, d);
	    result = TCL_ERROR;
	    goto done;
	}
	iResult = (long) d;
    }
    PUSH_OBJECT(Tcl_NewLongObj(iResult));

    /*
     * Reflect the change to stackTop back in eePtr.
     */

    done:
    Tcl_DecrRefCount(valuePtr);
    DECACHE_STACK_INFO();
    return result;
}

static int
ExprRoundFunc(interp, eePtr, clientData)
    Tcl_Interp *interp;		/* The interpreter in which to execute the
				 * function. */
    ExecEnv *eePtr;		/* Points to the environment for executing
				 * the function. */
    ClientData clientData;	/* Ignored. */
{
    StackItem *stackPtr;        /* Cached evaluation stack base pointer. */
    register int stackTop;	/* Cached top index of evaluation stack. */
    Tcl_Obj *valuePtr;
    Tcl_ObjType *tPtr;
    long i = 0;			/* Initialized to avoid compiler warning. */
    long iResult;
    double d, temp;
    int result = TCL_OK;

    /*
     * Set stackPtr and stackTop from eePtr.
     */
    
    CACHE_STACK_INFO();

    /*
     * Pop the argument from the evaluation stack.
     */

    valuePtr = POP_OBJECT();
    tPtr = valuePtr->typePtr;
    
    if (tPtr == &tclIntType) {
	i = valuePtr->internalRep.longValue;
    } else if (tPtr == &tclDoubleType) {
	d = valuePtr->internalRep.doubleValue;
    } else {			/* FAILS IF STRING REP HAS NULLS */
	char *s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
	
	if (TclLooksLikeInt(s)) {
	    result = Tcl_GetLongFromObj((Tcl_Interp *) NULL, valuePtr, &i);
	} else {
	    result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL, valuePtr, &d);
	}
	if (result != TCL_OK) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendToObj(Tcl_GetObjResult(interp),
		    "argument to math function didn't have numeric value", -1);
	    goto done;
	}
	tPtr = valuePtr->typePtr;
    }

    /*
     * Push a Tcl object with the result.
     */
    
    if (tPtr == &tclIntType) {
	iResult = i;
    } else {
	if (d < 0.0) {
	    if (d <= (((double) (long) LONG_MIN) - 0.5)) {
		tooLarge:
		Tcl_ResetResult(interp);
		Tcl_AppendToObj(Tcl_GetObjResult(interp),
		        "integer value too large to represent", -1);
		Tcl_SetErrorCode(interp, "ARITH", "IOVERFLOW",
			"integer value too large to represent",
			(char *) NULL);
		result = TCL_ERROR;
		goto done;
	    }
	    temp = (long) (d - 0.5);
	} else {
	    if (d >= (((double) LONG_MAX + 0.5))) {
		goto tooLarge;
	    }
	    temp = (long) (d + 0.5);
	}
	if (IS_NAN(temp) || IS_INF(temp)) {
	    TclExprFloatError(interp, temp);
	    result = TCL_ERROR;
	    goto done;
	}
	iResult = (long) temp;
    }
    PUSH_OBJECT(Tcl_NewLongObj(iResult));

    /*
     * Reflect the change to stackTop back in eePtr.
     */

    done:
    Tcl_DecrRefCount(valuePtr);
    DECACHE_STACK_INFO();
    return result;
}

/*
 *----------------------------------------------------------------------
 *
 * ExprCallMathFunc --
 *
 *	This procedure is invoked to call a non-builtin math function
 *	during the execution of an expression. 
 *
 * Results:
 *	TCL_OK is returned if all went well and the function's value
 *	was computed successfully. If an error occurred, TCL_ERROR
 *	is returned and an error message is left in the interpreter's
 *	result.	After a successful return this procedure pushes a Tcl object
 *	holding the result. 
 *
 * Side effects:
 *	None, unless the called math function has side effects.
 *
 *----------------------------------------------------------------------
 */

static int
ExprCallMathFunc(interp, eePtr, objc, objv)
    Tcl_Interp *interp;		/* The interpreter in which to execute the
				 * function. */
    ExecEnv *eePtr;		/* Points to the environment for executing
				 * the function. */
    int objc;			/* Number of arguments. The function name is
				 * the 0-th argument. */
    Tcl_Obj **objv;		/* The array of arguments. The function name
				 * is objv[0]. */
{
    Interp *iPtr = (Interp *) interp;
    StackItem *stackPtr;        /* Cached evaluation stack base pointer. */
    register int stackTop;	/* Cached top index of evaluation stack. */
    char *funcName;
    Tcl_HashEntry *hPtr;
    MathFunc *mathFuncPtr;	/* Information about math function. */
    Tcl_Value args[MAX_MATH_ARGS]; /* Arguments for function call. */
    Tcl_Value funcResult;	/* Result of function call as Tcl_Value. */
    register Tcl_Obj *valuePtr;
    Tcl_ObjType *tPtr;
    long i;
    double d;
    int j, k, result;
    
    Tcl_ResetResult(interp);
    
    /*
     * Set stackPtr and stackTop from eePtr.
     */
    
    CACHE_STACK_INFO();

    /*
     * Look up the MathFunc record for the function.
     * THIS FAILS IF THE OBJECT'S STRING REP CONTAINS NULLS.
     */

    funcName = Tcl_GetStringFromObj(objv[0], (int *) NULL);
    hPtr = Tcl_FindHashEntry(&iPtr->mathFuncTable, funcName);
    if (hPtr == NULL) {
	Tcl_AppendStringsToObj(Tcl_GetObjResult(interp),
		"unknown math function \"", funcName, "\"", (char *) NULL);
	result = TCL_ERROR;
	goto done;
    }
    mathFuncPtr = (MathFunc *) Tcl_GetHashValue(hPtr);
    if (mathFuncPtr->numArgs != (objc-1)) {
	panic("ExprCallMathFunc: expected number of args %d != actual number %d",
	        mathFuncPtr->numArgs, objc);
	result = TCL_ERROR;
	goto done;
    }

    /*
     * Collect the arguments for the function, if there are any, into the
     * array "args". Note that args[0] will have the Tcl_Value that
     * corresponds to objv[1].
     */

    for (j = 1, k = 0;  j < objc;  j++, k++) {
	valuePtr = objv[j];
	tPtr = valuePtr->typePtr;
	
	if (tPtr == &tclIntType) {
	    i = valuePtr->internalRep.longValue;
	} else if (tPtr == &tclDoubleType) {
	    d = valuePtr->internalRep.doubleValue;
	} else {
	    /*
	     * Try to convert to int first then double.
	     * FAILS IF STRING REP HAS NULLS.
	     */
	    
	    char *s = Tcl_GetStringFromObj(valuePtr, (int *) NULL);
	    
	    if (TclLooksLikeInt(s)) {
		result = Tcl_GetLongFromObj((Tcl_Interp *) NULL, valuePtr, &i);
	    } else {
		result = Tcl_GetDoubleFromObj((Tcl_Interp *) NULL,
			valuePtr, &d);
	    }
	    if (result != TCL_OK) {
		Tcl_AppendToObj(Tcl_GetObjResult(interp),
			"argument to math function didn't have numeric value", -1);
		goto done;
	    }
	    tPtr = valuePtr->typePtr;
	}

	/*
	 * Copy the object's numeric value to the argument record,
	 * converting it if necessary. 
	 */
	
	if (tPtr == &tclIntType) {
	    if (mathFuncPtr->argTypes[k] == TCL_DOUBLE) {
		args[k].type = TCL_DOUBLE;
		args[k].doubleValue = i;
	    } else {
		args[k].type = TCL_INT;
		args[k].intValue = i;
	    }
	} else {
	    if (mathFuncPtr->argTypes[k] == TCL_INT) {
		args[k].type = TCL_INT;
		args[k].intValue = (long) d;
	    } else {
		args[k].type = TCL_DOUBLE;
		args[k].doubleValue = d;
	    }
	}
    }

    /*
     * Invoke the function and copy its result back into valuePtr.
     */

    tcl_MathInProgress++;
    result = (*mathFuncPtr->proc)(mathFuncPtr->clientData, interp, args,
	    &funcResult);
    tcl_MathInProgress--;
    if (result != TCL_OK) {
	goto done;
    }

    /*
     * Pop the objc top stack elements and decrement their ref counts.
     */
		
    i = (stackTop - (objc-1));
    while (i <= stackTop) {
	valuePtr = stackPtr[i].o;
	Tcl_DecrRefCount(valuePtr);
	i++;
    }
    stackTop -= objc;
    
    /*
     * Push the call's object result.
     */
    
    if (funcResult.type == TCL_INT) {
	PUSH_OBJECT(Tcl_NewLongObj(funcResult.intValue));
    } else {
	d = funcResult.doubleValue;
	if (IS_NAN(d) || IS_INF(d)) {
	    TclExprFloatError(interp, d);
	    result = TCL_ERROR;
	    goto done;
	}
	PUSH_OBJECT(Tcl_NewDoubleObj(d));
    }

    /*
     * Reflect the change to stackTop back in eePtr.
     */

    done:
    DECACHE_STACK_INFO();
    return result;
}

/*
 *----------------------------------------------------------------------
 *
 * TclExprFloatError --
 *
 *	This procedure is called when an error occurs during a
 *	floating-point operation. It reads errno and sets
 *	interp->objResultPtr accordingly.
 *
 * Results:
 *	interp->objResultPtr is set to hold an error message.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

void
TclExprFloatError(interp, value)
    Tcl_Interp *interp;		/* Where to store error message. */
    double value;		/* Value returned after error;  used to
				 * distinguish underflows from overflows. */
{
    char *s;

    Tcl_ResetResult(interp);
    if ((errno == EDOM) || (value != value)) {
	s = "domain error: argument not in valid range";
	Tcl_AppendToObj(Tcl_GetObjResult(interp), s, -1);
	Tcl_SetErrorCode(interp, "ARITH", "DOMAIN", s, (char *) NULL);
    } else if ((errno == ERANGE) || IS_INF(value)) {
	if (value == 0.0) {
	    s = "floating-point value too small to represent";
	    Tcl_AppendToObj(Tcl_GetObjResult(interp), s, -1);
	    Tcl_SetErrorCode(interp, "ARITH", "UNDERFLOW", s, (char *) NULL);
	} else {
	    s = "floating-point value too large to represent";
	    Tcl_AppendToObj(Tcl_GetObjResult(interp), s, -1);
	    Tcl_SetErrorCode(interp, "ARITH", "OVERFLOW", s, (char *) NULL);
	}
    } else {			/* FAILS IF STRING REP CONTAINS NULLS */
	char msg[100];
	
	sprintf(msg, "unknown floating-point error, errno = %d", errno);
	Tcl_AppendToObj(Tcl_GetObjResult(interp), msg, -1);
	Tcl_SetErrorCode(interp, "ARITH", "UNKNOWN", msg, (char *) NULL);
    }
}

/*
 *----------------------------------------------------------------------
 *
 * Tcl_GetCommandFromObj --
 *
 *      Returns the command specified by the name in a Tcl_Obj.
 *
 * Results:
 *	Returns a token for the command if it is found. Otherwise, if it
 *	can't be found or there is an error, returns NULL.
 *
 * Side effects:
 *      May update the internal representation for the object, caching
 *      the command reference so that the next time this procedure is
 *	called with the same object, the command can be found quickly.
 *
 *----------------------------------------------------------------------
 */

Tcl_Command
Tcl_GetCommandFromObj(interp, objPtr)
    Tcl_Interp *interp;		/* The interpreter in which to resolve the
				 * command and to report errors. */
    register Tcl_Obj *objPtr;	/* The object containing the command's
				 * name. If the name starts with "::", will
				 * be looked up in global namespace. Else,
				 * looked up first in the current namespace
				 * if contextNsPtr is NULL, then in global
				 * namespace. */
{
    Interp *iPtr = (Interp *) interp;
    register ResolvedCmdName *resPtr;
    register Command *cmdPtr;
    Namespace *currNsPtr;
    int result;

    /*
     * Get the internal representation, converting to a command type if
     * needed. The internal representation is a ResolvedCmdName that points
     * to the actual command.
     */
    
    if (objPtr->typePtr != &tclCmdNameType) {
        result = tclCmdNameType.setFromAnyProc(interp, objPtr);
        if (result != TCL_OK) {
            return (Tcl_Command) NULL;
        }
    }
    resPtr = (ResolvedCmdName *) objPtr->internalRep.otherValuePtr;

    /*
     * Get the current namespace.
     */
    
    if (iPtr->varFramePtr != NULL) {
	currNsPtr = iPtr->varFramePtr->nsPtr;
    } else {
	currNsPtr = iPtr->globalNsPtr;
    }

    /*
     * Check the context namespace and the namespace epoch of the resolved
     * symbol to make sure that it is fresh. If not, then force another
     * conversion to the command type, to discard the old rep and create a
     * new one. Note that we verify that the namespace id of the context
     * namespace is the same as the one we cached; this insures that the
     * namespace wasn't deleted and a new one created at the same address
     * with the same command epoch.
     */
    
    cmdPtr = NULL;
    if ((resPtr != NULL)
	    && (resPtr->refNsPtr == currNsPtr)
	    && (resPtr->refNsId == currNsPtr->nsId)
	    && (resPtr->refNsCmdEpoch == currNsPtr->cmdRefEpoch)) {
        cmdPtr = resPtr->cmdPtr;
        if (cmdPtr->cmdEpoch != resPtr->cmdEpoch) {
            cmdPtr = NULL;
        }
    }

    if (cmdPtr == NULL) {
        result = tclCmdNameType.setFromAnyProc(interp, objPtr);
        if (result != TCL_OK) {
            return (Tcl_Command) NULL;
        }
        resPtr = (ResolvedCmdName *) objPtr->internalRep.otherValuePtr;
        if (resPtr != NULL) {
            cmdPtr = resPtr->cmdPtr;
        }
    }

    if (cmdPtr == NULL) {
	return (Tcl_Command) NULL;
    }
    return (Tcl_Command) cmdPtr;
}

/*
 *----------------------------------------------------------------------
 *
 * FreeCmdNameInternalRep --
 *
 *	Frees the resources associated with a cmdName object's internal
 *	representation.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Decrements the ref count of any cached ResolvedCmdName structure
 *	pointed to by the cmdName's internal representation. If this is 
 *	the last use of the ResolvedCmdName, it is freed. This in turn
 *	decrements the ref count of the Command structure pointed to by 
 *	the ResolvedSymbol, which may free the Command structure.
 *
 *----------------------------------------------------------------------
 */

static void
FreeCmdNameInternalRep(objPtr)
    register Tcl_Obj *objPtr;	/* CmdName object with internal
				 * representation to free. */
{
    register ResolvedCmdName *resPtr =
	(ResolvedCmdName *) objPtr->internalRep.otherValuePtr;

    if (resPtr != NULL) {
	/*
	 * Decrement the reference count of the ResolvedCmdName structure.
	 * If there are no more uses, free the ResolvedCmdName structure.
	 */
    
        resPtr->refCount--;
        if (resPtr->refCount == 0) {
            /*
	     * Now free the cached command, unless it is still in its
             * hash table or if there are other references to it
             * from other cmdName objects.
	     */
	    
            Command *cmdPtr = resPtr->cmdPtr;
            TclCleanupCommand(cmdPtr);
            ckfree((char *) resPtr);
        }
    }
}

/*
 *----------------------------------------------------------------------
 *
 * DupCmdNameInternalRep --
 *
 *	Initialize the internal representation of an cmdName Tcl_Obj to a
 *	copy of the internal representation of an existing cmdName object. 
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	"copyPtr"s internal rep is set to point to the ResolvedCmdName
 *	structure corresponding to "srcPtr"s internal rep. Increments the
 *	ref count of the ResolvedCmdName structure pointed to by the
 *	cmdName's internal representation.
 *
 *----------------------------------------------------------------------
 */

static void
DupCmdNameInternalRep(srcPtr, copyPtr)
    Tcl_Obj *srcPtr;		/* Object with internal rep to copy. */
    register Tcl_Obj *copyPtr;	/* Object with internal rep to set. */
{
    register ResolvedCmdName *resPtr =
        (ResolvedCmdName *) srcPtr->internalRep.otherValuePtr;

    copyPtr->internalRep.twoPtrValue.ptr1 = (VOID *) resPtr;
    copyPtr->internalRep.twoPtrValue.ptr2 = NULL;
    if (resPtr != NULL) {
        resPtr->refCount++;
    }
    copyPtr->typePtr = &tclCmdNameType;
}

/*
 *----------------------------------------------------------------------
 *
 * SetCmdNameFromAny --
 *
 *	Generate an cmdName internal form for the Tcl object "objPtr".
 *
 * Results:
 *	The return value is a standard Tcl result. The conversion always
 *	succeeds and TCL_OK is returned.
 *
 * Side effects:
 *	A pointer to a ResolvedCmdName structure that holds a cached pointer
 *	to the command with a name that matches objPtr's string rep is
 *	stored as objPtr's internal representation. This ResolvedCmdName
 *	pointer will be NULL if no matching command was found. The ref count
 *	of the cached Command's structure (if any) is also incremented.
 *
 *----------------------------------------------------------------------
 */

static int
SetCmdNameFromAny(interp, objPtr)
    Tcl_Interp *interp;		/* Used for error reporting if not NULL. */
    register Tcl_Obj *objPtr;	/* The object to convert. */
{
    Interp *iPtr = (Interp *) interp;
    char *name;
    Tcl_Command cmd;
    register Command *cmdPtr;
    Namespace *currNsPtr;
    register ResolvedCmdName *resPtr;

    /*
     * Get "objPtr"s string representation. Make it up-to-date if necessary.
     */

    name = objPtr->bytes;
    if (name == NULL) {
	name = Tcl_GetStringFromObj(objPtr, (int *) NULL);
    }

    /*
     * Find the Command structure, if any, that describes the command called
     * "name". Build a ResolvedCmdName that holds a cached pointer to this
     * Command, and bump the reference count in the referenced Command
     * structure. A Command structure will not be deleted as long as it is
     * referenced from a CmdName object.
     */

    cmd = Tcl_FindCommand(interp, name, (Tcl_Namespace *) NULL,
	    /*flags*/ 0);
    cmdPtr = (Command *) cmd;
    if (cmdPtr != NULL) {
	/*
	 * Get the current namespace.
	 */
	
	if (iPtr->varFramePtr != NULL) {
	    currNsPtr = iPtr->varFramePtr->nsPtr;
	} else {
	    currNsPtr = iPtr->globalNsPtr;
	}
	
	cmdPtr->refCount++;
        resPtr = (ResolvedCmdName *) ckalloc(sizeof(ResolvedCmdName));
        resPtr->cmdPtr        = cmdPtr;
        resPtr->refNsPtr      = currNsPtr;
        resPtr->refNsId       = currNsPtr->nsId;
        resPtr->refNsCmdEpoch = currNsPtr->cmdRefEpoch;
        resPtr->cmdEpoch      = cmdPtr->cmdEpoch;
        resPtr->refCount      = 1;
    } else {
	resPtr = NULL;	/* no command named "name" was found */
    }

    /*
     * Free the old internalRep before setting the new one. We do this as
     * late as possible to allow the conversion code, in particular
     * GetStringFromObj, to use that old internalRep. If no Command
     * structure was found, leave NULL as the cached value.
     */

    if ((objPtr->typePtr != NULL)
	    && (objPtr->typePtr->freeIntRepProc != NULL)) {
	objPtr->typePtr->freeIntRepProc(objPtr);
    }
    
    objPtr->internalRep.twoPtrValue.ptr1 = (VOID *) resPtr;
    objPtr->internalRep.twoPtrValue.ptr2 = NULL;
    objPtr->typePtr = &tclCmdNameType;
    return TCL_OK;
}

/*
 *----------------------------------------------------------------------
 *
 * UpdateStringOfCmdName --
 *
 *	Update the string representation for an cmdName object.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Generates a panic. 
 *
 *----------------------------------------------------------------------
 */

static void
UpdateStringOfCmdName(objPtr)
    Tcl_Obj *objPtr;		/* CmdName obj to update string rep. */
{
    /*
     * This procedure is never invoked since the internal representation of
     * a cmdName object is never modified.
     */

    panic("UpdateStringOfCmdName should never be invoked");
}
