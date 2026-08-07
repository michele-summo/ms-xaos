/*/////////////////////////////////////////////////////////////////////////////////////
// project : sFFe ( SegFault (or Segmentation Fault :) ) formula evalutaor )
// author  : Mateusz Malczak ( mateusz@malczak.info )
// wpage   :
///////////////////////////////////////////////////////////////////////////////////////
// possible config definitions
//   general
//	SFFE_DOUBLE - real math parser
//	SFFE_COMPLEX - complex math parser
//	SFFE_DEVEL - print extra info to stdout
//	SFFE_DIRECT_FPTR - use direct function pointers (!!!) omits payload
//	SFFE_DLL - Windows DLL
//
//   complex numbers (for SFFE_COMPLEX)
//	SFFE_CMPLX_GSL - uses GSL complex number routines
//	SFFE_CMPLX_ASM - uses my asm complex unit (compile it with NASM)
/////////////////////////////////////////////////////////////////////////////////////*/

#ifndef SFFE_H
#define SFFE_H
#include <stdlib.h>

#define SFFE_DIRECT_FPTR 1

#ifdef SFFE_REAL
#define SFFE_DOUBLE 1
#endif

/* --- */
/*TODO long double needed*/
#ifdef SFFE_CMPLX_ASM
#define SFFE_COMPLEX 1
typedef struct cmpx__ {
    double r, i;
} cmplx;
#define sfNumber cmplx
#elif SFFE_CMPLX_GSL
#define SFFE_COMPLEX 1
#include <gsl/gsl_complex.h>
typedef gsl_complex cmplx;
#define sfNumber gsl_complex
#elif SFFE_DOUBLE
#define sfNumber double
#endif

/* Size of the buffer a caller must provide in sffe.errormsg. Error messages
 * embed user-supplied formula text, so they are truncated to fit. */
#define SFFE_ERRORMSG_SIZE 200

/* parcnt value marking a function that takes any number of arguments. The
 * count it was actually called with is in argc, like for any other call. */
#define SFFE_VARIADIC ((unsigned char)-1)

enum sffe_error {
    MemoryError,
    UnbalancedBrackets,
    UnknownFunction,
    InvalidNumber,
    UnknownVariable,
    InvalidOperators,
    StackError,
    InvalidParameters,
    EmptyFormula,
};

typedef enum { sfvar_type_ptr, sfvar_type_managed_ptr } sfvartype;

/* One value in a compiled formula: either a leaf (a number, constant or
 * variable, with no operands) or the result of an operation, which points
 * directly at the values it consumes.
 *
 * Operands are listed right to left, so args[0] is the last one written in the
 * formula -- the order the sfaramN macros below expect. Holding them by
 * pointer means an operation can be evaluated wherever it sits in the program,
 * without the evaluator having to thread a stack through them. */
typedef struct sfargument__ {
    struct sfargument__ **args; /* NULL for a leaf */
    unsigned char argc;
    sfvartype type;
    sfNumber *value;
} sfarg;

/* sffe function prototype, parameters order is right-to-left (cdecl) */
typedef sfarg *(*sffptr)(sfarg *const a);

/* constats eval functions */
typedef void (*cfptr)(sfNumber *cnst);

/* Picks which argument of a lazy call is the live one, given how many it has.
 * Must not depend on the arguments themselves: it is consulted before any of
 * them has been evaluated. */
typedef unsigned int (*sfselptr)(unsigned int argc);

/* function type structure */
typedef struct {
    sffptr fptr;
    unsigned char parcnt;
    char name[20];
    /* When set, the call is compiled so that only the selected argument is
     * evaluated. Defaulted, so the table entries that do not want it can keep
     * listing just {function, arity, name}. */
    sfselptr sel = NULL;
} sffunction;

/* Compiled form of a lazy call: the operations of argument k are the range
 * [bounds[k], bounds[k + 1]) of the program, and bounds[nblocks] is where
 * execution resumes once the selected one has run. */
typedef struct sflazy__ {
    sfselptr select;
    unsigned int nblocks;
    unsigned int *bounds;
} sflazy;

/* basic sffe 'stack' operation ( function + result slot ) */
typedef struct {
    sfarg *arg;
#ifdef SFFE_DIRECT_FPTR
    sffptr fnc;
#else
    sffunction *fnc;
#endif
    /* Set only on the dispatch pseudo-operation standing in front of a lazy
     * call's arguments; it computes nothing itself. */
    sflazy *lazy;
} sfopr;

typedef struct {
    char *name;
    sfvartype type;
    sfNumber *value;
} sfvariable;

typedef struct sfcontext__ {
    unsigned int funcsCount; /* number of default / user functions */
    sffunction *functions;

    unsigned int constsCount;
    cfptr *constants;
} sffe_context;

/* SFFE main structure */
typedef struct {
    /*public*/
    const char *expression; /* parsed expression (read-only) */
    char *errormsg;         /* parser errors (read-only) */
    sfNumber *result;       /* evaluation result (read-only) */

    /* protected/private */
    unsigned int argCount; /* number of arguments in use */
    sfarg *args;

    unsigned int oprCount; /* number of operations in use */
    sfopr *oprs;

    unsigned int lazyCount; /* lazy calls in the expression, 0 for most */
    sflazy **lazy;

    unsigned int varCount; /* number of used variables */
    sfvariable *variables;

    unsigned int userfCount; /* number of user functions */
    sffunction *userf;
} sffe;

#define SFFE sffe
#define sffeparser sffe
#define sfparser sffe
#define SFFEPARSER sffe

/* 'stack' slot value */
#define sfvalue(p) (*((p)->value))

/* function parameters, right to left: sfaram1 is the last one written */
#define sfaram1(p) ((p)->args[0])
#define sfaram2(p) ((p)->args[1])
#define sfaram3(p) ((p)->args[2])
#define sfaram4(p) ((p)->args[3])
#define sfaram5(p) ((p)->args[4])
#define sfparamN(p, N) ((p)->args[(N) - 1])

/* create formula evaluator structure */
sffe *sffe_alloc(void);
/* free fe structure */
void sffe_free(sffe **parser);

/* parse expression 'expression' and store result in 'parser' struct, error (if
 * any) returned */
int sffe_parse(sffe **parser, const char *expression);

/* evaluate function and return evaluation result */
sfNumber sffe_eval(sffe *const parser);

/* user function with name 'vname', with 'parcnt' parameters and
 * defined with function pointed by 'funptr'*/
void *sffe_regfunc(sffe **parser, const char *vname, unsigned int parcnt,
                   sffptr funptr);

/* get already registered variable pointer, NULL if variable was not registered
 */
sfvariable *sffe_var(sffe *const parser, const char *name);

/* single variable 'vptrs' identified by name 'vchars' */
// void *sffe_regvar(sffe ** parser, sfNumber * vptrs, char vchars);
sfvariable *sffe_regvar(sffe **parser, sfNumber *vptrs, const char *name);

/* multiple variables */
void sffe_regvars(sffe **parser, unsigned int cN, sfNumber **vptrs,
                  char *const *names);

// sffunction *sffe_function_alloc(char *name, sffptr function_pointer, unsigned
// char paramsCount, void *payload);

// void sffe_function_free(sffunction* function);

/* set 'vptrs' as 'vchars' variable  */
sfNumber *sffe_setvar(sffe **parser, sfNumber vptrs, const char *name);

#ifdef SFFE_CMPLX_ASM
#include "sffe_cmplx_asm.h"
#elif SFFE_CMPLX_GSL
#include "sffe_cmplx_gsl.h"
#elif SFFE_DOUBLE
#include "sffe_real.h"
#endif

#endif
