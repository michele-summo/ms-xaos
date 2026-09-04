/*/////////////////////////////////////////////////////////////////////////////////////
// project : sFFe ( SegFault (or Segmentation Fault :) ) formula evalutor )
// author  : Mateusz Malczak ( mateusz@malczak.info )
// wpage   : malczak.info
/////////////////////////////////////////////////////////////////////////////////////*/

#include <cstdlib>
#include <cstdio>
#include <csetjmp>
#include <cctype>
#include <cstring>

#ifdef DEBUG
#define SFFE_DEVEL
#endif

#ifdef SFFE_DEVEL
#include <time.h>
#endif

#include "sffe.h"
#include "misc-f.h"
#include "i18n.h"

#ifdef SFFE_CMPLX_ASM
#include "sffe_cmplx_asm.h"
#elif SFFE_CMPLX_GSL
#include "sffe_cmplx_gsl.h"
#elif SFFE_DOUBLE
#include "sffe_real.h"
#endif

#ifdef SFFE_COMPLEX
#define sfset(arg, val)                                                        \
(arg)->args = NULL;                                                        \
    (arg)->argc = 0;                                                           \
    (arg)->omitted = false;                                                    \
    (arg)->value = (sfNumber *)malloc(sizeof(sfNumber));                       \
    if ((arg)->value) {                                                        \
        (arg)->type = sfvar_type_managed_ptr;                                  \
        cmplxset(*((arg)->value), (val), 0);                                   \
}
#else
#define sfset(arg, val)                                                        \
(arg)->args = NULL;                                                        \
    (arg)->argc = 0;                                                           \
    (arg)->omitted = false;                                                    \
    (arg)->value = (sfNumber *)malloc(sizeof(sfNumber));                       \
    if ((arg)->value) {                                                        \
        (arg)->type = sfvar_type_managed_ptr;                                  \
        *((arg)->value) = (val);                                               \
    }
#endif

/** utils */

void sffe_error_message(int errorcode, char *context, char *errormessage)
{
    /* The context is user-supplied formula text of unbounded length, so every
     * message is written with a bounded snprintf. */
    if (!context) {
        context = (char *)"";
    }

    switch (errorcode) {
    case MemoryError:
        snprintf(errormessage, SFFE_ERRORMSG_SIZE, "%s",
                 TR("Message", "Out of memory"));
        break;
    case UnbalancedBrackets:
        snprintf(errormessage, SFFE_ERRORMSG_SIZE,
                 TR("Message", "Unbalanced parentheses"), context);
        break;
    case UnknownFunction:
        snprintf(errormessage, SFFE_ERRORMSG_SIZE,
                 TR("Message", "Unknown function: %s"), context);
        break;
    case InvalidNumber:
        snprintf(errormessage, SFFE_ERRORMSG_SIZE,
                 TR("Message", "Invalid number: %s"), context);
        break;
    case UnknownVariable:
        snprintf(errormessage, SFFE_ERRORMSG_SIZE,
                 TR("Message", "Unknown variable: %s"), context);
        break;
    case InvalidOperators:
        snprintf(errormessage, SFFE_ERRORMSG_SIZE,
                 TR("Message", "Invalid operator: %s"), context);
        break;
    case StackError:
        snprintf(errormessage, SFFE_ERRORMSG_SIZE,
                 TR("Message", "Internal error occurred in formula: %s"),
                 context);
        break;
    case InvalidParameters:
        snprintf(errormessage, SFFE_ERRORMSG_SIZE,
                 TR("Message", "Function has incorrect parameter count: %s"),
                 context);
        break;
    case EmptyFormula:
        snprintf(errormessage, SFFE_ERRORMSG_SIZE,
                 TR("Message", "Empty formula"), context);
        break;
    }
}

void sffe_setup_error(sffe *parser, enum sffe_error errorcode, char *context)
{
    /* try to store error message */
    if (parser->errormsg) {
        sffe_error_message(errorcode, context, parser->errormsg);
    }
}

void sf_strdup(char **out, const char *in)
{
    size_t name_len = strlen(in);
    char *dup = (char *)malloc(name_len + 1);
    if (dup) {
        for (size_t i = 0; i < name_len; i += 1) {
            dup[i] = (char)tolower((int)in[i]);
        }
        dup[name_len] = '\0';
    }
    *out = dup;
}

unsigned char sf_priority(char *chr)
{
    switch (*chr) {
    case 'f':
        return 0x60;
    case '^':
        return 0x40;
    /* prefix sign: looser than '^' so that -z^2 is -(z^2), tighter than '*'
     * so that -z*c is (-z)*c */
    case 'u':
        return 0x30;
    case '/':
    case '*':
        return 0x20;
    case '+':
    case '-':
        return 0x00;
    default:
        return 0x80;
    }
}

/** debug functions */

#ifdef SFFE_DEVEL
/* in debug mode report errors on stdout */
void sffe_print_error(enum sffe_error errorcode, char *context)
{
    char errormessage[256];
    sffe_error_message(errorcode, context, errormessage);
    printf("Parser error : %s", errormessage);
}
#endif

/************************* custom functions */
/* all variables used in this section are defined depanding on complex number
 * implementation */
sffunction *sffe_function(char *fn, size_t len)
{
    /* sffnctscount - defined in sffe_cmplx_* file */
    for (unsigned char idx = sffnctsfirst; idx < sffnctscount; idx += 1) {
        if (!strncmp(fn, sfcmplxfunc[idx].name, len) &&
            strlen(sfcmplxfunc[idx].name) == len) {
            return (sffunction *)(sfcmplxfunc + idx);
        }
    }
    return NULL;
}

sffunction *sffe_operator(char op)
{
    switch (op) {
    case '^':
        return (sffunction *)sfcmplxfunc;
    case '+':
        return (sffunction *)sfcmplxfunc + 1;
    case '-':
        return (sffunction *)sfcmplxfunc + 2;
    case '*':
        return (sffunction *)sfcmplxfunc + 3;
    case '/':
        return (sffunction *)sfcmplxfunc + 4;
    default:
        break;
    }
    return NULL;
}

/* Prefix sign, the last of the operator entries. Unary '+' has no effect, so
 * it needs no operation at all. */
sffunction *sffe_unary_operator(char op)
{
    if (op == '-') {
        return (sffunction *)sfcmplxfunc + (sffnctsfirst - 1);
    }
    return NULL;
}

void *sffe_const(char *fn, size_t len, void *ptr)
{
    for (unsigned char idx = 0; idx < sfvarscount; idx += 1) {
        if (!strncmp(fn, sfcnames[idx], len) && strlen(sfcnames[idx]) == len) {
            sfcvals[idx]((sfNumber *)ptr);
            return ptr;
        }
    }
    return NULL;
}

/************************* custom function */

/* Starts compiling a lazy call. The parser owns every sflazy it creates, so
 * that abandoning a parse half way still releases them. */
sflazy *sffe_lazy_new(sffe *parser, sfselptr sel)
{
    sflazy **grown = (sflazy **)realloc(
        parser->lazy, (parser->lazyCount + 1) * sizeof(sflazy *));
    if (!grown) {
        return NULL;
    }
    parser->lazy = grown;

    sflazy *lz = (sflazy *)malloc(sizeof(sflazy));
    if (!lz) {
        return NULL;
    }
    lz->select = sel;
    lz->nblocks = 0;
    lz->bounds = NULL;
    lz->probe = NULL;
    lz->source = NULL;
    parser->lazy[parser->lazyCount++] = lz;
    return lz;
}

/* Records where one argument of a lazy call ends and the next begins. While
 * the call is open nblocks counts the boundaries stored; sffe_lazy_close turns
 * it into the number of blocks. */
bool sffe_lazy_bound(sflazy *lz, unsigned int op)
{
    unsigned int *grown = (unsigned int *)realloc(
        lz->bounds, (lz->nblocks + 1) * sizeof(unsigned int));
    if (!grown) {
        return false;
    }
    lz->bounds = grown;
    lz->bounds[lz->nblocks] = op;
    lz->nblocks += 1;
    return true;
}

/* Records which block a choice really runs. A block whose argument was left
 * empty runs the one before it, so that "ifiter(f, , , g)" spends three passes
 * on f: the empty places are extra passes of the argument they follow, and
 * where they lead is settled here rather than on every pass. */
bool sffe_lazy_source(sflazy *lz, unsigned int idx, bool omitted)
{
    unsigned int *grown = (unsigned int *)realloc(
        lz->source, (idx + 1) * sizeof(unsigned int));
    if (!grown) {
        return false;
    }
    lz->source = grown;
    lz->source[idx] = (omitted && idx) ? lz->source[idx - 1] : idx;
    return true;
}

void sffe_lazy_close(sflazy *lz)
{
    /* n arguments leave n + 1 boundaries behind */
    if (lz->nblocks) {
        lz->nblocks -= 1;
    }
}

sffe *sffe_alloc(void)
{
    sffe *rp = (sffe *)malloc(sizeof(sffe));
    if (!rp) {
        return NULL;
    }

    memset(rp, 0, sizeof(sffe));
    return rp;
}

void sffe_clear(sffe **parser)
{
    sffe *p = *parser;
    unsigned int i = 0;
    for (; i < p->argCount; i++) {
        if (p->args[i].type == sfvar_type_managed_ptr) {
            free(p->args[i].value);
        }
        free(p->args[i].args);
    }

    if (p->args) {
        free(p->args);
    }

    if (p->expression) {
        free((char *)p->expression);
    }

    if (p->oprs) {
        free(p->oprs);
    }

    if (p->lazy) {
        for (i = 0; i < p->lazyCount; i++) {
            if (p->lazy[i]) {
                free(p->lazy[i]->bounds);
                free(p->lazy[i]->source);
                free(p->lazy[i]);
            }
        }
        free(p->lazy);
        p->lazy = NULL;
    }
    p->lazyCount = 0;

    p->expression = NULL;
    p->args = NULL;
    p->oprs = NULL;
    p->result = NULL;
    /* The counts describe the arrays just released. Leaving them set makes the
     * next sffe_clear walk a freed NULL args array, which is exactly what
     * happens when a failed parse is followed by sffe_free or a re-parse. */
    p->argCount = 0;
    p->oprCount = 0;
}

void sffe_free(sffe **parser)
{
    sffe *p = *parser;

    sffe_clear(parser);

    if (p->userf) {
        free(p->userf);
    }

    if (p->varCount) {
        unsigned int i = 0;
        for (; i < p->varCount; i++) {
            if (p->variables[i].type == sfvar_type_managed_ptr) {
                free(p->variables[i].value);
            }
        }

        free(p->variables);
    }

    free(*parser);
    parser = NULL;
}

// avg time: 0.250767773750267
// avg time: 0.252280894999276

#ifdef SFFE_DIRECT_FPTR
#define sffe_call(op) (op)->fnc((op)->arg)
#else
#define sffe_call(op) (op)->fnc->fptr((op)->arg)
#endif

/* Runs operations [from, to). A dispatch operation runs just one of the
 * argument ranges behind it and then skips the rest, which is what makes a
 * lazy call cost one branch instead of all of them. Recurses so that lazy
 * calls nest. */
static void sffe_run(sfopr *const oprs, unsigned int from, unsigned int to)
{
    unsigned int pc = from;
    while (pc < to) {
        sfopr *op = oprs + pc;
        if (op->lazy) {
            sflazy *lz = op->lazy;
            unsigned int nsel = lz->nblocks;
            const sfNumber *probe = NULL;
            if (lz->probe) {
                /* the last argument is a value the selector needs rather than
                 * a branch: run it first, then choose among the others */
                nsel -= 1;
                sffe_run(oprs, lz->bounds[nsel], lz->bounds[nsel + 1]);
                probe = lz->probe->value;
            }
            unsigned int k = lz->select(nsel, probe);
            /* an empty argument runs the block of the one it follows */
            k = lz->source[k];
            sffe_run(oprs, lz->bounds[k], lz->bounds[k + 1]);
            pc = lz->bounds[lz->nblocks];
        } else {
            sffe_call(op);
            pc += 1;
        }
    }
}

/* Runs the operations in the order the parser emitted them. Each one already
 * knows where its operands live, so evaluating it is just a call: there is no
 * stack to maintain and no bookkeeping between operations. */
sfNumber sffe_eval(sffe *const parser)
{
    if (parser->lazyCount) {
        sffe_run(parser->oprs, 0, parser->oprCount);
    } else {
        /* Nothing to branch on, so do not pay for the test per operation. */
        sfopr *optr = parser->oprs;
        sfopr *optrl = optr + parser->oprCount;
        for (; optr != optrl; optr += 1) {
            sffe_call(optr);
        }
    }
    return *(parser->result);
}

sfvariable *sffe_var(sffe *const parser, const char *name)
{
    if (parser->varCount) {
        sfvariable *var = parser->variables;
        sfvariable *lvar = parser->variables + parser->varCount;
        while (var < lvar) {
            if (!strcasecmp(var->name, name)) {
                return var;
            }
            var += 1;
        }
    }

    return NULL;
}

sfvariable *sffe_regvar(sffe **parser, sfNumber *vptrs, const char *name)
{
    sffe *parser_ = *parser;
    sfvariable *var = sffe_var(parser_, name);
    if (var) {
        return var;
    }

    int vars_cnt = parser_->varCount + 1;
    parser_->variables = (sfvariable *)realloc(parser_->variables,
                                                vars_cnt * sizeof(sfvariable));
    if (!parser_->variables) {
        return NULL;
    }

    var = parser_->variables + parser_->varCount;
    if (vptrs == NULL) {
        vptrs = (sfNumber *)malloc(sizeof(sfNumber));
        memset(vptrs, 0, sizeof(sfNumber));
        var->type = sfvar_type_managed_ptr;
    } else {
        var->type = sfvar_type_ptr;
    }
    var->value = vptrs;
    sf_strdup(&var->name, name);

    parser_->varCount += 1;
    return var;
}

void sffe_regvars(sffe **parser, unsigned int cN, sfNumber **vptrs,
                  char *const *names)
{
    while (cN > 0) {
        cN -= 1;
        sffe_regvar(parser, (vptrs ? vptrs[cN] : NULL), names[cN]);
    }
}

sfNumber *sffe_setvar(sffe **parser, sfNumber vptrs, const char *name)
{
    sfNumber *value;
    sffe *parser_ = *parser;
    sfvariable *var = sffe_var(parser_, name);
    if (!var) {
        var = sffe_regvar(parser, NULL, name);
    }

    value = var->value;
    memcpy(value, &vptrs, sizeof(sfNumber));
    return value;
}

void *sffe_regfunc(sffe **parser, const char *vname, unsigned int parcnt,
                   sffptr funptr)
{
    sffe *parser_ = *parser;
    sffunction *sff;
    unsigned short i;

    parser_->userf = (sffunction *)realloc(
        parser_->userf, (parser_->userfCount + 1) * sizeof(sffunction));
    if (!parser_->userf) {
        return NULL;
    }

    sff = parser_->userf + parser_->userfCount;

    for (i = 0; i < strlen(vname); i++)
        sff->name[i] = tolower(vname[i]);
    sff->name[i] = 0;

    sff->parcnt = parcnt;
    sff->fptr = funptr;
    sff->sel = NULL; /* user functions get their arguments evaluated eagerly */
    /* and every argument written out: there is no table entry to say
     * otherwise, and the memory this came from says nothing at all */
    sff->probe = false;
    sff->optfrom = 0;

    parser_->userfCount += 1;
    return (void *)sff;
}

sfNumber *sffe_variable(sffe *const p, char *fname, size_t len)
{
    char fn[len + 1];
    fn[len] = '\0';
    memcpy(fn, fname, len);
    sfvariable *var = sffe_var(p, fn);
    if (var) {
        return var->value;
    }
    /* Nothing registered answers to it; the host may still know it. */
    if (p->resolve) {
        return p->resolve(p, fn);
    }
    return NULL;
}

/* Walks the compiled program rather than the text: the text would also match a
 * name that merely contains this one, and a variable registered and left
 * unnamed by the formula would look used. An operation holds its operands by
 * pointer, so a leaf that is read appears as an operand of something. */
int sffe_reads(sffe *const parser, const sfNumber *value)
{
    if (parser == NULL || value == NULL) {
        return 0;
    }
    /* A formula that is one value and no operations at all: there is nothing
     * to walk, and the value it hands back is the leaf itself. */
    if (parser->oprCount == 0) {
        return parser->result == value;
    }
    for (unsigned int i = 0; i < parser->oprCount; i += 1) {
        sfarg *op = parser->oprs[i].arg;
        if (op == NULL) {
            continue; /* a lazy dispatch, which holds no operands */
        }
        for (unsigned char k = 0; k < op->argc; k += 1) {
            if (op->args[k]->value == value) {
                return 1;
            }
        }
    }
    return 0;
}

sffunction *userfunction(const sffe *const p, char *fname, size_t len)
{
    unsigned char idx;
    for (idx = 0; idx < p->userfCount; idx += 1) {
        const char *name = p->userf[idx].name;
        if (strlen(name) == len && (!strncmp(fname, name, len))) {
            return (sffunction *)(p->userf + idx);
        }
    }
    return NULL;
}

#ifdef SFFE_COMPLEX
/* parse complex number in format
 * { [-+]ddd[.dddd[e[+-]ddd]] ; [-+]ddd[.dddd[e[+-]ddd]] }  */
char sffe_docmplx(char **str, sfarg **arg)
{

    char *ch1;
    ch1 = *str;
    number_t re = xstrtonum(ch1, str);
    if (ch1 == *str)
        return 1;
    if (*(*str)++ != ',')
        return 2;
    ch1 = *str;
    number_t im = xstrtonum(ch1, str);
    if (ch1 == *str)
        return 1;
    if (*(*str)++ != '}')
        return 2;

    cmplxset(*(*arg)->value, re, im);
    return 0;
}
#endif

char sffe_doname(char **str)
{
    do {
        *str += 1;
    } while (isalnum(**str) || **str == '_');

    if (**str == '(') {
        return 2; /* ( start of parameters  */
    }

    return 1;
}

int sffe_parse(sffe **parser, const char *expression)
{
    /**************variables */

    sffe *_parser;
    struct _operator { // @todo replace with stack/list
#ifdef SFFE_DEVEL
        char c; /* used in debug build to store operator character */
#endif
        const char *name;   /* function/operator name, for error messages */
        unsigned char type; /* store priority of the operator 'f' */
        unsigned char args; /* number of parameters */
        bool variadic;      /* takes any number of arguments */
        unsigned char seen; /* variadic only: arguments counted so far */
        sfselptr sel;       /* non-NULL: arguments are evaluated lazily */
        bool probe;         /* the last argument is read by the selector */
        sflazy *lazy;       /* the call being compiled, owned by the parser */
        unsigned char optfrom; /* first argument that may be left out, 1 based */
        unsigned int argi;     /* argument being read, counting from zero */
#ifdef SFFE_DIRECT_FPTR
        sffptr fnc;
#else
        sffunction *fnc;
#endif
    };

    struct __expression {
        struct _operator *stck;    /* operators on stack */
        unsigned int size;         /* number of items on stack */
        struct __expression *prev; /* previous stack */
    } * _tmp_exp, *_expression;

    sffunction **_functions; /* hold all functions used in expression in left -
                                to - right order */
    sffunction **_function;  /* currently expected function from 'fnctbl' */
    sfarg *_argument;

    char *tokens; /*tokenized form : (f(n)+f(n))*f(n)-n (f-func, n-num,const) */

    char *ech;
    char *ch1, *ch2;

    unsigned int ui1;
    unsigned char token;

    enum sffe_error err;
    const char *errctx; /* text quoted in the error message; NULL means ch1 */

    /* Number of values waiting to be consumed while phase 3 walks the tokens.
     * Every operator and function pops its operands off this count and pushes
     * its result back, so a formula is only well formed if each of them finds
     * enough operands and exactly one value is left at the end. Without this
     * the parser happily emits a call stack that sffe_eval walks off. */
    unsigned int _depth;

    /* Values produced so far and not yet consumed, innermost last. _depth is
     * its height. Holds slots of _parser->args, which is not reallocated once
     * phase 3 has started. */
    sfarg **_vstack;

    /* Next free slot of _parser->args to hand out as an operation result. The
     * leaves occupy the front of the array, one reserved slot per operation
     * follows. */
    unsigned int _result_slot;

    /* Lazy calls seen while tokenizing. Each needs one dispatch operation on
     * top of the operations its arguments compile to. */
    unsigned int lazy_calls;
    /* Room reserved for the argument places a call may have left empty, each
     * of which needs a node of its own, and the next of those nodes to hand
     * out. Counted before phase 3, which cannot grow the array it holds
     * pointers into. */
    unsigned int empty_places, empty_slot;

    /* Number of entries in _functions, which is one per function and operator
     * the tokenizer found. Not the same as the final _parser->oprCount: that
     * one counts the dispatch operations too. Only the SFFE_DEVEL tracing
     * reads it back. */
    [[maybe_unused]] unsigned int function_count;

    /**************used defines */

#define append_token(chr)                                                      \
    tokens = (char *)realloc(tokens, ui1 + 2);                                 \
        tokens[ui1++] = chr;                                                       \
        ch2 = tokens + ui1 - 1;                                                    \
        token = chr;                                                               \
        tokens[ui1] = '\0';

#define set_error(errno)                                                       \
    {                                                                          \
            err = errno;                                                           \
            break;                                                                 \
    }

/* Same, but names what the message should quote. In phase 3 ch1 points into
 * the tokenized form ("fn(n)"), which is meaningless to the user. */
#define set_error_at(errno, ctx)                                               \
    {                                                                          \
            err = errno;                                                           \
            errctx = ctx;                                                          \
            break;                                                                 \
    }

/* Emits one operation: takes its operands off the value stack, wires them into
 * a fresh result slot and pushes that back. Leaves err set instead of breaking
 * out, so the loops that flush the stack have to test !err themselves. */
#define pop_expression()                                                       \
    {                                                                          \
            struct _operator *_op = _expression->stck + (--_expression->size);     \
            unsigned char _nargs = _op->variadic ? _op->seen : _op->args;          \
            if (_depth < _nargs) {                                                 \
                err = ((_op->type & 0xE0) == 0x60) ? InvalidParameters              \
                                                   : InvalidOperators;             \
                errctx = _op->name;                                                \
            } else {                                                               \
                sfarg *_res = _parser->args + _result_slot;                        \
                /* "f()" takes nothing, and malloc(0) may or may not hand      \
                 * back a pointer; either way there is nothing to hold */      \
                _res->args =                                                       \
                    _nargs ? (sfarg **)malloc(_nargs * sizeof(sfarg *)) : NULL;    \
                if (_nargs && !_res->args) {                                       \
                    err = MemoryError;                                             \
                } else {                                                           \
                    _depth -= _nargs;                                              \
                    _res->argc = _nargs;                                           \
                    /* operands are right to left for the sfaramN macros */     \
                    for (unsigned char _k = 0; _k < _nargs; _k += 1) {             \
                        _res->args[_k] = _vstack[_depth + _nargs - 1 - _k];        \
                    }                                                              \
                    _res->value = (sfNumber *)malloc(sizeof(sfNumber));            \
                    if (!_res->value) {                                            \
                        err = MemoryError;                                         \
                    } else {                                                       \
                        _res->type = sfvar_type_managed_ptr;                       \
                        /* defined until the operation writes it */             \
                        cmplxset(*(_res->value), 0, 0);                            \
                        /* the value the selector reads, now that the operands \
                         * are known; args[0] is the last one written */       \
                        if (_op->lazy && _op->probe) {                             \
                            _op->lazy->probe = _res->args[0];                      \
                        }                                                          \
                        _vstack[_depth] = _res;                                    \
                        _depth += 1;                                               \
                        _result_slot += 1;                                         \
                        _parser->oprs[ui1].arg = _res;                             \
                        _parser->oprs[ui1].fnc = _op->fnc;                         \
                        _parser->oprs[ui1].lazy = NULL;                            \
                        ui1 += 1;                                                  \
                    }                                                              \
                }                                                                  \
            }                                                                      \
    }

/* An empty place is only a place the function said it could fill in itself.
 * Anywhere else it is an argument that was simply not given. */
#define check_empty(op)                                                        \
    if (!(op)->optfrom || (op)->argi + 1 < (op)->optfrom) {                    \
            set_error_at(InvalidParameters, (op)->name);                           \
    }

/* Stands an argument the call left empty on the value stack. It is a leaf like
 * any other, so the arity counting and the operand wiring need know nothing
 * about it; only the callee, reading the mark, does. The slot comes from the
 * places counted after phase 1, phase 3 not being free to grow an array whose
 * entries are already held by pointer.
 *
 * For a lazy call an empty argument means the argument before it, so it takes
 * that one's value rather than owning one -- and, not owning it, must not be
 * the node that frees it. For everyone else it is a zero, which is what poly
 * wants and what a function reading its own defaults never looks at. */
#define push_empty(op)                                                         \
    if (empty_slot >= _parser->argCount) {                                     \
            /* the reservation is an upper bound, so this cannot happen; and   \
             * saying so costs one test where being wrong costs the array */   \
            set_error_at(StackError, _parser->expression);                         \
    } else {                                                                   \
            sfarg *_hole = _parser->args + empty_slot;                             \
            empty_slot += 1;                                                       \
            _hole->args = NULL;                                                    \
            _hole->argc = 0;                                                       \
            _hole->omitted = true;                                                 \
            if ((op)->lazy && (op)->argi) {                                        \
                _hole->type = sfvar_type_ptr;                                      \
                _hole->value = _vstack[_depth - 1]->value;                         \
            } else {                                                               \
                _hole->value = (sfNumber *)malloc(sizeof(sfNumber));               \
                if (!_hole->value) {                                               \
                    err = MemoryError;                                             \
                } else {                                                           \
                    _hole->type = sfvar_type_managed_ptr;                          \
                    cmplxset(*(_hole->value), 0, 0);                               \
                }                                                                  \
            }                                                                      \
            _vstack[_depth] = _hole;                                               \
            _depth += 1;                                                           \
    }

#define max(a, b) ((a > b) ? a : b)

#ifdef SFFE_DEVEL
    clock_t begin = clock();
    printf("parse - BEGIN\n");
#endif

    /**************** code */
    _functions = NULL;
    errctx = NULL;
    _depth = 0;
    _vstack = NULL;
    _result_slot = 0;
    lazy_calls = 0;
    empty_places = 0;
    function_count = 0;
    tokens = (char *)malloc(1);
    err = MemoryError;
    _parser = *parser;

    /* Take both copies before touching the parser: callers legitimately pass
     * parser->expression straight back in, and sffe_clear would free it from
     * under us. 'work' is chewed up by phases 1 and 2; 'original' is what the
     * caller gets to read back afterwards. */
    char *work = NULL;
    char *original = NULL;
    sf_strdup(&work, expression);
    sf_strdup(&original, expression);
    if (!work || !original) {
        free(work);
        free(original);
        free(tokens);
        sffe_setup_error(_parser, MemoryError, NULL);
        return MemoryError;
    }

    /* clear all internal structures */
    if (_parser->expression) {
        sffe_clear(parser);
    }

    _parser->oprCount = 0;
    _parser->argCount = 0;

    ech = work;
    _parser->expression = work;

#ifdef SFFE_DEVEL
    printf(
        "\n|-----------------------------------------\n+ > %s[%d] - parsing\n|-----------------------------------------\n",
        __FILE__, __LINE__);
    printf("| input (len.=%tu): |%s|\n", strlen(_parser->expression),
           _parser->expression);
#endif

    /*! PHASE 1 !!!!!!!!! remove spaces, count brackets, change separators
     * ';' -> ',' and '[' ']' -> '{' '}'.
     *
     * Signs are deliberately left untouched. Phase 2 decides for each '+'/'-'
     * whether it is a prefix or an infix operator, so runs like "--" need no
     * rewriting here, and a prefix sign no longer has to be padded into a
     * subtraction (the "-x" -> "0-x" trick this phase used to do). */
    ch1 = NULL;
    ui1 = 0; /*brackets */
    ch2 = ech;

    /* skip leading spaces */
    while (isspace(*ech)) {
        ech += 1;
    }

    /*handle brackets and change ';'->',', '['->'{', ']'->'}' */
    while (*ech) {
        switch (*ech) {
        case '[':
            *ech = '{';
            break;
        case '(':
            ui1 += 1;
            break;
        case ']':
            *ech = '}';
            break;
        case ')':
            ui1 -= 1;
            break;
        case ';':
            *ech = ',';
            break;
        }

        /* compact in place: spaces are dropped, so ch2 trails ech */
        *ch2 = (char)tolower((int)*ech);
        ch2 += 1;

        /*skip spaces */
        do {
            ech += 1;
        } while (isspace(*ech));
    }
    *ch2 = '\0';

    if (ui1 && !err) {
        err = UnbalancedBrackets;
    }

#ifdef SFFE_DEVEL
    printf("| check (len.=%tu): |%s|\n", strlen(_parser->expression),
           _parser->expression);
#endif

    if (strlen(_parser->expression) == 0)
        err = EmptyFormula;

    /* Room for the places a call may have left empty -- the nothing between
     * two separators in "f(z, ,5)" -- each of which phase 3 gives a node of
     * its own, out of an array it is not free to grow.
     *
     * Counting the places one could be rather than the ones there are. A place
     * is empty when phase 2 produces no token between two separators, which is
     * not quite the same as there being nothing between them in the text: a
     * unary plus produces no token either. Reproducing that decision here
     * would be reproducing phase 2; taking one separator or closing bracket to
     * be one possible place cannot be short, since each of them takes at most
     * one, and costs a handful of nodes nobody asks for. */
    for (const char *scan = _parser->expression; *scan; scan++)
        if (*scan == ',' || *scan == ')')
            empty_places += 1;

    /*! PHASE 2 !!!!!!!! tokenize expression, lexical analysis (need
     * optimizations) */
    *tokens = '\0';
    ch2 = NULL;
    ui1 = 0;
    ch1 = NULL; /*string starting position */
    ech = (char *)_parser->expression;
    token = '('; /* in case of leading '-' */

    while (*ech && !err) {
        ch1 = ech;

        if (isalpha(*ech)) {
            switch (sffe_doname(&ech)) {
            case 1: /* const or variable */
                _parser->args = (sfarg *)realloc(
                    _parser->args, (_parser->argCount + 1) * sizeof(sfarg));

                if (!_parser->args) {
                    set_error(MemoryError);
                }

                _argument = _parser->args + (_parser->argCount++);
                _argument->type = sfvar_type_ptr;
                _argument->args = NULL; /* a leaf has no operands */
                _argument->argc = 0;
                _argument->omitted = false;
                _argument->value = (sfNumber *)sffe_variable(
                    _parser, ch1, (size_t)(ech - ch1));

                if (!_argument->value) {
                    sfset(_argument, 10.0); //? temporary const value
                    if (_argument->value) {
                        if (!sffe_const(ch1, (size_t)(ech - ch1),
                                        _argument->value)) {
                            *ech = 0; // terminate string after this symbol
                            set_error( UnknownVariable);
                        }
                    } else {
                        set_error(MemoryError);
                    }
                }

                token = 'n';
                break;

            case 2: /* function */
                _functions = (sffunction **)realloc(
                    _functions,
                    (_parser->oprCount + 1) * sizeof(sffunction *));

                if (!_functions) {
                    set_error(MemoryError);
                }

                _function = _functions + (_parser->oprCount++);
                *_function = NULL;

                if (_parser->userfCount) {
                    /*is it user defined function */
                    *_function = (sffunction *)(void *)userfunction(
                        _parser, ch1, (size_t)(ech - ch1));
                }

                if (!*_function) {
                    /*if not, is it build in function */
                    *_function = (sffunction *)(void *)sffe_function(
                        ch1, (size_t)(ech - ch1));
                }

                /* if not, or if the table lists the name with no implementation
                 * behind it, -> ERROR. Calling one of those used to jump
                 * through a NULL pointer. */
                if (!*_function || !(*_function)->fptr) {
                    *ech = 0; // terminate string after function name
                    set_error(UnknownFunction);
                }

                /* a lazy call also needs its dispatch operation */
                if ((*_function)->sel) {
                    lazy_calls += 1;
                }

                token = 'f';
                break;
            }
            /* is it a real number? */
        } else if (isdigit(*ech)) {

            /* numbers (this part can be optimized) */
            ch1 = ech; /* st = 1;  */

            number_t value = xstrtonum(ch1, &ech);
            if (ch1 == ech) {
                set_error(InvalidNumber);
            }

            /*epx */
            _parser->args = (sfarg *)realloc(
                _parser->args, (++_parser->argCount) * sizeof(sfarg));

            if (!_parser->args) {
                set_error(MemoryError);
            }

            _argument = _parser->args + _parser->argCount - 1;
            sfset(_argument, value);

            /*epx */
            token = 'n';
        }
#ifdef SFFE_COMPLEX
        /* if not, it can be complex number */
        else if (*ech == '{') {
            ech += 1;
            _parser->args = (sfarg *)realloc(
                _parser->args, (++_parser->argCount) * sizeof(sfarg));

            if (!_parser->args) {
                set_error(MemoryError);
            }

            _argument = _parser->args + _parser->argCount - 1;
            sfset(_argument, 0);

            if (sffe_docmplx(&ech, &_argument)) {
                set_error(InvalidNumber);
            }

            token = 'n';
        }
#endif
        /* a '+' or '-' with no left operand in sight is a prefix sign */
        else if (strchr("+-", (int)*ech) && strchr("(,+-*/^u", (int)token)) {
            ch1 = ech;
            ech += 1;

            /* unary plus is the identity, so it needs no operation at all */
            if (*ch1 == '+') {
                continue;
            }

            _functions = (sffunction **)realloc(
                _functions, (++_parser->oprCount) * sizeof(sffunction *));

            if (!_functions) {
                set_error(MemoryError);
            }

            _functions[_parser->oprCount - 1] = sffe_unary_operator(*ch1);
            token = 'u';
        }
        /* if not, we have operator */
        else {
            if (*ech != '(' && *ech != ')' && *ech != ',') {
                sffunction *function = sffe_operator(*ech);

                if (function) {
                    _functions = (sffunction **)realloc(
                        _functions,
                        (++_parser->oprCount) * sizeof(sffunction *));

                    if (!_functions) {
                        set_error(MemoryError);
                    }

                    _functions[_parser->oprCount - 1] = function;
                } else {
                    *(ech + 1) = 0; // terminate string after operator
                    set_error(InvalidOperators);
                }
            }

            ch1 = ech;
            token = *ech;
            ech += 1;
        }

        /* no error and already has any opcodes - check for skipped
         * multiplication. Handle nf, n(, )(, )f, )n, fn */
        if (!err && ui1 > 0) {
            if (token == 'f' || token == 'n' || token == '(') // last token
            {
                if (*ch2 == 'n' || *ch2 == ')') // last-1 token
                {
                    sffunction *oprptr = sffe_operator('*');
                    _functions = (sffunction **)realloc(
                        _functions,
                        (++_parser->oprCount) * sizeof(sffunction *));

                    if (!_functions) {
                        set_error(MemoryError);
                    }

                    /* if last token was function inject multiplication before
                     * it */
                    if (token == 'f') {
                        _functions[_parser->oprCount - 1] =
                            _functions[_parser->oprCount - 2];
                        _functions[_parser->oprCount - 2] =
                            (sffunction *)oprptr;
                    } else {
                        _functions[_parser->oprCount - 1] =
                            (sffunction *)oprptr;
                    }

                    // inject multiplication
                    unsigned char tmp = token;
                    append_token('*');
                    token = tmp;
                }
            }
        }

        append_token(token);
    }

    ech = tokens;

#ifdef SFFE_DEVEL
    printf(
        "| compiled expr.: |%s|\n| operations: %d\n| numbers,vars: %d\n| stack not.: ",
        tokens, _parser->oprCount, _parser->argCount);
#endif

    /*! PRE PHASE 3 !!!!! no operations in expression = single numeric value */
    if (!_parser->oprCount && _parser->argCount == 1) {
        _parser->oprs = (sfopr *)malloc(_parser->argCount * sizeof(sfopr));
        _parser->oprs[0].arg = (sfarg *)_parser->args;
        _parser->oprs[0].fnc = NULL;
        _parser->result = (sfNumber *)_parser->args->value;
    } else
        /*! PHASE 3 !!!!! create sffe 'stack' notation ]:-> */
        /* lots of memory operations are done here but no memory leaks should
           occur */
        if (!err) {
            /* One result slot per operation, after the leaves. From here on
             * _parser->args must not move: the operands are held by pointer. */
            _result_slot = _parser->argCount;
            function_count = _parser->oprCount;
            /* the empty places get their nodes after the result slots, where
             * the memset below leaves them saying they own nothing */
            empty_slot = _parser->argCount + _parser->oprCount;
            ui1 = _parser->argCount + _parser->oprCount + empty_places;
            _parser->args = (sfarg *)realloc(_parser->args, ui1 * sizeof(sfarg));
            /* Zero the slots we just added: sffe_clear frees whatever they say
             * they own, so an uninitialised one becomes a free() of garbage. */
            memset(_parser->args + _parser->argCount, 0,
                   (ui1 - _parser->argCount) * sizeof(sfarg));
            _parser->argCount = ui1;
            /* room for the dispatch operations too; they produce no value and
             * so have no result slot above */
            _parser->oprs = (sfopr *)malloc((_parser->oprCount + lazy_calls) *
                                            sizeof(sfopr));
            /* An operation consumes at least one value and leaves one, so the
             * stack never grows past the number of slots. */
            _vstack = (sfarg **)malloc((ui1 + 1) * sizeof(sfarg *));
            if (!_parser->args || !_parser->oprs || !_vstack) {
                err = MemoryError;
            }
            ch1 = NULL; /* number */
            _depth = 0;
            /* leaves are consumed in the order phase 2 appended them */
            _argument = _parser->args;

            /* stacks ( stores operations and controls parameters count inside of
         * brackts blocks ) */
            _expression =
                (struct __expression *)malloc(sizeof(struct __expression));
            _expression->size = 0; /* 0-stack is empty, but ready to write (one slot
                                  allocated), >0-number of element on stack */
            _expression->stck =
                (struct _operator *)malloc(sizeof(struct _operator));
            _expression->prev = NULL;
            memset(_expression->stck, 0, sizeof(struct _operator));

            ui1 = 0; /* used in defines */
            _function = _functions;

            while (*ech && !err) {
                switch (*ech) {
                    /*  O */
                case '+':
                case '-':
                case '*':
                case '/':
                case '^': {
#ifdef SFFE_DEVEL
                    if (ch1) {
                        printf("%c", *ch1);
                    }
#endif

                    unsigned char type = sf_priority(ech);
                    /* there is an operator on stack */
                    if (_expression->size) {
                        /* remove all operators with higher, or equal priority
                         */
                        while (type <=
                               _expression->stck[_expression->size - 1].type) {
                            pop_expression();
#ifdef SFFE_DEVEL
                            printf("%c",
                                   _expression->stck[_expression->size].c);
#endif

                            if (err || _expression->size == 0) {
                                break;
                            }
                        }
                        if (err) {
                            break;
                        }

                        _expression->stck = (struct _operator *)realloc(
                            _expression->stck,
                            (_expression->size + 1) * sizeof(struct _operator));
                    }

                    sffunction *function = *_function;

#ifdef SFFE_DEVEL
                    struct _operator *opstck =
                        &_expression->stck[_expression->size];
                    opstck->c = *ech;
#endif

                    /* store operator priority */
                    _expression->stck[_expression->size].type = type;
                    /* every infix operator in the table is binary */
                    _expression->stck[_expression->size].args = 2;
                    _expression->stck[_expression->size].variadic = false;
                    _expression->stck[_expression->size].sel = NULL;
                    _expression->stck[_expression->size].probe = false;
                    _expression->stck[_expression->size].lazy = NULL;
                    _expression->stck[_expression->size].optfrom = 0;
                    _expression->stck[_expression->size].argi = 0;
                    _expression->stck[_expression->size].name = function->name;

                    /* get function pointer */
#ifdef SFFE_DIRECT_FPTR
                    _expression->stck[_expression->size].fnc = function->fptr;
#else
                    _expression->stck[_expression->size].fnc = function;
#endif

                    _expression->size += 1;

                    _function += 1;
                    ch1 = NULL;
                } break;
                    /* F  */
                case 'f': {
                    _expression->stck = (struct _operator *)realloc(
                        _expression->stck,
                        (_expression->size + 1) * sizeof(struct _operator));

                    sffunction *function = *_function;

                    struct _operator *opstck =
                        &_expression->stck[_expression->size];
#ifdef SFFE_DEVEL
                    opstck->c = 'f';
#endif

                    bool anyargs = (function->parcnt == SFFE_VARIADIC);
                    unsigned char parcnt =
                        anyargs ? 0 : (function->parcnt & 0x1F);
                    /* mark operator as a function, and store number of
                     * available parameters. A variadic call records 0, which
                     * makes the fixed-arity checks below inert; its real
                     * count is tallied in 'seen' as commas are consumed. */
                    opstck->type = 0x60 | parcnt;
                    opstck->args = parcnt;
                    opstck->variadic = anyargs;
                    opstck->seen = 1;
                    opstck->name = function->name;
                    /* the '(' that follows starts the lazy bookkeeping */
                    opstck->sel = function->sel;
                    opstck->probe = function->probe;
                    opstck->lazy = NULL;
                    opstck->optfrom = function->optfrom;
                    opstck->argi = 0;

                    /* get function pointer */
#ifdef SFFE_DIRECT_FPTR
                    _expression->stck[_expression->size].fnc = function->fptr;
#else
                    _expression->stck[_expression->size].fnc = function;
#endif

                    _expression->size += 1;

                    _function += 1;
                    ch1 = NULL;

                    // consume ()
                    //                    if(!parcnt)
                    //                    {
                    //                        ech += 2;
                    //                    }

                } break; // skip to ( ???
                    /* u - prefix sign */
                case 'u': {
                    _expression->stck = (struct _operator *)realloc(
                        _expression->stck,
                        (_expression->size + 1) * sizeof(struct _operator));

                    struct _operator *opstck =
                        &_expression->stck[_expression->size];
#ifdef SFFE_DEVEL
                    opstck->c = 'u';
#endif

                    /* A prefix operator takes what comes after it, so nothing
                     * already on the stack can be resolved yet: unlike an infix
                     * operator it pops nothing before pushing itself. That is
                     * also what makes "--z" stack instead of cancelling out. */
                    opstck->type = sf_priority(ech);
                    opstck->args = 1;
                    opstck->variadic = false;
                    opstck->sel = NULL;
                    opstck->probe = false;
                    opstck->lazy = NULL;
                    opstck->optfrom = 0;
                    opstck->argi = 0;
                    opstck->name = (*_function)->name;

#ifdef SFFE_DIRECT_FPTR
                    opstck->fnc = (*_function)->fptr;
#else
                    opstck->fnc = *_function;
#endif

                    _expression->size += 1;

                    _function += 1;
                    ch1 = NULL;
                } break;
                    /* (  */
                case '(': {
                    /* If this parenthesis opens the arguments of a lazy call,
                     * reserve the dispatch operation that will stand in front
                     * of them, and start recording where each argument's code
                     * begins. */
                    if (_expression->size) {
                        struct _operator *call =
                            &_expression->stck[_expression->size - 1];
                        if (call->sel && !call->lazy) {
                            call->lazy = sffe_lazy_new(_parser, call->sel);
                            if (!call->lazy) {
                                set_error(MemoryError);
                            }
                            _parser->oprs[ui1].arg = NULL;
                            _parser->oprs[ui1].fnc = NULL;
                            _parser->oprs[ui1].lazy = call->lazy;
                            ui1 += 1;
                            if (!sffe_lazy_bound(call->lazy, ui1)) {
                                set_error(MemoryError);
                            }
                        }
                    }

                    /* store current stack */
                    _tmp_exp = (struct __expression *)malloc(
                        sizeof(struct __expression));
                    _tmp_exp->prev = _expression;
                    _expression = _tmp_exp;
                    _expression->size = 0;
                    _expression->stck =
                        (struct _operator *)malloc(sizeof(struct _operator));

#ifdef SFFE_DEVEL
                    _expression->stck[0].c = '_';
#endif

                    token = 0;
                } break;
                    /*  ; */
                case ',': {
#ifdef SFFE_DEVEL
                    if (ch1) {
                        printf("%c", *ch1);
                    }
#endif
                    ch1 = NULL;

                    /* if there is something on stack, flush if we need to read
                     * next parameter */
                    while (_expression->size && !err) {
                        pop_expression();
#ifdef SFFE_DEVEL
                        printf("%c", _expression->stck[_expression->size].c);
#endif
                    }
                    if (err) {
                        break;
                    }

                    struct __expression *pstack = _expression->prev;

                    /* A comma only separates arguments of a call. Outside one
                     * there is no enclosing stack to look at, and reading it
                     * anyway is how "z,c" used to crash the parser. */
                    if (!pstack || !pstack->size ||
                        (pstack->stck[pstack->size - 1].type & 0xE0) != 0x60) {
                        set_error_at(InvalidOperators, ",");
                    }

                    struct _operator *opstck =
                        &pstack->stck[pstack->size -
                                      1]; // here is last function before
                        // opening new op stack

                    /* Nothing at all stands between this separator and the
                     * one, or the bracket, before it: the argument was left
                     * empty for the function to fill in. */
                    bool hole =
                        (ech > tokens) && (ech[-1] == '(' || ech[-1] == ',');
                    if (hole) {
                        check_empty(opstck);
                        push_empty(opstck);
                        if (err) {
                            break;
                        }
                    }

                    /* one argument of a lazy call has just been compiled */
                    if (opstck->lazy &&
                        (!sffe_lazy_bound(opstck->lazy, ui1) ||
                         !sffe_lazy_source(opstck->lazy, opstck->argi, hole))) {
                        set_error(MemoryError);
                    }
                    opstck->argi += 1;

                    if (opstck->variadic) {
                        /* no declared arity to check against, just tally one
                         * more argument. The only ceiling is what 'seen' can
                         * hold. */
                        if (opstck->seen == 0xFF) {
                            set_error_at(InvalidParameters, opstck->name);
                        }
                        opstck->seen += 1;
                        break;
                    }

                    /* wrong number of parameters */
                    if ((opstck->type & 0x1f) == 1) {
                        set_error_at(InvalidParameters, opstck->name);
                    }

                    /* reduce a number of allowed parameters */
                    opstck->type = 0x60 | max(0, (opstck->type & 0x1f) - 1);
                } break;
                    /* )  */
                case ')': {
#ifdef SFFE_DEVEL
                    if (ch1) {
                        printf("%c", *ch1);
                    }
#endif
                    ch1 = NULL;

                    /* if there is something on stack, flush it we need to read
                     * next parameter */
                    while (_expression->size && !err) {
                        pop_expression();
#ifdef SFFE_DEVEL
                        printf("%c", _expression->stck[_expression->size].c);
#endif
                    }
                    if (err) {
                        break;
                    }

                    /* no stack available = stack overrelesed */
                    if (!_expression->prev) {
                        set_error(StackError);
                    }

                    _tmp_exp = _expression;
                    _expression = _tmp_exp->prev;

                    /* destroy block stack */
                    free(_tmp_exp->stck);
                    free(_tmp_exp);

                    /* parser was reading function, at the top of current stack
                     * is a function. identified by '*.t==3' */
                    if (_expression->size) {
                        struct _operator *opstck =
                            &_expression
                                 ->stck[_expression->size -
                                        1]; // here is last function before
                            // opening new op stack
                        if ((opstck->type & 0xE0) == 0x60) {

                            /* the last argument was left empty, exactly as one
                             * between two separators is */
                            bool hole = (ech > tokens) && (ech[-1] == ',');
                            if (hole) {
                                check_empty(opstck);
                                push_empty(opstck);
                                if (err) {
                                    break;
                                }
                            }

                            /* "f()": nothing was written at all, which only a
                             * function whose every argument has a default can
                             * mean. It is a call of no arguments, not a call of
                             * one that was left out, so no place is taken for
                             * it and none was counted. */
                            if (ech > tokens && ech[-1] == '(') {
                                if (opstck->optfrom != 1) {
                                    set_error_at(InvalidParameters,
                                                 opstck->name);
                                }
                                opstck->type = 0x60;
                                opstck->args = 0;
                                opstck->seen = 0;
                            }

                            /* wrong number of parameters */
                            if ((opstck->type & 0x1f) > 1) {
                                set_error_at(InvalidParameters, opstck->name);
                            }

                            /* Close the lazy call: the last argument ends here
                             * and this is also where the dispatch resumes, on
                             * the call's own operation emitted just below. */
                            if (!err && opstck->lazy) {
                                if (!sffe_lazy_bound(opstck->lazy, ui1) ||
                                    !sffe_lazy_source(opstck->lazy,
                                                      opstck->argi, hole)) {
                                    set_error(MemoryError);
                                }
                                sffe_lazy_close(opstck->lazy);
                            }

                            if (!err) {
                                pop_expression();
#ifdef SFFE_DEVEL
                                printf("%c",
                                       _expression->stck[_expression->size].c);
#endif
                                if (_expression->size) {
                                    _expression->stck =
                                        (struct _operator *)realloc(
                                            _expression->stck,
                                            (_expression->size) *
                                                sizeof(struct _operator));
                                }
                            }
                        }
                    }

                } break;
                    /* n */
                case 'n':
                    ch1 = ech;
                    /* leaves were appended by phase 2 in this same order */
                    _vstack[_depth] = _argument;
                    _depth += 1;
                    _argument += 1;
                    break;
                }
                ech += 1;
            }

            if (!err) {

#ifdef SFFE_DEVEL
                if (ch1) {
                    printf("%c", *ch1);
                }
#endif

                /*clean up _expression */
                while (_expression) {
                    while (_expression->size && !err) {
                        pop_expression();
#ifdef SFFE_DEVEL
                        printf("%c", _expression->stck[_expression->size].c);
#endif
                    }

                    free(_expression->stck);
                    _tmp_exp = _expression->prev;
                    free(_expression);
                    _expression = _tmp_exp;
                }

                /* Every operand must have been consumed by exactly one
                 * operation, leaving just the result behind. Anything else
                 * means the call stack we built does not match the arguments
                 * we have, and sffe_eval would read past them. */
                if (!err && _depth != 1) {
                    err = _depth ? InvalidOperators : EmptyFormula;
                    errctx = _parser->expression;
                }

                /* ui1 counted the dispatch operations too, so it is the real
                 * length of the program the evaluator has to run. */
                if (!err) {
                    _parser->oprCount = ui1;
                }

#ifdef SFFE_DEVEL
                printf("\n| numbers: ");
                for (ui1 = 0; ui1 < _parser->argCount; ui1 += 1) {
                    if ((_parser->args + ui1)->value) {
#ifdef SFFE_COMPLEX
                        printf(" %g%+gI", real((*(_parser->args + ui1)->value)),
                               imag((*(_parser->args + ui1)->value)));
#else
                        printf(" %g", (*(_parser->args + ui1)->value));
#endif
                    } else {
                        printf(" [_]");
                    }
                }

                /* _functions holds what the tokenizer found, so it is indexed
                 * by function_count and not by oprCount, which also counts the
                 * dispatch operations the parser inserted. */
                printf("\n| functions fnctbl:");
                for (ui1 = 0; ui1 < function_count; ui1 += 1) {
                    printf(" 0x%.6X [%s]", (int)(size_t)_functions[ui1]->fptr,
                           _functions[ui1]->name);
                }

                printf("\n| functions used ptrs:");
                for (ui1 = 0; ui1 < _parser->oprCount; ui1 += 1) {
                    if (_parser->oprs[ui1].lazy) {
                        printf(" [lazy dispatch]");
                    } else {
                        printf(" 0x%.6X", (int)(size_t)_parser->oprs[ui1].fnc);
                    }
                }

                double time_spent = (double)(clock() - begin) / CLOCKS_PER_SEC;
                printf("\n| compiled in  %f s", time_spent);

                printf(
                    "\n|-----------------------------------------\n+ < %s[%d] - parsing\n|-----------------------------------------\n",
                    __FILE__, __LINE__);

#endif
            } else {
                /* prevent memory leaks */

                /* clean up stack */
                while (_expression) {
                    free(_expression->stck);
                    _tmp_exp = _expression->prev;
                    free(_expression);
                    _expression = _tmp_exp;
                }
            }

            /* set up evaluation result pointer (result is stored in last operation
         * return). Only the successful path has a fully built operation array;
         * on error the tail of it is still uninitialised, and with no
         * operations at all there is no last one to read. */
            if (!err) {
                if (!_parser->oprCount) {
                    err = EmptyFormula;
                } else {
                    _parser->result =
                        (sfNumber *)(_parser->oprs + _parser->oprCount - 1)
                            ->arg->value;

                    if (!_parser->result)
                        err = MemoryError;
                }
            }
        }

    if (err) {
        /* ch1 only points at readable source text while phase 2 is running;
         * phase 3 walks the tokenized form and names the culprit via errctx. */
        char *context = (char *)(errctx ? errctx : ch1);
#ifdef SFFE_DEVEL
        sffe_print_error(err, context);
#endif
        sffe_setup_error(_parser, err, context);
        sffe_clear(&_parser);
    }

    /*undefine defines */
#undef priority
#undef sfpopstack
#undef insertfnc
#undef code
#undef errset
#undef max

    /* free lookup tables */
    free(tokens);
    free(_functions);
    free(_vstack);

#ifdef SFFE_DEVEL
    printf("\nparse - END\n");
#endif

    /* Hand back the formula as it was given to us. Phases 1 and 2 normalise
     * and cut up the working copy, which is not what a caller reading
     * parser->expression (to display or to save it) wants to see. On the error
     * path sffe_clear has already released the working copy. */
    if (_parser->expression) {
        free((char *)_parser->expression);
    }
    _parser->expression = original;

    return err;
}

#undef sfset
#undef sfvar
