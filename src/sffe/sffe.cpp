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
(arg)->value = (sfNumber *)malloc(sizeof(sfNumber));                       \
    if ((arg)->value) {                                                        \
        (arg)->type = sfvar_type_managed_ptr;                                  \
        cmplxset(*((arg)->value), (val), 0);                                   \
}
#else
#define sfset(arg, val)                                                        \
(arg)->value = (sfNumber *)malloc(sizeof(sfNumber));                       \
    if ((arg)->value) {                                                        \
        (arg)->type = sfvar_type_managed_ptr;                                  \
        *((arg)->value) = (val);                                               \
}
#endif

/** utils */

void sffe_error_message(int errorcode, char *context, char *errormessage)
{
    switch (errorcode) {
    case MemoryError:
        sprintf(errormessage, "%s", TR("Message", "Out of memory"));
        break;
    case UnbalancedBrackets:
        sprintf(errormessage, TR("Message", "Unbalanced parentheses"),
                context);
        break;
    case UnknownFunction:
        sprintf(errormessage, TR("Message", "Unknown function: %s"),
                context);
        break;
    case InvalidNumber:
        sprintf(errormessage, TR("Message", "Invalid number: %s"), context);
        break;
    case UnknownVariable:
        sprintf(errormessage, TR("Message", "Unknown variable: %s"),
                context);
        break;
    case InvalidOperators:
        sprintf(errormessage, TR("Message", "Invalid operator: %s"),
                context);
        break;
    case StackError:
        sprintf(errormessage,
                TR("Message", "Internal error occurred in formula: %s"),
                context);
        break;
    case InvalidParameters:
        sprintf(errormessage,
                TR("Message", "Function has incorrect parameter count: %s"),
                context);
        break;
    case EmptyFormula:
        sprintf(errormessage, TR("Message", "Empty formula"), context);
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
    for (unsigned char idx = 5; idx < sffnctscount; idx += 1) {
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
    p->expression = NULL;
    p->args = NULL;
    p->oprs = NULL;
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

sfNumber sffe_eval(sffe *const parser)
{
    sfopr *optr = parser->oprs;
    sfopr *optro = optr;
    sfopr *optrl = optr + parser->oprCount;
    for (; optr != optrl; optr += 1, optro += 1) {
        optro->arg->parg = optro->arg - 1;
        sfarg *arg = optr->arg;
#ifdef SFFE_DIRECT_FPTR
        arg->parg = optr->fnc(arg)->parg;
#else
        arg->parg = optr->fnc->fptr(arg)->parg;
#endif
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
    return NULL;
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
        unsigned char type; /* store priority of the operator 'f' */
        unsigned char args; /* number of parameters */
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
    sfarg *_argument, *_arg_itr;

    char *tokens; /*tokenized form : (f(n)+f(n))*f(n)-n (f-func, n-num,const) */

    char *ech;
    char *ch1, *ch2, *buf, *bufp;

    bool num_mode;

    unsigned int ui1, buflen;
    unsigned char token;

    enum sffe_error err;

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

#define insert_fnc_slot()                                                      \
    for (_arg_itr = _parser->args + _parser->argCount - 1;                     \
                                                                               _arg_itr > _argument; _arg_itr -= 1) {                                \
            *_arg_itr = *(_arg_itr - 1);                                           \
    }                                                                          \
        sfset(_arg_itr, -1.0);

#define pop_expression()                                                       \
    {                                                                          \
            _expression->size -= 1;                                                \
            insert_fnc_slot();                                                     \
            _parser->oprs[ui1].arg = (sfarg *)_argument;                           \
            _parser->oprs[ui1].fnc = _expression->stck[_expression->size].fnc;     \
            ui1 += 1;                                                              \
            _argument += 1;                                                        \
    }

#define max(a, b) ((a > b) ? a : b)

#ifdef SFFE_DEVEL
    clock_t begin = clock();
    printf("parse - BEGIN\n");
#endif

    /**************** code */
    _functions = NULL;
    ech = (char *)expression;
    tokens = (char *)malloc(1);
    err = MemoryError;
    _parser = *parser;

    /* clear all internal structures */
    if (_parser->expression) {
        sffe_clear(parser);
    }

    _parser->oprCount = 0;
    _parser->argCount = 0;

    ech = (char *)malloc(strlen(expression) + 1);
    sf_strdup(&ech, expression);
    _parser->expression = ech;

#ifdef SFFE_DEVEL
    printf(
        "\n|-----------------------------------------\n+ > %s[%d] - parsing\n|-----------------------------------------\n",
        __FILE__, __LINE__);
    printf("| input (len.=%tu): |%s|\n", strlen(_parser->expression),
           _parser->expression);
#endif

    /*! PHASE 1 !!!!!!!!! remove spaces, count brackets, change decimal
     * separators ',' to '.', remove multiple operators eg. ++--++1 -> 1, -+++2
     * -> -2 */
    ch1 = NULL;
    ui1 = 0; /*brackets */
    ch2 = ech;
    buflen = strlen(_parser->expression) * 2;

    buf = new char[buflen + 1];
    *buf = '\0';
    bufp = buf;

    /* skip leading spaces */
    while (isspace(*ech)) {
        ech += 1;
    }

    num_mode = false;
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

        *ch2 = (char)tolower((int)*ech);

        if (*ech == '{') {
            num_mode = true;
        } else if (*ech == '}') {
            num_mode = false;
        }
        if (!num_mode && strchr("+-", (int)*ech)) {
            if (!*buf) {
                *bufp = '0';
            } else if(strchr("(,", (int)*bufp)) {
                bufp++;
                *bufp = '0';
            }
        }

        /*fix multiple arithm operators */
        if (*bufp && strchr("+-/*^", (int)*ech) && strchr("+-/*^", (int)*bufp)) {
            if (*bufp == '-' && *ech == '-') {
                //*ch1 = '+';
                *bufp = '+';
            } else if (*bufp == '-' && *ech == '+') {
                //*ch1 = '-';
                *bufp = '-';
            } else if (*bufp == '+' && *ech == '-') {
                //*ch1 = '-';
                *bufp = '-';
            } else if (*bufp == *ech) {
                //*ch1 = *ech;
                *bufp = *ech;
            } else if (strchr("+-", (int)*ech)) {
                //ch1 = ++ch2;
                bufp++;
                *bufp = *ch2;
            } else if (*bufp != *ech) {
                err = InvalidOperators;
                break;
            }
        } else {

            //ch1 = ch2;
            if (*buf) {
                bufp++;
            }
            *bufp = *ch2;
            ch2 ++;
        }

        /*skip spaces */
        do {
            ech += 1;
        } while (isspace(*ech));
    }
    bufp++;
    *bufp = '\0';
    *ch2 = '\0';

    /* _parser->expression = (char *)realloc((char *)_parser->expression,
                                          strlen(_parser->expression) + 1); */
    const char* _tempexpr = _parser->expression;
    _parser->expression = buflen > 0 ? (char *)realloc((char *)buf, strlen(buf) + 1) : new char{'\0'};

    if (ui1 && !err) {
        err = UnbalancedBrackets;
    }

#ifdef SFFE_DEVEL
    printf("| check (len.=%tu): |%s|\n", strlen(_parser->expression),
           _parser->expression);
#endif

    if (strlen(_parser->expression) == 0)
        err = EmptyFormula;

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

                /* if not -> ERROR */
                if (!*_function) {
                    *ech = 0; // terminate string after function name
                    set_error(UnknownFunction);
                }

                token = 'f';
                break;
            }
            /* is it a real number? */
        } else if (isdigit(*ech) ||
                   (strchr("/*^(", (int)token) && strchr("+-", *ech))) {

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
            /* add value slots for uses operators/functions */
            ui1 = _parser->argCount + _parser->oprCount;
            _parser->args = (sfarg *)realloc(_parser->args, ui1 * sizeof(sfarg));
            memset(_parser->args + _parser->argCount, 0,
                   _parser->oprCount * sizeof(sfarg));
            _parser->argCount = ui1;
            _argument = _parser->args;
            _parser->oprs = (sfopr *)malloc(_parser->oprCount * sizeof(sfopr));
            ch1 = NULL; /* number */

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
                    if (ch1) {
#ifdef SFFE_DEVEL
                        printf("%c", *ch1);
#endif
                        _argument += 1;
                    }

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

                            if (_expression->size == 0) {
                                break;
                            }
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

                    unsigned char parcnt = function->parcnt & 0x1F;
                    /* mark operator as a function, and store number of
                     * available parameters (0 - unlimited) */
                    opstck->type = 0x60 | parcnt;
                    opstck->args = parcnt;

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
                    /* (  */
                case '(': {
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
                    /* check if anything has been read !!! */
                    if (ch1) {
#ifdef SFFE_DEVEL
                        printf("%c", *ch1);
#endif
                        _argument += 1;
                        ch1 = NULL;
                    }

                    /* if there is something on stack, flush if we need to read
                     * next parameter */
                    while (_expression->size) {
                        pop_expression();
#ifdef SFFE_DEVEL
                        printf("%c", _expression->stck[_expression->size].c);
#endif
                    }

                    struct __expression *pstack = _expression->prev;
                    struct _operator *opstck =
                        &pstack->stck[pstack->size -
                                      1]; // here is last function before
                        // opening new op stack

                    /* wrong number of parameters */
                    if ((opstck->type & 0x1f) == 1) {
                        set_error(InvalidParameters);
                    }

                    /* reduce a number of allowed parameters */
                    opstck->type = 0x60 | max(0, (opstck->type & 0x1f) - 1);
                } break;
                    /* )  */
                case ')': {
                    if (ch1) {
#ifdef SFFE_DEVEL
                        printf("%c", *ch1);
#endif
                        _argument += 1;
                    }
                    ch1 = NULL;

                    /* if there is something on stack, flush it we need to read
                     * next parameter */
                    while (_expression->size) {
                        pop_expression();
#ifdef SFFE_DEVEL
                        printf("%c", _expression->stck[_expression->size].c);
#endif
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

                            /* wrong number of parameters */
                            if ((opstck->type & 0x1f) > 1) {
                                set_error(InvalidParameters);
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
                    break;
                }
                ech += 1;
            }

            if (!err) {

                if (ch1) {
#ifdef SFFE_DEVEL
                    printf("%c", *ch1);
#endif
                    _argument += 1;
                }

                /*clean up _expression */
                while (_expression) {
                    while (_expression->size) {
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

                /* set up formula call stack */
                (_parser->args)->parg = NULL;

                for (ui1 = 1; ui1 < _parser->argCount; ui1 += 1) {
                    (_parser->args + ui1)->parg = (_parser->args + ui1 - 1);
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

                printf("\n| functions fnctbl:");
                for (ui1 = 0; ui1 < _parser->oprCount; ui1 += 1) {
                    printf(" 0x%.6X [%s]", (int)(size_t)_functions[ui1]->fptr,
                           _functions[ui1]->name);
                }

                printf("\n| functions used ptrs:");
                for (ui1 = 0; ui1 < _parser->oprCount; ui1 += 1) {
                    printf(" 0x%.6X", (int)(size_t)_parser->oprs[ui1].fnc);
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
         * return) */
            _parser->result =
                (sfNumber *)(_parser->oprs + _parser->oprCount - 1)->arg->value;

            if (!_parser->result)
                err = MemoryError;
        }

    if (err) {
#ifdef SFFE_DEVEL
        sffe_print_error(err, ch1);
#endif
        sffe_setup_error(_parser, err, ch1);
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

#ifdef SFFE_DEVEL
    printf("\nparse - END\n");
#endif
    //TODO
    delete _parser->expression;
    _parser->expression = _tempexpr;
    return err;
}

#undef sfset
#undef sfvar
