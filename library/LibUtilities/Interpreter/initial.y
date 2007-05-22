%{
/*****************************************************************************
 * Synopsis
 * --------
 * initial.y: yacc code for a simple function interpreter, which also allows
 * lookup for named tokens (all are double precision).
 *
 * Modelled on "hoc3" in Chapter 8 of "The UNIX Programming Environment"
 * by Kernighan & Pike, Prentice-Hall 1984.  A hash-table data structure is
 * used in place of the linked list of Symbols which they employed.
 * Hashing code from Chapter 6 of "The C Programming Language", 2nd Edn,
 * by Kernighan & Ritchie, Prentice-Hall 1988.
 *
 * Summary
 * -------
 * void    yy_initialize (void);
 * void    yy_help       (void);
 * void    yy_show       (void);
 * integer yy_dump       (char*, const integer);
 *
 * double  yy_interpret  (const char*);
 *
 * void    yy_vec_init   (const char*, const char*);
 * void    yy_vec_interp (const integer, ...);
 *
 * Notes
 * -----
 * 1. yy_initialize must be called before other routines will work.
 * 2. yy_help prints a summary of available functions on stdout.
 * 3. yy_show prints a summary of installed variables on stdout.
 * 4. yy_dump is similar to yy_show, but prints into a string of given length.
 * 5. yy_interpret is the central routine.  It can be used:
 *    (a), to install a named symbol in internal tables: e.g. "name = value";
 *    (b), to retrieve the value of an installed symbol: e.g. "name";
 *    (c), to evaluate a function string: e.g. "cos(x)*exp(-t)".
 * 6. yy_vec_init is used to set up the interpreter for vector evaluation.
 * 7. yy_vec_interp subsequently used for "vectorized" calls to yy_interpret.
 *
 * Operators
 * ---------
 * Unary:      -
 * Binary:     -, +, *, /, ^ (exponentiation), ~ (atan2), & (hypot), % (fmod)
 * Functions:  sin,  cos,  tan,  abs, floor, ceil, int, heav (Heaviside),
 *             asin, acos, atan, log, log10, exp,  sqrt,
 *             sinh, cosh, tanh, asinh, acosh, atanh,
 *             erf, erfc, lgamma,
 *             j0, j1, y0, y1, jn.
 * Procedures: jn, yn, rad, ang, rejn, imjn, jacobi, womcos, womsin.
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>

#define HASHSIZE 199
#define HASHSEED 31
#define VEC_MAX  32

typedef double (*PFD)( ); /* NB: no arguments -- non-ANSI (on purpose). */

typedef struct symbol 
{
  char* name;
  short type;

  union 
  {
      double val;		/* -- If VAR.   */
      PFD    ptr;		/* -- If BLTIN. */
  } u;

  struct symbol* next;
} Symbol;

static double
  Heavi  (double), 
  Jn     (double,double),
  Yn     (double,double),
  Rad    (double,double),
  Ang    (double,double),
  ReJn   (double,double,double),
  ImJn   (double,double,double),
  Jacobi (double,double,double,double),
  Womcos (double,double,double,double,double),
  Womsin (double,double,double,double,double);

static unsigned hash     (const char*);
static Symbol*  lookup   (const char*);
static Symbol*  install  (const char*, const integer, const double);
static void*    emalloc  (const size_t);

       int      yyparse (void);
static int      yylex   (void);
static void     yyerror (char*);

static double  value;
static Symbol* hashtab[HASHSIZE];
static char    func_string[STR_MAX], *cur_string;
static integer nvec = 0;
static Symbol* vs[VEC_MAX];
extern int     errno;

static struct {			    /* -- Built-in functions. */
  char* name;
  short narg;
  PFD   func;
} builtin[] = {
  "cos"   ,  1, cos     ,
  "sin"   ,  1, sin     ,
  "tan"   ,  1, tan     ,
  "exp"   ,  1, exp     ,
  "sinh"  ,  1, sinh    ,
  "cosh"  ,  1, cosh    ,
  "tanh"  ,  1, tanh    ,
  "erf"   ,  1, erf     ,
  "erfc"  ,  1, erfc    ,
  "int"   ,  1, rint    ,
  "abs"   ,  1, fabs    ,
  "floor" ,  1, floor   ,
  "ceil"  ,  1, ceil    ,
  "acos"  ,  1, acos    ,
  "asin"  ,  1, asin    ,
  "atan"  ,  1, atan    ,
  "acosh" ,  1, acosh   ,
  "asinh" ,  1, asinh   ,
  "atanh" ,  1, atanh   ,
  "log"   ,  1, log     ,
  "log10" ,  1, log10   ,
  "sqrt"  ,  1, sqrt    ,
  "heav"  ,  1, Heavi   ,
  "j0"    ,  1, j0      ,
  "j1"    ,  1, j1      ,
  "y0"    ,  1, y0      ,
  "y1"    ,  1, y1      ,
  "lgamma",  1, lgamma  ,

  "jn"    ,  2, Jn      ,
  "yn"    ,  2, Yn      ,
  "rad"   ,  2, Rad     ,
  "ang"   ,  2, Ang     ,

  "rejn"  ,  3, ReJn    ,
  "imjn"  ,  3, ImJn    ,
  
  "jacobi",  4, Jacobi  ,

  "womcos",  5, Womcos  ,
  "womsin",  5, Womsin  ,

  NULL, NULL, NULL
};

#include "LibUtilities/Interpreter/defaults.h"

%}
/* -- Yacc grammar follows: */
%union {			/* -- yacc stack type      */
  double  val;			/* -- actual value         */
  Symbol* sym;			/* -- symbol table pointer */
}
%token <val>   NUMBER
%token <sym>   VAR BLTIN_UNARY BLTIN_BINARY BLTIN_TERNARY BLTIN_QUATERNARY BLTIN_QUINTERNARY UNDEF
%type  <val>   expr asgn
%right '='
%left  '+' '-'
%left  '*' '/'
%left  UNARYMINUS
%right '^' '~' '&'		/* -- exponentiation, atan2, hypot */
%%
list:     /* nothing */
        | list '\n'
        | list asgn '\n'
        | list expr '\n'     { value = $2; }
        ;
asgn:     VAR '=' expr       { $$=$1->u.val=$3; $1->type = VAR; }
        ;
expr:     NUMBER
        | VAR                { if ($1->type == UNDEF) {
				 message ("yyparse: undefined variable ",
					  $1->name, WARNING);
			       }
			       $$ = $1->u.val;
			     }
        | asgn
        | BLTIN_UNARY       '(' expr ')'
          { $$ = (*($1->u.ptr))($3); }
        | BLTIN_BINARY      '(' expr ',' expr ')'
          { $$ = (*($1->u.ptr))($3,$5); }
        | BLTIN_TERNARY     '(' expr ',' expr ',' expr ')'
          { $$ = (*($1->u.ptr))($3,$5,$7); }
        | BLTIN_QUATERNARY  '(' expr ',' expr ',' expr ',' expr ')'
          { $$ = (*($1->u.ptr))($3,$5,$7,$9); }
        | BLTIN_QUINTERNARY '(' expr ',' expr ',' expr ',' expr ',' expr ')'
          { $$ = (*($1->u.ptr))($3,$5,$7,$9,$11); }
        | expr '+' expr      { $$ = $1 + $3; }
        | expr '-' expr      { $$ = $1 - $3; }
        | expr '*' expr      { $$ = $1 * $3; }
        | expr '/' expr      { if ($3 == 0.0) 
				 message ("yyparse",
					  "division by zero", ERROR);
			       $$ = $1 / $3;
			     }
        | expr '^' expr      { $$ = pow   ($1, $3); }
        | expr '&' expr      { $$ = hypot ($1, $3); }
        | expr '~' expr      { $$ = atan2 ($1, $3); }
        | expr '%' expr      { $$ = fmod  ($1, $3); }
        | '(' expr ')'       { $$ = $2; }
        | '-' expr %prec UNARYMINUS { $$ = -$2; }
        ;
%%


void yy_initialize (void)
/* ------------------------------------------------------------------------- *
 * Load lookup tables and symbol table with default values.
 *
 * This routine should be called at start of run-time.
 * ------------------------------------------------------------------------- */
{
  static   integer initialized = 0;
  register integer i;
  register Symbol* s;

  if (!initialized) {
    for (i = 0; consts[i].name; i++)
      install (consts[i].name, VAR, consts[i].cval);

    for (i = 0; builtin[i].name; i++)
      switch (builtin[i].narg) {
      case 1:
	s = install (builtin[i].name, BLTIN_UNARY, 0.0);
	s -> u.ptr = builtin[i].func;
	break;
      case 2:
	s = install (builtin[i].name, BLTIN_BINARY, 0.0);
	s -> u.ptr = builtin[i].func;
	break;
      case 3:
	s = install (builtin[i].name, BLTIN_TERNARY, 0.0);
	s -> u.ptr = builtin[i].func;
	break;
      case 4:
	s = install (builtin[i].name, BLTIN_QUATERNARY, 0.0);
	s -> u.ptr = builtin[i].func;
	break;
      case 5:
	s = install (builtin[i].name, BLTIN_QUINTERNARY, 0.0);
	s -> u.ptr = builtin[i].func;
	break;
      default:
	message ("yy_initialize", "never happen", ERROR);
	break;
      }

    initialized = 1;
  }
}


double yy_interpret (const char* s)
/* ------------------------------------------------------------------------- *
 * Given a string, interpret it as a function using yacc-generated yyparse.
 * ------------------------------------------------------------------------- */
{
  if (strlen (s) > STR_MAX)
    message ("yy_interpret: too many characters passed:\n", s, ERROR);
  
  strcat (strcpy (func_string, s), "\n");
  
  cur_string = func_string;
  
  yyparse ();
  return value;
}


void yy_vec_init (const char* names,
		  const char* fn   )
/* ------------------------------------------------------------------------- *
 * Set up the vector parser.
 *
 * names contains a list of variable names  e.g. "x y z",
 * fn    contains a function for evaluation e.g. "sin(x)*cos(y)*exp(z)".
 *
 * Valid separator characters in name are space, tab, comma, (semi-)colon.
 * Function string can contain previously-defined symbols (e.g. PI).
 * ------------------------------------------------------------------------- */
{
  char    routine   [] = "vecInit()";
  char    separator [] = " ,:;\t";
  char    tmp       [STR_MAX];
  char*   p;
  Symbol* s;

  if (strlen (fn) == 0)
    message (routine, "empty function string", ERROR);
  else if (strlen (fn) > STR_MAX)
    message (routine, "too many characters in function string", ERROR);

  nvec = 0;
  strcpy (tmp, names);

  p = strtok (tmp, separator);
  do {
    if (nvec++ > VEC_MAX) message (routine, "too many variables", ERROR); 
    vs[nvec-1] = (s = lookup (p)) ? s : install (p, VAR, 0.0);
  } while (p = strtok (NULL, separator));

  
  strcat (strcpy (func_string, fn), "\n");
}


void yy_vec_interp (const int ntot, ...)
/* ------------------------------------------------------------------------- *
 * Vector parser.  Following ntot there should be passed a number of
 * pointers to double (vectors), of which there should be in number the
 * number of variables named previously to vecInit, plus one: the result of
 * continually re-parsing the string "fn" is placed in the last vector, for
 * a total of ntot parsings.
 *
 * To follow on from the previous example, four vectors would be passed,
 * i.e.  vecInterp(ntot, x, y, z, u); the result fn(x,y,z) is placed in u.
 * ------------------------------------------------------------------------- */
{
  char             routine[] = "yy_vec_interp";
  register int i, n;
  double*          x[VEC_MAX];
  double*          fx = NULL;
  va_list          ap;
  
  va_start (ap, ntot);
  for (i = 0; i < nvec; i++) {
    x[i] = NULL;
    if (!(x[i] = va_arg (ap, double*)))
	message (routine, "not enough vectors passed..1", ERROR);
  }
  if (!(fx = va_arg (ap, double*)))
    message (routine, "not enough vectors passed..2", ERROR);
  va_end (ap);

  for (n = 0; n < ntot; n++) {
    cur_string = func_string;
    for (i = 0; i < nvec; i++) vs[i]->u.val = x[i][n];
    yyparse ();
    fx[n] = value;
  }
}


void yy_help (void)
/* ------------------------------------------------------------------------- *
 * Print details of callable functions to stderr.
 * ------------------------------------------------------------------------- */
{
  fprintf 
    (stderr, 
     "Unary:      -\n"
     "Binary:     -, +, *, /, ^ (exponentiation), "
     "~ (atan2), & (hypot), %% (fmod)\n"
     "Functions:  sin,  cos,  tan,  abs, floor, ceil, int, heav (Heaviside),\n"
     "            asin, acos, atan, log, log10, exp,  sqrt,\n"
     "            sinh, cosh, tanh, asinh, acosh, atanh,\n"
     "            erf, erfc, lgamma,\n"
     "            j0, j1, y0, y1, jn\n"
     "Procedures: jn, yn, rad, ang, rejn, imjn, jacobi, womcos, womsin\n");
}


void yy_show (void)
/* ------------------------------------------------------------------------- *
 * Print details of installed variables to stderr.
 * ------------------------------------------------------------------------- */
{
  register int i;
  register Symbol* sp;

  for (i = 0; i < HASHSIZE; i++)
    for (sp = hashtab[i]; sp; sp = sp -> next)
      if (sp -> type == VAR) 
	fprintf (stderr, "%-15s = %-.17g\n", sp -> name, sp -> u.val);
}


int yy_dump (char*         str,
		 const int max)
/* ------------------------------------------------------------------------- *
 * Load description of internal variables into string, to length max.
 * If string overflows, return 0, else 1.
 * ------------------------------------------------------------------------- */
{
  register int i, n = 0;
  register Symbol* sp;
  char             buf[FILENAME_MAX];

  for (i = 0; i < HASHSIZE; i++)
    for (sp = hashtab[i]; sp; sp = sp -> next)
      if (sp -> type == VAR) {
	sprintf (buf, "%-15s = %-.17g\n", sp -> name, sp -> u.val);
	if ((n += strlen (buf)) > max - 2)
	  return 0;
	else
	  strcat (str, buf);
      }

  return 1;
}


static int yylex (void)
/* ------------------------------------------------------------------------- *
 * Lexical analysis routine called by yyparse, using string loaded by
 * yy_interpret.
 * ------------------------------------------------------------------------- */
{
  register int c;

  while ((c = *cur_string++) == ' ' || c == '\t');

  if (c == EOF) return 0;

  if (c == '.' || isdigit (c)) {
    yylval.val = strtod (--cur_string, &cur_string);
    return NUMBER;
  }

  if (isalpha (c)) {
    register Symbol* s;
    char             sbuf[STR_MAX];
    register char*   p = sbuf;
    do {
      *p++ = c;
    } while ((c = *cur_string++) != EOF && (isalnum (c) || c == '_'));
    cur_string--;
    *p = '\0';
    if ((s = lookup (sbuf)) == NULL) s = install (sbuf, UNDEF, 0.0);
    yylval.sym = s;
    return (s -> type == UNDEF) ? VAR : s -> type;
  }

  return c;
}


static void yyerror (char *s)
/* ------------------------------------------------------------------------- *
 * Handler for yyparse syntax errors.
 * ------------------------------------------------------------------------- */
{
  message ("yyparse", s, WARNING);
}


static unsigned hash (const char* s)
/* ------------------------------------------------------------------------- *
 * Generate hash table index.
 * ------------------------------------------------------------------------- */
{
  register unsigned hashval;

  for (hashval = 0; *s != '\0'; s++) hashval = *s + HASHSEED * hashval;
  
  return hashval % HASHSIZE;
}


static Symbol* lookup (const char* s)
/* ------------------------------------------------------------------------- *
 * Find s in symbol hashtable.
 * ------------------------------------------------------------------------- */
{
  register Symbol* sp;
  
  for (sp = hashtab[hash (s)]; sp; sp = sp->next)
    if (strcmp (s, sp->name) == 0) return sp;

  return NULL;
}


static Symbol* install (const char*   s,
			const int t,
			const double  d)
/* ------------------------------------------------------------------------- *
 * Install s in symbol hashtable.
 * ------------------------------------------------------------------------- */
{
  register Symbol*  sp;
  register unsigned hashval;

  if (!(sp = lookup (s))) {	/* -- Not found, install in hashtab. */
    sp = (Symbol *) emalloc (sizeof (Symbol));
    if (sp == NULL || (sp -> name = strdup (s)) == NULL) return NULL;
    hashval          = hash (s);
    sp -> next       = hashtab[hashval];
    hashtab[hashval] = sp;
  }

  sp -> type  = t;
  sp -> u.val = d;
  
  return sp;
}


static void *emalloc (const size_t n)
/* ------------------------------------------------------------------------- *
 * Check return from malloc.
 * ------------------------------------------------------------------------- */
{
  void* p;

  if (!(p = (void *) malloc (n))) message ("emalloc", "out of memory", ERROR);
 
  return p;
}

static double Heavi (double x) { return (x >= 0.0) ? 1.0 : 0.0; }

static double Rad (double x, double y) { return hypot (x, y);  }
static double Ang (double x, double y) { return atan2 (y, x);  }
static double Jn  (double i, double x) { return jn((int)i, x); }
static double Yn  (double i, double x) { return yn((int)i, x); }

static double Jacobi (double z, double n, double alpha, double beta)
/* ------------------------------------------------------------------------- *
 * Return value of the n_th order Jacobi polynomial
 *   P^(alpha,beta)_n(z) alpha > -1, beta > -1 at z.
 * ------------------------------------------------------------------------- */
{
  const double     apb = alpha + beta;
  register int i,k;
  double           a1,a2,a3,a4;
  double           poly, polyn1, polyn2;
  
  polyn2 = 1.0;
  polyn1 = 0.5*(alpha - beta + (alpha + beta + 2.0)*z);
  
  for (k = 2; k <= n; ++k){
    a1 =  2.0*k*(k + apb)*(2.0*k + apb - 2.0);
    a2 = (2.0*k + apb - 1.0)*(alpha*alpha - beta*beta);
    a3 = (2.0*k + apb - 2.0)*(2.0*k + apb - 1.0)*(2.0*k + apb);
    a4 =  2.0*(k + alpha - 1.0)*(k + beta - 1.0)*(2.0*k + apb);
    
    a2 /= a1;
    a3 /= a1;
    a4 /= a1;
    
    poly   = (a2 + a3*z)*polyn1 - a4*polyn2;
    polyn2 = polyn1;
    polyn1 = poly  ;
  }

  return poly;
}

/* -- Complex Bessel function, ex netlib, and functions that use it. */

void F77NAME(zbesj) (const double*, const double*, const double*, 
		     const int*, const int*, double*, double*, 
		     int*, int*);

void zbesj (const double *x, const double *y, const double ord, 
	    const int Kode, const int n, double *ReJ, 
	    double *ImJ, int* nz, int* ierr) {
  int N = n;
  int K = Kode;
  double  order = ord;

  F77NAME(zbesj) (x, y, &order, &K, &N, ReJ, ImJ, nz, ierr);
  if (*ierr) message ("initial.y", "zbesj", ERROR);
}

static double ReJn (double n, double x,  double y)
{
  int nz, ierr;
  double  rej, imj;
  
  zbesj (&x,&y,n,1,1,&rej,&imj,&nz,&ierr);
  return rej;
}

static double ImJn (double n, double x, double y)
{
  int nz, ierr;
  double  rej, imj;
  
  zbesj (&x,&y,n,1,1,&rej,&imj,&nz,&ierr);
  return imj;
}

static double Womersley (double A,
			 double B,
			 double r,
			 double R, 
			 double mu, 
			 double wnum,
			 double t)
/* ------------------------------------------------------------------------- *
 * Calculate the Womersley solution at r for a pipe of radius R and
 * wave number wnum.  The solution is assumed to be set so that the
 * spatial mean of the flow satisfies
 *   u_avg(r) = A cos (wnum t) + B sin(wnum t) 
 * ------------------------------------------------------------------------- */
{
  double x,y;

  if (r > R) message ("Womersley", "r > R", ERROR);

  if (wnum == 0) /* Return Poiseuille flow with mean of 1. */
    return 2*(1-r*r/R/R);
  else {
    int ierr,nz;
    double  cr,ci,J0r,J0i,rej,imj,re,im,fac;
    double  isqrt2 = 1.0/sqrt(2.0);
    static  double R_str, wnum_str,mu_str;
    static  double Jr,Ji,alpha,j0r,j0i, isqrt;

    /* For case of repeated calls with same parameters, 
     * store those independent of r.
     */

    if ((R != R_str) || (wnum != wnum_str) || (mu != mu_str)) {
      double retmp[2],imtmp[2];
      alpha = R*sqrt(wnum/mu);

      re  = -alpha*isqrt2;
      im  =  alpha*isqrt2;
      zbesj(&re,&im,0,1,2,retmp,imtmp,&nz,&ierr);
      j0r = retmp[0]; j0i = imtmp[0];
      rej = retmp[1]; imj = imtmp[1];

      fac = 1/(j0r*j0r+j0i*j0i);
      Jr = 1+2*fac/alpha*((rej*j0r+imj*j0i)*isqrt2-(imj*j0r - rej*j0i)*isqrt2);
      Ji = 2*fac/alpha*((rej*j0r+imj*j0i)*isqrt2 + (imj*j0r - rej*j0i)*isqrt2);

      R_str = R; wnum_str = wnum; mu_str = mu;
    }

    /* setup cr, ci from pre-stored value of Jr & Ji */

    fac = 1/(Jr*Jr + Ji*Ji);
    cr  =  (A*Jr - B*Ji)*fac;
    ci  = -(A*Ji + B*Jr)*fac;
    
    /* setup J0r, J0i */

    re  = -alpha*isqrt2*r/R;
    im  =  alpha*isqrt2*r/R;
    zbesj(&re,&im,0,1,1,&rej,&imj,&nz,&ierr);
    fac = 1/(j0r*j0r+j0i*j0i);
    J0r = 1-fac*(rej*j0r+imj*j0i);
    J0i = -fac*(imj*j0r-rej*j0i);

    return (cr*J0r - ci*J0i)*cos(wnum*t) - (ci*J0r + cr*J0i)*sin(wnum*t);
  }
}

static double Womsin (double r, double R, double mu, double wnum, double t) {
  if (wnum == 0) return 0.0; /* No sin term for zeroth mode. */
  
  return Womersley (0.0, 1.0, r, R, mu, wnum, t);
}

static double Womcos (double r, double R, double mu, double wnum, double t) {
  return Womersley (1.0, 0.0, r, R, mu, wnum, t);
}

/* ------------------------------------------------------------------------- *
 * Error checking is done for cases where we can have illegal input
 * values (typically negative), otherwise we accept exception returns.
 *
 * 1/6/2001: After much trouble with underflow errors, I've turned these off.
 * We just accept whatever comes back (typically a nan).
 *
 * example: static double Log (double x) { return errcheck (log (x), "log"); }
 * ------------------------------------------------------------------------- */


static double errcheck (const double d,
			const char*  s)
/* ------------------------------------------------------------------------- *
 * Check result of math library call. We don't use this anymore!
 * ------------------------------------------------------------------------- */
{
  if (errno == EDOM) {
    errno = 0;
    message ("errcheck: argument out of domain in call to", s, ERROR);
  } else if (errno == ERANGE) {
    errno = 0;
    message ("errcheck: result out of range in call to", s, ERROR);
  }

  return d;
}

#undef HASHSIZE
#undef HASHSEED
#undef VEC_MAX
