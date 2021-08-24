#include "Polylib.h"
#include <cstdio>
#include <cmath>

using namespace std;
using namespace Polylib;

/* --------------------------------------------------------------------
   To compile:
   g++ -g -c Polylib.cpp -I ./ -I ../../
   g++ -g -c Polylib_test.cpp -I ./ -I ../../
   g++ -g -o polytest Polylib_test.o Polylib.o -lm
 * --------------------------------------------------------------------*/

/* -------------------------------------------------------------------
   This is a routine to test the integration, differentiation and
   interpolation routines in the polylib.c.

   First, it performs the integral

      /1      alpha   beta  alpha,beta
     |   (1-x)   (1+x)     P (x)       dx  = 0
     /-1                    n

   for all   -0.5 <= alpha <= 5   (increments of 0.5)
             -0.5 <= beta  <= 5   (increments of 0.5)

   using np points where
          NPLOWER <= np <= NPUPPER
                2     <= n  <= 2*np - delta

   delta = 1 (gauss), 2(radau), 3(lobatto).
   The integral is evaluated and if it is larger that EPS then the
   value of alpha,beta,np,n and the integral is printed to the screen.

   After every alpha value the statement
       "finished checking all beta values for alpha = #"
   is printed

   The routine then evaluates the derivate of

          d   n      n-1
      -- x  = n x
      dx

   for all   -0.5 <= alpha <= 5   (increments of 0.5)
             -0.5 <= beta  <= 5   (increments of 0.5)

   using np points where
          NPLOWER <= np <= NPUPPER
                2     <= n  <= np - 1

   The error is check in a pointwise sense and if it is larger than
   EPS then the value of alpha,beta,np,n and the error is printed to
   the screen. After every alpha value the statement
       "finished checking all beta values for alpha = #"
   is printed

   Finally the routine  evaluates the interpolation of

             n      n
        z  to  x

   where z are the quadrature zeros and x are the equispaced points

                  2*i
        x    =   -----   - 1.0    (0 <= i <= np-1)
     i       (np-1)


   for all   -0.5 <= alpha <= 5   (increments of 0.5)
             -0.5 <= beta  <= 5   (increments of 0.5)

   using np points where
          NPLOWER <= np <= NPUPPER
                2     <= n  <= np - 1

   The error is check in a pointwise sense and if it is larger than
   EPS then the value of alpha,beta,np,n and the error is printed to
   the screen. After every alpha value the statement
      "finished checking all beta values for alpha = #"
   is printed

   The above checks are performed for all the Gauss, Gauss-Radau and
   Gauss-Lobatto points. If you want to disable any routine then set
      GAUSS_INT, GAUSS_RADAU_INT, GAUSS_LOBATTO_INT = 0
   for the integration rouintes
      GAUSS_DIFF,GAUSS_RADAU_DIFF, GAUSS_LOBATTO_DIFF = 0
   for the differentiation routines
      GAUSS_INTERP,GAUSS_RADAU_INTERP, GAUSS_LOBATTO_INTERP = 0
   for the interpolation routines.
   ------------------------------------------------------------------*/

#define NPLOWER  5
#define NPUPPER 15
#define EPS  1e-12

#define GAUSS_INT            1
#define GAUSS_RADAUM_INT     1
#define GAUSS_RADAUP_INT     1
#define GAUSS_LOBATTO_INT    1
#define GAUSS_DIFF           1
#define GAUSS_RADAUM_DIFF    1
#define GAUSS_RADAUP_DIFF    1
#define GAUSS_LOBATTO_DIFF   1
#define GAUSS_INTERP         1
#define GAUSS_RADAUM_INTERP  1
#define GAUSS_RADAUP_INTERP  1
#define GAUSS_LOBATTO_INTERP 1

/* local routines */
double    ddot (int, double *, int, double *, int);
double   *dvector (int, int);
int      test_gamma_fraction();

int main(){
  int np,n,i;
  double *z,*w,*p,sum=0,alpha,beta,*d,*dt;

  z  = dvector(0,NPUPPER-1);
  w  = dvector(0,NPUPPER-1);
  p  = dvector(0,NPUPPER-1);

  d  = dvector(0,NPUPPER*NPUPPER-1);
  dt = dvector(0,NPUPPER*NPUPPER-1);

#if GAUSS_INT
  /* Gauss Integration */
  printf("Begin checking Gauss integration\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){

      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwgj(z,w,np,alpha,beta);
    for(n = 2; n < 2*np-1; ++n){
      jacobfd(np,z,p,NULL,n,alpha,beta);
      sum = ddot(np,w,1,p,1);
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d integal was %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }

      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
      }
  printf("Finished checking Gauss Integration\n");
#endif

#if GAUSS_RADAUM_INT
  /* Gauss Radau Integration */
  printf("Begin checking Gauss Radau Integration\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){
      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwgrjm(z,w,np,alpha,beta);
    for(n = 2; n < 2*np-2; ++n){
      jacobfd(np,z,p,NULL,n,alpha,beta);
      sum = ddot(np,w,1,p,1);
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d integal was %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }

      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Radau (z=-1) Integration\n");
#endif


#if GAUSS_RADAUP_INT
  /* Gauss Radau Integration */
  printf("Begin checking Gauss Radau Integration\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){
      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwgrjp(z,w,np,alpha,beta);
    for(n = 2; n < 2*np-2; ++n){
      jacobfd(np,z,p,NULL,n,alpha,beta);
      sum = ddot(np,w,1,p,1);
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d integal was %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }

      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Radau (z=1) Integration\n");
#endif

#if GAUSS_LOBATTO_INT
  /* Gauss Lobatto Integration */
  printf("Begin checking Gauss Lobatto integration\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){

      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwglj(z,w,np,alpha,beta);
    for(n = 2; n < 2*np-3; ++n){
      jacobfd(np,z,p,NULL,n,alpha,beta);
      sum = ddot(np,w,1,p,1);
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d integal was %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }

      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Lobatto Integration\n");
#endif

#if GAUSS_DIFF
  printf("Begin checking differentiation through Gauss points\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){

      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwgj(z,w,np,alpha,beta);
    for(n = 2; n < np-1; ++n){
      Dgj(d,z,np,alpha,beta);

      for(i = 0; i < np; ++i) p[i] = pow(z[i],n);
      sum = 0;
      for(i = 0; i < np; ++i)
        sum += fabs(ddot(np,d+i*np,1,p,1) - n*pow(z[i],n-1));
      sum /= np;
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d difference %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }
      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Jacobi differentiation\n");
#endif

#if GAUSS_RADAUM_DIFF
  printf("Begin checking differentiation through Gauss Radau points\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){

      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwgrjm(z,w,np,alpha,beta);
    for(n = 2; n < np-1; ++n){
      Dgrjm(d,z,np,alpha,beta);
      for(i = 0; i < np; ++i) p[i] = pow(z[i],n);
      sum = 0;
      for(i = 0; i < np; ++i)
        sum += fabs(ddot(np,d+i*np,1,p,1) - n*pow(z[i],n-1));
      sum /= np;
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d difference %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }
      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Radau (z=-1) differentiation\n");
#endif

#if GAUSS_RADAUP_DIFF
  printf("Begin checking differentiation through Gauss Radau (z=1) points\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){

      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwgrjp(z,w,np,alpha,beta);
    for(n = 2; n < np-1; ++n){
      Dgrjp(d,z,np,alpha,beta);
      for(i = 0; i < np; ++i) p[i] = pow(z[i],n);
      sum = 0;
      for(i = 0; i < np; ++i)
        sum += fabs(ddot(np,d+i*np,1,p,1) - n*pow(z[i],n-1));
      sum /= np;
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d difference %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }
      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Radau (z=1) differentiation\n");
#endif

#if GAUSS_LOBATTO_DIFF
  printf("Begin checking differentiation through Gauss Lobatto points\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){

      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwglj(z,w,np,alpha,beta);
    for(n = 2; n < np-1; ++n){
      Dglj(d,z,np,alpha,beta);
      for(i = 0; i < np; ++i) p[i] = pow(z[i],n);
      sum = 0;
      for(i = 0; i < np; ++i)
        sum += fabs(ddot(np,d+i*np,1,p,1) - n*pow(z[i],n-1));
      sum /= np;
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d difference %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }
      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Lobatto differentiation\n");
#endif

  /* check interpolation routines */
#if GAUSS_INTERP
  printf("Begin checking interpolation through Gauss points\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){

      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwgj(z,w,np,alpha,beta);
    for(n = 2; n < np-1; ++n){
      for(i = 0; i < np; ++i) {
        w[i] = 2.0*i/(double)(np-1)-1.0;
        p[i] = pow(z[i],n);
      }
      Imgj(d,z,w,np,np,alpha,beta);
      sum = 0;
      for(i = 0; i < np; ++i)
        sum += fabs(ddot(np,d+i*np,1,p,1) - pow(w[i],n));
      sum /= np;
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d difference %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }
      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Jacobi interpolation\n");
#endif

#if GAUSS_RADAUM_INTERP
  printf("Begin checking Interpolation through Gauss Radau (z=-1) points\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){

      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwgrjm(z,w,np,alpha,beta);
    for(n = 2; n < np-1; ++n){
      for(i = 0; i < np; ++i) {
        w[i] = 2.0*i/(double)(np-1)-1.0;
        p[i] = pow(z[i],n);
      }
      Imgrjm(d,z,w,np,np,alpha,beta);
      sum = 0;
      for(i = 0; i < np; ++i)
        sum += fabs(ddot(np,d+i*np,1,p,1) - pow(w[i],n));
      sum /= np;
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d difference %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }
      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Radua Jacobi (z=-1) interpolation\n");
#endif
#if GAUSS_RADAUP_INTERP
  printf("Begin checking Interpolation through Gauss Radau (z=1) points\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){

      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwgrjp(z,w,np,alpha,beta);
    for(n = 2; n < np-1; ++n){
      for(i = 0; i < np; ++i) {
        w[i] = 2.0*i/(double)(np-1)-1.0;
        p[i] = pow(z[i],n);
      }
      Imgrjp(d,z,w,np,np,alpha,beta);
      sum = 0;
      for(i = 0; i < np; ++i)
        sum += fabs(ddot(np,d+i*np,1,p,1) - pow(w[i],n));
      sum /= np;
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d difference %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }
      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Radau (z=1) interpolation\n");
#endif

#if GAUSS_LOBATTO_INTERP
  printf("Begin checking Interpolation through Gauss Lobatto points\n");
  alpha = -0.5;
  while(alpha <= 5.0){
    beta = -0.5;
    while(beta <= 5.0){

      for(np = NPLOWER; np <= NPUPPER; ++np){
    zwglj(z,w,np,alpha,beta);
    for(n = 2; n < np-1; ++n){
      for(i = 0; i < np; ++i) {
        w[i] = 2.0*i/(double)(np-1)-1.0;
        p[i] = pow(z[i],n);
      }
      Imglj(d,z,w,np,np,alpha,beta);
      sum = 0;
      for(i = 0; i < np; ++i)
        sum += fabs(ddot(np,d+i*np,1,p,1) - pow(w[i],n));
      sum /= np;
      if(fabs(sum)>EPS)
        printf("alpha = %lf, beta = %lf, np = %d, n = %d difference %lg\n"
           ,alpha,beta,np,n,sum);
    }
      }
      beta += 0.5;
    }
    printf("finished checking all beta values for alpha = %lf\n",alpha);
    alpha += 0.5;
  }
  printf("Finished checking Gauss Lobatto interploation\n");
#endif
  test_gamma_fraction();


  free(z); free(w); free(p); free(d); free(dt);
  return 0;
}

double ddot (int n, double *x, int incx, double *y, int incy)
{
  register double sum = 0.;

  while (n--) {
    sum += (*x) * (*y);
    x   += incx;
    y   += incy;
  }
  return sum;
}


double *dvector(int nl,int nh)
{
  double *v;

  v = (double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  return v-nl;
}

int test_gamma_fraction(){
  double a = 362880.;
  double c = Polylib::gammaFracGammaF(10, 0., 0, 1.);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 10.,1., c/a-1.);
  c = Polylib::gammaFracGammaF(1, 0., 12, -2.);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 1.,10., c*a-1.);

  a = 30.;
  c = Polylib::gammaFracGammaF(4, 3., 5, 0.);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 5.,7., c/a-1.);
  c = Polylib::gammaFracGammaF(5, 0., 7, 0.);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 7.,5., c*a-1.);

  a = 21651.09375;
  c = Polylib::gammaFracGammaF(7, 3.5, 5, 0.5);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 10.5,5.5, c/a-1.);
  c = Polylib::gammaFracGammaF(5, 0.5, 10, 0.5);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 5.5,10.5, c*a-1.);

  a = 97429.921875;
  c = Polylib::gammaFracGammaF(10, 0.5, 5, -0.5);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 10.5,4.5, c/a-1.);
  c = Polylib::gammaFracGammaF(5, -0.5, 11, -0.5);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 4.5,10.5, c*a-1.);

  a = 2279.0625;
  c = Polylib::gammaFracGammaF(10, -0.5, 5, 0.5);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 9.5,5.5, c/a-1.);
  c = Polylib::gammaFracGammaF(5, 0.5, 10, -0.5);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 5.5,9.5, c*a-1.);

  a = 639383.8623046875;
  c = Polylib::gammaFracGammaF(10, 0.5, 0, 0.5);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 10.5,0.5, c/a-1.);
  c = Polylib::gammaFracGammaF(0, 0.5, 10, 0.5);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 0.5,10.5, c*a-1.);

  a = 5967498288235982848.;
  c = Polylib::gammaFracGammaF(100, 0., 90, 0.5);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 100.,90.5, c/a-1.);
  c = Polylib::gammaFracGammaF(90, 0.5, 100, 0.);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 90.5,100., c*a-1.);

  a = 5618641603777298169856.;
  c = Polylib::gammaFracGammaF(200, 0., 191, -0.5);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 200.,190.5, c/a-1.);
  c = Polylib::gammaFracGammaF(190, 0.5, 200, 0.);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 190.5,200., c*a-1.);

  a = 77396694214720029196288.;
  c = Polylib::gammaFracGammaF(200, 0., 190, 0.);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 200.,190., c/a-1.);
  c = Polylib::gammaFracGammaF(190, 0., 200, 0.);
  printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", 190.,200., c*a-1.);

  int Q=2;
  for(double alpha=-0.5; alpha<=150.; alpha+=0.5)
  {
      for(double beta=-0.5; beta<=150.; beta+=0.5)
      {
          a = Polylib::gammaF(Q+alpha)/Polylib::gammaF(Q+beta);
          c = Polylib::gammaFracGammaF(Q, alpha, Q, beta) / a - 1.;
          if(fabs(c)>5.e-15)
          {
              printf("alpha = %7.2lf, beta = %7.2lf, difference %9.2e\n", Q+alpha,Q+beta, c);
          }
      }
  }

  return 0;
}