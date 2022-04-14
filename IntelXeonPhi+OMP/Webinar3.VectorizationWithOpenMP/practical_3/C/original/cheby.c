/*
CWA - Feb 2014
--------------
Generate some (xbar, fr) values of a function, evaluate the
coefficients of a Chebyshev polynomial which interpolates
the generated values exactly at the given points (fe), evaluate
the polynomial at the given points and then find the Euclidean
norm of the difference between the original fr and the evaluated
fr.

Serial implementation.

*/

#include <stdlib.h>
#include <stdio.h>



#define DEBUG 0
#define USE_DOUBLES 1



#include "cheby.h"
#include "eps.h"



void find_Tcoeffs(real_t *, real_t *, int);
void evaluate_Tpoly(real_t *, real_t *, int, real_t *, int);
void gen_data(real_t *, real_t *, int);



int main(int argc, char **argv)
{
   int i, j;
   int n1;
   real_t *xbar, *fr, *coeffs, *fe, var, dif, macheps, nrm2, tol;
   TIMER_T t1, t2;

   for (i=1; i<argc; i++)
   {

      n1 = atoi(argv[i]);
      if (n1<2) { fprintf(stderr,"Error, n1 must be > 1.\n"); continue; }

      /* Allocate space for the generated function and 
         x-values */
      xbar = (real_t *)malloc(sizeof(real_t)*n1);
      fr = (real_t *)malloc(sizeof(real_t)*n1);

      /* Generate some data from a function */
      gen_data(fr, xbar, n1);

      printf("Input %d, size = %d\n",i,n1);

      DEBUGPRINT("Input data (x, f(x)):\n");
      for (j=0; DEBUG && j<n1; j++)
         DEBUGPRINT("%d %f %f\n",j,xbar[j],fr[j]);

      /* Allocate space for the evaluated function and 
         the coefficients */
      fe = (real_t *)malloc(sizeof(real_t)*n1);
      coeffs = (real_t *)malloc(sizeof(real_t)*n1);

      /* Find the Chebyshev coefficients */
      t1 = TIME();
      find_Tcoeffs(coeffs, fr, n1);
      t2 = TIME();
      printf("Time to find coefficients = %f secs.\n",TIMER_RES(t1,t2));

      DEBUGPRINT("Chebyshev coefficients:\n");
      for (j=0; DEBUG && j<n1; j++)
         DEBUGPRINT("%e\n",coeffs[j]);

      /* Evaluate using full poly - should be exact */
      t1 = TIME();
      evaluate_Tpoly(fe, xbar, n1, coeffs, n1);
      t2 = TIME();
      printf("Time to evaluate poly = %f secs.\n",TIMER_RES(t1,t2));

      DEBUGPRINT("Original and evaluated values:\n");
      var = 0.0;
#pragma omp parallel for schedule(static) reduction(+:var)
      for (j=0; j<n1; j++)
      {
         DEBUGPRINT("%f %f\n",fr[j],fe[j]);
         dif = fr[j]-fe[j];
         var = var + dif*dif;
      }

      nrm2 = SQRT(var);
      EPS(&var,&macheps);
      tol = (real_t)(n1)*(real_t)(n1)*macheps;
      if (nrm2<=tol)
         printf("Result is OK (Euclidean norm: %e)\n",nrm2);
      else
         printf("Error:\n    Euclidean norm %e greater than tolerance %e\n",
                 nrm2,tol);

      free(xbar);
      free(fr);
      free(fe);
      free(coeffs);
   }
   return 0;
}



void evaluate_Tpoly(real_t *fe, real_t *xbar, int nx, real_t *coeffs, int nplus1)
{
   int i, n, k;
   real_t factor, dk, bk;

   n = nplus1 - 1;
  
#pragma omp parallel for schedule(static) 
   for (i=0; i<nx; i++)
   {
      factor = 2.0*(1.0 + xbar[i]);
      dk = 0.0;
      bk = 0.0;
      for (k=n; k>0; k--)
      {
         dk = coeffs[k] - dk + factor*bk;
         bk = dk - bk;
      }
      fe[i] = 0.5*coeffs[0] - dk + 0.5*factor*bk;
   }
}



void find_Tcoeffs(real_t *coeffs, real_t *fr, int nplus1)
{
   int i, k;
   int n, nless1;
   real_t fln, piby2n, f0, halffn, fli, factor, bk, dk;
   if (nplus1<2) return;
   if (nplus1==2)
   {
      coeffs[0] = fr[0] + fr[1];
      coeffs[1] = 0.5*(fr[0] - fr[1]);
      return;
   }
   n = nplus1 - 1;
   fln = (real_t)n;
   nless1 = n - 1;
   piby2n = 0.5*PI/fln;
   f0 = fr[0];
   halffn = 0.5*fr[n];

#pragma omp parallel for schedule(static)
   for (i=0; i<nplus1; i++)
   {
      fli = (real_t)i;
      factor = SIN(piby2n*fli);
      factor = 4.0*factor*factor;
      dk = halffn;
      bk = halffn;
      for (k=nless1; k>0; k--)
      {
         dk = fr[k] + dk - factor*bk;
         bk = bk + dk;
      }
      coeffs[i] = (f0+2.0*dk-factor*bk)/fln;
   }
   coeffs[n] = 0.5*coeffs[n];
}



void gen_data(real_t *f, real_t *x, int n1)
{
   int n, i;

   n = n1 - 1;

#pragma omp parallel for schedule(static)
   for (i=0; i<n1; i++)
   {
      x[i] = COS(PI*(real_t)i/(real_t)n);
      f[i] = EXP(x[i]);
   }
}
