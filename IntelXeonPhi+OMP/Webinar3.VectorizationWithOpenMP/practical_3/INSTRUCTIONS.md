# Vectorizing Code with OpenMP


This exercise demonstrates how code modification may be necessary to assist the compiler
in using efficient vector instructions.

This directory contains an OpenMP version of the Chebyshev program used in Exercise 2. 
As in the previous exercise: make a new directory as a clean copy of the original
before making changes so that you have a good point of comparison.


## Part 1: Building with Different Options

The Makefile is set up to compile an executable at –O3. 
Verify that when compiled at –O2 the program is slower. 
Use a problem size of 500000 (i.e. 10x that used in Practical 2). 

Inspect the compiler output when using the `–qopt-report=5 -qopt-report-phase=vec` flags. 

The file of interest here will be cheby.optrpt, which will give compiler feedback
on what could be vectorized and what could not be vectorized. For example from cheby.c
the code

```C
182 #pragma omp parallel for schedule(static)
183    for (i=0; i<n1; i++)
184    {
185       x[i] = COS(PI*(real_t)i/(real_t)n);
186       f[i] = EXP(x[i]);
187    }
```

was automatically vectorized by the compiler, and it provided the following output in
cheby.optrpt

```
212 LOOP BEGIN at cheby.c(183,4)
213 <Multiversioned v1>
214    remark #15389: vectorization support: reference x has unaligned access   [ cheby.c(185,7) ]
215    remark #15389: vectorization support: reference x has unaligned access   [ cheby.c(186,14) ]
216    remark #15389: vectorization support: reference f has unaligned access   [ cheby.c(186,7) ]
217    remark #15381: vectorization support: unaligned access used inside loop body
218    remark #15305: vectorization support: vector length 8
219    remark #15309: vectorization support: normalized vectorization overhead 0.089
220    remark #15300: LOOP WAS VECTORIZED
221    remark #15442: entire loop may be executed in remainder
222    remark #15450: unmasked unaligned unit stride loads: 1
223    remark #15451: unmasked unaligned unit stride stores: 2
224    remark #15475: --- begin vector loop cost summary ---
225    remark #15476: scalar loop cost: 245
226    remark #15477: vector loop cost: 23.870
227    remark #15478: estimated potential speedup: 9.590
228    remark #15482: vectorized math library calls: 2
229    remark #15487: type converts: 1
230    remark #15488: --- end vector loop cost summary ---
231 LOOP END
```


Note that not all loops are automatically vectorized, for example one of the 
most expensive loops in cheby.c:

```C
156 #pragma omp parallel for schedule(static)
157    for (i=0; i<nplus1; i++)
158    {
159       fli = (real_t)i;
160       factor = SIN(piby2n*fli);
161       factor = 4.0*factor*factor;
162       dk = halffn;
163       bk = halffn;
164       for (k=nless1; k>0; k--)
165       {
166          dk = fr[k] + dk - factor*bk;
167          bk = bk + dk;
168       }
169       coeffs[i] = (f0+2.0*dk-factor*bk)/fln;
170    }
```

fails to vectorize. The Intel compiler tries to give us some information we can use to potentially
reformulate this loop in a way that can be vectorized:


```
 99 LOOP BEGIN at cheby.c(157,4) inlined into cheby.c(70,7)
 100    remark #15344: loop was not vectorized: vector dependence prevents vectorization
 101    remark #15346: vector dependence: assumed OUTPUT dependence between k line 164 and k line 164
 102    remark #15346: vector dependence: assumed OUTPUT dependence between k line 164 and k line 164
 103
 104    LOOP BEGIN at cheby.c(164,7) inlined into cheby.c(70,7)
 105       remark #15344: loop was not vectorized: vector dependence prevents vectorization
 106       remark #15346: vector dependence: assumed FLOW dependence between  line 167 and  line 166
 107       remark #15346: vector dependence: assumed ANTI dependence between  line 166 and  line 167
 108    LOOP END
 109 LOOP END
```


Here the Intel compiler (correctly) detects a dependency which makes vectorization incorrect.
We can verify this by forcing vectorization with an OpenMP directive and seeing that the resulting output is
incorrect (computed error is well above threshold).

## Part 2: Enabling vectorization

The Intel compiler failed to vectorize here because as the code currently is, it would be incorrect to do so.
We can transform the code into an equivalent form which has no dependencies. The code transformation
I suggest here is known as "unroll-and-jam." The unroll-and-jam strategy takes a nested loop sequence such as


```C
for(int i=0;i<n;i++){
  for(int j=0;j<m;j++){
    //Loop carried dependencies here
    //May not vectorize!
  }
}
```


and changes it into 


```C
for(int i=0;i<n;i+=BLOCK){
  int iend=MIN(i+BLOCK,n); //In case BLOCK does not divide n.
  for(int j=0;j<m;j++){
    //Dependencies here,
    //May not vectorize!
    #pragma omp simd
    for(int ii=i;ii<iend;ii++){
      //but possibly no dependencies in
      //this new inner loop.
    }
  }
}
```


I will give an outline of how this can be done for cheby.c on the earlier example, and leave as an exercise for you to apply this code transformation
wherever Intel compiler could not vectorize, by following the vectorization report.


Remember to measure performance to make sure you actually got a benefit! You should see at least an 8X speedup on Knights Landing.




### Outline for one loop

Recall that the loop 

```C
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
```

would not vectorize because of the dependency in the inner loop. We can apply unroll-and-jam here to yield a
*new* inner loop which contains no dependencies. An outline of how this can be achieved follows:


```C
#define MIN(x,y) ( (x)<(y)?(x):(y) )
#define BLOCK 32 
#pragma omp parallel for schedule(static)
   for (i=0; i<nplus1; i+=BLOCK) // Increment by block parameter
                                // This is the "unroll" step.
   {

   //This cleanly handles the case
   //where BLOCK does not
   //divide nplus1 evenly
   int iend=MIN(i+BLOCK,nplus1);


// These variables will need to become
// arrays so that each vector lane
// can safely modify its own copy
//      fli = (real_t)i;
//      factor = SIN(piby2n*fli);
//      factor = 4.0*factor*factor;
//      dk = halffn;
//      bk = halffn;
  real_t flis[BLOCK];
  real_t factors[BLOCK];
  real_t dks[BLOCK];
  real_t bks[BLOCK];


//Initializing these now becomes a loop
#pragma omp simd
      for(int ii=i;ii<iend;ii++){
        flis[ii-i]=(real_t)ii;
        factors[ii-i]=4.0*SIN(piby2n*fli)*SIN(piby2n*fli);
        dks[ii-i]=halffn;
        bks[ii-i]=halffn;
      }

      for (k=nless1; k>0; k--)
      {

        //Now use the new inner loop inside
        //the old inner loop.
        //This is the "jam" part of unroll-and-jam
        //Because we "jam" the outer loop in the inner loop.

        #pragma omp simd
        for(int ii=i;ii<iend;ii++){
          //Don't forget to use our new arrays
          // dk = fr[k] + dk - factor*bk;
          // bk = bk + dk;
          dks[ii-i]=fr[k]+dks[ii-i]-factors[ii-i]*bks[ii-i];
          bks[ii-i]=bks[ii-i]+dks[ii-i];
         }
      }
      for(int ii=i;ii<iend;ii++){
        coeffs[ii] = (f0+2.0*dks[ii-i]-factors[ii-i]*bks[ii-i])/fln;
      }
   }
```


