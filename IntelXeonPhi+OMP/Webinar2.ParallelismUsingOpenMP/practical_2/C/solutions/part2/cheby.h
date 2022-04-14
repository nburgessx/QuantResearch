/* Pre-processor stuff that needn't be changed */

#include <time.h>
#include <omp.h>
#include <math.h>

#define PI 3.14159265358979323846



/* Handle Debugging */
#define DEBUGPRINT(...) if (DEBUG) printf(__VA_ARGS__)



/* Handle Timing */
#if _OPENMP
#define TIMER_T double
#else
#define TIMER_T clock_t
#endif

#if _OPENMP
#define HAVE_OMP 1
#else
#define HAVE_OMP 0
#endif

#define TIME() HAVE_OMP ? omp_get_wtime() : clock()
#define TIMER_RES(t1,t2) HAVE_OMP ? t2-t1 : (float)(t2-t1)/CLOCKS_PER_SEC



/* Handle Precision */
#if USE_DOUBLES

typedef double real_t;
#define SIN sin
#define COS cos
#define EXP exp
#define SQRT sqrt
#define EPS eps

#else

typedef float real_t;
#define SIN sinf
#define COS cosf
#define EXP expf
#define SQRT sqrtf
#define EPS epsf

#endif /* end if USE_DOUBLES */


