#include "coeffs.h"
#include "string.h"

#define DIM 13
#include <string.h>

void fill_coeffs_(float* in){ memcpy(in,coeffs,DIM*DIM*sizeof(float)); }
