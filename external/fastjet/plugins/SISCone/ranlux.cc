// file: ranlux.xpp
#include "ranlux.h"
#include <stdlib.h>
#include <stdio.h>

/* This is a lagged fibonacci generator with skipping developed by Luescher.
   The sequence is a series of 24-bit integers, x_n, 

   x_n = d_n + b_n

   where d_n = x_{n-10} - x_{n-24} - c_{n-1}, b_n = 0 if d_n >= 0 and
   b_n = 2^24 if d_n < 0, c_n = 0 if d_n >= 0 and c_n = 1 if d_n < 0,
   where after 24 samples a group of p integers are "skipped", to
   reduce correlations. By default p = 199, but can be increased to
   365.

   The period of the generator is around 10^171. 

   From: M. Luescher, "A portable high-quality random number generator
   for lattice field theory calculations", Computer Physics
   Communications, 79 (1994) 100-110.

   Available on the net as hep-lat/9309020 at http://xxx.lanl.gov/

   See also,

   F. James, "RANLUX: A Fortran implementation of the high-quality
   pseudo-random number generator of Luscher", Computer Physics
   Communications, 79 (1994) 111-114

   Kenneth G. Hamilton, F. James, "Acceleration of RANLUX", Computer
   Physics Communications, 101 (1997) 241-248

   Kenneth G. Hamilton, "Assembler RANLUX for PCs", Computer Physics
   Communications, 101 (1997) 249-253  */

namespace siscone{

static const unsigned long int mask_lo = 0x00ffffffUL;  // 2^24 - 1
static const unsigned long int mask_hi = ~0x00ffffffUL;
static const unsigned long int two24 = 16777216;        // 2^24


// internal generator structure
//------------------------------
typedef struct {
  unsigned int i;
  unsigned int j;
  unsigned int n;
  unsigned int skip;
  unsigned int carry;
  unsigned long int u[24];
} ranlux_state_t;


// internal generator state
//--------------------------
ranlux_state_t local_ranlux_state;


// incrementation of the generator state
//---------------------------------------
static inline unsigned long int increment_state(){
  unsigned int i = local_ranlux_state.i;
  unsigned int j = local_ranlux_state.j;
  long int delta = local_ranlux_state.u[j] - local_ranlux_state.u[i] 
    - local_ranlux_state.carry;

  if (delta & mask_hi){
    local_ranlux_state.carry = 1;
    delta &= mask_lo;
  } else {
    local_ranlux_state.carry = 0;
  }

  local_ranlux_state.u[i] = delta;
  
  if (i==0)
    i = 23;
  else
    i--;

  local_ranlux_state.i = i;

  if (j == 0)
    j = 23;
  else
    j--;

  local_ranlux_state.j = j;

  return delta;
}


// set generator state
//---------------------
static void ranlux_set(unsigned long int s){
  int i;
  long int seed;
  
  if (s==0)
    s = 314159265;      /* default seed is 314159265 */
  
  seed = s;
  
  /* This is the initialization algorithm of F. James, widely in use
     for RANLUX. */

  for (i=0;i<24;i++){
    unsigned long int k = seed/53668;
    seed = 40014*(seed-k*53668)-k*12211;
    if (seed<0){
      seed += 2147483563;
    }
    local_ranlux_state.u[i] = seed%two24;
  }

  local_ranlux_state.i = 23;
  local_ranlux_state.j = 9;
  local_ranlux_state.n = 0;
  local_ranlux_state.skip = 389-24; // 389 => best decorrelation

  if (local_ranlux_state.u[23]&mask_hi){
    local_ranlux_state.carry = 1;
  } else {
    local_ranlux_state.carry = 0;
  }
}


// generator initialization
//--------------------------
void ranlux_init(){
  // seed the generator
  ranlux_set(0);
}


// get random number
//-------------------
unsigned long int ranlux_get(){
  const unsigned int skip = local_ranlux_state.skip;
  unsigned long int r = increment_state();
  
  local_ranlux_state.n++;

  if (local_ranlux_state.n == 24){
    unsigned int i;
    local_ranlux_state.n = 0;
    for (i = 0; i < skip; i++)
      increment_state();
  }

  return r;
}

// print generator state
//-----------------------
void ranlux_print_state(){
  size_t i;
  unsigned char *p = (unsigned char *) (&local_ranlux_state);
  const size_t n = sizeof (ranlux_state_t);

  for (i=0;i<n;i++){
    /* FIXME: we're assuming that a char is 8 bits */
    printf("%.2x", *(p+i));
  }
}

}
