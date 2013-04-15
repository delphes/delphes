//! \file ranlux.h

#ifndef __RANLUX_H__
#define __RANLUX_H__

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

/// initialize 'ranlux' generator
void ranlux_init();

/// generate random value (24 bits)
unsigned long int ranlux_get();

/// save state of the generator
void ranlux_print_state();

}
#endif
