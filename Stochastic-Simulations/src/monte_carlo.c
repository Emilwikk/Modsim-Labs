#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"
extern int nnsize;


int monte_carlo_sweep(Par *par, double *pos) {
  int i, j, d, naccept = 0;
  double r2, ediff;
  double newpos[D], dist[D];

#ifdef FAST
  check_neighbor_list(par, pos);
#endif

						 
  return naccept;
}
