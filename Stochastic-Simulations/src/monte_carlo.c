#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"
extern int nnsize;


int monte_carlo_sweep(Par *par, double *pos) {
  int i, j, d, naccept = 0;
  double r2, ediff, delta;
  double newpos[D], dist[D];

  for (i = 0; i < par->n; i++){
    for (d = 0; d < D; d++){
      delta = par->b * dran_sign();
      newpos[d] = pos[i * D + d] + delta;
    }
  for (j = 0; d < par->n; i++){
    if (j == i){
      ediff = ediff;
    } else {
      double dist_new = dist2(par->L, newpos, &pos[j * D], dist);
      double dist_old = dist2(par->L, &pos[i * D], &pos[j * D], dist);
      ediff += pair_interaction(dist_new)-pair_interaction(dist_old); 
    }
    double alpha = fmin(1,exp(-ediff/par->T));
    double xi = dran();
    if (xi < alpha) {
      for (d = 0; d < D; d++){
      pos[i * D + d] = newpos[d];
      naccept++;
      }
    }
    
  }
}

#ifdef FAST
  check_neighbor_list(par, pos);
#endif


  return naccept;
}
