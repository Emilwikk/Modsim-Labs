#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"
extern int nnsize;


int monte_carlo_sweep(Par *par, double *pos) {
  int i, j, d, naccept = 0;
  double ediff, change;
  double newpos[D], dist[D];

  for (i = 0; i < par->n; i++){
    ediff = 0.0;
    for (d = 0; d < D; d++){
      change = par->b * dran_sign();
      // Apply periodic boundary condition and update store new position
      newpos[d] = inrange(pos[i * D + d] + change, par->L[d]);
    }
    for (j = 0; j < par->n; j++){
        if (j == i){
            continue;
        } 

        // Calc (short for calculate) distances
        double dist2_new = dist2(par->L, newpos, &pos[j * D], dist);
        double dist2_old = dist2(par->L, &pos[i * D], &pos[j * D], dist);

        // Interactions happen only if particle distace is smaller than cutoff distance
        if(dist2_new < CUT*CUT){
            ediff += pair_interaction(dist2_new);                
        }
        if(dist2_old < CUT*CUT){
            ediff -= pair_interaction(dist2_old); 
        }
        
        
    }
    double alpha = fmin(1,exp(-ediff/par->T));
    double xi = dran();
    if (xi < alpha) {
        // Update position of particle i if condition is met
        for (d = 0; d < D; d++){
            pos[i * D + d] = newpos[d];
        }
        // Count successful moves
        naccept++;
    }
  }

  #ifdef FAST
    check_neighbor_list(par, pos);
  #endif


  return naccept;
}
