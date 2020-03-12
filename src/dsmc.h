#ifndef DSMC_H
#define DSMC_H

#include "fmatrix.h"
#include "configuration.h"

void run_dsmc(fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & vmesh, fmatrix & dens_n, configuration & config, string output_name);

#endif