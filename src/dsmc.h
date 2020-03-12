#ifndef DSMC_H
#define DSMC_H

#include "fmatrix.h"
#include "configuration.h"
#include "fields.h"

void run_dsmc(mesh_set & mesh, fmatrix & dens_n, configuration & config, string output_name);

#endif