#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

#include "fmatrix.h"
#include "configuration.h"

void load_cross_sections(configuration & config);
void delete_cross_sections_arrays();

// e-I cross sections
extern fmatrix elastic_cs;
extern fmatrix ionization_cs;
extern fmatrix* excitation_cs;

// I+-I cross sections
extern fmatrix isotropic_cs;
extern fmatrix backscattering_cs;

#endif // !CROSS_SECTIONS_H


