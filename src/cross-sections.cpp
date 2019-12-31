
#include "cross-sections.h"

#include <string>

#include "fmatrix.h"
#include "input-output.h"
#include "config.h"


// e-He cross sections
fmatrix elastic_cs;
fmatrix ionization_cs;
fmatrix excitation_1_cs;
fmatrix excitation_2_cs;
fmatrix* excitation_cs;

// He+-He cross sections
fmatrix isotropic_cs;
fmatrix backscattering_cs;

void load_cross_sections(){
    // e-I cross sections
    elastic_cs = load_csv(CROSS_SECTIONS_PATH + ELS_PATH);
    ionization_cs = load_csv(CROSS_SECTIONS_PATH + IZ_PATH);
    excitation_cs = new fmatrix[N_EXC];
    for(int i = 0; i < N_EXC; i++) excitation_cs[i] = load_csv(CROSS_SECTIONS_PATH + E_EXC_PATH[i]);
    
    // I+-I cross sections
    isotropic_cs = load_csv(CROSS_SECTIONS_PATH + ISO_PATH);
    backscattering_cs = load_csv(CROSS_SECTIONS_PATH + BS_PATH);
}

void delete_cross_sections_arrays(){
    delete[] excitation_cs;
    delete[] E_EXC_PATH;
    delete[] E_EXC;
}