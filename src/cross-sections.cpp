
#include "cross-sections.h"

#include <string>

#include "fmatrix.h"
#include "input-output.h"
#include "configuration.h"


// e-He cross sections
fmatrix elastic_cs;
fmatrix ionization_cs;
fmatrix* excitation_cs;

// He+-He cross sections
fmatrix isotropic_cs;
fmatrix backscattering_cs;

void load_cross_sections(configuration & config){
    // e-I cross sections

    tmatrix<string> exc_path = config.ss("ugas/exc_path");
    int n_exc = config.i("p/n_exc");


    elastic_cs = load_csv(config.s("project/cross_sections_path") + config.s("ugas/els_path"));
    
    ionization_cs = load_csv(config.s("project/cross_sections_path") + config.s("ugas/iz_path"));
    
    excitation_cs = new fmatrix[n_exc];
    for(int i = 0; i < n_exc; i++) excitation_cs[i] = load_csv(config.s("project/cross_sections_path") + exc_path.val[i]);
    
    // I+-I cross sections
    isotropic_cs = load_csv(config.s("project/cross_sections_path") + config.s("ugas/iso_path"));
    backscattering_cs = load_csv(config.s("project/cross_sections_path") + config.s("ugas/bs_path"));
}

void delete_cross_sections_arrays(){
    delete[] excitation_cs;
    // delete[] E_EXC_PATH;
    // delete[] E_EXC;
}