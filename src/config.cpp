
#include "config.h"
#include "INIReader.h"

#include <string>
#include <cmath>
#include <iostream>

using namespace std;

//[project]
 string OUTPUT_PATH;
 string CROSS_SECTIONS_PATH;

//; ----------------------------- Constants -------------------------------------
//[physical]
 double Q                      ;//Elementary charge [C]
 double EPS_0                  ;//Vacuum permittivity [F * m^-1]
 double PI                     ;//Pi []
 double K_BOLTZ                ;//Boltzmann constant [J * K^-1]

//; ----------------------------- Domain ----------------------------------------
//[geometry]
 double L_X                    ;//Lenght of simulation domain in x [m]
 double L_Y                    ;//Lenght of simulation domain in y [m]
 int N_MESH_X                  ;//Number of grid points in the x direction (n_benchmark + 1) []
 int N_MESH_Y                  ;//Number of grid points in the y direction []
 int N_THRUSTER                ;//Number of nodes that represent the thruster []
 double A_X                    ;//Grid distorsion coefficient in the x direction []
 double A_Y                    ;//Grid distorsion coefficient in the y direction []

//[time]
 double DT                     ;//Time step [s]
 int N_STEPS                   ;//Number of simulation steps to be executed []
 int N_STEPS_DSMC              ;//Number of DSMC simulation steps to be executed []
 int N_AVERAGE                 ;//Number of average steps []
 int N_AVERAGE_DSMC            ;//Number of DSMC average steps []
 int K_SUB                  ;//Ion subcycling factor []
 int K_SUB_DSMC             ;//DSMC subcycling factor []

//[boundary-conditions]
 double VOLT_0                 ;//Voltage on the thruster [V]
 double VOLT_1                 ;//Voltage on the chamber walls [V]

//; ----------------------------- Species ---------------------------------------

//[thruster]
 double MASS_FLOW_RATE        ;//Mass flow rate of thruster [sccm]
 double ETA_PROPELLANT        ;//Propellant utilization efficiency []

//[particles]
 double N_FACTOR              ;//Macro-particle weight factor (i.e. how many particles one macro-particle represents) []
 double N_FACTOR_DSMC         ;//Macro-particle weight factor for neutrals []
 int N_MAX_PARTICLES          ;//Maximum number of particles in the array []

//[neutral]
 string EXPM_NEUTRAL          ;//Expansion model of neutral flow: 'dsmc' or 'constant'
 double T_NEUTRAL             ;//Temperature of neutral gas [eV]
 double N_NEUTRAL             ;//Density of neutral gas [m^-3]

//[electrons]
 double M_EL                  ;//Electron mass [kg]
 double Q_EL                  ;//Electron charge [C]
 double T_EL                  ;//Initial electron temperature [eV]
 double I_EL                  ;//Electron injection current [A/m]

//[ions]
 string GAS_NAME              ;//Name of gas: 'helium' or 'xenon'
 double Q_I                   ;//Ion charge [kg]
 double T_I                   ;//Initial ion temperature [eV]
 double I_I                   ;//Ion injection current [A/m]
 double MACH_I                ;//Injection mach number []

 double M_I                   ;//Ion mass [kg]
 double* E_EXC                ;//First excitation energy [eV]
 double E_IZ                  ;//Ionization energy [eV]
 int N_EXC                    ;//Number of excited states considered (1 or 2)

 string ELS_PATH;             ;//Path to elastic cross section
 string* E_EXC_PATH;          ;//Paths to excitation cross sections
 string IZ_PATH               ;//Path ionization cross section
 string BS_PATH               ;//Path backscattering cross section
 string ISO_PATH              ;//Path isotropic cross section

// Calculated parameters
 double DX;
 double DY;
 double VOLT_0_NORM;
 double VOLT_1_NORM;
 double GAMMA;
 double ALPHA;
 double N_INJ_N;
 double N_INJ_I;
 double N_INJ_EL;
 double VD_I;


void load_config_file(string filename){

    
    INIReader reader(filename);

    // Project configuration
    OUTPUT_PATH = reader.Get("project", "OUTPUT_PATH", "");
    CROSS_SECTIONS_PATH = reader.Get("project", "CROSS_SECTIONS_PATH", "");

    // Physical constants
    Q = reader.GetReal("physical", "Q", -1);
    EPS_0 = reader.GetReal("physical", "EPS_0", -1);
    PI = reader.GetReal("physical", "PI", -1);
    K_BOLTZ = reader.GetReal("physical", "K_BOLTZ", -1);

    // Geometry
    L_X = reader.GetReal("geometry", "L_X", -1);
    L_Y = reader.GetReal("geometry", "L_Y", -1);
    N_MESH_X = (int) reader.GetReal("geometry", "N_MESH_X", -1);
    N_MESH_Y = (int) reader.GetReal("geometry", "N_MESH_Y", -1);
    N_THRUSTER = (int) reader.GetReal("geometry", "N_THRUSTER", -1);
    A_X = reader.GetReal("geometry", "A_X", -1);
    A_Y = reader.GetReal("geometry", "A_Y", -1);

    // Time
    DT = reader.GetReal("time", "DT", -1); 
    N_STEPS = (int) reader.GetReal("time", "N_STEPS", -1);
    N_AVERAGE = (int) reader.GetReal("time", "N_AVERAGE", -1);
    K_SUB =  (int) reader.GetReal("time", "K_SUB", -1);
    N_STEPS_DSMC = (int) reader.GetReal("time", "N_STEPS_DSMC", -1);
    N_AVERAGE_DSMC = (int) reader.GetReal("time", "N_AVERAGE_DSMC", -1);
    K_SUB_DSMC =  (int) reader.GetReal("time", "K_SUB_DSMC", -1);

    // Boundary conditions
    VOLT_0 = reader.GetReal("boundary-conditions", "VOLT_0", -1);
    VOLT_1 = reader.GetReal("boundary-conditions", "VOLT_1", -1);   
    
    // Plasma
    N_FACTOR = reader.GetReal("particles", "N_FACTOR", -1);
    N_MAX_PARTICLES = (int) reader.GetReal("particles", "N_MAX_PARTICLES", -1);

    
    // ---------------------- species ---------------------------

    // neutral
    T_NEUTRAL = reader.GetReal("neutral", "T_NEUTRAL", -1);
    N_NEUTRAL = reader.GetReal("neutral", "N_NEUTRAL", -1);
    EXPM_NEUTRAL = reader.Get("neutral", "EXPM_NEUTRAL", "");
    MASS_FLOW_RATE = reader.GetReal("thruster", "MASS_FLOW_RATE", -1);
    ETA_PROPELLANT = reader.GetReal("thruster", "ETA_PROPELLANT", -1);
    N_FACTOR_DSMC = reader.GetReal("particles", "N_FACTOR_DSMC", -1);

    // Particle 1 - Electrons
    M_EL = reader.GetReal("electrons", "M_EL", -1);
    Q_EL = reader.GetReal("electrons", "Q_EL", -1);
    T_EL = reader.GetReal("electrons", "T_EL", -1);
    I_EL = reader.GetReal("electrons", "I_EL", -1);

    // Particle 2 - Ions
    GAS_NAME = reader.Get("ions", "GAS_NAME", "");
    Q_I = reader.GetReal("ions", "Q_I", -1);
    T_I = reader.GetReal("ions", "T_I", -1);
    I_I = reader.GetReal("ions", "I_I", -1);
    // J_I = ETA_PROPELLANT * 4.477962e17 * MASS_FLOW_RATE * Q;
    MACH_I = reader.GetReal("ions", "MACH_I", -1);
    
    M_I = reader.GetReal(GAS_NAME, "M_I", -1);
    N_EXC = (int) reader.GetReal(GAS_NAME, "N_EXC", -1);
    E_IZ = reader.GetReal(GAS_NAME, "E_IZ", -1);
    E_EXC = new double[N_EXC];

    ELS_PATH = reader.Get(GAS_NAME, "ELS_PATH", "");
    E_EXC_PATH = new string[N_EXC];
    IZ_PATH = reader.Get(GAS_NAME, "IZ_PATH", "");
    BS_PATH = reader.Get(GAS_NAME, "BS_PATH", "");
    ISO_PATH = reader.Get(GAS_NAME, "ISO_PATH", "");

    for(int i=0; i<N_EXC; i++){
        E_EXC[i] = reader.GetReal(GAS_NAME, "E_EXC" + to_string(i + 1), -1);
        E_EXC_PATH[i] =  reader.Get(GAS_NAME, "EXC" + to_string(i + 1) + "_PATH", "");
    }

    // Calculated parameters
    DX = L_X / ((double) N_MESH_X - 1);
    DY = L_Y / ((double) N_MESH_Y - 1);
    VOLT_0_NORM = VOLT_0 * Q * pow(DT, 2) / (M_EL * pow(DX, 2));
    VOLT_1_NORM = VOLT_1 * Q * pow(DT, 2) / (M_EL * pow(DX, 2));
    
    VD_I = MACH_I * sqrt(T_EL * Q / M_I);
    N_INJ_EL = I_EL * DT / (Q * N_FACTOR);
    N_INJ_I  = I_I  * DT / (Q * N_FACTOR);
    
    N_INJ_N  = (1 - ETA_PROPELLANT) * (4.477962e17 * MASS_FLOW_RATE) * DT / (N_FACTOR_DSMC);
    
    GAMMA = 2 * N_FACTOR * pow(Q, 2) * pow(DT, 2) / (M_EL * EPS_0 * pow(DX, 2));
    ALPHA = K_SUB * M_EL / M_I;

}
