
#include "config.h"
#include "INIReader.h"

#include <string>
#include <cmath>
#include <iostream>

using namespace std;

// Project configuration
string OUTPUT_PATH;
string CROSS_SECTIONS_PATH;

// Physical constants
double Q;
double EPS_0;
double PI;
double K_BOLTZ;

// Solver parameters
int N_MESH_X;
int N_MESH_Y;
double A_X;
double A_Y;
int N_MAX_PARTICLES;
int N_STEPS;
int N_AVERAGE;
int K_SUB;

// Physical parameters
double L_X;
double L_Y;
double FREQ;
double VOLT_0;
double VOLT_1;      
double T_NEUTRAL;
double N_NEUTRAL;
double PLASMA_DENSITY;
int N_THRUSTER;

// Particle 1 - Electrons
double M_EL;
double Q_EL;
double T_EL;
double J_EL;
double N_INJ_EL;

// Particle 2 - Ions
string GAS_NAME;
double Q_I;
double MACH_I;
double VD_I;
double T_I;
double J_I;
double N_INJ_I;

// Neutrals
string EXPM_NEUTRAL;
double MASS_FLOW_RATE;
double ETA_PROPELLANT;
double N_FACTOR_NEUTRAL;
double N_INJ_N;

double M_I;
int N_EXC;
double E_IZ;
double* E_EXC;

string ELS_PATH;
string* E_EXC_PATH;
string IZ_PATH;
string BS_PATH;
string ISO_PATH;

// Operational parameters
double DT;
double DX;
double DY;
double VOLT_0_NORM;
double VOLT_1_NORM;
double N_FACTOR;
double GAMMA;
double ALPHA;
double SCCM_2_KGS;

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

    // Boundary conditions
    VOLT_0 = reader.GetReal("boundary-conditions", "VOLT_0", -1);
    VOLT_1 = reader.GetReal("boundary-conditions", "VOLT_1", -1);   
    
    // Plasma
    PLASMA_DENSITY = reader.GetReal("plasma", "PLASMA_DENSITY", -1);
    N_FACTOR = reader.GetReal("particles", "N_FACTOR", -1);
    N_MAX_PARTICLES = (int) reader.GetReal("particles", "N_MAX_PARTICLES", -1);

    
    // ---------------------- species ---------------------------

    // neutral
    T_NEUTRAL = reader.GetReal("neutral", "T_NEUTRAL", -1);
    N_NEUTRAL = reader.GetReal("neutral", "N_NEUTRAL", -1);
    EXPM_NEUTRAL = reader.Get("neutral", "EXPM_NEUTRAL", "");
    MASS_FLOW_RATE = reader.GetReal("thruster", "MASS_FLOW_RATE", -1);
    ETA_PROPELLANT = reader.GetReal("thruster", "ETA_PROPELLANT", -1);
    N_FACTOR_NEUTRAL = reader.GetReal("particles", "N_FACTOR_NEUTRAL", -1);

    // Particle 1 - Electrons
    M_EL = reader.GetReal("electrons", "M_EL", -1);
    Q_EL = reader.GetReal("electrons", "Q_EL", -1);
    T_EL = reader.GetReal("electrons", "T_EL", -1);
    J_EL = reader.GetReal("electrons", "J_EL", -1);

    // Particle 2 - Ions
    GAS_NAME = reader.Get("ions", "GAS_NAME", "");
    Q_I = reader.GetReal("ions", "Q_I", -1);
    T_I = reader.GetReal("ions", "T_I", -1);
    J_I = reader.GetReal("ions", "J_I", -1);
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
    N_INJ_EL = J_EL * DT / (Q * N_FACTOR);
    N_INJ_I  = J_I  * DT / (Q * N_FACTOR);
    
//    SCCM_2_KGS = 4.477962e17; // Set this to 1/M_I if the MFR is given in kg/s
    N_INJ_N  = (1 - ETA_PROPELLANT) * (4.477962e17 * MASS_FLOW_RATE) * DT / (N_FACTOR_NEUTRAL);
    
    GAMMA = 2 * N_FACTOR * pow(Q, 2) * pow(DT, 2) / (M_EL * EPS_0 * pow(DX, 2));
    ALPHA = K_SUB * M_EL / M_I;

}
