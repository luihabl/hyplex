
#include "config.h"
#include "INIReader.h"

#include <string>
#include <cmath>

using namespace std;

// INIReader reader("config.ini");

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
int PARTICLE_PER_CELL;
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

// Particle 1 - Electrons
double M_EL;
double Q_EL;
double T_EL;

// Particle 2 - Ions
string GAS_NAME;
double Q_I;
double T_I;
  
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
double N_TOTAL;
double N_FACTOR;
double GAMMA;
double ALPHA;

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

    // Solver parameters
    N_MESH_X = (int) reader.GetReal("solver", "N_MESH_X", -1);
    N_MESH_Y = (int) reader.GetReal("solver", "N_MESH_Y", -1);
    A_X = reader.GetReal("solver", "A_X", -1);
    A_Y = reader.GetReal("solver", "A_Y", -1);
    N_MAX_PARTICLES = (int) reader.GetReal("solver", "N_MAX_PARTICLES", -1);
    PARTICLE_PER_CELL = (int) reader.GetReal("solver", "PARTICLE_PER_CELL", -1);
    N_STEPS = (int) reader.GetReal("solver", "N_STEPS", -1);
    N_AVERAGE = (int) reader.GetReal("solver", "N_AVERAGE", -1);
    K_SUB =  (int) reader.GetReal("solver", "K_SUB", -1);

    // Physical parameters
    L_X = reader.GetReal("global", "L_X", -1);
    double L_Y = reader.GetReal("global", "L_Y", -1);
    // L_Y = 0.5 * L_X * (N_MESH_Y - 1) / (N_MESH_X - 1);
    FREQ = reader.GetReal("global", "FREQ", -1);
    VOLT_0 = reader.GetReal("global", "VOLT_0", -1);
    VOLT_1 = reader.GetReal("global", "VOLT_1", -1);      
    T_NEUTRAL = reader.GetReal("global", "T_NEUTRAL", -1) * K_BOLTZ / Q;
    N_NEUTRAL = reader.GetReal("global", "N_NEUTRAL", -1);
    PLASMA_DENSITY = reader.GetReal("global", "PLASMA_DENSITY", -1);

    // Particle 1 - Electrons
    M_EL = reader.GetReal("electrons", "M_EL", -1);
    Q_EL = reader.GetReal("electrons", "Q_EL", -1);
    T_EL = reader.GetReal("electrons", "T_EL", -1) * K_BOLTZ / Q;

    // Particle 2 - Ions
    GAS_NAME = reader.Get("ions", "GAS_NAME", "");
    Q_I = reader.GetReal("ions", "Q_I", -1);
    T_I = reader.GetReal("ions", "T_I", -1) * K_BOLTZ / Q;
    
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

    // Operational parameters
    DT = 1.0 / (reader.GetReal("global", "DT_DENOMINATOR", -1) * FREQ);  // [benchmark-parameter]
    DX =  L_X / ((double) N_MESH_X - 1);
    DY =  L_Y / ((double) N_MESH_Y - 1);
    VOLT_0_NORM = VOLT_0 * Q * pow(DT, 2) / (M_EL * pow(DX, 2));
    VOLT_1_NORM = VOLT_1 * Q * pow(DT, 2) / (M_EL * pow(DX, 2));
    N_TOTAL = PARTICLE_PER_CELL * (N_MESH_X - 1) * (N_MESH_Y - 1);
    N_FACTOR = PLASMA_DENSITY * L_X * L_Y / N_TOTAL;
    GAMMA = 2 * N_FACTOR * pow(Q, 2) * pow(DT, 2) / (M_EL * EPS_0 * pow(DX, 2));
    ALPHA = K_SUB * M_EL / M_I;
}
