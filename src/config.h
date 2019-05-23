#ifndef CONFIG_H
#define CONFIG_H

/*
	Declaration of configuration global variables. Values are set in config.cpp
*/

#include <string>
using namespace std;

// Project configuration
extern string OUTPUT_PATH;
extern string CROSS_SECTIONS_PATH;

// Physical constants
extern double Q;
extern double EPS_0;
extern double PI;
extern double K_BOLTZ;

// Solver parameters
extern int N_MESH_X;          
extern int N_MESH_Y;           
extern double A_X;
extern double A_Y;
extern int N_MAX_PARTICLES;
extern int PARTICLE_PER_CELL;     
extern int N_STEPS;             
extern int N_AVERAGE;           
extern int K_SUB;

// Physical parameters
extern double L_X;
extern double L_Y; 
extern double FREQ;
extern double VOLT_0;
extern double VOLT_1;         
extern double T_NEUTRAL;
extern double N_NEUTRAL;      
extern double PLASMA_DENSITY;
extern int N_THRUSTER; 

// Particle 1 - Electrons
extern double M_EL;
extern double Q_EL;
extern double T_EL;

// Particle 2 - Ions
extern string GAS_NAME;
extern double M_I;
extern double Q_I;
extern double MACH_I;
extern double VD_I;
extern double T_I;  
extern double E_IZ;
extern int N_EXC;
extern double* E_EXC;

extern string ELS_PATH;
extern string IZ_PATH;
extern string BS_PATH;
extern string ISO_PATH;
extern string* E_EXC_PATH;

// Operational parameters
extern double DT; 
extern double DX;
extern double DY;
extern double VOLT_0_NORM;
extern double VOLT_1_NORM;
extern double N_TOTAL;
extern double N_FACTOR;
extern double GAMMA;
extern double ALPHA;

void load_config_file(string filename);
#endif
