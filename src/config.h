#ifndef CONFIG_H
#define CONFIG_H

/*
	Declaration of configuration global variables. Values are set in config.cpp
*/

#include <string>
using namespace std;

//; ----------------------------- General ---------------------------------------
//[project]
extern string OUTPUT_PATH;
extern string INPUT_PATH;
extern string CROSS_SECTIONS_PATH;

//[simulation]
extern string INITIAL_STATE;

//; ----------------------------- Constants -------------------------------------
//[physical]
extern double Q                      ;//Elementary charge [C]
extern double EPS_0                  ;//Vacuum permittivity [F * m^-1]
extern double PI                     ;//Pi []
extern double K_BOLTZ                ;//Boltzmann constant [J * K^-1]

//; ----------------------------- Domain ----------------------------------------
//[geometry]
extern double L_X                    ;//Lenght of simulation domain in x [m]
extern double L_Y                    ;//Lenght of simulation domain in y [m]
extern int N_MESH_X                  ;//Number of grid points in the x direction (n_benchmark + 1) []
extern int N_MESH_Y                  ;//Number of grid points in the y direction []
extern int N_THRUSTER                ;//Number of nodes that represent the thruster []
extern double A_X                    ;//Grid distorsion coefficient in the x direction []
extern double A_Y                    ;//Grid distorsion coefficient in the y direction []

//[time]
extern double DT                     ;//Time step [s]
extern double FREQ					 ;//Frequency of the signal [Hz]
extern int N_STEPS                   ;//Number of simulation steps to be executed []
extern int N_STEPS_DSMC              ;//Number of DSMC simulation steps to be executed []
extern int N_AVERAGE                 ;//Number of average steps []
extern int N_AVERAGE_DSMC            ;//Number of DSMC average steps []
extern int K_SUB                     ;//Ion subcycling factor []
extern int K_SUB_DSMC                ;//DSMC subcycling factor []

//[boundary-conditions]
extern double VOLT_0                 ;//Voltage on the thruster [V]
extern double VOLT_1                 ;//Voltage on the chamber walls [V]
extern double V_SB 				     ;//Self-bias voltage [V]
extern double V_RF 					 ;//Amplitude of the RF signal [V]
extern double C_CAP					 ;//Capacitance of thruster connection [F]
extern string OB_TYPE                ;//Type of outer boundary condition: neumann, dirichlet

//; ----------------------------- Species ---------------------------------------

//[thruster]
extern double MASS_FLOW_RATE        ;//Mass flow rate of thruster [sccm]
extern double ETA_PROPELLANT        ;//Propellant utilization efficiency []

//[particles]
extern double N_FACTOR              ;//Macro-particle weight factor (i.e. how many particles one macro-particle represents) []
extern double N_FACTOR_DSMC         ;//Macro-particle weight factor for neutrals []
extern int N_MAX_PARTICLES          ;//Maximum number of particles in the array []

//[neutral]
extern string EXPM_NEUTRAL          ;//Expansion model of neutral flow: 'dsmc' or 'constant'
extern string MCC_COLL              ;//Toggle MCC collisions: 'on' or 'off'
extern double T_NEUTRAL             ;//Temperature of neutral gas [eV]
extern double N_NEUTRAL             ;//Density of neutral gas [m^-3]

//[electrons]
extern string INJ_MODEL             ;//Electron injection model: 'constant', 'pulsed', 'balanced' or 'square'
extern double SQR_DUTY_CYCLE        ;//Duty cycle of square pulsed injection []
extern double M_EL                  ;//Electron mass [kg]
extern double Q_EL                  ;//Electron charge [C]
extern double T_EL                  ;//Initial electron temperature [eV]
extern double I_EL                  ;//Electron injection current [A/m]
extern double V_DRIFT_EL            ;//Electron injection drift velocity [m/s]

//[ions]
extern string GAS_NAME              ;//Name of gas: 'helium' or 'xenon'
extern double Q_I                   ;//Ion charge [kg]
extern double T_I                   ;//Initial ion temperature [eV]
extern double I_I                   ;//Ion injection current [A/m]
extern double V_DRIFT_I             ;//Ion injection drift velocity [m/s]

extern double M_I                   ;//Ion mass [kg]
extern double* E_EXC                ;//First excitation energy [eV]
extern double E_IZ                  ;//Ionization energy [eV]
extern int N_EXC                    ;//Number of excited states considered (1 or 2)

extern string ELS_PATH;             ;//Path to elastic cross section
extern string* E_EXC_PATH;          ;//Paths to excitation cross sections
extern string IZ_PATH               ;//Path ionization cross section
extern string BS_PATH               ;//Path backscattering cross section
extern string ISO_PATH              ;//Path isotropic cross section

// Calculated parameters
extern double DX;
extern double DY;
extern double OMEGA_I;
extern int RF_PERIOD_I;
extern double VOLT_0_NORM;
extern double VOLT_1_NORM;
extern double GAMMA;
extern double ALPHA;
extern double K_PHI;
extern double K_Q;
extern double N_INJ_N;
extern double N_INJ_I;
extern double N_INJ_EL;
extern double K_INJ_EL;

void load_config_file(string filename);
#endif
