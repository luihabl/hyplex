
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>

#include "mpi.h"
#include "arg-parser.h"
#include "simulation.h"

#include "fields.h"
#include "fmatrix.h"
#include "fmath.h"
#include "mcc.h"

#include "particles.h"
#include "particles-in-mesh.h"
#include "rsolver.h"
#include "input-output.h"
#include "num-tools.h"
#include "dsmc.h"
#include "state-info.h"

#include "configuration.h"
#include "diagnostics.h"
#include "clock.h"

#define DEFAULT_CONFIG_PATH "input/config/config.yaml"

using namespace std;
using namespace std::chrono;

// ----------------------------- Main function --------------------------------
int main(int argc, char* argv[])
{
    MPI_Init(&argc,&argv);
    
    mpi_info mpi;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi.size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
    mpi.size_double = (double) mpi.size;

    // ------------------- Loading configuration ------------------------------

    argparser arg(argc, argv);
    string config_path = arg.get("config", DEFAULT_CONFIG_PATH);
    configuration config(arg.get("config", DEFAULT_CONFIG_PATH));
    config.set_job_name(arg.get("name", ""));

    // ------------------- Starting simulation --------------------------------

    simulation sim(mpi, config);
    sim.gen_neutral_field();
    sim.load_state();
    sim.run();

    MPI_Finalize();
	return 0;
}
