#include <string>

#include "mpi.h"
#include "arg-parser.h"
#include "simulation.h"
#include "configuration.h"

#define DEFAULT_CONFIG_PATH "input/config/config.yaml"

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
    std::string config_path = arg.get("config", DEFAULT_CONFIG_PATH);
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
