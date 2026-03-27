//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//   Stubs for unit tests for suzuki_trotter_mpi()
//------------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>
#include <vector>

// Vampire headers
#include "vmpi.hpp"

namespace vmpi{

    int num_core_atoms = 1;
    int num_bdry_atoms = 0;

    double ComputeTime = 0.0;
    double WaitTime = 0.0;
    double TotalComputeTime = 0.0;
    double TotalWaitTime = 0.0;
    double SwapTimer(double local, double& wait_time){
        (void)wait_time;
        return local;
    }

    void mpi_init_halo_swap(){}
    void mpi_complete_halo_swap(){}
    void mpi_init_halo_swap_coords(){}
    void mpi_complete_halo_swap_coords(){}
    void barrier(){}

}
