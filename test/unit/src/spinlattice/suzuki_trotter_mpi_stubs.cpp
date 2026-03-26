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
#include "atoms.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "spinlattice/internal.hpp"

namespace sld{
   
    void suzuki_trotter_parallel_init(std::vector<double> &x,
                        std::vector<double> &y,
                        std::vector<double> &z,
                        double min_dim[3],
                        double max_dim[3]){

        // Initialises 8 octants if required
        if(sld::internal::c_octants.empty()){
            sld::internal::c_octants.resize(8);
            sld::internal::b_octants.resize(8);
        }

        // Clear existing octant data
        for(int i = 0; i < 8; ++i){
            sld::internal::c_octants[i].clear();
            sld::internal::b_octants[i].clear();
        }

        // For testing: puts atom 0 in first octant
        if(atoms::num_atoms > 0){
            sld::internal::c_octants[0].push_back(0);
        }

        suzuki_trotter_parallel_initialized = true;
    }

} // end namespace sld

namespace vmpi{

    int num_core_atoms = 1;
    int num_bdry_atoms = 0;

    void mpi_init_halo_swap(){}
    void mpi_complete_halo_swap(){}
    void mpi_init_halo_swap_coords(){}
    void mpi_complete_halo_swap_coords(){}
    void barrier(){}

}
