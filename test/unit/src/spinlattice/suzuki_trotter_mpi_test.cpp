//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "spinlattice/internal.hpp"
#include "../unit_tests.hpp"

namespace ut{
    namespace spinlattice{

        void sld_setup_mpi(){
            atoms::num_atoms = 1;
            atoms::type_array.assign(1, 0);
            atoms::x_spin_array.assign(1, 1.0);
            atoms::y_spin_array.assign(1, 0.0);
            atoms::z_spin_array.assign(1, 0.0);
            atoms::x_coord_array.assign(1, 0.0);
            atoms::y_coord_array.assign(1, 0.0);
            atoms::z_coord_array.assign(1, 0.0);
            atoms::x_velo_array.assign(1, 0.0);
            atoms::y_velo_array.assign(1, 0.0);
            atoms::z_velo_array.assign(1, 0.0);
            atoms::neighbour_list_start_index.assign(1, 0);
            atoms::neighbour_list_end_index.assign(1, 0);
            atoms::neighbour_list_array.assign(1, 0);

            mp::dt_SI = 1.0e-15;
            mp::dt = 1.0e-3;
            sim::temperature = 300.0;
            sim::time = 0;
            sim::equilibration_time = 1;
            sim::noise_index = 0.0;
            sim::M_decimation = 1;
            sim::atom_idx_z.assign(1, 0.0);
            sim::coarse_noise_field.assign(1, 0.0);

            sld::internal::mp.resize(1);
            sld::internal::x0_coord_array = atoms::x_coord_array;
            sld::internal::y0_coord_array = atoms::y_coord_array;
            sld::internal::z0_coord_array = atoms::z_coord_array;
            sld::internal::forces_array_x.assign(1, 0.0);
            sld::internal::forces_array_y.assign(1, 0.0);
            sld::internal::forces_array_z.assign(1, 0.0);
            sld::internal::fields_array_x.assign(1, 0.0);
            sld::internal::fields_array_y.assign(1, 0.0);
            sld::internal::fields_array_z.assign(1, 0.0);
            sld::internal::potential_eng.assign(1, 0.0);
            sld::internal::sumJ.assign(1, 0.0);
            sld::internal::sumC.assign(1, 0.0);

            sld::internal::x_coord_storage_array.assign(1, 0.0);
            sld::internal::y_coord_storage_array.assign(1, 0.0);
            sld::internal::z_coord_storage_array.assign(1, 0.0);
            sld::internal::c_octants.assign(8, std::vector<int>());
            sld::internal::b_octants.assign(8, std::vector<int>());

            sld::internal::classical_noise_first_call = true;
            sld::internal::spin_noise_generation_first_call = true;
            sld::internal::debug_banner_printed = false;
            sld::suzuki_trotter_parallel_initialized = false;
        }

        int suzuki_trotter_mpi_gaussian_noise_test(){
            sld::internal::use_llgq_thermostat = false;
            sld_setup_mpi();

            std::ostringstream captured;
            std::streambuf* original_buffer = std::cout.rdbuf(captured.rdbuf());

            sld::suzuki_trotter_step_parallel(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::type_array);

            std::cout.rdbuf(original_buffer);
            std::string expected_string = "DEBUG active: entering sld::suzuki_trotter_step_parallel().";
            std::string expected_string_2 = "Adding lattice noise using classical noise.";

            if (string_error(captured.str(), expected_string, "DEBUG not active, please switch on in makefile")) return 1;
            if (string_error(captured.str(), expected_string_2, "Classical noise not called")) return 1;

            return 0;
        }

        int suzuki_trotter_mpi_llgq_noise_test(){
            sld::internal::use_llgq_thermostat = true;
            sld::internal::initialise_noise = false;
            sld_setup_mpi();

            std::ostringstream captured;
            std::streambuf* original_buffer = std::cout.rdbuf(captured.rdbuf());

            sld::suzuki_trotter_step_parallel(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::type_array);

            std::cout.rdbuf(original_buffer);
            std::string expected_string = "DEBUG active: entering sld::suzuki_trotter_step_parallel().";
            std::string expected_string_2 = "'Use_LLGQ_Thermostat' enabled, initialising quantum noise...";

            if (string_error(captured.str(), expected_string, "Quantum noise not initialised")) return 1;
            if (string_error(captured.str(), expected_string_2, "Quantum noise not added to lattice noise")) return 1;

         return 0;
      }

    } // end namespace spinlattice

    int spinlattice_mpi_tests(const bool verbose){
        if(verbose) std::cout << "Testing spin-lattice MPI module" << std::endl;

        int error_count = 0;
        error_count += spinlattice::suzuki_trotter_mpi_gaussian_noise_test();
        error_count += spinlattice::suzuki_trotter_mpi_llgq_noise_test();

        if(verbose) std::cout << "================================" << std::endl;
        if(error_count == 0) std::cout << " spin-lattice MPI : PASS " << std::endl;
        else std::cout << " spin-lattice MPI : FAIL " << error_count << std::endl;
        if(verbose) std::cout << "================================" << std::endl;

        return error_count;
    }

} // end namespace ut
