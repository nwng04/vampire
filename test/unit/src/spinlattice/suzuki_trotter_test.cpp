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

namespace ut{
   namespace spinlattice{

      int bool_error(const bool value, const bool expected_value, const std::string& name){
         if(value == expected_value) return 0;

         std::cout << "FAIL: " << name << " expected " << expected_value << " but got " << value << std::endl;
         return 1;
      }

      int string_error(const std::string& string, const std::string& phrase_to_find, const std::string& name){
         if(string.find(phrase_to_find) != std::string::npos) return 0;

         std::cout << "FAIL: " << name << " expected output containing: " << phrase_to_find << std::endl;
         return 1;
      }

      void sld_setup(){
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
         sld::internal::classical_noise_first_call = true;
         sld::internal::spin_noise_generation_first_call = true;
      }

      int suzuki_trotter_gaussian_lattice_noise_called_test(){
         sld::internal::use_llgq_thermostat = false;
         sld_setup();

         std::ostringstream captured;
         std::streambuf* original_buffer = std::cout.rdbuf(captured.rdbuf());

         int suzuki_trotter_output = sld::suzuki_trotter();

         std::cout.rdbuf(original_buffer);
         std::string expected_string = "DEBUG active: entering sld::suzuki_trotter().";
         std::string expected_string_2 = "Adding lattice noise using classical noise.";

         if (bool_error(suzuki_trotter_output == EXIT_SUCCESS, true, "sld::suzuki_trotter return code")) return 1;
         if (string_error(captured.str(), expected_string, "DEBUG not active, please switch on in makefile")) return 1;
         if (string_error(captured.str(), expected_string_2, "Classical noise not called")) return 1;

         return 0;
      }

      int suzuki_trotter_LLGQ_noise_called_test(){
         sld::internal::use_llgq_thermostat = true;
         sld_setup();

         std::ostringstream captured;
         std::streambuf* original_buffer = std::cout.rdbuf(captured.rdbuf());

         int suzuki_trotter_output = sld::suzuki_trotter();

         std::cout.rdbuf(original_buffer);
         std::string expected_string = "Use_LLGQ_Thermostat' enabled, initialising quantum noise...";
         std::string expected_string_2 = "Adding lattice noise using LLGQ thermostat (quantum noise).";

         if (bool_error(suzuki_trotter_output == EXIT_SUCCESS, true, "sld::suzuki_trotter return code")) return 1;
         if (string_error(captured.str(), expected_string, "Quantum noise not initialised")) return 1;
         if (string_error(captured.str(), expected_string_2, "Quantum noise not added to lattice noise")) return 1;

         return 0;
      }

   } // end namespace spinlattice

   int spinlattice_tests(const bool verbose){
      if(verbose) std::cout << "Testing spin-lattice module" << std::endl;

      int error_count = 0;
      error_count += spinlattice::suzuki_trotter_gaussian_lattice_noise_called_test();
      error_count += spinlattice::suzuki_trotter_LLGQ_noise_called_test();

      if(verbose) std::cout << "================================" << std::endl;
      if(error_count == 0) std::cout << " spin-lattice : PASS " << std::endl;
      else std::cout << " spin-lattice : FAIL " << error_count << std::endl;
      if(verbose) std::cout << "================================" << std::endl;

      return error_count;
   }

} // end namespace ut
