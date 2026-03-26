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
// Stubs for unit tests for suzuki_trotter()

#include <cstdint>
#include <vector>

#include "material.hpp"
#include "sim.hpp"
#include "sld.hpp"

namespace mp {
   double dt_SI = 1.0e-15;
   double dt = 1.0e-3;
   std::vector<materials_t> material(1);
}

namespace sim {

   uint64_t time = 0;
   uint64_t equilibration_time = 0;
   double temperature = 300.0;
   int hamiltonian_simulation_flags[10] = {0};

   std::vector<double> coarse_noise_field;
   double noise_index = 0.0;
   int M_decimation = 1;
   std::vector<double> atom_idx_z;

   double estimate_cutoff_omega_cdf(double, double) {
      return 1.0;
   }

   void assign_unique_indices(int) {}

   void calculate_random_fields(int, int, double, int, double, int) {}

   double get_noise(const std::vector<double>&, double, int, size_t) {
      return 0.0;
   }
}

namespace mtrandom {
   double gaussian() {
      return 0.0;
   }
}

namespace sld {

   void compute_fields(const int,
                     const int,
                     const std::vector<int>&,
                     const std::vector<int>&,
                     const std::vector<int>&,
                     const std::vector<int>&,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     std::vector<double>&,
                     std::vector<double>&,
                     std::vector<double>&,
                     std::vector<double>&,
                     std::vector<double>&,
                     std::vector<double>&) {}

   void compute_forces(const int,
                     const int,
                     const std::vector<int>&,
                     const std::vector<int>&,
                     const std::vector<int>&,
                     const std::vector<int>&,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     const std::vector<double>&,
                     std::vector<double>&,
                     std::vector<double>&,
                     std::vector<double>&,
                     std::vector<double>&) {}

   double compute_effective_J(const int, const int, std::vector<double>&) {
      return 0.0;
   }

   double compute_effective_C(const int, const int, std::vector<double>&) {
      return 0.0;
   }
}
