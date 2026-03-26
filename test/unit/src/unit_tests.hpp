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

// include header for test functions
#pragma once

namespace ut{

   // simple struct specifying modules to test
   struct module_t {
      bool utility = false;
      bool spinlattice = false;
      bool sld_interface = false;
   };

   // module level functions
   int utility_tests(const bool verbose);
   int spinlattice_tests(const bool verbose);
   //int spinlattice_mpi_tests(const bool verbose);
   int sld_interface_tests(const bool verbose);

   // verification methods
   int bool_error(const bool value, const bool expected_value, const std::string& name);
   int string_error(const std::string& string, const std::string& phrase_to_find, const std::string& name);

}
