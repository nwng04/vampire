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
    namespace sld{
        int interface_quantum_thermostat_test(){
            int ec = 0;
            std::string prefix = "spin-lattice";
            std::string word = "quantum-thermostat";
            std::string value = "true";
            bool expected_value = true;

            bool matched = ::sld::match_input_parameter(prefix, word, value, "", 0);

            if(!matched) return 1;
            if (::sld::internal::use_llgq_thermostat != expected_value) return 1;

            else return 0;
        }
    } // end namespace sld

    int sld_interface_tests(const bool verbose){
        if(verbose) std::cout << "Testing sld interface module" << std::endl;

        int error_count = 0;
        
        error_count += sld::interface_quantum_thermostat_test();

        if(verbose) std::cout << "================================" << std::endl;
        if(error_count == 0) std::cout << " sld interface : PASS " << std::endl;
        else std::cout << " sld interface : FAIL " << error_count << std::endl;
        if(verbose) std::cout << "================================" << std::endl;

        return error_count;
    }
    
} // end namespace ut