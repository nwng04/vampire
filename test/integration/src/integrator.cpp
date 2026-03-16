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

// module headers
#include "internal.hpp"
#include "sld.hpp"

//------------------------------------------------------------------------------
// Test to verify consistent time evolution for different integrators
//     Note: this is a regression test rather than absolute as some tests
//           are stochastic in nature
//------------------------------------------------------------------------------
bool integrator_test(const std::string integrator, const std::string dir, const double expected_values[], const int filerow, const std::string executable){

   // get root directory
   std::string path = std::filesystem::current_path();

   std::cout << "Testing integrator for " << dir << "   \t: " << std::flush;

   // change directory
   if( !vt::chdir(path+"/data/"+dir) ) return false;

   // run vampire
   int vmp = vt::system(executable);
   if( vmp != 0){
      std::cerr << "Error running vampire. Returning as failed test." << std::endl;
      // Ensure subsequent tests run from the integration root directory.
      vt::chdir(path);
      return false;
   }

   std::string line;

   // open output file
   std::ifstream ifile;
   ifile.open("output");
   if(!ifile.is_open()){
      std::cout << "FAIL | output file not found" << std::endl;
      vt::chdir(path);
      return false;
   }

   // Read a deterministic row from output.
   for(int i=0; i<filerow; i++) getline(ifile, line);

   std::stringstream liness(line);
   double v1 = 0.0;
   double vx = 0.0;
   double vy = 0.0;
   double vz = 0.0;
   double vm = 0.0;

   liness >> v1 >> vx >> vy >> vz >> vm;

   // cleanup
   //vt::system("rm -f output log");

   // return to parent directory
   if( !vt::chdir(path) ) return false;

   // now test value obtained from code
   bool ok = false;

   const double ratiox = vx/expected_values[1];
   const double ratioy = vy/expected_values[2];
   const double ratioz = vz/expected_values[3];
   ok = ratiox >0.99999 && ratiox < 1.00001 && ratioy >0.99999 && ratioy < 1.00001 && ratioz >0.99999 && ratioz < 1.00001;

   if(integrator == "suzuki-trotter"){
      const double tol_t = 1.0e-16;
      const double tol_m = 1.0e-16;

      ok = ok &&
           std::abs(v1 - expected_values[0]) <= tol_t &&
           std::abs(vm - expected_values[4]) <= tol_m;
   }

   if(ok){
      std::cout << "OK" << std::endl;
      return true;
   }
   else{
      std::cout << "FAIL | expected: " << expected_values[0] << "\t" << expected_values[1] << "\t" << expected_values[2] << "\t" << expected_values[3] << "\t" << expected_values[4] << "\tobtained: " << v1 << "\t" << vx << "\t" << vy << "\t" << vz << "\t" << vm << "\t" << line << std::endl;
      return false;
   }

}
