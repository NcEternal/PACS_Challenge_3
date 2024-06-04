#pragma once
#ifndef HH_Input_Output_HH
#define HH_Input_Output_HH

#include "Parameters.hpp"

namespace laplace {

/* Funtion to read the parameters of the problem from a file */
parameters read(const std::string& filename);

/* Overload of the previous function that lets only rank 0 read the file */
parameters read(const std::string& filename, int rank);

/* Function to print the parameters of the problem on the std::cout */
void print_data(const parameters& params);

/* Overload of the previous function that lets only rank 0 write to std::cout */
void print_data(const parameters& params, int rank);

}

#endif